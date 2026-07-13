
# 1. prepare SNP files --------------------------------------------------------------------------------------------

# 1-1 prepare SNP files for cnv inference
prepare_files_for_cnvInfer = function(
        index_snp_file, r2_cutoff = 0.5, sample_ethic = "EUR",
        output_dir = "./", sample_name = "",
        vcf_file = "data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf",
        chromosome_annotation = "chr"
) {
        # 1. generate ld snp info df from index snp file
        # 2. add wgs info
        # Input: 1. index snp file ; 2. the ethnicity of the sample (not the SNPs)

        ldsnp_info_list = get_ldsnp_info(
                index_snp_file = index_snp_file,
                population = sample_ethic,
                output_dir = output_dir,
                r2_cutoff = r2_cutoff,
                for_cnv_call = T
        )

        ldsnp_info_addVcf_list = get_vcf_info_main(
                snp_info_file = ldsnp_info_list$output_file,
                vcf_file = vcf_file,
                output_dir = output_dir,
                sample_name = sample_name,
                chromosome_annotation = chromosome_annotation
        )

        return (ldsnp_info_addVcf_list$snp_info_addVcf_df)

}

# 1-2. add locus information
add_locus_info = function(snp_info_addVcf_df) {
        # Aim: to add locus information (index SNPs in high LD are considered as one locus)

        # generate {indexSNP object} for each index SNPs
        index_snps = unique(as.character(snp_info_addVcf_df$query_snp))
        snp_list = vector("list", length(index_snps))
        names(snp_list) = index_snps

        for (i in 1:length(index_snps)) {
                snp_info_addVcf_df_sel_indexSnp = snp_info_addVcf_df[snp_info_addVcf_df$query_snp == index_snps[i], ]
                snp_list[[i]] = new("indexSNP", rsID = index_snps[i],
                                    ldSNPs = as.character(snp_info_addVcf_df[snp_info_addVcf_df$rsID %in% snp_info_addVcf_df$query_snp &
                                                                        snp_info_addVcf_df$query_snp == index_snps[i], "rsID"]),
                                    ldInfo = snp_info_addVcf_df[snp_info_addVcf_df$rsID %in% snp_info_addVcf_df$query_snp &
                                                                        snp_info_addVcf_df$query_snp == index_snps[i], "D."])
        }
        # calculate locus and locus D'
        index_snp_mat = matrix(nrow = length(snp_list), ncol = length(snp_list))

        for (i in 1:length(snp_list)) {
                for (j in 1:length(snp_list)) {
                        index_snp_mat[i, j] = ldRelation(snp_list[[i]], snp_list[[j]])
                }
        }

        index_snp_info_df = data.frame(query_snp = index_snps)
        rownames(index_snp_info_df) = index_snps
        index_snp_info_df$locus = NA
        index_snp_info_df$locus_D = NA

        for (i in 1:nrow(index_snp_info_df)) {
                locus_snp_idx = which(!is.na(index_snp_mat[i, ]))[1]
                index_snp_info_df$locus[i] = index_snps[locus_snp_idx]
                index_snp_info_df$locus_D[i] = index_snp_mat[i, locus_snp_idx]
        }

        # add locus information
        snp_info_addVcf_df$locus = index_snp_info_df[as.character(snp_info_addVcf_df$query_snp), "locus"]
        snp_info_addVcf_df$locus_D = index_snp_info_df[as.character(snp_info_addVcf_df$query_snp), "locus_D"]

        return (list(snp_info_addVcf_df = snp_info_addVcf_df,
                     index_snp_info_df = index_snp_info_df))
}


# 1-3. modify vcf information for further analysis
mod_vcf_info = function(snp_info_addVcf_df) {
        # Aim: to modify vcf df to make it more suitable for further analysis

        snp_info_addVcf_df$allele_1_count = NA
        snp_info_addVcf_df$allele_2_count = NA

        # 1. change column names
        names(snp_info_addVcf_df)[grepl("ref_count", names(snp_info_addVcf_df))] = "ref_count"
        names(snp_info_addVcf_df)[grepl("alt_count", names(snp_info_addVcf_df))] = "alt_count"
        # head(snp_info_addVcf_df)

        # 2. sort df by query snp / locus
        if ("locus" %in% names(snp_info_addVcf_df)) {
                snp_info_addVcf_df = dplyr::arrange(snp_info_addVcf_df, locus)
                D_prime_sign = as.numeric(as.character(snp_info_addVcf_df$D.)) * snp_info_addVcf_df$locus_D > 0
        } else {
                snp_info_addVcf_df = dplyr::arrange(snp_info_addVcf_df, query_snp)
                D_prime_sign = as.numeric(as.character(snp_info_addVcf_df$D.)) > 0
        }

        # 3. modify vcf information by D'
        snp_info_addVcf_df$allele_1_count[D_prime_sign] = snp_info_addVcf_df$ref_count[D_prime_sign]
        snp_info_addVcf_df$allele_2_count[D_prime_sign] = snp_info_addVcf_df$alt_count[D_prime_sign]

        snp_info_addVcf_df$allele_1_count[!D_prime_sign] = snp_info_addVcf_df$alt_count[!D_prime_sign]
        snp_info_addVcf_df$allele_2_count[!D_prime_sign] = snp_info_addVcf_df$ref_count[!D_prime_sign]

        return(snp_info_addVcf_df)
}


# 2. get loci cnv raw ---------------------------------------------------------------------------------------------
get_het_locus_summary_df = function(snp_info_addVcf_df, sample_ethic = "EUR", r2_cutoff = 0.5,
                                    distance_threshold = 50, min_ldsnp_num = 3, read_count_cutoff = 500,
                                    index_snp_info_df) {
        # To calculate CNV number for each locus
        # step 1: collect distritbution total SNPs for each locus

        # calculate het-SNP count by LD
        snp_info_addVcf_df_narm = snp_info_addVcf_df[complete.cases(snp_info_addVcf_df) &
                                                             snp_info_addVcf_df$r2 >= r2_cutoff, ]
        snp_info_addVcf_df_narm_mk = mark_close_snp(snp_info_addVcf_df_narm, distance_threshold = distance_threshold)
        snp_info_addVcf_df_narm_rm_close_snp = dplyr::filter(snp_info_addVcf_df_narm_mk, for_cnv_calculation == "Y")
        het_snp_summary_df_locus = dplyr::summarise(
                dplyr::group_by(snp_info_addVcf_df_narm_rm_close_snp, locus),
                sum_allele_1_count = sum(allele_1_count),
                sum_allele_2_count = sum(allele_2_count),
                .groups = "drop"
        )
        # filter low count het SNPs
        sel_high_count_locus = table(snp_info_addVcf_df_narm_rm_close_snp$locus) > min_ldsnp_num # an input_SNP must have >= 3 het SNPs in LD
        sel_locus = names(table(snp_info_addVcf_df_narm_rm_close_snp$locus)[sel_high_count_locus])
        het_snp_summary_df_locus = dplyr::filter(
                het_snp_summary_df_locus,
                sum_allele_1_count + sum_allele_2_count > read_count_cutoff & locus %in% sel_locus
        )

        # add SNPs in the same locus
        het_snp_summary_df = add_locus_snp(het_snp_summary_df_locus, index_snp_info_df)
        return (het_snp_summary_df)
}

# 2-1. Mark SNPs that are too close
mark_close_snp = function(snp_info_addVcf_df, distance_threshold = 50) {
        snp_info_addVcf_df$for_cnv_calculation = NA
        for (i in 1: nrow(snp_info_addVcf_df)) {
                if (i == 1)  {
                        snp_info_addVcf_df$for_cnv_calculation[i] = "Y"
                } else if (snp_info_addVcf_df$chr[i] == snp_info_addVcf_df$chr[i-1] &
                           abs(snp_info_addVcf_df$pos[i] - snp_info_addVcf_df$pos[i-1]) < distance_threshold ) {
                        snp_info_addVcf_df$for_cnv_calculation[i] = "N"
                } else {
                        snp_info_addVcf_df$for_cnv_calculation[i] = "Y"
                }
        }

        return (snp_info_addVcf_df)
}

# 2-2. add SNPs in the same locus
add_locus_snp = function(het_snp_summary_df_locus, index_snp_info_df) {
        index_snp_info_df$sum_allele_2_count =  index_snp_info_df$sum_allele_1_count = NA
        for (i in 1:nrow(index_snp_info_df)) {
                locus_i = index_snp_info_df$locus[i]
                if (locus_i %in% as.character(het_snp_summary_df_locus$locus)) {
                        j = which(as.character(het_snp_summary_df_locus$locus) == locus_i)
                        if (index_snp_info_df$locus_D[i] > 0) {
                                index_snp_info_df$sum_allele_1_count[i] = het_snp_summary_df_locus$sum_allele_1_count[j]
                                index_snp_info_df$sum_allele_2_count[i] = het_snp_summary_df_locus$sum_allele_2_count[j]
                        } else {
                                index_snp_info_df$sum_allele_1_count[i] = het_snp_summary_df_locus$sum_allele_2_count[j]
                                index_snp_info_df$sum_allele_2_count[i] = het_snp_summary_df_locus$sum_allele_1_count[j]
                        }
                }
        }

        het_snp_summary_df_querySNP = index_snp_info_df[complete.cases(index_snp_info_df), ]
        rownames(het_snp_summary_df_querySNP) = 1:nrow(het_snp_summary_df_querySNP)

        return (het_snp_summary_df_querySNP)
}

# get CNV info ---------------------------------------------------------------------------------------------------

get_loci_cnv_byEthic = function(
        index_snp_file, cnv_param_list,
        sample_name = "", sample_ethic = "EUR",
        output_dir = "./", chromosome_annotation = "chr"
) {
# Aim: to get loci cnv information of a particular ethic group
# Input: index snp file, vcf file, ethic
# Output: het_snp_summary_df_ethic (summarized info of loci cnv based on a particular ethic group)

        # 1. process snp data
        ## 1) get ld snp table; 2) add vcf information
        snp_info_addVcf_df_raw = prepare_files_for_cnvInfer(
                index_snp_file = index_snp_file,
                r2_cutoff = cnv_param_list$r2_cutoff,
                sample_name = sample_name,
                sample_ethic = sample_ethic,
                output_dir = output_dir,
                vcf_file = cnv_param_list$vcf_file_for_cnv,
                chromosome_annotation = chromosome_annotation
        )
        ## modify vcf information for further analysis
        snp_info_addVcf_list = add_locus_info(snp_info_addVcf_df_raw)
        snp_info_addVcf_df = mod_vcf_info(snp_info_addVcf_list$snp_info_addVcf_df)

        # 2. get loci cnv by calculating summing up r2 SNP distribution
        het_snp_summary_df_r2 = get_het_locus_summary_df(
                snp_info_addVcf_df = snp_info_addVcf_df,
                sample_ethic = sample_ethic,
                r2_cutoff = cnv_param_list$r2_cutoff,
                distance_threshold = cnv_param_list$distance_threshold,
                min_ldsnp_num = cnv_param_list$min_ldsnp_num,
                read_count_cutoff = cnv_param_list$read_count_cutoff,
                index_snp_info_df = snp_info_addVcf_list$index_snp_info_df
        )

        return(het_snp_summary_df_r2)

}

