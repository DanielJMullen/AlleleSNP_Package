# 02-get genotype by querying vcf files
## ==================================================================================
## Function: get_vcf_info_main
###  Aim: To get the genotypes of input SNPs by vcf files
###  INPUT: input snp file (.csv)
###  OUTPUT: input snp file with information of genotypes of vcf files
## ==================================================================================

# outline ---------------------------------------------------------------------------------------------------------

# # 1. read snp file into data frame
# # snp_list should contain both the df and GRange object
#
# # 2. read peak files into a list of GRange object
#
# # 3. add overlapping info to original snp table

# 0. load packages ------------------------------------------------------------------------------------------------
# [deprecated]


# toolbox: get short names

get_short_fileName = function(fileName_fullPath) {
        str_vec = strsplit(fileName_fullPath, "/")[[1]]
        str_vec[length(str_vec)]
}

# 1. generate output file name ------------------------------------------------------------------------------------


gen_output_file_vcfInfo = function(snp_info_file, output_dir = "./", sample_name = "") {
        snp_batch_id = gsub("\\.[^.]*$", "", basename(snp_info_file))
        sample_text = ifelse(sample_name == "", "", paste0("_", sample_name))
        if (grepl("\\.tsv$", snp_info_file, ignore.case = T)) {
                output_file = paste0(output_dir, "/", snp_batch_id, sample_text, "_genotypeInfo.tsv")
        } else {
                output_file = paste0(output_dir, "/", snp_batch_id, sample_text, "_genotypeInfo.csv")
        }

        return(output_file)
}


# 2. read vcf files into a list -----------------------------------------------------------------------------------


gen_vcf_list = function(vcf_dir = NA, vcf_file = NA, snp_info_gr) {
        if (!is.na(vcf_dir)) {
                vcf_files = list.files(
                        vcf_dir,
                        pattern = "\\.vcf(\\.gz)?$",
                        full.names = T
                )
                vcf_files_short = list.files(
                        vcf_dir,
                        pattern = "\\.vcf(\\.gz)?$",
                        full.names = F
                )
        }
        if (!is.na(vcf_file)) {
                vcf_files = vcf_file
                vcf_files_short = get_short_fileName(vcf_file)
        }
        if (!exists("vcf_files") || length(vcf_files) == 0) {
                stop("[ERROR] No .vcf or .vcf.gz files found")
        }
        vcf_list = vector("list", length(vcf_files))
        snp_param = VariantAnnotation::ScanVcfParam(which = snp_info_gr)
        for (i in seq_along(vcf_files)) {
                cat("reading", vcf_files[i], "\n")
                vcf_list[[i]] = read_vcfFile(vcf_files[i], snp_param)
        }
        names(vcf_list) = vcf_files_short
        return(vcf_list)
}


read_vcfFile = function(vcf_file, param) {
        if (grepl("\\.vcf\\.gz$", vcf_file, ignore.case = T)) {
                compressVcf = vcf_file
        } else {
                compressVcf = Rsamtools::bgzip(vcf_file, tempfile())
        }
        idx = Rsamtools::indexTabix(compressVcf, "vcf")
        tab = Rsamtools::TabixFile(compressVcf, idx)
        vcf = VariantAnnotation::readVcf(tab, "hg19", param)
        return(vcf)
}


.process_vcf_file = function(vcf_file, vcf_file_short, snp_param) {
        cat("reading", vcf_file, "\n")

        vcf_data = read_vcfFile(vcf_file, snp_param)
        vcf_software = get_vcf_software(vcf_file_short)
        vcf_data_het = sel_hetSnps_vcf(vcf_data)
        vcf_data_df = process_vcf_data(vcf_data_het, vcf_software)

        return(vcf_data_df)
}


# 3. process vcf list data -------------------------------------------------------------------------------------
process_vcf_list = function(vcf_list) {
# To process vcf list
# Output: to get a processed df for each vcf elemnt in the list
# - row: SNPs
# - col: seqnames start end width strand paramRangeID REF ALT QUAL FILTER ref_count alt_count
# helper functions: get_vcf_software, sel_hetSnps_vcf, process_vcf_data

        vcf_df_list = replicate(length(vcf_list), list())
        names(vcf_df_list) = names(vcf_list)

        for (i in 1:length(vcf_list)) {
                vcf_data_i = vcf_list[[i]]
                # get the software that generates the vcf file
                vcf_software_i = get_vcf_software (names(vcf_list)[[i]])
                # subset vcf with het genotypes
                vcf_data_i_het = sel_hetSnps_vcf(vcf_data_i)
                # process vcf data to data.frame format
                vcf_data_i_df = process_vcf_data(vcf_data_i_het, vcf_software_i)
                # assign vcf_df_list with processed data.frame
                vcf_df_list[[i]] = vcf_data_i_df
        }

        # return to our df
        return (vcf_df_list)
}


# 3-1 helper function
get_vcf_software = function(vcf_fileName) {
# Aim: to get the software that was used to generate the vcf file
        if (grepl("[Gg][Aa][Tt][Kk]", vcf_fileName)) {
                return ("GATK")
        }
        if (grepl("[Ss][Aa][Mm][Tt][Oo][Oo][Ll][Ss]", vcf_fileName)) {
                return ("Samtools")
        }
        if (grepl("[Bb][Ii][Ss][Ss][Nn][Pp]", vcf_fileName)) {
                return ("BisSNP")
        }
        else {
                stop("Please add the software name you used in generating vcf file to your vcf file name.
                     for example: A549.GATK.vcf")
        }
}

# 3-2 helper function
sel_hetSnps_vcf = function(vcf_data){
# Aim: to select heterozygous SNPs from vcf data

        genotype = VariantAnnotation::geno(vcf_data)$GT
        vcf_data_het = vcf_data[genotype[,1] == "0/1"]

        return(vcf_data_het)
}

# 3-3 helper function
process_vcf_data = function(vcf_data, vcf_software = "GATK") {
# Aim: to transform vcf data to data.frame format
# Input: vcf_data, (vcf annotation format)
# Output: vcf data frame, with genotype, counts of ref and alt alleles added
# helper function: extract_allelic_dist

        # calculate ref and alt allele counts
        vcf_allelic_dist = extract_allelic_dist(vcf_data, vcf_software)

        ref_count_vec = sapply(vcf_allelic_dist, function(x) {
                as.numeric(strsplit(x, ",")[[1]][1])
        })
        alt_count_vec = sapply(vcf_allelic_dist, function(x) {
                as.numeric(strsplit(x, ",")[[1]][2])
        })

        # add ref and alt allele counts information
        vcf_data_gr = SummarizedExperiment::rowRanges(vcf_data)
        names(vcf_data_gr) = 1:length(vcf_data_gr)
        # vcf_data_df = as.data.frame(vcf_data_gr) # This might not work in R package. Don't know why
        vcf_data_df = data.frame(seqnames = vcf_data_gr@seqnames, vcf_data_gr@ranges, vcf_data_gr@elementMetadata)
        vcf_data_df$ref_count = ref_count_vec
        vcf_data_df$alt_count = alt_count_vec

        # modify the data format (from DNAString to character)
        vcf_data_df$ALT = sapply(vcf_data_df$ALT, function(x) as.character(unlist(x)))
        # vcf_data_df$ALT = sapply(vcf_data_df$ALT, function(x) toString(x))

        return(vcf_data_df)
}

# 3-3-1 helper function
extract_allelic_dist = function(vcf_data, vcf_software = "GATK"){
# Aim: to extract allelic distribution in vcf file generated by various softwares

        if (vcf_software == "GATK") {
                AD = VariantAnnotation::geno(vcf_data)$AD
                AD_1 = gsub("c\\(","",paste(AD))
                AD_2 = gsub("\\)","", AD_1)
                AD_3 = gsub(":",", ",AD_2)
                AD_4 = gsub(", ", ",", AD_3)
                return(AD_4)
        }

        if (vcf_software == "Samtools") {
                DP4 = VariantAnnotation::info(vcf_data)$DP4
                AD = sapply(DP4, function(x) paste0(x[1] + x[2], ",", x[3] + x[4]))
                return(AD)
        }

        if (vcf_software == "BisSNP") {
                DP4 = VariantAnnotation::geno(vcf_data)$DP4[, , c(1, 2, 3, 4)]
                AD = apply(DP4, 1, function(x){paste0(x[1]+x[2],',',x[3]+x[4])})
                return(AD)
        }
}

# 4. to add genotype information to snp info data frame -----------------------------------------------------------

add_vcf_info = function(snp_info_df, vcf_df_list) {
# Aim: to get genotype information added

        # add snp genotype information
        for (i in 1:length(vcf_df_list)) {
                vcf_df_i = vcf_df_list[[i]]
                vcf_file_i = names(vcf_df_list)[i]
                # select SNPs with het genotypes
                sel_rows = as.character(snp_info_df$rsID) %in% as.character(vcf_df_i$paramRangeID)

                # add vcf information
                snp_info_df[, vcf_file_i] = sel_rows
                snp_info_df[sel_rows, paste0(vcf_file_i, "_ref_count")] = vcf_df_i$ref_count
                snp_info_df[sel_rows, paste0(vcf_file_i, "_alt_count")] = vcf_df_i$alt_count
        }

        return(snp_info_df)
}



# main function ---------------------------------------------------------------------------------------------------

get_vcf_info_main = function(
        snp_info_file,
        vcf_dir = NA,
        vcf_file = NA,
        output_dir = "./",
        sample_name = "",
        output_file = NA,
        chromosome_annotation = "chr",
        n_cores = 1
) {

        # 0. load packages
        cat("get vcf information for SNPs ... \n")

        # 0. generate output file name
        if (is.na(output_file)) {
                output_file = gen_output_file_vcfInfo(snp_info_file = snp_info_file,
                                                      output_dir = output_dir,
                                                      sample_name = sample_name)
        }
        cat("    output file name:", output_file, '\n')
        # 1. read snp file into data frame
        snp_info_list = read_inputSNP_file(
                snp_info_file = snp_info_file,
                chromosome_annotation = chromosome_annotation
        )

        # 2. read and process vcf data
        if (n_cores == 1) {
                vcf_list = gen_vcf_list(vcf_dir = vcf_dir,
                                        vcf_file = vcf_file,
                                        snp_info_gr = snp_info_list$snp_info_gr)

                # 3. process vcf data
                vcf_df_list = process_vcf_list(vcf_list)
        } else {
                if (!is.na(vcf_dir)) {
                        vcf_files = list.files(
                                vcf_dir,
                                pattern = "\\.vcf(\\.gz)?$",
                                full.names = T
                        )
                        vcf_files_short = list.files(
                                vcf_dir,
                                pattern = "\\.vcf(\\.gz)?$",
                                full.names = F
                        )
                }
                if (!is.na(vcf_file)) {
                        vcf_files = vcf_file
                        vcf_files_short = get_short_fileName(vcf_file)
                }
                if (!exists("vcf_files") || length(vcf_files) == 0) {
                        stop("[ERROR] No .vcf or .vcf.gz files found")
                }

                snp_param = VariantAnnotation::ScanVcfParam(
                        which = snp_info_list$snp_info_gr
                )

                vcf_df_list = parallel::mclapply(
                        seq_along(vcf_files),
                        function(i) {
                                .process_vcf_file(
                                        vcf_file = vcf_files[i],
                                        vcf_file_short = vcf_files_short[i],
                                        snp_param = snp_param
                                )
                        },
                        mc.cores = min(n_cores, length(vcf_files))
                )
                names(vcf_df_list) = vcf_files_short
        }

        # 4. add overlapping info to original snp table
        snp_info_addVcf_df = add_vcf_info (snp_info_df = snp_info_list$snp_info_df,
                                           vcf_df_list = vcf_df_list)

        if (!identical(output_file, F)) {
                if (grepl("\\.tsv$", output_file, ignore.case = T)) {
                        write.tsv0(snp_info_addVcf_df, output_file)
                } else {
                        write.csv0(snp_info_addVcf_df, output_file)
                }
        }

        cat("vcf information added ... \n")
        return(list(snp_info_addVcf_df = snp_info_addVcf_df, output_file = output_file))
}

