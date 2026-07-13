# F2.5-get_encodeCnv_info
## ==================================================================================
## Function: get_encodeCnv_info_main
###  Aim: To get the cnv of input SNPs by vcf/encode cnv data files
###  INPUT: input snp file (.csv), vcf/encode cnv data
###  OUTPUT: input snp file with information of cnv
## ==================================================================================

# helper functions ------------------------------------------------------------------------------------------------

# generate snp-cnv table
# Aim: to generate snp-cnv table: rows - snps; cols - encode samples;
gen_snp_cnv_table = function(snp_info_df, snp_info_gr) {

        cat("Add ENCODE cnv information...\n")
        # encode_cnv_gr_list = readRDS(cnv_grList_file)
        utils::data("encode_cnv_gr_list", envir = environment())
        snp_cnv_table = data.frame(matrix(nrow = length(snp_info_gr), ncol = length(encode_cnv_gr_list)))
        colnames(snp_cnv_table) = names(encode_cnv_gr_list)

        for (i in 1:length(encode_cnv_gr_list)) {
                overlaps_df_i = as.data.frame(GenomicRanges::findOverlaps(snp_info_gr, encode_cnv_gr_list[[i]]))
                cnv_snp_info_i = encode_cnv_gr_list[[i]][overlaps_df_i$subjectHits]$name
                snp_cnv_table[overlaps_df_i$queryHits, i] = cnv_snp_info_i
        }

        # merge data frame
        snp_cnv_df = data.frame(snp_info_df, snp_cnv_table)

        return(snp_cnv_df)
}


# main function ---------------------------------------------------------------------------------------------------

get_encodeCnv_info_main = function(
        snp_info_file,
        output_dir = "./",
        sample_name = "",
        output_file = NA,
        chromosome_annotation = "chr"
) {

        # 1. read snp info file
        snp_info_list = read_inputSNP_file(
                snp_info_file = snp_info_file,
                chromosome_annotation = chromosome_annotation
        )

        # 2. generate snp-cnv table (row: snps, col: all encode cnv samples)
        snp_info_addCnv_df = gen_snp_cnv_table(snp_info_df = snp_info_list$snp_info_df,
                                               snp_info_gr = snp_info_list$snp_info_gr)
        # head(snp_cnv_df)

        # 3. write down snp-cnv table
        if (is.na(output_file)) {
                output_file = gen_output_file_cnvInfo(snp_info_file = snp_info_file,
                                                      output_dir = output_dir,
                                                      sample_name = sample_name)
        }
        cat("    output file name:", output_file, '\n')

        if (!identical(output_file, F)) {
                if (grepl("\\.tsv$", output_file, ignore.case = T)) {
                        write.tsv0(snp_info_addCnv_df, output_file)
                } else {
                        write.csv0(snp_info_addCnv_df, output_file)
                }
        }
        cat("cnv information added ... \n")

        return(list(snp_info_addCnv_df = snp_info_addCnv_df, output_file = output_file))
}
