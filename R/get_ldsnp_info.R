# F1.1-get ld snp information
## ==================================================================================
## Function: get_ldsnp_info_main
###  Aim: To get HaploReg files from an input SNP file which contains a few SNPs and its corresponding risky population
###  INPUT: index SNP list (a .csv file), with the 1st col the rsID and the 2nd col Population
###  OUTPUT: information table of index + LD SNP from HaploReg
## ==================================================================================
# ^ OUTDATED, CHANGE LATER


# main function ---------------------------------------------------------------------------------------------------


# Get HaploReg output for idx SNPs and assoc LD SNPs
# Uses Daniel's httr:POST() workflow
get_ldsnp_info = function(
        index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv",
        population = NA,
        r2_cutoff = 0.5,
        for_cnv_call = F,
        output_dir = "./",
        output_file = NA,
        failed_snp_file = NA,
        haploreg_url = "https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php",
        haploreg_timeout = 180,
        sleep_range = c(3, 6)
) {

        # fmt output
        cat("[INFO] Getting high-LD SNP info from index SNPs ... \n")
        if (!dir.exists(output_dir)) {
                dir.create(output_dir, recursive = T)
        }
        if (.is_na_scalar(output_file)) {
                output_file = .gen_output_file_ldsnp(
                        index_snp_file = index_snp_file,
                        population = population,
                        r2_cutoff = r2_cutoff,
                        for_cnv_call = for_cnv_call,
                        output_dir = output_dir
                )
        }
        cat("[INFO] Output file:", output_file, "\n")
        if (.is_na_scalar(failed_snp_file)) {
                if (identical(output_file, F)) {
                        failed_snp_file = F
                } else {
                        failed_snp_file = sub("\\.[^.]+$", "_haploreg_failed_snps.txt", output_file)
                }
        }
        if (!identical(failed_snp_file, F)) {
                cat("[INFO] Failed HaploReg SNP file:", failed_snp_file, "\n")
        }

        # read input SNP file, normalize to 2 cols: rsID, population
        index_snp_df = .read_index_snp_file(index_snp_file = index_snp_file, population = population)
        is_rsID = !is.na(index_snp_df$rsID) & substr(index_snp_df$rsID, 1, 2) == "rs"
        non_rsID_SNPs = index_snp_df[!is_rsID, "rsID"]
        index_snp_df = index_snp_df[is_rsID, ]
        if (nrow(index_snp_df) == 0) {
                warning("[WARN] No rsID-formatted SNPs found in index_snp_file")
                ldsnp_info_df = .empty_haploreg_result()
                if (!identical(output_file, F)) {
                        .write_ld_output(ldsnp_info_df, output_file)
                }
                if (!identical(failed_snp_file, F)) {
                        .write_failed_haploreg_snps(failed_snps = character(), failed_snp_file = failed_snp_file)
                }
                return(list(
                        ldsnp_info_df = ldsnp_info_df,
                        output_file = output_file,
                        failed_snp_file = failed_snp_file,
                        failed_haploreg_snp_IDs = character(),
                        non_rsID_SNPs = non_rsID_SNPs
                ))
        }

        # query HaploReg
        # bad responses return empty df (for now)
        haploreg_list = mapply(
                FUN = .get_haploreg_snp,
                snp_id = index_snp_df$rsID,
                population_string = index_snp_df$population,
                MoreArgs = list(
                        r2_cutoff = r2_cutoff,
                        haploreg_url = haploreg_url,
                        haploreg_timeout = haploreg_timeout,
                        sleep_range = sleep_range
                ),
                SIMPLIFY = F
        )
        failed_haploreg_snp_IDs = unique(unlist(
                lapply(seq_along(haploreg_list), function(i) {
                        failed = isTRUE(attr(haploreg_list[[i]], "haploreg_failed"))
                        if (failed) {
                                index_snp_df$rsID[i]
                        } else {
                                character()
                        }
                }),
                use.names = F
        ))
        failed_haploreg_snp_IDs = failed_haploreg_snp_IDs[!is.na(failed_haploreg_snp_IDs) & failed_haploreg_snp_IDs != ""]
        if (!identical(failed_snp_file, F)) {
                .write_failed_haploreg_snps(failed_snps = failed_haploreg_snp_IDs, failed_snp_file = failed_snp_file)
        }

        if (length(haploreg_list) == 0 || all(vapply(haploreg_list, nrow, numeric(1)) == 0)) {
                warning("[WARN] (0) LD SNPs returned by HaploReg.")
                ldsnp_info_df = .empty_haploreg_result()
        } else {
                haploreg_df_raw = do.call(rbind, haploreg_list)
                haploreg_df_raw[] = lapply(haploreg_df_raw, function(x) {
                        if (is.character(x)) {
                                gsub("[\r\n\t]", "", x)
                        } else {
                                x
                        }
                })
                if (for_cnv_call) {
                        haploreg_df = haploreg_df_raw
                } else {
                        haploreg_df = .combine_duplicate_ld_snps(haploreg_df_raw, index_snps = index_snp_df$rsID)
                }

                # remove SNPs that didn't produce usable HaploReg LD vals
                non_haploreg_snp_IDs = haploreg_df[is.na(haploreg_df$r2) | haploreg_df$r2 == "" | grepl("NA", haploreg_df$r2), "rsID"]
                haploreg_df = haploreg_df[!(is.na(haploreg_df$r2) | haploreg_df$r2 == "" | grepl("NA", haploreg_df$r2)),]

                # remove indels and SNPs with multiple alternate alleles.
                indel_or_multiple_alt_snp_IDs = haploreg_df[is.na(haploreg_df$ref) | is.na(haploreg_df$alt) | nchar(as.character(haploreg_df$ref)) > 1 | nchar(as.character(haploreg_df$alt)) > 1, "rsID"]
                haploreg_df = haploreg_df[!(is.na(haploreg_df$ref) | is.na(haploreg_df$alt) | nchar(as.character(haploreg_df$ref)) > 1 | nchar(as.character(haploreg_df$alt)) > 1),]

                # remove SNPs with no chr or pos info
                no_position_noted_snp_IDs = haploreg_df[ is.na(haploreg_df$chr) | is.na(haploreg_df$pos_hg38) | haploreg_df$chr == "" | haploreg_df$pos_hg38 == "", "rsID"]
                haploreg_df = haploreg_df[!(is.na(haploreg_df$chr) | is.na(haploreg_df$pos_hg38) | haploreg_df$chr == "" | haploreg_df$pos_hg38 == ""),]

                # reformat df
                names(haploreg_df)[names(haploreg_df) == "pos_hg38"] = "pos"
                leading_cols = c("rsID", "chr", "pos", "ref", "alt")
                remaining_cols = setdiff(names(haploreg_df), leading_cols)
                ldsnp_info_df = haploreg_df[, c(leading_cols, remaining_cols), drop = F]
                chr_sort = suppressWarnings(as.numeric(gsub("chr", "", ldsnp_info_df$chr)))
                if (!all(is.na(chr_sort))) {
                        ldsnp_info_df = ldsnp_info_df[order(chr_sort, as.numeric(ldsnp_info_df$pos)), ]
                }
        }

        if (!exists("non_haploreg_snp_IDs")) non_haploreg_snp_IDs = character()
        if (!exists("indel_or_multiple_alt_snp_IDs")) indel_or_multiple_alt_snp_IDs = character()
        if (!exists("no_position_noted_snp_IDs")) no_position_noted_snp_IDs = character()
        if (!identical(output_file, F)) { # if output_file == F, do not write down file
                .write_ld_output(ldsnp_info_df, output_file)
        }
        cat("[INFO] High-LD SNPs added ... \n")
        return(list(
                ldsnp_info_df = ldsnp_info_df,
                output_file = output_file,
                failed_snp_file = failed_snp_file,
                failed_haploreg_snp_IDs = unique(failed_haploreg_snp_IDs),
                non_rsID_SNPs = unique(non_rsID_SNPs),
                non_haploreg_snp_IDs = unique(non_haploreg_snp_IDs),
                indel_or_multiple_alt_snp_IDs = unique(indel_or_multiple_alt_snp_IDs),
                no_position_noted_snp_IDs = unique(no_position_noted_snp_IDs)
        ))
}


# helper functions -----------------------------------------------------------------------------------------------


.is_na_scalar = function(x) {
        length(x) == 1 && is.na(x)
}


.gen_output_file_ldsnp = function(
        index_snp_file,
        population,
        r2_cutoff,
        for_cnv_call = F,
        output_dir = "./"
) {
        snp_batch_id = gsub("\\.[^.]*$", "", basename(index_snp_file))
        if (!.is_na_scalar(population)) {
                population_text = paste(population, collapse = "_")
                if (for_cnv_call) {
                        output_file = paste0(output_dir, "/", snp_batch_id, "_cnv_", population_text, "_", r2_cutoff, ".csv")
                } else {
                        output_file = paste0(output_dir, "/", snp_batch_id, "_", population_text, "_", r2_cutoff, ".csv")
                }
        } else {
                output_file = paste0(output_dir, "/", snp_batch_id, "_riskPop_", r2_cutoff, ".csv")
        }
        return(output_file)
}


.read_index_snp_file = function(index_snp_file, population = NA) {
        sep = ifelse(grepl("\\.tsv$", index_snp_file, ignore.case = T), "\t", ",")
        index_snp_raw = utils::read.table(
                index_snp_file,
                header = F,
                sep = sep,
                quote = "",
                comment.char = "",
                stringsAsFactors = F,
                check.names = F,
                fill = T
        )
        first_row = tolower(as.character(unlist(index_snp_raw[1, ])))
        known_header_terms = c(
                "snp", "rsid", "id", "population", "paper",
                "population..asian.eur.", "populations_identifying"
        )
        if (any(first_row %in% known_header_terms)) {
                names(index_snp_raw) = make.names(as.character(unlist(index_snp_raw[1, ])), unique = T)
                index_snp_raw = index_snp_raw[-1, , drop = F]
                rownames(index_snp_raw) = NULL
        } else {
                names(index_snp_raw) = paste0("V", seq_len(ncol(index_snp_raw)))
        }

        if ("SNP" %in% names(index_snp_raw)) {
                snp_col = "SNP"
        } else if ("rsID" %in% names(index_snp_raw)) {
                snp_col = "rsID"
        } else if ("rsid" %in% names(index_snp_raw)) {
                snp_col = "rsid"
        } else if ("id" %in% names(index_snp_raw)) {
                snp_col = "id"
        } else {
                snp_col = "V1"
        }

        if (!.is_na_scalar(population)) {
                index_snp_df = data.frame(
                        rsID = unique(as.character(index_snp_raw[[snp_col]])),
                        population = paste(population, collapse = ","),
                        stringsAsFactors = F
                )
                index_snp_df$population = .normalize_population_string(index_snp_df$population)
                return(index_snp_df)
        }

        pop_col_candidates = c("population", "Population", "Population..Asian.Eur.", "populations_identifying", "V2")
        pop_col = pop_col_candidates[pop_col_candidates %in% names(index_snp_raw)][1]
        if (is.na(pop_col)) {
                # if DNE, fallback to 4 major HR pops from Daniel's scripts
                index_snp_df = data.frame(
                        rsID = unique(as.character(index_snp_raw[[snp_col]])),
                        population = "AFR,AMR,ASN,EUR",
                        stringsAsFactors = F
                )
        } else if (snp_col == "SNP") {
                unique_snps = unique(as.character(index_snp_raw[[snp_col]]))
                population_list = lapply(unique_snps, function(snp_i) {
                        snp_info = index_snp_raw[index_snp_raw[[snp_col]] == snp_i, , drop = F]
                        raw_population_string = paste(snp_info[[pop_col]], collapse = ",")
                        .normalize_population_string(raw_population_string)
                })
                index_snp_df = data.frame(
                        rsID = unique_snps,
                        population = unlist(population_list),
                        stringsAsFactors = F
                )
        } else {
                index_snp_df = data.frame(
                        rsID = as.character(index_snp_raw[[snp_col]]),
                        population = as.character(index_snp_raw[[pop_col]]),
                        stringsAsFactors = F
                )
                index_snp_df$population = .normalize_population_string(index_snp_df$population)
                index_snp_df = unique(index_snp_df)
        }

        return(index_snp_df)
}


.normalize_population_string = function(population_string) {
        population_vec = unlist(strsplit(paste(population_string, collapse = ","), split = ","))
        population_vec = trimws(population_vec)
        population_vec = population_vec[population_vec != ""]
        population_vec = gsub(" ", "", population_vec)
        population_vec = toupper(population_vec)
        population_vec[population_vec %in% c("ASIAN", "ASN/EAS", "EAS", "SAS")] = "ASN"
        population_vec[population_vec %in% c("EUR", "EURO", "EUROPEAN")] = "EUR"
        population_vec[population_vec %in% c("AFR", "AFRICAN")] = "AFR"
        population_vec[population_vec %in% c("AMR", "LATINO", "LATINAMERICAN")] = "AMR"
        population_vec = sort(unique(population_vec[population_vec %in% c("AFR", "AMR", "ASN", "EUR")]))
        if (length(population_vec) == 0) {
                population_vec = c("AFR", "AMR", "ASN", "EUR")
        }
        return(paste(population_vec, collapse = ","))
}


.get_haploreg_snp = function(
        snp_id,
        population_string,
        r2_cutoff,
        haploreg_url,
        haploreg_timeout,
        sleep_range
) {
        populations = unlist(strsplit(as.character(population_string), split = ","))
        populations = trimws(populations)
        populations = populations[!is.na(populations) & populations != ""]
        snp_population_results_list = lapply(populations, function(pop_i) {
                .get_haploreg_population(
                        snp_id = snp_id,
                        population = pop_i,
                        r2_cutoff = r2_cutoff,
                        haploreg_url = haploreg_url,
                        haploreg_timeout = haploreg_timeout
                )
        })
        population_failures = vapply(
                snp_population_results_list,
                function(result_i) {isTRUE(attr(result_i, "haploreg_failed"))},
                logical(1)
        )
        population_failure_reasons = vapply(
                snp_population_results_list,
                function(result_i) {
                        reason = attr(result_i, "failure_reason")
                        if (is.null(reason) || length(reason) == 0 || is.na(reason)) {
                                ""
                        } else {
                                as.character(reason)
                        }
                },
                character(1)
        )

        if (length(snp_population_results_list) == 0 || all(vapply(snp_population_results_list, nrow, numeric(1)) == 0)) {
                snp_results_df = .empty_haploreg_result()
        } else {
                snp_results_df = do.call(rbind, snp_population_results_list)
        }
        attr(snp_results_df, "haploreg_failed") = any(population_failures)
        if (any(population_failures)) {
                attr(snp_results_df, "failed_populations") = populations[population_failures]
                attr(snp_results_df, "failure_reasons") = population_failure_reasons[population_failures]
        } else {
                attr(snp_results_df, "failed_populations") = character()
                attr(snp_results_df, "failure_reasons") = character()
        }

        # anti rate limit from Daniel's script
        if (!is.null(sleep_range) && length(sleep_range) == 2 && all(sleep_range >= 0)) {
                Sys.sleep(stats::runif(1, min(sleep_range), max(sleep_range)))
        }
        return(snp_results_df)
}


.get_haploreg_population = function(
        snp_id,
        population,
        r2_cutoff,
        haploreg_url,
        haploreg_timeout
) {
        payload = list(query = snp_id, ldThresh = r2_cutoff, ldPop = population, submit = "Submit", output = "text")
        response_error = NULL
        response = tryCatch(
                {httr::POST(haploreg_url, body = payload, encode = "form", httr::timeout(haploreg_timeout))},
                error = function(e) {response_error <<- conditionMessage(e); NULL}
        )
        if (is.null(response)) {
                return(.failed_haploreg_result(reason = paste0("request error: ", response_error)))
        }
        if (httr::http_error(response)) {
                return(.failed_haploreg_result(reason = paste0("HTTP ", httr::status_code(response))))
        }
        content_error = NULL
        content = tryCatch(
                {httr::content(response, "text", encoding = "UTF-8")},
                error = function(e) {content_error <<- conditionMessage(e); ""}
        )
        if (!is.null(content_error)) {
                return(.failed_haploreg_result(reason = paste0("content error: ", content_error)))
        }
        if (length(content) == 0 || all(is.na(content)) || !nzchar(paste(content, collapse = ""))) {
                return(.failed_haploreg_result(reason = "empty response"))
        }
        content = paste(content, collapse = "\n")

        if (grepl("<html|Gateway Time-out|504|server didn't respond", content, ignore.case = T)) {
                return(.failed_haploreg_result(reason = "HTML or gateway timeout response"))
        }

        content_lines = strsplit(content, split = "\n")[[1]]
        content_lines = gsub("\r", "", content_lines)
        content_lines = content_lines[content_lines != ""]
        if (length(content_lines) == 0) {
                return(.failed_haploreg_result(reason = "response contained no header"))
        }
        header = trimws(strsplit(content_lines[1], split = "\t")[[1]])
        header_lower = tolower(header)
        if (!any(header_lower %in% c("rsid", "rs_id")) || !any(header_lower %in% c("chr", "chrom", "chromosome"))) {
                return(.failed_haploreg_result(reason = "unexpected response header"))
        }

        parse_error = NULL
        results_table = tryCatch(
                {utils::read.delim(text = paste(content_lines, collapse = "\n"),
                        header = T,
                        sep = "\t",
                        stringsAsFactors = F,
                        check.names = F,
                        fill = T
                )},
                error = function(e) {parse_error <<- conditionMessage(e); NULL}
        )
        if (is.null(results_table)) {
                return(.failed_haploreg_result(reason = paste0("table parsing error: ", parse_error)))
        }

        # valid Haploreg response with no rows
        if (nrow(results_table) == 0) {
                return(.successful_empty_haploreg_result())
        }
        results_table = .standardize_haploreg_columns(
                results_table = results_table,
                snp_id = snp_id,
                population = population
        )
        if (nrow(results_table) == 0) {
                return(.successful_empty_haploreg_result())
        }
        r2_numeric = suppressWarnings(
                as.numeric(results_table$r2)
        )
        results_table = results_table[!is.na(r2_numeric) & r2_numeric >= r2_cutoff,]
        if (nrow(results_table) == 0) {
                return(.successful_empty_haploreg_result())
        }
        attr(results_table, "haploreg_failed") = F
        attr(results_table, "failure_reason") = NA_character_
        return(results_table)
}


.empty_haploreg_result = function() {
        data.frame(
                rsID = character(),
                chr = character(),
                pos_hg38 = character(),
                ref = character(),
                alt = character(),
                r2 = character(),
                D. = character(),
                is_query_snp = character(),
                query_snp = character(),
                population = character(),
                reference_population = character(),
                query_snp_rsid = character(),
                stringsAsFactors = F
        )
}


.successful_empty_haploreg_result = function() {
        result = .empty_haploreg_result()
        attr(result, "haploreg_failed") = F
        attr(result, "failure_reason") = NA_character_
        return(result)
}


.failed_haploreg_result = function(reason) {
        result = .empty_haploreg_result()
        attr(result, "haploreg_failed") = T
        attr(result, "failure_reason") = as.character(reason)
        return(result)
}


.write_failed_haploreg_snps = function(failed_snps, failed_snp_file) {
        failed_snps = unique(as.character(failed_snps))
        failed_snps = failed_snps[!is.na(failed_snps) & nzchar(failed_snps)]
        failed_dir = dirname(failed_snp_file)
        if (!dir.exists(failed_dir)) {
                dir.create(failed_dir, recursive = T)
        }
        writeLines(failed_snps, con = failed_snp_file)
        cat("[INFO] Failed HaploReg SNPs written:", length(failed_snps), "\n")
        invisible(failed_snp_file)
}


.standardize_haploreg_columns = function(results_table, snp_id, population) {
        names(results_table) = trimws(names(results_table))

        # changing col name to be R compatible
        if ("D'" %in% names(results_table) && !("D." %in% names(results_table))) {
                names(results_table)[names(results_table) == "D'"] = "D."
        }
        
        required_cols = names(.empty_haploreg_result())
        for (col in required_cols) {
                if (!(col %in% names(results_table))) {
                        results_table[[col]] = NA
                }
        }

        results_table$chr = gsub("Array", "", as.character(results_table$chr))
        results_table$chr = trimws(results_table$chr)
        results_table$rsID = as.character(results_table$rsID)
        results_table$pos_hg38 = as.character(results_table$pos_hg38)
        results_table$ref = as.character(results_table$ref)
        results_table$alt = as.character(results_table$alt)
        results_table$r2 = as.character(results_table$r2)
        results_table$D. = as.character(results_table$D.)

        # add idx SNP/pop info
        results_table$query_snp = snp_id
        results_table$query_snp_rsid = snp_id
        results_table$population = population
        results_table$reference_population = population

        # is row idx SNP
        results_table$is_query_snp = as.character(
                as.numeric(results_table$rsID == snp_id)
        )

        # remove blank/bad rows
        results_table = results_table[
                !(is.na(results_table$rsID) | results_table$rsID == ""),
        ]

        return(results_table[, required_cols, drop = F])
}


.combine_duplicate_ld_snps = function(haploreg_df, index_snps) {
        if (nrow(haploreg_df) == 0) {
                return(.empty_haploreg_result())
        }
        unique_snps = unique(as.character(haploreg_df$rsID))
        haploreg_df_proc = data.frame()
        for (snp in unique_snps) {
                snp_entries = haploreg_df[haploreg_df$rsID == snp, , drop = F]
                r2_numeric = suppressWarnings(as.numeric(snp_entries$r2))
                ord = order(r2_numeric, decreasing = T, na.last = T)
                snp_entries = snp_entries[ord, , drop = F]
                base = snp_entries[1, , drop = F]
                base$r2 = paste(snp_entries$r2, collapse = "|")
                base$D. = paste(snp_entries$D., collapse = "|")
                base$is_query_snp = as.character(as.numeric(snp %in% index_snps))
                base$query_snp = paste(unique(snp_entries$query_snp), collapse = "|")
                base$query_snp_rsid = base$query_snp
                base$population = paste(unique(snp_entries$population), collapse = "|")
                base$reference_population = paste(unique(snp_entries$reference_population), collapse = "|")
                haploreg_df_proc = rbind(haploreg_df_proc, base)
        }
        return(haploreg_df_proc)
}


.write_ld_output = function(ldsnp_info_df, output_file) {
        if (grepl("\\.tsv$", output_file, ignore.case = T)) {
                write.tsv0(ldsnp_info_df, output_file)
        } else {
                write.csv0(ldsnp_info_df, output_file)
        }
}


# backwards compatibility ------------------------------------------------------------------------------------
# not sure if this is needed, but just in case, leaving it here for now


get_ldsnp_info_main = function(
        index_snp_file = "./data/input_snps/LUC_Index_SNPs_20160607.csv",
        population = NA,
        r2_cutoff = 0.5,
        for_cnv_call = F,
        output_dir = "./",
        output_file = NA,
        failed_snp_file = NA
) {
        get_ldsnp_info(
                index_snp_file = index_snp_file,
                population = population,
                r2_cutoff = r2_cutoff,
                for_cnv_call = for_cnv_call,
                output_dir = output_dir,
                output_file = output_file,
                failed_snp_file = failed_snp_file
        )
}