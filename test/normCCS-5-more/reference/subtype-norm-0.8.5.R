subtypeNorm <- function (df, n_norm_subtype = 2, norm_roc_method = "Raw", norm_rank_marker = "RR", 
    numCores = NULL, verbose = TRUE, norm_search_method = c("bruteforce", 
        "dp")) 
{
    norm_search_method <- match.arg(norm_search_method)
    norm_validate_input(df, n_norm_subtype, norm_roc_method, 
        numCores)
    if (length(norm_rank_marker) != 1L || is.na(norm_rank_marker) || 
        !norm_rank_marker %in% c("RR", "NRR", "RR+NRR", "NRR+RR")) {
        stop("norm_rank_marker must be one of 'RR', 'NRR', 'RR+NRR', 'NRR+RR'.")
    }
    search_started <- proc.time()[["elapsed"]]
    if (T) {
        df2 <- ddply(df, c("Cohort", "Subtype", "response"), 
            dplyr::summarize, nSample = length(SampleIDs), tumor_type = unique(tumor_type))
        df3 <- ddply(df2, c("Cohort", "Subtype", "tumor_type"), 
            dplyr::summarize, size = sum(nSample, na.rm = TRUE), 
            nResponse = ifelse(length(nSample[response %in% c(1, 
                "1")]) == 0, 0, nSample[response %in% c(1, "1")]), 
            response_rate = nResponse/size)
        df4 <- ddply(df3, c("Cohort"), dplyr::reframe, Subtype = "All", 
            size = sum(size, na.rm = TRUE), nResponse = sum(nResponse, 
                na.rm = TRUE), response_rate = nResponse/size)
        df6 <- NULL
        for (i in 1:nrow(df3)) {
            df3.i <- df3[i, ]
            df3.i_all_response_rate <- df4[df4$Cohort %in% df3.i$Cohort & 
                df4$Subtype == "All", "response_rate"]
            df3.i$normalized_response_rate <- (df3.i$response_rate - 
                df3.i_all_response_rate)/df3.i_all_response_rate
            df6 <- rbind(df6, df3.i)
        }
        df_cutoff <- df6 %>% ddply(c("Subtype"), dplyr::summarize, 
            medianNRR = median(normalized_response_rate, na.rm = T), 
            medianRR = median(response_rate, na.rm = T))
        if (norm_rank_marker == "RR+NRR") {
            df_cutoff <- arrange(df_cutoff, desc(medianRR), desc(medianNRR))
        }
        else if (norm_rank_marker == "NRR+RR") {
            df_cutoff <- arrange(df_cutoff, desc(medianNRR), 
                desc(medianRR))
        }
        else if (norm_rank_marker == "RR") {
            df_cutoff <- arrange(df_cutoff, desc(medianRR))
        }
        else if (norm_rank_marker == "NRR") {
            df_cutoff <- arrange(df_cutoff, desc(medianNRR))
        }
        else {
            stop("Not available rank marker. Please use one of \"RR\" and \"NRR\"!")
        }
        if (n_norm_subtype > nrow(df_cutoff)) {
            stop("n_norm_subtype must not exceed the number of observed subtypes.")
        }
        if (norm_search_method == "dp") {
            result <- subtype_norm_dp(df, df_cutoff, n_norm_subtype, 
                norm_roc_method, search_started, verbose)
            result$SearchDiagnostics$rank_marker <- norm_rank_marker
            return(result)
        }
        n_partitions <- choose(nrow(df_cutoff) - 1, n_norm_subtype - 
            1)
        estimated_result_bytes <- n_partitions * (length(unique(df$Cohort)) + 
            1) * 192
        if (!is.finite(n_partitions) || n_partitions > 1e+06 || 
            estimated_result_bytes > 2 * 1024^3) {
            stop("Brute-force search exceeds the partition/result-memory budget; use norm_search_method = 'dp'.")
        }
        combinations <- combn(seq_len(nrow(df_cutoff) - 1), n_norm_subtype - 
            1)
    }
    if (T) {
        if (is.null(numCores)) {
            numCores <- detectCores() - 1
        }
        registerDoParallel(cores = numCores)
        on.exit(stopImplicitCluster(), add = TRUE)
        if (verbose) 
            LuckyVerbose("subtypeNorm: Parallel processing - Find best combination...", 
                levels = 2)
        if (ncol(combinations) >= 10 * numCores) {
            XL <- cut_vector(1:ncol(combinations), nsplit = numCores)
        }
        else {
            XL <- cut_vector(1:ncol(combinations), nsplit = 4)
        }
        system.time(df_perform <- foreach(i = 1:length(XL), .combine = "rbind", 
            .packages = c("plyr", "dplyr", "dplyr", "luckyBase", 
                "pROC")) %dopar% {
            batch_results <- lapply(XL[[i]], function(k) {
                x <- combinations[, k]
                medianNRR_category <- rep(NA, nrow(df_cutoff))
                for (j in 0:length(x)) {
                  lower_limit_j <- ifelse(j == 0, 1, (x[j] + 
                    1))
                  upper_limit_j <- ifelse(j == length(x), nrow(df_cutoff), 
                    x[j + 1])
                  medianNRR_category[lower_limit_j:upper_limit_j] <- n_norm_subtype - 
                    j
                }
                df_cutoff$medianNRR_category <- medianNRR_category - 
                  1
                df_cutoff_summary <- summarize(group_by(df_cutoff, 
                  medianNRR_category), annotation = paste0(Subtype, 
                  collapse = ","))
                df$normSubtype <- convert(df$Subtype, "Subtype", 
                  "medianNRR_category", df_cutoff) %>% as.numeric()
                info_str <- paste(paste(df_cutoff_summary$medianNRR_category, 
                  df_cutoff_summary$annotation, sep = "="), collapse = "; ")
                df_perform_i1 <- cbind(Cohort = "All", Info = info_str, 
                  nSample = nrow(df), get_perform_markers(df$response, 
                    df$normSubtype, n_norm_subtype = n_norm_subtype, 
                    norm_roc_method = norm_roc_method), stringsAsFactors = FALSE)
                df_perform_i2 <- df %>% ddply(., c("Cohort"), 
                  dplyr::summarize, Info = info_str, nSample = length(Cohort), 
                  get_perform_markers(response, normSubtype, 
                    n_norm_subtype = n_norm_subtype, norm_roc_method = norm_roc_method))
                rbind(df_perform_i1, as.data.frame(df_perform_i2))
            })
            do.call("rbind", batch_results)
        })
        stopImplicitCluster()
    }
    if (verbose) 
        LuckyVerbose("subtypeNorm: Summary...", levels = 2)
    system.time(df_perform_summary <- df_perform %>% group_by(Info) %>% 
        dplyr::summarize(medianROCAUC = median(ROCAUC, na.rm = TRUE)) %>% 
        ungroup() %>% arrange(desc(medianROCAUC)))
    if (verbose) 
        LuckyVerbose("subtypeNorm: Done!", levels = 2)
    return(list(Raw = df_perform, Summary = df_perform_summary))
}
