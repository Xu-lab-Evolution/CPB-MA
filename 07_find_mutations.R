###############################################################################
## De-novo mutation filtering pipeline from a joint-called VCF
###############################################################################
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(data.table)
  library(vcfR)
})

## --- Configuration -----------------------------------------------------------
setwd("/lustre/project/m2_jgu-evolpest/xus/CheckCPB/11_jointDenovo")
VCF_FILE <- "isec/0000.vcf"
DEPTH_FILE     <- "/lustre/miifs01/project/m2_jgu-evolpest/xus/CheckCPB/03_mpileup/MutationAccumulation-Stats_update.csv"
PED_FILE       <- "/lustre/project/m2_jgu-evolpest/xus/CheckCPB/01_Variantscalling/pedigree.ped"
CALLABLE_SITES_FILE <- "Callable_sites.tsv"

SEX_CHR <- "chr6"   # set your sex chromosome label

## Filtering knobs
ENFORCE_POP_LEAK_CHECK <- TRUE
PARENT_ALT_MAX_COUNT   <- 1
PARENT_ALT_MAX_FRAC    <- 0
CHILD_BALANCE_P_MIN    <- 0.05


## --- Extract ALT and TOTAL depths from AD matrix ---------------------
get_depths_vcfR <- function(ad_matrix, sample_name) {
  sample_col <- ad_matrix[, sample_name, drop = TRUE]
  sample_col <- as.character(sample_col)
  sample_col[is.na(sample_col) | sample_col == ""] <- "."
  split_ad <- strsplit(sample_col, ",", fixed = TRUE)
  to_int <- function(x) suppressWarnings(as.integer(x))
  alt <- vapply(split_ad, function(x) {
    if (length(x) < 2 || is.na(x[1]) || x[1] == ".") return(0L)
    v <- to_int(x[2]); if (is.na(v)) 0L else v
  }, integer(1))
  tot <- vapply(split_ad, function(x) {
    if (length(x) == 0 || is.na(x[1]) || x[1] == ".") return(0L)
    v <- to_int(x); v[is.na(v)] <- 0L; sum(v)
  }, integer(1))
  list(alt = alt, tot = tot)
}

## --- Trio-level filtering (vcfR) --------------------------------------------
find_trio_denovo <- function(vcf_obj, FatherID, MotherID, ChildID, ChildSex, depth_data) {
  trio_ids <- c(FatherID, MotherID, ChildID)
  trio_mean_depths <- depth_data$Mean_depth[match(trio_ids, depth_data$SampleID)]
  if (any(is.na(trio_mean_depths))) {
    warning("Samples not found in depth file for trio: ", paste(trio_ids, collapse = ", "), ". Skipping depth filter.")
    vcf_filtered_depth <- vcf_obj
  } else {
    trio_avg_depth <- mean(trio_mean_depths, na.rm = TRUE)
    dp_mat_trio <- extract.gt(vcf_obj, element = "DP", as.numeric = TRUE)[, trio_ids, drop = FALSE]
    keep_indices <- apply(dp_mat_trio, 1, function(dp_row) {
      all(dp_row >= trio_avg_depth / 2 & dp_row <= trio_avg_depth * 2, na.rm = TRUE)
    })
    keep_indices[is.na(keep_indices)] <- FALSE
    vcf_filtered_depth <- vcf_obj[keep_indices, ]
  }
  if (nrow(vcf_filtered_depth@gt) == 0) return(NULL)
  gt <- extract.gt(vcf_filtered_depth, element = "GT")
  ad <- extract.gt(vcf_filtered_depth, element = "AD")
  parent_ad_f <- get_depths_vcfR(ad, FatherID)
  parent_ad_m <- get_depths_vcfR(ad, MotherID)
  parent_ok_vec <- function(alt, tot) { (alt <= PARENT_ALT_MAX_COUNT) | (ifelse(tot > 0, alt/tot, 0) <= PARENT_ALT_MAX_FRAC) }
  parent_ad_ok <- parent_ok_vec(parent_ad_f$alt, parent_ad_f$tot) & parent_ok_vec(parent_ad_m$alt, parent_ad_m$tot)
  is_autosomal <- vcf_filtered_depth@fix[, "CHROM"] != SEX_CHR
  autosomal_base <- (gt[, ChildID] %in% c("0/1", "0|1")) & (gt[, FatherID] %in% c("0/0", "0|0", "0")) & (gt[, MotherID] %in% c("0/0", "0|0", "0"))
  child_depths <- get_depths_vcfR(ad, ChildID)
  prob <- 2 * pbinom(pmin(child_depths$alt, child_depths$tot - child_depths$alt), size = child_depths$tot, prob = 0.5)
  prob[is.na(prob) | child_depths$tot == 0] <- 0
  if (ENFORCE_POP_LEAK_CHECK) {
    other_samples <- setdiff(colnames(vcf_filtered_depth@gt), c("FORMAT", FatherID, MotherID, ChildID))
    pop_alt_sum <- integer(nrow(ad))
    if (length(other_samples) > 0) { for (s in other_samples) pop_alt_sum <- pop_alt_sum + get_depths_vcfR(ad, s)$alt }
    no_leaks_in_pop <- (pop_alt_sum <= 3)
  } else { no_leaks_in_pop <- rep(TRUE, nrow(ad)) }
  keep_autosomal <- is_autosomal & autosomal_base & (prob >= CHILD_BALANCE_P_MIN) & parent_ad_ok & no_leaks_in_pop
  is_sex_chr <- vcf_filtered_depth@fix[, "CHROM"] == SEX_CHR
  keep_sex_chr <- rep(FALSE, nrow(vcf_filtered_depth@gt))
  if (any(is_sex_chr)) {
    if (ChildSex == 2) {
      sex_base <- (gt[, ChildID] %in% c("0/1", "0|1")) & (gt[, FatherID] %in% c("0/0", "0|0", "0")) & (gt[, MotherID] %in% c("0/0", "0|0", "0"))
      keep_sex_chr <- is_sex_chr & sex_base & (prob >= CHILD_BALANCE_P_MIN) & parent_ad_ok & no_leaks_in_pop
    } else {
      sex_base <- (gt[, ChildID] %in% c("1/1", "1|1", "1")) & (gt[, FatherID] %in% c("0/0", "0|0", "0")) & (gt[, MotherID] %in% c("0/0", "0|0", "0"))
      keep_sex_chr <- is_sex_chr & sex_base & parent_ad_ok & no_leaks_in_pop
    }
  }
  final_keep <- keep_autosomal | keep_sex_chr
  dn_vcf <- vcf_filtered_depth[which(final_keep), ]
  if (nrow(dn_vcf@gt) == 0) return(NULL)
  tbl <- data.table(CHROM=dn_vcf@fix[,"CHROM"], POS=as.integer(dn_vcf@fix[,"POS"]), REF=dn_vcf@fix[,"REF"], ALT=dn_vcf@fix[,"ALT"],
                    GT_Child=extract.gt(dn_vcf,"GT")[,ChildID], GT_Father=extract.gt(dn_vcf,"GT")[,FatherID], GT_Mother=extract.gt(dn_vcf,"GT")[,MotherID],
                    DP_Child=extract.gt(dn_vcf,"DP",as.numeric=TRUE)[,ChildID], DP_Father=extract.gt(dn_vcf,"DP",as.numeric=TRUE)[,FatherID], DP_Mother=extract.gt(dn_vcf,"DP",as.numeric=TRUE)[,MotherID])
  return(tbl)
}

###############################################################################
## Calculate Mutation Rate with 95% CI
###############################################################################
calculate_and_summarize_mutation_rate <- function(denovo_table, callable_sites_table) {
  message("-----------------------------------------------------------")
  message("Calculating mutation rates...")

  # Count de novo mutations for each offspring
  mutation_counts <- denovo_table[, .N, by = .(`Offspring sample ID`)]
  setnames(mutation_counts, "N", "DeNovoCount")

  # Ensure the callable sites table has the correct column name for merging
  if ("Child" %in% names(callable_sites_table)) {
    setnames(callable_sites_table, "Child", "Offspring sample ID")
  }

  # Merge mutation counts with callable sites data
  merged_data <- merge(callable_sites_table, mutation_counts, by = "Offspring sample ID", all.x = TRUE)
  merged_data[is.na(DeNovoCount), DeNovoCount := 0] # Set count to 0 for trios with no mutations

  # Calculate the mutation rate
  merged_data[, MutationRate := DeNovoCount / (2 * Callable_bp)]

  # Calculate 95% CI for the mutation rate using exact Poisson confidence intervals
  # Lower bound for count
  ci_count_lower <- 0.5 * qchisq(0.025, 2 * merged_data$DeNovoCount)
  ci_count_lower[merged_data$DeNovoCount == 0] <- 0 # Set lower bound to 0 if count is 0
  
  # Upper bound for count
  ci_count_upper <- 0.5 * qchisq(0.975, 2 * (merged_data$DeNovoCount + 1))

  # Convert count CIs to rate CIs
  merged_data[, CI_95_lower := ci_count_lower / (2 * Callable_bp)]
  merged_data[, CI_95_upper := ci_count_upper / (2 * Callable_bp)]

  message("Mutation rate calculation complete.")
  return(merged_data)
}

###############################################################################
## Main
###############################################################################
message("Loading pedigree and depth data...")
ped_data_raw <- unique(fread(PED_FILE, header = FALSE))
standard_ped_cols <- c("Family", "ChildID", "FatherID", "MotherID", "Sex", "Phenotype")
setnames(ped_data_raw, names(ped_data_raw)[1:ncol(ped_data_raw)], standard_ped_cols[1:ncol(ped_data_raw)])
ped_data <- ped_data_raw
setnames(ped_data, c("ChildID", "FatherID", "MotherID"), c("Child", "Father", "Mother"))

depth_data <- fread(DEPTH_FILE, sep = ";")
if (!"SampleID" %in% names(depth_data) && "Sample" %in% names(depth_data)) {
  setnames(depth_data, "Sample", "SampleID")
}

message("Loading callable sites data...")
callable_sites_data <- fread(CALLABLE_SITES_FILE)
message(sprintf("Loaded callable sites for %d trios.", nrow(callable_sites_data)))

message("Loading large joint VCF file. This may take a while...")
full_vcf <- read.vcfR(VCF_FILE, verbose = FALSE)
message("VCF file loaded successfully.")

message("Filtering to keep only bi-allelic variants...")
is_biallelic <- is.biallelic(full_vcf); is_biallelic[is.na(is_biallelic)] <- FALSE
full_vcf <- full_vcf[is_biallelic, ]
message(format(nrow(full_vcf@gt), big.mark = ","), " bi-allelic variants remaining for analysis.")

master_list <- list()
message("Starting de novo analysis for ", nrow(ped_data), " trios...")
for (i in seq_len(nrow(ped_data))) {
  trio <- ped_data[i, ]
  message(sprintf("Analyzing trio %d/%d: Child %s (Father %s, Mother %s)", i, nrow(ped_data), trio$Child, trio$Father, trio$Mother))
  if (!all(c(trio$Father, trio$Mother, trio$Child) %in% colnames(full_vcf@gt))) {
    message("  -> Skipping trio: One or more samples not found in VCF.")
    next
  }
  dn_table <- find_trio_denovo(full_vcf, trio$Father, trio$Mother, trio$Child, trio$Sex, depth_data)
  if (!is.null(dn_table) && nrow(dn_table) > 0) {
    dn_table[, `:=`(`Father sample ID`=trio$Father, `Mother sample ID`=trio$Mother, `Offspring sample ID`=trio$Child, `Crossing Family ID`=trio$Family)]
    master_list[[length(master_list) + 1]] <- dn_table
    message("  -> Found ", nrow(dn_table), " potential de novo mutations.")
  } else { message("  -> No de novo mutations found.") }
}

if (length(master_list) > 0) {
  message("Combining results from all trios...")
  master_table <- rbindlist(master_list, use.names = TRUE, fill = TRUE)
  snp_table   <- master_table[nchar(REF) == 1 & nchar(ALT) == 1]
  indel_table <- master_table[nchar(REF) != 1 | nchar(ALT) != 1]
  message(sprintf("Total de novo variants found: %d SNPs, %d Indels.", nrow(snp_table), nrow(indel_table)))
  format_output_table <- function(dt) { dt[, .(CHR=CHROM, POS, REF, ALT, `Father Genotype`=GT_Father, `Mother Genotype`=GT_Mother, `Offspring Genotype`=GT_Child, `Father sample ID`, `Mother sample ID`, `Offspring sample ID`, `Crossing Family ID`, `Depth of Offspring`=DP_Child, `Depth of Father`=DP_Father, `Depth of M
other`=DP_Mother)] }
  if (nrow(snp_table) > 0) { fwrite(format_output_table(snp_table), "all_denovo_snps_from_joint_vcf.tsv", sep = "\t", quote = FALSE); message("✓ SNP summary written to: all_denovo_snps_from_joint_vcf.tsv") }
  if (nrow(indel_table) > 0) { fwrite(format_output_table(indel_table), "all_denovo_indels_from_joint_vcf.tsv", sep = "\t", quote = FALSE); message("✓ Indel summary written to: all_denovo_indels_from_joint_vcf.tsv") }
} else {
  message("No de novo mutations were found across any trios.")
  master_table <- data.table()
}


## --- Mutation Rate Calculation and Output ---
mutation_rate_results <- calculate_and_summarize_mutation_rate(master_table, callable_sites_data)

message("--- Mutation Rates per Trio ---")
print(mutation_rate_results)
fwrite(mutation_rate_results, "mutation_rates_per_trio.tsv", sep = "\t", quote = FALSE)
message("✓ Per-trio mutation rate summary written to: mutation_rates_per_trio.tsv")

## --- Global Mutation Rate Calculation ---
total_mutations <- sum(mutation_rate_results$DeNovoCount)
total_callable_bp <- sum(mutation_rate_results$Callable_bp)
global_rate <- total_mutations / (2 * total_callable_bp)

# CI for the total count
ci_count_lower_global <- 0.5 * qchisq(0.025, 2 * total_mutations)
ci_count_upper_global <- 0.5 * qchisq(0.975, 2 * (total_mutations + 1))
if (total_mutations == 0) ci_count_lower_global <- 0

# CI for the global rate
ci_rate_lower_global <- ci_count_lower_global / (2 * total_callable_bp)
ci_rate_upper_global <- ci_count_upper_global / (2 * total_callable_bp)

global_summary <- data.table(
  Description = "All Families Combined",
  TotalDeNovoCount = total_mutations,
  TotalCallable_bp = total_callable_bp,
  MutationRate = global_rate,
  CI_95_lower = ci_rate_lower_global,
  CI_95_upper = ci_rate_upper_global
)

message("\n--- Global Mutation Rate Summary ---")
print(global_summary)
fwrite(global_summary, "mutation_rate_global_summary.tsv", sep = "\t", quote = FALSE)
message("✓ Global mutation rate summary written to: mutation_rate_global_summary.tsv")
message("-----------------------------------------------------------")
