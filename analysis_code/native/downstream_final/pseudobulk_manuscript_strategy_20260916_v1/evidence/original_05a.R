suppressPackageStartupMessages(library(data.table))

TEMPLATE <- "/Users/taki/fsdownload/pseudobulk_target_discovery_v1/sample_disease_map_template.tsv"
CROSSWALK <- "/Users/taki/Desktop/Codex任务1/sample_metadata_crosswalk_v1.tsv"
OUTPUT <- "/Users/taki/Desktop/Codex任务1/sample_disease_map.tsv"
AUDIT <- "/Users/taki/Desktop/Codex任务1/sample_disease_map_audit.tsv"

template <- fread(TEMPLATE)
crosswalk <- fread(CROSSWALK)

# Collapse technical rows to the biological-sample level. Study overrides below
# are authoritative for the paired primary datasets.
bio <- crosswalk[, .(
  curated_disease_labels = paste(sort(unique(
    disease_standardized[!is.na(disease_standardized) & disease_standardized != ""]
  )), collapse = ";")
), by = biological_sample_uid]

map <- merge(template, bio, by = "biological_sample_uid", all.x = TRUE)
map[, condition_for_pseudobulk := tolower(as.character(condition_v2))]
map[, disease_context := NULL]
map[, disease_context := NA_character_]

# Normalize curated disease labels for descriptive/non-primary records.
label <- tolower(map$curated_disease_labels)
map[grepl("atherosclerosis|coronary artery disease", label),
    disease_context := "Atherosclerosis"]
map[grepl("aortic aneurysm", label),
    disease_context := "Aortic_aneurysm"]
map[grepl("aortic dissection", label),
    disease_context := "Aortic_dissection"]
map[grepl("pulmonary arterial hypertension", label) &
      !grepl("chronic thromboembolic", label),
    disease_context := "Pulmonary_arterial_hypertension"]
map[grepl("chronic thromboembolic pulmonary hypertension", label),
    disease_context := "CTEPH"]
map[grepl("calcific aortic valve|non-fibrocalcific aortic valve", label),
    disease_context := "Calcific_aortic_valve_disease"]

# Source-verified study-level corrections for labels that were blank or encoded
# as control descriptors in the original metadata.
map[study_canonical == "GSE166676", disease_context := "Aortic_aneurysm"]
map[study_canonical %in% c("GSE213740", "GSE222318", "GSE254132"),
    disease_context := "Aortic_dissection"]
map[study_canonical == "GSE220774",
    disease_context := "Calcific_aortic_valve_disease"]
map[study_canonical == "PRJNA562645", `:=`(
  disease_context = "Calcific_aortic_valve_disease",
  condition_for_pseudobulk = fifelse(
    grepl("calcific aortic valve disease", curated_disease_labels, ignore.case = TRUE),
    "disease", "healthy"
  )
)]

# Primary pseudobulk is limited to studies with at least three independent
# donors per condition and an internal healthy comparison. AAA and CAVD remain
# mapped for future sensitivity work but are not used in the primary model
# because their matched control donor counts are insufficient.
primary_studies <- c(
  "GSE159677", "GSE248610",                         # atherosclerosis
  "GSE213740", "GSE222318", "GSE254132",          # aortic dissection
  "GSE228644"                                       # PAH
)
map[, include_primary := study_canonical %in% primary_studies]
map[, mapping_note := fifelse(
  include_primary,
  "Primary: internally controlled study with adequate donor replication",
  "Excluded from primary disease pseudobulk; retained only as mapped metadata"
)]

if (map[include_primary == TRUE, any(is.na(disease_context) | disease_context == "")]) {
  stop("A primary sample lacks disease_context")
}
if (map[include_primary == TRUE, any(!condition_for_pseudobulk %in% c("healthy", "disease"))]) {
  stop("A primary sample has invalid condition_for_pseudobulk")
}

audit <- unique(map[, .(
  biological_sample_uid, donor_uid_v2, study_canonical,
  disease_context, condition_for_pseudobulk, include_primary
)])[, .(
  n_biological_samples = .N,
  n_independent_donors = uniqueN(donor_uid_v2)
), by = .(study_canonical, disease_context, condition_for_pseudobulk, include_primary)]
setorder(audit, -include_primary, disease_context, study_canonical, condition_for_pseudobulk)

keep <- names(template)
keep <- setdiff(keep, c("disease_context", "include_primary", "mapping_note"))
keep <- unique(c(keep, "condition_for_pseudobulk", "disease_context", "include_primary", "mapping_note"))
fwrite(map[, ..keep], OUTPUT, sep = "\t", na = "")
fwrite(audit, AUDIT, sep = "\t", na = "")

cat("Wrote:", OUTPUT, "\n")
cat("Audit:", AUDIT, "\n\n")
print(audit[include_primary == TRUE])
