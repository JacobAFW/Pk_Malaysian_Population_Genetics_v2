# Stage B-2 helper: attach sample metadata to Pk.fam and write Pk.dups.
# Lifted from scripts/gadi/Analyses.Rmd L495-579 (paths rewritten to data/gadi/).
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(readxl); library(janitor)
})
source(here::here("scripts/_setup.R"))

PLINK_DIR <- file.path(DATA_PROC, "plink")
fam_path  <- file.path(PLINK_DIR, "Pk.fam")
stopifnot(file.exists(fam_path))

metadata <- read_csv(file.path(PROJ, "data/gadi/full_dataset_011223.csv"),
                     show_col_types = FALSE) |>
  rename(Sample = studycode) |>
  add_column(Cluster = NA) |>
  mutate(sevemal = ifelse(grepl("Um", sevemal), "Um",
                          ifelse(sevemal == "1", "Sm",
                                 ifelse(sevemal == "0", "Um", sevemal)))) |>
  mutate(Sample = str_to_upper(Sample))

original_metadata <- read_xlsx(file.path(PROJ, "data/gadi/PK_Sabah_Sample_naming_indexes.xlsx")) |>
  select(1:3) |>
  mutate(subjectid = ifelse(str_length(subjectid) == 1, paste0("00", subjectid),
                            ifelse(str_length(subjectid) == 2, paste0("0", subjectid), subjectid))) |>
  unite(Sample2, c("group", "subjectid"), sep = "") |>
  rename(Sample = sampleid)

archived_data <- read_csv(file.path(PROJ, "data/gadi/Pk_clusters_metadata.csv"),
                          show_col_types = FALSE) |>
  select(1, 3, 5) |>
  rename(district = area) |>
  mutate(Cluster = str_remove(Group, "-Pk")) |>
  select(-Group) |>
  rbind(read_csv(file.path(PROJ, "data/gadi/Pk_clusters_peninsular_metadata.csv"),
                 show_col_types = FALSE) |>
    select(2, 8) |>
    rename(Sample = ENA_accession_no_ES, district = Location) |>
    add_column(Cluster = "Peninsular")) |>
  add_column(enroldate = NA, age = NA, sex = NA, ethnicity = NA, village = NA,
             occupation = NA, Hb = NA, platelets = NA, parasitemia = NA, fever = NA,
             respiratory_rate = NA, oxygen_saturation = NA, systolic_BP = NA,
             diastolic_BP = NA, G6PD = NA, pregnant = NA, creatinine = NA, study = NA,
             bicarbonate = NA, PCR = NA, sevemal = NA, bilirubin = NA, glucose = NA,
             bleeding_severity = NA)

metadata <- metadata |>
  rbind(archived_data) |>
  select(Sample, sevemal) |>
  # de-dup before joining to .fam — multi-mapping rows inflate fam row count
  # (without this, NKU2_5 etc. get duplicated by the left_join and break plink)
  distinct(Sample, .keep_all = TRUE)

Plink_fam_raw <- read.table(fam_path, sep = " ", stringsAsFactors = FALSE)

Plink_fam <- Plink_fam_raw |>
  filter(!grepl("PK_SB_DNA", V1)) |>
  mutate(Sample = str_remove(V1, "_DK.*")) |>
  mutate(Sample = str_replace(Sample, "_", "-")) |>
  left_join(metadata, by = "Sample") |>
  select(-Sample) |>
  rbind(
    Plink_fam_raw |>
      filter(grepl("PK_SB_DNA", V1)) |>
      mutate(Sample = str_remove(V1, "_DK.*")) |>
      left_join(original_metadata, by = "Sample") |>
      mutate(Sample = Sample2) |>
      left_join(metadata, by = "Sample") |>
      select(-c(Sample, Sample2))
  )

write.table(Plink_fam, fam_path, quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)

Plink_fam |> select(-(2:6)) |> write_csv(file.path(PLINK_DIR, "Pk.csv"))

# Duplicates
dups <- read.table(fam_path, sep = " ", stringsAsFactors = FALSE) |>
  mutate(Sample = str_remove(V1, "_DK.*")) |>
  mutate(Sample = str_replace(Sample, "_", "-")) |>
  group_by(Sample) |>
  filter(n() > 1) |>
  ungroup() |>
  select(1:6)

write.table(dups, file.path(PLINK_DIR, "Pk.dups"),
            quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
message(sprintf("[10b] fam rows = %d  | dup rows = %d  | unique Samples = %d",
                nrow(Plink_fam), nrow(dups), n_distinct(str_replace(str_remove(Plink_fam$V1,"_DK.*"),"_","-"))))
