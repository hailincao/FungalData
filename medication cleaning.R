library(dplyr)
library(tibble)
library(writexl)
library(tidyr)
library(stringr)

entryS <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/Report 9-Volunteer Medical History - Report 9-Volunteer Medical History.csv")
mapping <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/Original Study Mapping - Sheet3.csv")

entryS <- entryS %>%
  left_join(mapping, by = c("Your.Biome.Health.App.ID" = "Biome.Health.App.ID")) %>%
  mutate(Your.Biome.Health.App.ID = STUDY.ID) %>%
  select(-STUDY.ID)

entryS_med <- entryS %>%
  select(
    Your.Biome.Health.App.ID,
    If.you.are.currently.taking.any.other.medications..please.list..1..the.medication..2..for.what.reason..3..how.long.you.ve.been.taking.the.medication.
  )

write_xlsx(entryS_med, "/Users/caoyang/Desktop/Tetel Lab/datasets/entryS_med.xlsx")


medications <- tibble(
  study_id = character(),
  medication_name = character(),
  dosage = character(),
  SSRI = character()
)

medications <- medications %>%
  add_row(study_id = "59", medication_name = "Adderall ER", dosage = "20mg", SSRI = "0") %>%
  add_row(study_id = "59", medication_name = "Pristiq", dosage = "25mg", SSRI = "0") 

medications <- medications %>%
  add_row(study_id = "51", medication_name = "Prozac", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "51", medication_name = "Celebrex", dosage = NA, SSRI = "0") 

medications <- medications %>%
  add_row(study_id = "27", medication_name = "Zoloft", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "27", medication_name = "Methylphenidate", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "11", medication_name = "Spirinolactone", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "11", medication_name = "Minocycline", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "11", medication_name = "Escitalopram", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "11", medication_name = "Aripiprazole", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "39", medication_name = "Budesonide", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "16", medication_name = "Ventolin", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "16", medication_name = "Advair", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "16", medication_name = "Alprazolam", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "16", medication_name = "Bupropion", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "16", medication_name = "Escitalopram", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "16", medication_name = "Olanzipine", dosage = NA, SSRI = "1")

medications <- medications %>%
  add_row(study_id = "28", medication_name = "Spironolactone", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "28", medication_name = "Sertraline", dosage = NA, SSRI = "1")

medications <- medications %>%
  add_row(study_id = "32", medication_name = "Amitriptyline", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "14", medication_name = "Zoloft", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "14", medication_name = "Propranolol", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "30", medication_name = "Trileptal", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "30", medication_name = "Sertraline", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "30", medication_name = "Atorvastatin", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "3", medication_name = "Flush and Defend by Uqora", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "18", medication_name = "Adderall", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "31", medication_name = "Bupropion", dosage = "300mg", SSRI = "0") %>%
  add_row(study_id = "31", medication_name = "Methylphenidate HCI", dosage = "10mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "33", medication_name = "Lexapro", dosage = "20mg", SSRI = "1") %>%
  add_row(study_id = "33", medication_name = "Adderall XR", dosage = "20mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "60", medication_name = "Escitalopram", dosage = "20mg", SSRI = "1") %>%
  add_row(study_id = "60", medication_name = "Melatonin", dosage = "5mg", SSRI = "0") %>%
  add_row(study_id = "60", medication_name = "Vitamin D", dosage = "25mcg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "1", medication_name = "Lamictal", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "1", medication_name = "Gabapentin", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "24", medication_name = "Prozac", dosage = NA, SSRI = "1")

medications <- medications %>%
  add_row(study_id = "52", medication_name = "Dextroamphetamine", dosage = "30mg 24hr Capsule", SSRI = "0") %>%
  add_row(study_id = "52", medication_name = "Dextroamphetamine", dosage = "7.5mg Tablet", SSRI = "0") %>%
  add_row(study_id = "52", medication_name = "Albuterol", dosage = "90mcg/actuation inhaler", SSRI = "0") %>%
  add_row(study_id = "52", medication_name = "Rizatriptan", dosage = "5mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "65", medication_name = "Focalin", dosage = "15mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "9", medication_name = "Fluticasone Propionate", dosage = "nasal spray and eye drops", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "50", medication_name = "Spironolactone", dosage = "100mg", SSRI = "0")

#the following data comes from Copy of prehandcoded_med_data
prehandcoded_med_data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/prehandcoded_med_data.csv")

coded_med_df <- prehandcoded_med_data %>%
  select(biome_id, other_medicine) %>%
  filter(!is.na(other_medicine), other_medicine != "") %>%
  mutate(
    other_medicine = str_replace_all(other_medicine, "[\r\n]", ",")  # Convert newlines to commas
  ) %>%
  separate_rows(other_medicine, sep = ",|;") %>%  # Split by comma or semicolon
  mutate(
    medication_name = str_trim(str_remove(other_medicine, "\\(.*?\\)")),  # Remove parentheses content
    notes = str_extract(other_medicine, "\\(.*?\\)") %>% str_remove_all("[()]")  # Optional: capture what was in parentheses
  ) %>%
  filter(medication_name != "") %>%
  select(biome_id, medication_name, notes)

coded_med_df <- coded_med_df %>%
  mutate(
    # Move all dosage/schedule/form info to notes (if any)
    notes = if_else(
      str_detect(medication_name, "\\d+\\s*(mg|mcg|g|IU|units|ml|mL|%)|XR|IR|SR|daily|night|morning|each"),
      str_extract(medication_name, ".*?(\\d+\\s*(mg|mcg|g|IU|units|ml|mL|%)?.*)"),
      notes
    ),
    # Extract clean medication name (first word(s) before dosage/form)
    medication_name = str_remove(medication_name, "(\\s+\\d+.*)|(\\s+XR|IR|SR|daily|once|each).*"),
    medication_name = str_trim(medication_name)
  ) 

medications <- medications %>%
  add_row(study_id = "14", medication_name = "Propranolol", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "14", medication_name = "Zoloft", dosage = NA, SSRI = "1") %>%
  add_row(study_id = "14", medication_name = "Low-Ogestrel", dosage = NA, SSRI = "0")


medications <- medications %>%
  add_row(study_id = "30", medication_name = "Orilissa", dosage = NA, SSRI = "0")


medications <- medications %>%
  add_row(study_id = "62", medication_name = "Alyacen", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "74", medication_name = "Bupropion", dosage = "200mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "73", medication_name = "Loestrin", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "73", medication_name = "Lexapro", dosage = "5mg", SSRI = "1") %>%
  add_row(study_id = "73", medication_name = "Zyrtech", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "74", medication_name = "Bupropion", dosage = "200mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "75", medication_name = "Sertraline", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "18", medication_name = "Adderall XR", dosage = "15mg", SSRI = "0") %>%
  add_row(study_id = "18", medication_name = "Adderall IR", dosage = "10mg", SSRI = "0")
  
medications <- medications %>%
  add_row(study_id = "66", medication_name = "Adderall", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "69", medication_name = "Thyroid medication", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "69", medication_name = "Levothyroxine", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "69", medication_name = "Multivitamins", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "69", medication_name = "Vitamin D", dosage = NA, SSRI = "0")


medications <- medications %>%
  add_row(study_id = "50", medication_name = "spironolactone", dosage = "100mg", SSRI = "0")

medications <- medications %>%
  add_row(study_id = "16", medication_name = "Lamotrigine", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "16", medication_name = "Albuterol", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "16", medication_name = "Alprazolam", dosage = NA, SSRI = "0")

medications <- medications %>%
  add_row(study_id = "3", medication_name = "Esterylla birth control", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "3", medication_name = "Uquora", dosage = NA, SSRI = "0") %>%
  add_row(study_id = "3", medication_name = "Seed probiotic", dosage = NA, SSRI = "0")

medications <- medications %>%
  arrange(study_id)

ssri_ids <- medications %>%
  filter(SSRI == "1") %>%
  pull(study_id) %>%
  unique()
#there are 11 of them


adderall_ids <- medications %>%
  filter(str_detect(medication_name, regex("Adderall", ignore_case = TRUE))) %>%
  distinct(study_id) %>%
  pull(study_id) %>%
  unique()










