

model = FBA4Greatmod.generation(fba_mat = "inst/data/Ec_core.mat")

###
diet = readr::read_delim("inst/data/diets/vmh/EU_average.tsv",
                  delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE)
diet$Reaction = gsub("\\[", replacement = "\\(", diet$Reaction)
diet$Reaction = gsub("\\]", replacement = "\\)", diet$Reaction)
diet$`Flux Value` = -1*(diet$`Flux Value`)/24

model = setDiet(model, dietf = diet)
model = setDiet.name(model,diet_name = "EU_average")

writeFBAfile(model,"Ec_core_EU_average")


