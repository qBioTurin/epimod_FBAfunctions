
setwd("~/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input")

library(readr)

load("./CDmodels/CD196HemeSink/CD196HemeSink.RData")

diets = list()
files <- list.files(path = "Diets/vmh/", pattern = "\\.tsv$")

for (i in 1:length(files)) {
  diet = read_delim(paste("Diets/vmh/", files[i], sep = ""), 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
  
  diet$Reaction = gsub("\\[", replacement = "\\(", diet$Reaction)
  diet$Reaction = gsub("\\]", replacement = "\\)", diet$Reaction)
  
  # flux in mmol/human*day
  
  diet$`Flux Value` = -1*(diet$`Flux Value`)/24
  
  # flux in mmol/human*h
  
  diets[[i]] = diet
  names(diets)[i] = unlist(strsplit(files[i], split='.', fixed=TRUE))[1]
  
}

diets[["Template"]] = diets[["Unhealthy"]]
diets[["Zero"]] = diets[["Unhealthy"]]
diets[["ZeroEssential"]] = diets[["Unhealthy"]]

# model.mat@react_id[which(model.mat@react_id %in% diets[["Template"]]$Reaction)]
diets[["Template"]]$`Flux Value`[which(diets[["Template"]]$Reaction %in% model.mat@react_id)] = 
  model.mat@lowbnd[which(model.mat@react_id %in% diets[["Template"]]$Reaction)]

diets[["Zero"]]$`Flux Value` = (diets[["Zero"]]$`Flux Value`)*1e-50

diets[["ZeroEssential"]]$`Flux Value`[which(diets[["Unhealthy"]]$Reaction %in% 
                                              c("EX_pro_L(e)", "EX_leu_L(e)", "EX_ile_L(e)", 
                                                "EX_val_L(e)", "EX_trp_L(e)", "EX_cys_L(e)"))] = rep(1e-50, 6)

save(diets, file = "Diets/diets.RData")
load("./Diets/diets.RData")
