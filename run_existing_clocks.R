source("mammalmethylclockR_fns.R")
source("age-transformations-zoller.R")

load("output/dolphin_data.Rdata")
# from https://github.com/jazoller96/mammalian-methyl-clocks/tree/main/MammalMethylClock-Package
# unzip, then R/sysdata.rds
load("data/clocks_metadatabase.rda")

betas = data.frame(x)
betas$DOLPHIN.ID = sample_info$DOLPHIN.ID
sample_info$Tissue = "Skin"
sample_info$SpeciesLatinName = "Tursiops truncatus"

yxs.list <- alignDatToInfo(info = sample_info, dat = betas, "DOLPHIN.ID", "DOLPHIN.ID")

ys <- yxs.list[[1]]
xs <- yxs.list[[2]]

ys$Tissue <- factor(ys$Tissue) # need to be factorized for package functions
ys$SpeciesLatinName <- factor(ys$SpeciesLatinName) # need 

OUTVAR <- "Age"
COLVAR <- "Tissue"
ys.output <- predictAge(xs, ys, tissue.names = c("Skin"), species.name = "Tursiops truncatus")
head(ys.output)

save(ys.output, file = "output/existing_clocks.Rdata")

## checks how many cpg sites in existing clocks are present in our data

tissue.names = "Skin"
species.name = "Tursiops truncatus"

rows_match <- intersect(
  unique(unlist(lapply(c(tissue.names, "ALL"), grep, clocks_metadatabase$Tissues))),
  unique(unlist(lapply(c(species.name, "ALL"), grep, clocks_metadatabase$SpeciesNames))))
clocks_metadatabase$ClockName[rows_match]

# Robeck2021 = BottlenoseSkinAge.LogLinear2
rowid = which(clocks_metadatabase$ClockName == "Coef.BottlenoseSkinAge.LogLinear2")
cpgs_Robeck2021 = clock_coefficients_list[[clocks_metadatabase$CoefficientsName[rowid]]]
cpgs_Robeck2021 = cpgs_Robeck2021[!is.na(cpgs_Robeck2021$Coef.BottlenoseSkinAge.LogLinear2), "var"]
cpgs_Robeck2021 = cpgs_Robeck2021[cpgs_Robeck2021 != "(Intercept)"]
length(cpgs_Robeck2021)
sum(!(cpgs_Robeck2021 %in% names(betas)))

# Unpub20xx = CetaceanTursiopsTruncatusSkinAge.LogLinear2
rowid = which(clocks_metadatabase$ClockName == "Coef.CetaceanTursiopsTruncatusSkinAge.LogLinear2")
cpgs_Zoller2025 = clock_coefficients_list[[clocks_metadatabase$CoefficientsName[rowid]]]
cpgs_Zoller2025 = cpgs_Zoller2025[!is.na(cpgs_Zoller2025$Coef.CetaceanTursiopsTruncatusSkinAge.LogLinear2), "var"]
cpgs_Zoller2025 = cpgs_Zoller2025[cpgs_Zoller2025 != "(Intercept)"]
cpgs_Zoller2025
length(cpgs_Zoller2025)
sum(!(cpgs_Zoller2025 %in% names(betas)))

# Universal = Universal3_Age.LogLinearRelAdult
rowid = which(clocks_metadatabase$ClockName == "Coef.Universal3_Age.LogLinearRelAdult")
cpgs_Universal = clock_coefficients_list[[clocks_metadatabase$CoefficientsName[rowid]]]
cpgs_Universal = cpgs_Universal[!is.na(cpgs_Universal$Coef.Universal3_Age.LogLinearRelAdult), "var"]
cpgs_Universal = cpgs_Universal[cpgs_Universal != "(Intercept)"]
cpgs_Universal
length(cpgs_Universal)
sum(!(cpgs_Universal %in% names(betas)))


