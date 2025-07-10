library(MammalMethylClock)
library(dplyr)

load("output/dolphin_data_10f.Rdata")

betas = data.frame(x)
betas$DOLPHIN.ID = sample_info$DOLPHIN.ID
sample_info$Tissue = "Skin"
sample_info$SpeciesLatinName = "Tursiops truncatus"

yxs.list <- alignDatToInfo(info = sample_info, dat = betas, "DOLPHIN.ID", "DOLPHIN.ID")

ys <- yxs.list[[1]]
xs <- yxs.list[[2]]

# tt = getClockDatabase()
# write.csv(tt, file = "output/ClockDatabase.csv", row.names = FALSE)

ys$Tissue <- factor(ys$Tissue) # need to be factorized for package functions
ys$SpeciesLatinName <- factor(ys$SpeciesLatinName) # need 

OUTVAR <- "Age"
COLVAR <- "Tissue"
ys.output <- predictAge(xs, ys, tissue.names = c("Skin"), species.name = "Tursiops truncatus")
head(ys.output)
DNAmAge.BottlenoseSkin
DNAmAge.BottlenoseSkinAge.LogLinear2

save(ys.output, file = "output/existing_clocks.Rdata")
