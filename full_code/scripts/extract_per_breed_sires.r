## Runs in main directory to work out which sires should be kept for various breed-based analyses
sire_models <- read.csv("files/sire_models.csv")
sires_to_keep_all      <- sire_models[,1:2]
sires_to_keep_jersey   <- sire_models[sire_models$sire_jersey_propn >= 15/16, 1:2]
sires_to_keep_holstein <- sire_models[sire_models$sire_jersey_propn <=  1/16, 1:2]

write.table(sires_to_keep_all, "files/sires_to_keep_all", row.names = F, col.names = F)
write.table(sires_to_keep_jersey, "files/sires_to_keep_jersey", row.names = F, col.names = F)
write.table(sires_to_keep_holstein, "files/sires_to_keep_holstein", row.names = F, col.names = F)