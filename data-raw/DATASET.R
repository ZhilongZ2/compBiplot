## code to prepare `ECU.MF` dataset goes here

ECU.MF.raw <- read.csv("data-raw/ECU-MF.csv")

colnames(ECU.MF.raw) <-  c("sex","age","leukocytes","neutrophilsP","lymphocytesP",
                           "monocytesP","eosinophilsP","basophilsP","neutrophils",
                           "lymphocytes","NLR","monocytes","eosinophils",
                           "basophils","redbloodcells","mcv","mch","mchc",
                           "rdwP","hemoglobin","hematocritP","platelets","mpv","Condition")

dim(ECU.MF.raw)
sapply(ECU.MF.raw, is.numeric)
ECU.MF.raw$sex <- as.factor(ECU.MF.raw$sex)

sum(is.na(ECU.MF.raw))

ECU.MF.counts.raw <- ECU.MF.raw[, c(9,10,12:14)]

diff_wbc <- ECU.MF.raw$leukocytes - rowSums(ECU.MF.counts.raw)
rdiff_wbc <- abs(diff_wbc)/rowSums(ECU.MF.counts.raw)

exclude_index <- which(rdiff_wbc > 0.1)

# Filtering
ECU.MF.counts <- ECU.MF.counts.raw[-exclude_index, ]
ECU.MF <- ECU.MF.raw[-exclude_index, ]
diff_wbc_filtered <- ECU.MF$leukocytes - rowSums(ECU.MF.counts)
rdiff_wbc_filtered <- abs(diff_wbc_filtered)/rowSums(ECU.MF.counts)

# Remove redundant variables and categorical variable sex for later logratio PCA
ECU.MF$Condition <- factor(ECU.MF$Condition)
ECU.MF <- ECU.MF[, c(-1,-3,-4:-8,-11)]

usethis::use_data(ECU.MF, overwrite = TRUE)
