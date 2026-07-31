library(luckyBase)
Plus.library(c('dplyr','digest'))

# Detect OS and set root_path accordingly
sysname <- Sys.info()[["sysname"]]
if (identical(sysname, "Windows")) {
  # Use the root of drive E on Windows
  root_path <- "E:/"
  root_path_sync <- "E:/Sync/"
  n_cores <- 16
} else if (identical(sysname, "Darwin")) {
  # Use the mounted volume path on macOS
  root_path <- "/Volumes/2T01/winE/"
  root_path_sync <- "/Volumes/share/docker/Sync/data/Sync/"
  n_cores <- max(1, parallel::detectCores(logical = TRUE) - 3)
} else {
  stop(sprintf("Unsupported OS: %s", sysname))
}

project_version <- "PADv20250720"
data_version <- "PADv20240810"
model_version <- "PADv20240911"
project.dir <- paste0(root_path, "iProjects/RCheck/GSClassifier/routine01/ccs/", project_version)

# Data (RNA-expression; not TSP)
data_all <- readRDS(paste0(root_path, "iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_", data_version, ".rds"))

# Model
resCCS <- readRDS(paste0(root_path_sync, "Project/", project_version, "/models/resCCS_", paramMD5, ".rds"))
resCCS@Repeat$model.dir <- paste0(root_path, "iProjects/RCheck/GSClassifier/routine01/ccs/", model_version)
print(resCCS@Data$scaller.performance[, 1:3])
print(length(unique(resCCS@Data$CCS)))

# Params of DR and clustering
# paramMD5 <- "5ff3a2de76e6cf902e765e8224f9cb66"
# record <- "0=15,19,22,25,26,31,32,36,50,56,66,68,72,81,82,85,89,90,92,95,99,102,104,105,116,122,142,143,147,151,152,153,154,155,158,163,167,168,172,175,177,191,197,206; 1=140,6,10,134,33,30,138,45,2,3,5,11; 2=8,54,65,51,27,28,21,18,13,34,17,107,121,75,78,144,156,162,24,42,12,55,29; 3=161,79,39,164,53,69,84,61,149,216,14,94,157,83,16,137,20,44,63,49,48,71,47,141,46,77,96,58"
# df_param <- readRDS(paste0(root_path, "iProjects/RCheck/GSClassifier/routine01/ccs/PADv20250629/fire-turning/dr&clustering/PADv20250629_Parameters_DR&Clustering.rds"))
# df_param[df_param$paramMD5 %in% paramMD5, ]
# record_version <- digest(record, algo = "md5") %>% substr(., 1, 8)

