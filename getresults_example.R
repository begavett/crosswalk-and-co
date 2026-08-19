## grab results files, combine, and save
library(tidyverse)

datadir <- "/home/homedir/"

d <- .4
n <- 5000
reps <- 1000
resultlist <- lapply(1:reps, function(rep)
{
    if (!file.exists(paste0(datadir, "results/withci_n", n, "_d", d*10, "/result", rep, ".RDS"))) {
        return()
    } else {
        readRDS(paste0(datadir, "results/withci_n", n, "_d", d*10, "/result", rep, ".RDS")) %>%
            mutate(rep = rep,
                   n_sample = n,
                   b1 = d)
    }
})
allres_ci <- do.call(rbind.data.frame, resultlist)

missnums <- which(sapply(resultlist, is.null))
print(paste0(missnums, collapse = ","), quote = F)

runsum <- do.call(rbind.data.frame, readRDS(paste0(datadir, "results/noci_n", n, "_d", d*10, "/run_summaryn", n, "_d", d*10, ".RDS")))


saveRDS(allres_ci, paste0(datadir, "results/allresn", n, "_d", d*10, "_ci_jun2026.RDS"))




