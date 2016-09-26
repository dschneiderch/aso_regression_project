library(tidyverse)
source("scripts/functions_modelfitting.R")

geoarea2 = "tuo"
ires = "500m"
fscasource = "aso"  #needs to be aso if ires < 500m
pathout = "output/splitsample-modeling/"
dir.create(pathout, recursive = TRUE)

allasomdls = readRDS(paste0(pathout, "asomdls_augment_", ires, ".rds"))
print("read rds done.")

allasomdls <- allasomdls %>% mutate(yr = substr(dte, 1, 4))

for (iyr in unique(allasomdls$yr)) {
  saveRDS(allasomdls %>% select(dte, yr, .id, asodte, phvaso_aug_glmmdl) %>% filter_(~yr ==
    iyr), paste0(pathout, "phvasomdls_augment_", iyr, "_", ires, ".rds"))
    # saveRDS(allasomdls %>% select(dte, yr, .id, asodte, phvaso_coef_glmmdl) %>%
    # filter_(~yr == iyr), paste0(pathout, 'phvasomdls_coef_', iyr, '_', ires,
    # '.rds'))
    saveRDS(allasomdls %>% select(dte, yr, .id, asodte, phvasofsca_aug_glmmdl) %>%
    filter_(~yr == iyr), paste0(pathout, "phvasofscamdls_augment_", iyr, "_",
    ires, ".rds"))
    # saveRDS(allasomdls %>% select(dte, yr, .id, asodte, phvasofsca_coef_glmmdl) %>%
    # filter_(~yr == iyr), paste0(pathout, 'phvasofscamdls_coef_', iyr, '_', ires,
    # '.rds'))
  }
  rm(allasomdls)

  print("yrly save done.")

  phvaso_aug = bind_rows(readRDS(paste0(pathout, "phvasomdls_augment_2013_", ires,
  ".rds")), readRDS(paste0(pathout, "phvasomdls_augment_2014_", ires, ".rds")),
  readRDS(paste0(pathout, "phvasomdls_augment_2015_", ires, ".rds")), readRDS(paste0(pathout,
    "phvasomdls_augment_2016_", ires, ".rds")))

    # phvaso_aug=allasomdls %>% select(-phvasofsca_aug_glmmdl)
    stats = phvaso_aug %>% filter(dte != asodte) %>% group_by(dte, asodte, yr) %>% unnest(phvaso_aug_glmmdl) %>%
    mutate(swe=ifelse(swe<0,NA,swe)) %>%
    summarise(rmse = rootmse(swe, swehat), pctrmse = rmse/mean(swe, na.rm = T) *
    100)

    # rm(phvaso_aug)

    # stats= allasomdls %>% filter(dte!=asodte) %>% group_by(dte,asodte,yr) %>%
    # unnest(phvaso_aug_glmmdl) %>% summarise(rmse=rootmse(swe,swehat),
    # pctrmse=rmse/mean(swe,na.rm=T)*100)

    readr::write_tsv(stats, path = paste0(pathout, "errors_allasodates_phvaso_", ires,
    ".txt"))


    stats %>% mutate(asoyr = substr(asodte, 1, 4)) %>% filter(yr != asoyr) %>% group_by(dte) %>%
    summarise(bestasodte = asodte[which.min(rmse)], bestrmse = min(rmse, na.rm = T),
    bestpctrmse = pctrmse[which.min(rmse)]) %>% readr::write_tsv(., path = paste0(pathout,
      "bestasodates_phvaso_", ires, ".txt"))

      print("phvaso done.")

    phvasofsca_aug = bind_rows(readRDS(paste0(pathout, "phvasofscamdls_augment_2013_",
      ires, ".rds")), readRDS(paste0(pathout, "phvasofscamdls_augment_2014_", ires,
      ".rds")), readRDS(paste0(pathout, "phvasofscamdls_augment_2015_", ires, ".rds")),
      readRDS(paste0(pathout, "phvasofscamdls_augment_2016_", ires, ".rds")))

      # phvasofsca_aug=allasomdls %>% select(-phvaso_aug_glmmdl)
      stats=
      phvasofsca_aug %>%
      filter(dte!=asodte) %>%
      group_by(dte,asodte,yr) %>%
      unnest(phvasofsca_aug_glmmdl) %>%
      mutate(swe=ifelse(swe<0,NA,swe)) %>%
      summarise(rmse=rootmse(swe,swehat),
      pctrmse=rmse/mean(swe,na.rm=T)*100)

      # rm(phvasofsca_aug)
      readr::write_tsv(stats, path = paste0(pathout, "errors_allasodates_phvasofsca_",
      ires, ".txt"))

      stats %>% mutate(asoyr = substr(asodte, 1, 4)) %>% filter(yr != asoyr) %>% group_by(dte) %>%
      summarise(bestasodte = asodte[which.min(rmse)], bestrmse = min(rmse, na.rm = T),
      bestpctrmse = pctrmse[which.min(rmse)]) %>% readr::write_tsv(., path = paste0(pathout,
        "bestasodates_phvasofsca_", ires, ".txt"))

        print("phvasofsca done.")
