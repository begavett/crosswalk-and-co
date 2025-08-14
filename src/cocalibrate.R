cocalibrate <- function(rg_dat, 
                        fg_dat, 
                        rg_items,
                        fg_items) {
  
  library(mirt)
  library(dplyr)
  library(magrittr)
  
  linking_items <- intersect(rg_items, fg_items)
  
  if(length(linking_items) == 0) stop ("No linking items")
  
  mirt_rg <- mirt(data = rg_dat %>%
                    select(all_of(rg_items)),
                  model = 1,
                  itemtype = "graded",
                  technical = list(NCYCLES = 5000))
  
  itembank_rg <- data.frame(coef(mirt_rg, simplify = TRUE)$items) %>%
    mutate(item = rownames(.))
  
  linking_itembank_rg <- itembank_rg %>%
    filter(item %in% linking_items)
  
  # multipleGroup estimates parameters using all data,
  # which I don't think is right
  # I think we want to estimate parameters in Group 1 first
  # Otherwise the estimated paramters in a MG model will be based on
  # both G1 and G2 data
  # mirt_mg <- multipleGroup(data = bind_rows(rg_dat, fg_dat) %>%
  #                 select(all_of(union(rg_items, fg_items))),
  #               model = 1,
  #               group = c(rg_dat$Group, fg_dat$Group),
  #               itemtype = "graded",
  #               invariance = c(linking_items, "free_means", "free_vars"), SE = TRUE)
  # 
  # coef(mirt_mg, simplify = TRUE)
  # fscores(mirt_mg, full.scores.SE = TRUE)
  
  
  mirt_fg0 <- mirt(data = fg_dat %>%
                     select(all_of(fg_items)),
                   model = 1,
                   itemtype = "graded",
                   pars = "values",
                   technical = list(NCYCLES = 5000))
  
  mirt_fg0_partable <- mirt_fg0
  
  for(i in linking_items) {
    
    i_pars <- linking_itembank_rg %>%
      t() %>%
      data.frame() %>%
      mutate(par = rownames(.)) %>%
      select(all_of(i), par) %>%
      filter(!is.na(!!sym(i))) %>%
      slice(1:(n()-1)) %>%
      pull(par)
    
    for(p in i_pars) {
      mirt_fg0_partable$value[mirt_fg0_partable$item == i & mirt_fg0_partable$name == p] <- linking_itembank_rg[i, p]
      mirt_fg0_partable$est[mirt_fg0_partable$item == i & mirt_fg0_partable$name == p] <- FALSE
    }
    
    mirt_fg0_partable$est[mirt_fg0_partable$name == "MEAN_1"] <- TRUE
    mirt_fg0_partable$est[mirt_fg0_partable$name == "COV_11"] <- TRUE
  }
  
  mirt_fg <- mirt(data = fg_dat %>%
                     select(all_of(fg_items)),
                   model = 1,
                   itemtype = "graded",
                   pars = mirt_fg0_partable,
                  technical = list(NCYCLES = 5000))
  
  coef(mirt_rg, simplify = TRUE)
  coef(mirt_fg, simplify = TRUE)
  
  fscores_rg <- data.frame(fscores(mirt_rg, full.scores.SE = TRUE)) %>%
    mutate(Group = "Group 1")
  fscores_fg <- data.frame(fscores(mirt_fg, full.scores.SE = TRUE)) %>%
    mutate(Group = "Group 2")
  
  return(list(mirt_rg = mirt_rg, 
              mirt_fg = mirt_fg, 
              fscores_rg = fscores_rg, 
              fscores_fg = fscores_fg,
              linking_items = linking_items))
}
