simCog <- function(itempars, itemnames, edu_prob_1 = .3, b0 = 0, b1 = .2, n_sample = 500, cogname = "Mem", groupname) {
  
  require(simDAG)
  require(dplyr)
  require(magrittr)
  
  itemtypes <- itempars %>%
    select(starts_with("d")) %>%
    rowwise() %>%
    mutate(nThresh = sum(!is.na(c_across(everything())))) %>%
    ungroup() %>%
    mutate(itemtype = ifelse(nThresh == 1, "dich", "graded")) %>%
    pull(itemtype)
  
  edu <- rbernoulli(n_sample, p = prop_high_edu, output = "numeric")
  var_edu <- prop_high_edu * (1 - prop_high_edu)
  sigma2 <- 1 - (b1^2 * var_edu)
  sigma <- sqrt(sigma2)
  
  thetas <- b0 + b1*edu + rnorm(n_sample, mean = 0, sd = sigma)
  thetas_df <- data.frame(theta = thetas, edu = edu)
  
  describe(thetas_df)


  lm(theta ~ edu, data = thetas_df) %>%
    summary()

  lm(theta ~ edu, data = thetas_df) %>%
    lm.beta::lm.beta() %>%
    summary()
  
  sim_item_resp <- mirt::simdata(a = data.matrix(itempars %>% pull(a1)),
                                      d = as.vector(itempars %>% dplyr::select(starts_with("d"))) %>%
                                        unlist() %>%
                                        as.numeric() %>%
                                        matrix(nrow = nrow(itempars %>% dplyr::select(starts_with("d"))), byrow = FALSE),
                                      itemtype = itemtypes,
                                      Theta = data.matrix(thetas)) %>%
    as_tibble() %>%
    set_names(ItemNames)
  
  cog_mirt <- mirt(data = sim_item_resp, model = paste0(cogname, ' = 1-', ncol(sim_item_resp)),
                      itemtype = "graded")
  
  sim_item_resp <- fscores(cog_mirt, full.scores.SE = TRUE) %>%
    data.frame() %>%
    set_names(paste0(cogname, c("_FS", "_FS_SE"))) %>%
    bind_cols(sim_item_resp)
  
  cog_fscores <- bind_cols(thetas_df, sim_item_resp) %>%
    mutate(Group = groupname)
  
  return(cog_fscores)
}
