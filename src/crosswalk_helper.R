getscore <- function(scendat, itemnms, scorenm)
{
  data <- scendat %>% select(all_of(itemnms)) 
  mirtout <- mirt(data, model = 1, itemtype = c("graded"),
                  technical = list(NCYCLES = 1000))
  score <- fscores(mirtout, full.scores.SE = FALSE)
  colnames(score) <- scorenm
  scendat <- cbind.data.frame(scendat, score)
  ## left_join(cbind.data.frame(score, ID = data$ID), by = "ID")
  return(scendat)
}

getcrosswalk <- function(groupdata, cog1, cog2)
{
  groupdata <- groupdata %>% select(all_of(c(cog1, cog2)), Dementia)
  crswkout <- crosswalk(cog1, cog2, condition_by = "Dementia",
                        condition_loop = TRUE, data = groupdata)
  ##control = list(nboot = 100, seed = 42, ncores = 4L)
}

