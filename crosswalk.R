### Run crosswalk package 
library(cogxwalkr)


getscore <- function(scendat, itemnms, scorenm)
{
    data <- scendat %>% select(itemnms) 
    mirtout <- mirt(data, model = 1, itemtype = c("graded"))
    score <- fscores(mirtout, full.scores.SE = FALSE)
    colnames(score) <- scorenm
    scendat <- cbind.data.frame(scendat, score)
        ## left_join(cbind.data.frame(score, ID = data$ID), by = "ID")
    return(scendat)
}

getcrosswalk <- function(groupdata, cog1, cog2)
{
   groupdata <- groupdata %>% select(cog1, cog2, "Dementia")
   crswkout <- crosswalk(cog1, cog2, condition_by = "Dementia",
                         condition_loop = TRUE, data = groupdata)
    ##control = list(nboot = 100, seed = 42, ncores = 4L)
}


crswkoutlist <- outdfs <- datasets <- vector("list", 4)
for(scen in 1:4)
{
    itemnmslist <- list(item_lists[[combos[s == scen & g == 1, id]]],
                        item_lists[[combos[s == scen & g == 2, id]]])
    groupdatalist <- list(scenario_data[[scen]] %>% filter(Group == "Group 1") %>%
                          getscore(itemnmslist[[1]], "score1") %>%
                          getscore(itemnmslist[[2]], "score2"),
                          scenario_data[[scen]] %>% filter(Group == "Group 2") %>%
                          getscore(itemnmslist[[1]], "score1") %>%
                          getscore(itemnmslist[[2]], "score2"))
    datasets[[scen]] <- groupdatalist
    crswkoutlist[[scen]] <- list(getcrosswalk(groupdatalist[[1]], "score2", "score1"),
                                 getcrosswalk(groupdatalist[[2]], "score1", "score2"))
    edufits <- list(lm(score2 ~ edu, groupdatalist[[2]]),
                    lm(score1 ~ edu, groupdatalist[[1]]))
    coefs = sapply(1:2, function(i)
    {
        coef(crswkoutlist[[scen]][[i]]$fit)*coef(edufits[[i]])["edu"]
    })
    outdfs[[scen]] <- cbind.data.frame(scenario = scen, crosswalk_to = c("Group 2", "Group 1"),
                                       coefs = coefs)
}
allout <- do.call(rbind.data.frame, outdfs)
