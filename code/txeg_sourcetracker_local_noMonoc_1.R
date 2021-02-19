
# load libraries and functions
library(dplyr); packageVersion("dplyr")
source("/home/briana.whitaker/txeg/SourceTracker.r")
#source("/home/briana.whitaker/txeg/SourceTracker_modSubFunctions.r")

# load map_st
load("/home/briana.whitaker/txeg/TXEG_SourceTracker_mapFile_noMonocot.RData")

# load feast_otu_sub
load("/home/briana.whitaker/txeg/TXEG_SourceTracker_otuTable_noMonocot.RData")

# print minDepth
minDepth <- min(rowSums(feast_otu_sub)); minDepth 


######
siteyears <- levels(map_st$SiteYear)
out_loc <- list()
for(i in (1:length(siteyears))) {
    require(dplyr)
    source('/home/briana.whitaker/txeg/SourceTracker.r')
    #source('/home/briana.whitaker/txeg/SourceTracker_modSubFunctions.r')
    minDepth <- min(rowSums(feast_otu_sub)) 
    name <- siteyears[i]  #SiteYear code
    #map file
    map <- map_st %>% filter(SiteYear == name)  # one site at a time
    rownames(map) <- map$SampleID
    train.ix <- which(map$SourceSink=='source')
    test.ix <- which(map$SourceSink=='sink')
    envs <- map$Env
    # samples as rows otu subset
    otu <- feast_otu_sub[map$SampleID,]
    # if we tune alphas (takes a really really long time)
    #tune.results <- tune.st.mod(otu[train.ix,], envs[train.ix], 
    #                ntrials = 20, rarefaction_depth =  minDepth)  #recommended to use 25
    # if we don't use tuning (this is what we finally settled on, after all that back-n-forth)
    alpha1 <- alpha2 <- 0.001
    # train SourceTracker object on training data
    st <- sourcetracker(otu[train.ix,], envs[train.ix], 
                           rarefaction_depth =  minDepth) #close to min
    # Estimate source proportions in test data
    results <- predict(st, otu[test.ix,], burnin=100, nrestarts=20,  #20 or 100
                      alpha1 = alpha1, alpha2 = alpha2,
                      rarefaction_depth =  minDepth, full.results = TRUE)
    obj <- results$full.results
    out_loc <- base::rbind(out_loc, cbind("props" = c(unname(results$proportions)), 
                  "props_sd" =  c(unname(results$proportions_sd)),
                  "sourceType" = results$train.envs,
                  "siteyear" = rep(name, length(results$train.envs)),
          "draws" = rep( dim(results$draws)[1], length(results$train.envs)), 
          "alpha1" = rep(alpha1, length(results$train.envs)),
          "alpha2" = rep(alpha2, length(results$train.envs)) ))
    save(obj, file = base::paste("/home/briana.whitaker/txeg/", name,
                                 "_asvList_local_iter1.RData", sep = "")) #saves 'obj' object
}

write.csv(out_loc, "/home/briana.whitaker/txeg/TXEG_SourceTracker_2013_20draws_Local_No1.csv")

