
# load libraries and functions
library(dplyr); packageVersion("dplyr")
source("/home/briana.whitaker/txeg/SourceTracker.r")
#source("/home/briana.whitaker/txeg/SourceTracker_modSubFunctions.r")

# load map_grad
load("/home/briana.whitaker/txeg/TXEG_SourceTracker_mapGradient_noMonocot.RData")

# load feast_otu_sub (same as local analysis)
load("/home/briana.whitaker/txeg/TXEG_SourceTracker_otuTable_noMonocot.RData")

# print minDepth
minDepth <- min(rowSums(feast_otu_sub)); minDepth 

######
siteyears_grad <- levels(map_grad$SiteYear)
out_grad <- list()
for(i in (1:length(siteyears_grad))) {
    require(dplyr)
    source("/home/briana.whitaker/txeg/SourceTracker.r")
    #source("/home/briana.whitaker/txeg/SourceTracker_modSubFunctions.r")
    minDepth <- min(rowSums(feast_otu_sub))
    name <- siteyears_grad[i]  #SiteYear code
    #map file
    map <- rbind(map_grad %>% filter(Env == name),  #subsets the Phal which has the site as its name (1-at-a-time)
                      map_grad %>% filter(SourceSink == "source")) #whole gradient sources
    rownames(map) <- map$SampleID
    map <- droplevels(map)  #drop excess sink levels here
    # samples as rows otu subset
    otu <- feast_otu_sub[rownames(feast_otu_sub) %in% map$SampleID,]
    # need to be in the same sample order!
    otu <- otu[order(rownames(otu)), ]
    map <- map[order(rownames(map)), ]
    #identical(rownames(otu), rownames(map)) #TRUE
    train.ix <- which(map$SourceSink=='source')
    test.ix <- which(map$SourceSink=='sink')
    envs <- map$Env
    # if we don't use tuning (this is what we finally settled on, after all that back-n-forth)
    alpha1 <- alpha2 <- 0.001
    # train SourceTracker object on training data
    st <- sourcetracker(otu[train.ix,], envs[train.ix],
                           rarefaction_depth =  minDepth)
    # Estimate source proportions in test data
    results <- predict(st, otu[test.ix,], nrestarts=20,  #20 or 100
                      alpha1 = alpha1, alpha2 = alpha2,
                      rarefaction_depth =  minDepth, full.results = TRUE)
    obj <- results$full.results
    out_grad <- base::rbind(out_grad, cbind("props" = c(unname(results$proportions)), 
                  "props_sd" =  c(unname(results$proportions_sd)),
                  "sourceType" = results$train.envs,
                  "siteyear" = rep(name, length(results$train.envs)),
          "draws" = rep( dim(results$draws)[1], length(results$train.envs)), 
          "alpha1" = rep(alpha1, length(results$train.envs)),
          "alpha2" = rep(alpha2, length(results$train.envs)) ))
    save(obj, file = base::paste("/home/briana.whitaker/txeg/", name,
                                 "_asvList_gradient_iter1.RData", sep = "")) #saves 'obj' object
}

write.csv(out_grad, "/home/briana.whitaker/txeg/TXEG_SourceTracker_2013_20draws_Gradient_No1.csv")


