
# load libraries and functions
library(dplyr); packageVersion("dplyr")
source("/home/briana.whitaker/txeg/SourceTracker.r")
#source("/home/briana.whitaker/txeg/SourceTracker_modSubFunctions.r")

# load map_grad (same as gradient analysis)
load("/home/briana.whitaker/txeg/TXEG_SourceTracker_mapGradient_noMonocot.RData")

# load feast_otu_sub (same as local analysis)
load("/home/briana.whitaker/txeg/TXEG_SourceTracker_otuTable_noMonocot.RData")

# load spatial distance file [spat_2013]
load("/home/briana.whitaker/txeg/TXEG_SiteDistances_2013_list.RData")

# print minDepth
minDepth <- min(rowSums(feast_otu_sub)); minDepth 



######
siteyears_reg <- levels(map_grad$SiteYear)
out_reg <- list()
for(i in (1:length(siteyears_reg))) {
    require(dplyr)
    source("/home/briana.whitaker/txeg/SourceTracker.r")
    #source("/home/briana.whitaker/txeg/SourceTracker_modSubFunctions.r")
    minDepth <- min(rowSums(feast_otu_sub)) 
    name <- siteyears_reg[i]  #SiteYear code
    #identify regional sites by 125km distance cutoff
    regions.df <- spat_2013 %>% filter(cols == name) %>% filter(dist <= 125) # filter sites w/in 125km
    regions <- regions.df$rows
    avgSpatDist <- base::mean(regions.df$dist)
    sdSpatDist <- stats::sd(regions.df$dist)
    #map file
    map.0 <- map_grad[map_grad$SiteYear %in% name | map_grad$SiteYear %in% regions,]
    map <- rbind( map.0 %>% filter(Env == name),  #subsets the Phal which has the site as its name (1-at-a-time)
                  map.0 %>% filter(SourceSink == "source")  ) #regional sources, but not sinks
    rownames(map) <- map$SampleID
    map <- droplevels(map)  #drop excess sink levels here
    # samples as rows otu subset
    otu <- feast_otu_sub[rownames(feast_otu_sub) %in% map$SampleID,]
    # need to be in the same sample order!
    otu <- otu[order(rownames(otu)), ]
    map <- map[order(rownames(map)), ]
    #identical(rownames(otu), rownames(map)) #debugging check, must be TRUE
    train.ix <- which(map$SourceSink=='source')
    test.ix <- which(map$SourceSink=='sink')
    envs <- map$Env
    # if we don't use tuning (this is what we finally settled on, after all that back-n-forth)
    alpha1 <- alpha2 <- 0.001
    
    # train SourceTracker object on training data
    st <- sourcetracker(otu[train.ix,], envs[train.ix], 
                           rarefaction_depth =  minDepth) #close to min
    # Estimate source proportions in test data
    results <- predict(st, otu[test.ix,], nrestarts=20,  #20 or 100
                      alpha1 = alpha1, alpha2 = alpha2,
                      rarefaction_depth =  minDepth, full.results = TRUE)
    obj <- results$full.results
    out_reg <- base::rbind(out_reg, cbind("props" = c(unname(results$proportions)), 
                  "props_sd" =  c(unname(results$proportions_sd)),
                  "sourceType" = results$train.envs,
                  "siteyear" = rep(name, length(results$train.envs)),
          "draws" = rep( dim(results$draws)[1], length(results$train.envs)), 
          "alpha1" = rep(alpha1, length(results$train.envs)),
          "alpha2" = rep(alpha2, length(results$train.envs)),
          "numInRegion" = rep(length(regions), length(results$train.envs)),
          "avgSpatRegion" = rep(avgSpatDist, length(results$train.envs)),
          "sdSpatRegion" = rep(sdSpatDist, length(results$train.envs)) ))
    save(obj, file = base::paste("/home/briana.whitaker/txeg/", name, 
                                 "_asvList_regional_iter1.RData", sep = "")) #saves 'obj' object
}

write.csv(out_reg, "/home/briana.whitaker/txeg/TXEG_SourceTracker_2013_20draws_Regional_No1.csv")


