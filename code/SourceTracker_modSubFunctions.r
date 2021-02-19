## Core Functions 
# sourcetracker
# predict.sourcetracker
# plot.sourcetracker.fit
## Internal Functions
# run.gibbs
# tune.st    ----->    MODIFIED
# eval.fit    ----->    MODIFIED
# rarefy
# plot.sourcetracker.pie, plot.sourcetracker.bar, plot.sourcetracker.dist, sourcetracker.error.bars, plot.eval
# save.mapping.file
# sortmatrix, jsdmatrix, jsd, kld


######## MODIFIED

# modified to allow different rarefaction depth and ntrials
eval.fit.mod <- function(otus, envs, ntrials, rarefaction_depth, 
             individual.samples=TRUE, verbosity=1, ...){
    train.envs <- sort(unique(envs))
    V <- length(train.envs)
    env.sizes <- table(envs)
    
    # make sure each pair of envs gets picked
    # build up all pairs of samples, each column is a pair
    # each source env gets to be first and second sample once
    # pairs <- expand.grid(1:V,1:V)
    # pairs <- pairs[pairs[,1]!=pairs[,2],]
    # make nreps pairs randomly
    pairs <- NULL
    for(i in 1:ntrials){
        pairs <- rbind(pairs, sample(V,size=2))
    }
    
    mixtures <- runif(ntrials)
    y <- matrix(0,nrow=ntrials, ncol=V+1)
    yhat <- matrix(0,nrow=ntrials, ncol=V+1)
    yhat.sd <- matrix(0,nrow=ntrials, ncol=V+1)
    colnames(y) <- c(as.character(train.envs),'Unknown')
    colnames(yhat) <- c(as.character(train.envs),'Unknown')
    newsamples <- NULL
    allenvs <- NULL
    for(i in 1:ntrials){
        env1 <- pairs[i,1]
        env2 <- pairs[i,2]
        allenvs <- rbind(allenvs, c(env1, env2))
        if(verbosity > 1){
           cat(sprintf('%d of %d: %.2f*%s + %.2f*%s: \n',i,ntrials,mixtures[i], train.envs[env1],1-mixtures[i], train.envs[env2]))
        } else if(verbosity > 0){
            cat('.')
        }

        # all indices of each environment
        env1.ix.all <- which(envs == train.envs[env1])
        env2.ix.all <- which(envs == train.envs[env2])
        
        if(individual.samples){
            # get one sample from each env
            # cast as list so that sample doesn't misinterpret a length-1 vector
            env1.ix <- sample(as.list(env1.ix.all),size=1)[[1]]
            env2.ix <- sample(as.list(env2.ix.all),size=1)[[1]]
            
            # train sourcetracker, hold out entire second env. and first env. sample
            # note: don't hold out first sample if that env has only one sample
            if(length(env1.ix.all) == 1){
                st <- sourcetracker(otus[-env2.ix.all,], envs[-env2.ix.all], rarefaction_depth=rarefaction_depth)
            } else {
                st <- sourcetracker(otus[-c(env1.ix,env2.ix.all),], envs[-c(env1.ix,env2.ix.all)], rarefaction_depth=rarefaction_depth)
            }
            
            # make fake sample, weighted mixture of two source samples
            s1 <- otus[env1.ix,]
            s2 <- otus[env2.ix,]
            
        } else {
            # train sourcetracker, hold out entire second env.
            st <- sourcetracker(otus[-env2.ix.all,], envs[-env2.ix.all], rarefaction_depth=rarefaction_depth)
            
            # make fake sample as mixture of _environment_ means
            s1 <- colSums(rarefy(otus[env1.ix.all,], maxdepth=rarefaction_depth))
            s2 <- colSums(rarefy(otus[env2.ix.all,], maxdepth= gv  ))
        }
        
        newsample <- mixtures[i] * s1/sum(s1) + (1-mixtures[i]) * s2/sum(s2)
        newsample <- round(100000 * newsample)
        newsample <- matrix(newsample, nrow=1)
        newsample <- rarefy(newsample,maxdepth=ceiling(sum(s1+s2)/2))
        newsamples <- rbind(newsamples, newsample)
        y[i,env1] <- mixtures[i]
        y[i,V+1] <- 1-mixtures[i]
        
        # test on fake sample
        results <- predict(st, newsample, rarefaction_depth=rarefaction_depth, verbosity=verbosity-1, ...)
        for(j in 1:ncol(results$proportions)){
            whichenv <- which(colnames(yhat) == colnames(results$proportions)[j])
            yhat[i,whichenv] <- results$proportions[,j]
            yhat.sd[i,whichenv] <- results$proportions_sd[,j]
        }
    }

    # calculate RMSE
    se <- as.numeric((y[,-V] - yhat[,-V])**2)
    mse <- mean(se)
    se.sem <- sd(se)/sqrt(length(se))
    rmse <- mse**.5
    rmse.sem <- se.sem**.5
    
    return(list(y=y,yhat=yhat,yhat.sd=yhat.sd,newsamples=newsamples, 
            env.pairs=allenvs, train.envs=train.envs, rmse=rmse, rmse.sem=rmse.sem))
}

#modified to test multiple values of alpha1, but not alpha2 or beta priors
## also inherits eval.fit.mod to allow for precise rarefaction depth
tune.st.mod <- function(otus, envs, ntrials, rarefaction_depth,
            individual.samples=TRUE, alpha1=10**(-3:0), alpha2=1e-1, 
            beta=10, verbosity=0, ...){
    results <- list()
    alphas <- expand.grid(rev(alpha1), alpha2)
    colnames(alphas) <- c('alpha1','alpha2')
    rmse <- numeric(nrow(alphas))
    rmse.sem <- numeric(nrow(alphas))
    for(i in 1:nrow(alphas)){
        cat(sprintf('Loop %d of %d, alpha1=%f, alpha2=%f ',i,nrow(alphas),alphas[i,1], alphas[i,2]))
        if(verbosity > 2) cat('\n')
        results[[i]] <- eval.fit.mod(otus, envs, ntrials, rarefaction_depth,
                                individual.samples=individual.samples,
                                alpha1=alphas[i,1], alpha2=alphas[i,2], beta=beta, verbosity=verbosity-1, ...)
        rmse[i] <- results[[i]]$rmse
        rmse.sem[i] <- results[[i]]$rmse.sem
        if(verbosity > 0) cat(sprintf('RMSE = %.3f +/- %.3f\n',rmse[i], rmse.sem[i]))
    }
    ## choose alpha as most conservative value of alpha2 (smallest)
    ## then most conservative value of alpha1 (largest)
    ## that gives pseudo-r2 within 1 sem of the max.
    # best.ix <- min(which(rmse <= min(rmse + rmse.sem)))
    ## Alternative: simply choose the lowest rmse
    best.ix <- which.min(rmse)
    best.rmse <- rmse[best.ix]
    best.alpha1 <- alphas[best.ix,1]
    best.alpha2 <- alphas[best.ix,2]
    return(list(alphas=alphas, rmse=rmse, rmse.sem=rmse.sem, best.rmse=best.rmse, 
                best.alpha1=best.alpha1, best.alpha2=best.alpha2, 
                results=results))
}