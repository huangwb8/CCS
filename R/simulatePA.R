

#' @title Simulate Probabilistic Array
#' @description Simulate Probabilistic Array
#' @param nTissue the amount of tissues
#' @param nCohort the amount of cohorts
#' @param nSample the amount of samples
#' @param nCluster the amount of clusters
#' @param mean the mean of Gaussian noise
#' @param sd the standard deviation of Gaussian noise
#' @param seed random seed
#' @return A list of sample data, which includes the real types and related Probabilistic Array
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' @export
simulatePA <- function(nTissue = 15,
                       nCohort = 10,
                       nSample = 100,
                       nCluster = 4,
                       mean = 0, sd = 1,
                       seed = 2023){

  # Test
  if(F){
    nTissue = 15
    nCluster = 4
    nCohort = 10
    nSample = 100
    seed = 2023
    whichOne=1
    mean = 0
    sd = 0.5
  }

  # create a vector
  simulateRA_one <- function(nCluster=4,
                             whichOne=1,
                             mean = 0,
                             sd = 0.5,
                             seed = 1){
    vt <- rep(0,nCluster)
    vt[whichOne] <- 1
    set.seed(seed); vt <- softmax(vt + rnorm(nCluster, mean = mean, sd = sd))
    return(vt)
  }

  # PA colname
  dat_colname <- NULL
  for(i in 1:nTissue){
    for(j in 1:nCohort){
      dat_colname <- c(dat_colname, paste0('t',i,'c',j))
    }
  }

  # create an array
  set.seed(seed); seed1 <- sample(1:100000, nTissue, replace = F)
  list_dat <- list()
  for(i in 1:length(seed1)){ # i=1
    set.seed(seed1[i]); seed2 <- sample(100001:200000, nCohort, replace = F)
    for(j in 1:length(seed2)){ # j=1
      set.seed(seed2[j]); seed3 <- sample(200001:300000, nSample, replace = F)
      for(k in 1:length(seed3)){ # k=1

        set.seed(seed3[k]); seed4 <- sample(300001:400000, 2, replace = F)
        set.seed(seed4[1]); whichOne <- sample(1:nCluster, 1, replace = F)
        set.seed(seed4[2]); seed5 <- sample(400001:1000000, nCohort * nTissue, replace = F)
        dat <- sapply(seed5,function(x)simulateRA_one(nCluster,whichOne = whichOne,mean = mean,sd = sd, seed = x))
        rownames(dat) <- paste('type', 1:nrow(dat), sep = '')
        colnames(dat) <- dat_colname
        list_dat[[paste0('Tissue',i,collapse = '')]][[paste0('t',i,'c',j)]][[paste0('Sample',k,collapse = '')]] <- list(expr = dat, subtype = whichOne)
      }
    }
  }

  # Output
  return(list_dat)
}
