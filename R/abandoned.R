
# Abandoned

if(F){

  if(nrow(d4) > 1000){
    set.seed(seeds[6]);
    numComplete <- NbClust2(
      data = d4[sample(1:nrow(d4), 1000),],
      distance = "euclidean",
      min.nc = min.nc,
      max.nc= max.nc,
      method = "ward.D2",
      index = "all",
      verbose = verbose
    )
    dis <- dist(d4, method = "euclidean")
    hc <- hclust(dis, method =  "ward.D2")
    y2 <- cutree(hc, length(unique(numComplete$Best.partition)))
  } else {
    set.seed(seeds[6]);
    numComplete <- NbClust2(
      data = d4,
      distance = "euclidean",
      min.nc = min.nc,
      max.nc= max.nc,
      method = "ward.D2",
      index = "all",
      verbose = verbose
    )
    y2 <- numComplete$Best.partition
  }
}

