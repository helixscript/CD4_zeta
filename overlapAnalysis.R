library(tidyverse)
library(GenomicRanges)

sites <- readRDS('sites.rds')

sites$start <- sites$position
sites$end <- sites$position
sites$siteGroup <- NA

expandSitesBy <- 3

g <- makeGRangesFromDataFrame(sites, keep.extra.columns = TRUE) + expandSitesBy
r <- GenomicRanges::reduce(g, with.revmap = TRUE)

siteGroup <- 1
invisible(lapply(r$revmap, function(x){
  if(length(x) > 1){
    sites[sites$posid %in% unique(g[x]$posid),]$siteGroup <<- siteGroup
    siteGroup <<- siteGroup + 1
  }
}))

o <- select(sites[! is.na(sites$siteGroup),], siteGroup, patient, dataSource, posid, internalSampleID, cellType, timePoint, reads, estAbund, relAbund, miseqid) %>% arrange(siteGroup) 

openxlsx::write.xlsx(o, 'overlappingSites.xlsx')
