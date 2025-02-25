library(tidyverse)
library(GenomicRanges)

# Merge data sets together -- mergedData has had 454 and Illumina data previously merged and illuminaOnlyData 
# includes data patient data for which only Illumina data has been collected.
#---------------------------------------------------------------------------------------------------------------------------
mergedData <- bind_rows(lapply(list.files('data', pattern = 'intSites.csv', recursive = TRUE, full.names = TRUE), read_csv))
mergedData$patient <- toupper(mergedData$patient)
mergedData <- select(mergedData, patient, seqnames, start, strand, dataSource, GTSP, cellType, timePoint, reads, estAbund, relAbund)
names(mergedData) <- c('patient', 'seqnames', 'position', 'strand', 'dataSource', 'internalSampleID', 'cellType', 'timePoint', 'reads', 'estAbund', 'relAbund')


illuminaOnlyData <- read_tsv('data/intSites.tsv', show_col_types = FALSE)
illuminaOnlyData$patient <- toupper(illuminaOnlyData$subject)
illuminaOnlyData <- subset(illuminaOnlyData, ! patient %in% mergedData$patient)
illuminaOnlyData$dataSource <- 'Illumina'
illuminaOnlyData <- select(illuminaOnlyData, patient, chromosome, position, strand, dataSource, internalSampleID, cellType, timePoint, reads, estAbund, relAbund)
names(illuminaOnlyData) <- c('patient', 'seqnames', 'position', 'strand', 'dataSource', 'internalSampleID', 'cellType', 'timePoint', 'reads', 'estAbund', 'relAbund')

sites <- bind_rows(mergedData, illuminaOnlyData)
sites$posid <- paste0(sites$seqnames, sites$strand, sites$position)



# Overlap analysis
#---------------------------------------------------------------------------------------------------------------------------
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

o <- select(sites[! is.na(sites$siteGroup),], siteGroup, patient, dataSource, posid, internalSampleID, cellType, timePoint, reads, estAbund, relAbund) %>% arrange(siteGroup) 

openxlsx::write.xlsx(o, 'overlappingSites.xlsx')
