library(tidyverse)
library(RMySQL)

# Merge data sets together -- mergedData has had 454 and Illumina data previously merged and illuminaOnlyData
# includes data patient data for which only Illumina data has been collected.
# ---------------------------------------------------------------------------------------------------------------------------
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

dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- dbGetQuery(dbConn, 'select sampleName, miseqid from samples where sampleName like "%GTSP%"')
intSitesamples$sampleName <- gsub('\\-\\d+$', '', intSitesamples$sampleName)
intSitesamples <- distinct(intSitesamples)
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))

intSitesamples <- bind_rows(lapply(split(intSitesamples, intSitesamples$sampleName), function(x){
  if(n_distinct(x$miseqid) > 1) x <- tibble(sampleName = x$sampleName[1], miseqid = paste0(sort(unique(x$miseqid)), collapse = ' / '))
  x
}))

sites <- left_join(sites, intSitesamples, by = c('internalSampleID' = 'sampleName'))
saveRDS(sites, 'sites.rds')
