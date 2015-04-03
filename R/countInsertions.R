## count number of insertions within exons or introns

read.count.tables <- function(file){
  # read in count table
  header = c("chr", "start", "end",
             "id", "mapq", "strand", 
             "cigar", "chr_feature", "start_feature",
             "end_feature", "feature", "strand_feature", 
             "ensembl_gene", "genesymbol", "ensembl_transcript", 
             "num")
  df <- read.table(file, sep="\t")
  colnames(df) <- header
  return(df)  
}

separate.silent.mutagenic.insertions <- function(df_insertions){
  # diff orientation --> muatgenic
  # same orientation --> silent
  mutagenic <- df_insertions[df_insertions$strand == df_insertions$strand_feature,]
  silent <- df_insertions[df_insertions$strand != df_insertions$strand_feature,]
  return(list(mutagenic, silent))
}

count.mutagenic.insertions <- function(count_exon, count_intron){
  mutagenic_insertions <- rbind(count_exon, count_intron)
  genes <- unique(mutagenic_insertions$genesymbol)
  mutagenic_insertions$genesymbol <- factor(mutagenic_insertions$genesymbol, levels=genes)
  count_table <- as.data.frame(table(mutagenic_insertions$genesymbol))
  # only count unique insertions
  total_insertions <- length(unique(mutagenic_insertions$id))
  count_table$no_insertions <- total_insertions - count_table$Freq
  return(count_table)
}

write2file <- function(count_table, filename){
  fn_out <- filename
  write.table(count_table, file=paste(fn_out),
              append=FALSE, quote=FALSE, 
              row.names=FALSE, col.names=TRUE, 
              dec=".", sep="\t")
}


