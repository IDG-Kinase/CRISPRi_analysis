calc_reads_percent_per_day <- function(alignments) {
  alignments = alignments %>%
    group_by(Day, clone) %>%
    summarise(total_reads = sum(NumReads)) %>%
    select(Day,clone,total_reads) %>%
    left_join(alignments) %>%
    mutate(Percentage_reads = NumReads/total_reads)
  
  return(alignments)
}

ratio_to_day_0 <- function(alignments) {
  alignments = alignments %>%
    filter(Day == 0) %>%
    rename(Percentage_reads_D0 = Percentage_reads) %>%
    ungroup() %>%
    select(Name,Percentage_reads_D0) %>%
    left_join(alignments) %>%
    mutate(D0_ratio = Percentage_reads/Percentage_reads_D0)
  
  return(alignments)
}