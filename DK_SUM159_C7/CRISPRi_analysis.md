SUM 159 Clone 7 CRISPRi Analysis
================

Taking a look at the CRISPRi screen for the dark kinases in SUM 159 cells.

Loading Data
============

``` r
alignments  = rbind(
  read.delim(here('DK_SUM159_C7/salmon_output/SUM159_C7_D0_S4_R1_001.fastq.gz/quant.sf')) %>%
    mutate(Day = 0),
  read.delim(here('DK_SUM159_C7/salmon_output/SUM159_C7_D14_S5_R1_001.fastq.gz/quant.sf')) %>%
    mutate(Day = 14),
  read.delim(here('DK_SUM159_C7/salmon_output/SUM159_C7_D28_S6_R1_001.fastq.gz/quant.sf')) %>%
    mutate(Day = 28)
) %>% mutate(gene_name = str_extract(Name,"[^_]*"), 
             gene_name_plus_id = str_extract(Name,"[^_]*_[^_]*_[^_]*"),
             cell_line = "SUM159",
             clone = 7)
```

Calculations and Conversions
============================

The number of reads observed for each day is variable, so instead of doing the downstream calculations on the read count, we will instead convert these into the percentage of reads from that day. We'll also pick out the percentage of reads from Day 0 and use that to calculate the shifts in reads between Day 0 and Day 14/28.

``` r
#Count the total reads per day, we will use this correct for the depth of sequencing coverage
alignments = calc_reads_percent_per_day(alignments)
```

    ## Joining, by = c("Day", "clone")

``` r
#Pick out the first day reads and determine the ratio of first day to other days
alignments = ratio_to_day_0(alignments)
```

    ## Joining, by = "Name"

``` r
#remove a few unused columns
alignments = alignments %>%
  select(-TPM,-Length,-EffectiveLength)

write_rds(alignments,here('DK_SUM159_C7/read_counts.rds'))
```

Non-target Analysis
===================

With the ratios calculated, my next step is to take a look at the non-targetting sequences. These sequences shouldn't be targetting any genetic locus and so can act as a negative control for the other targeting sequences.

``` r
non_target_counts = alignments %>% filter(gene_name == "non-targeting")

ggplot(non_target_counts, aes(x=Day,y=Percentage_reads, color=gene_name_plus_id)) + 
  geom_line() + theme_berginski() + labs(x="Days",y="Percentage of Total Reads",color="Target ID")
```

![](CRISPRi_analysis_files/figure-markdown_github/non-targeting-percentage-1.png)

The raw percentage counts look good, all roughly stable without any changes. Let's look at the ratios with day 0.

``` r
ggplot(non_target_counts, aes(x=Day,y=D0_ratio, color=gene_name_plus_id)) + 
  geom_line() + theme_berginski() + labs(x="Days",y="Day 0 Ratio",color="Target ID")
```

![](CRISPRi_analysis_files/figure-markdown_github/non-targeting-ratio-1.png)

``` r
ggplot(non_target_counts, aes(x=as.factor(Day),y=D0_ratio)) + 
  geom_boxplot() + theme_berginski() + labs(x="Days",y="Day 0 Ratio")
```

![](CRISPRi_analysis_files/figure-markdown_github/non-targeting-ratio-2.png)

One last point, let's look at the confidence intervals for the ratio values. I'll use the bootstrap confidence interval.

``` r
non_target_conf = non_target_counts %>% 
  #Toss out day one as the ratios are meaningless there
  filter(Day != 0) %>% 
  group_by(Day) %>% 
  summarise(boot.low = boot.ci(boot(D0_ratio, function(data,index) {mean(data[index])}, 10000), type="bca", conf=0.95)$bca[4],
            boot.high = boot.ci(boot(D0_ratio, function(data,index) {mean(data[index])}, 10000), type="bca", conf=0.95)$bca[5],
            quantile.low = boot.ci(boot(D0_ratio, function(data,index) {quantile(data[index],c(0.05))}, 10000), type="bca", conf=0.0001)$bca[4])

non_target_conf %>%
  kable()
```

<table>
<thead>
<tr>
<th style="text-align:right;">
Day
</th>
<th style="text-align:right;">
boot.low
</th>
<th style="text-align:right;">
boot.high
</th>
<th style="text-align:right;">
quantile.low
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.7495567
</td>
<td style="text-align:right;">
0.9803475
</td>
<td style="text-align:right;">
0.4758460
</td>
</tr>
<tr>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
0.8454633
</td>
<td style="text-align:right;">
1.0907758
</td>
<td style="text-align:right;">
0.5718818
</td>
</tr>
</tbody>
</table>
OK, this looks good, the confidence intervals include 1. Let's move onto the rest of the data.

Gene Ratio Values
=================

As a first pass, let's take a look at the distribution of day 0 ratios for both day 14 and 28.

``` r
#remove the non-targeting alignments and day zero ratios
alignments_gene = alignments %>% filter(gene_name != "non-targeting" & Day != 0)
ggplot(alignments_gene, aes(x=D0_ratio, y = stat(density), color=as.factor(Day))) + 
  geom_freqpoly(bins=50) + theme_berginski() +
  labs(x="Ratio with Day 0",y="Density") +
  xlim(c(0,5))
```

    ## Warning: Removed 5 rows containing non-finite values (stat_bin).

    ## Warning: Removed 4 rows containing missing values (geom_path).

![](CRISPRi_analysis_files/figure-markdown_github/Day0-dist-1.png)

Looks like the spread increases by day 28, which should be expected due to each of the sequences having more time to act on the growth phenotype. Is there are a relationship between the ratios at day 14 and 28?

``` r
alignments_day_comp = alignments_gene %>% select(Name,Day,D0_ratio,gene_name_plus_id) %>% spread(Day,D0_ratio,sep="_")
ggplot(alignments_day_comp, aes(x=Day_14,y=Day_28)) + 
  geom_point(alpha=0.1) + geom_density2d() + geom_smooth(method="lm") +
  theme_berginski() +
  labs(x ="Day 14/Day 0", y="Day 28/Day 0")
```

![](CRISPRi_analysis_files/figure-markdown_github/Day14-vs-28-1.png)

I'm not sure how correlated I expected this to be, but the 0.76 value seems a bit low, but maybe that's just par for this type of data.

Determining the Screening Hits
==============================

There is no fixed method for determining which genes are "hits" given the ratios observed. This difficulty is compounded by the fact that some of these guide sequences simply won't be effective at reducing gene expression. Thus, we can't differtiate between an ineffective sequence and an effective sequence that doesn't affect growth. Despite that, we can still look to see how many sequences associated with a given gene show up below a given ratio threshold, as a first attempt I've used a 2-fold decrease.

``` r
alignment_summary = alignments %>%
  filter(Day != 0, Day != 14) %>%
  group_by(gene_name,Day) %>%
  summarise(number_below_2fold = sum(D0_ratio < 0.5, na.rm=T),
            number_above_2fold = sum(D0_ratio > 2, na.rm=T),
            average_ratio = mean(D0_ratio),
            p_val = tidy(wilcox.test(D0_ratio,non_target_counts$D0_ratio[non_target_counts$Day == 28]))$p.value)
```

    ## Warning in wilcox.test.default(D0_ratio,
    ## non_target_counts$D0_ratio[non_target_counts$Day == : cannot compute exact
    ## p-value with ties

Now, let's filter this down to just those with at least two sequences below 2-fold decrease.

``` r
dropout_hits = alignment_summary %>% 
  filter(number_below_2fold > 1) %>% 
  arrange(desc(number_below_2fold)) %>% 
  select(gene_name,number_below_2fold)

write_rds(dropout_hits,here('DK_SUM159_C7/dropout_hits.rds'))

dropout_hits %>%
  kable(col.names = c("Gene Name","Number Sequences 2-Fold Decrease"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
Gene Name
</th>
<th style="text-align:right;">
Number Sequences 2-Fold Decrease
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PI4KA
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
PKMYT1
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
RIOK2
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
PRPF4B
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
TP53RK
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
CDK12
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
NEK1
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
RIOK1
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
PIP5K1A
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
PIK3C2B
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
SBK3
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
TLK2
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
CAMKK1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
CAMKK2
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
CDK10
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
CLK3
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
HIPK1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
HIPK4
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
ICK
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRK1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
MAP3K15
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
MAST3
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
NEK8
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
PAK7
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
PIP4K2C
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
SCYL1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
STK3
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
STK31
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
STK32B
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
TTBK1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
VRK2
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
WNK3
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>
.

Genes with guides with 2-fold Increase in Adundance
===================================================

``` r
dropin_hits = alignment_summary %>% 
  filter(number_above_2fold > 1) %>% 
  arrange(desc(number_above_2fold)) %>% 
  select(gene_name,number_above_2fold)

dropin_hits %>%
  kable(col.names = c("Gene Name","Number Sequences 2-Fold Increase"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
Gene Name
</th>
<th style="text-align:right;">
Number Sequences 2-Fold Increase
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PSKH1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
STK32C
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>
.

Visualizing the Screening Hits
==============================

    ## Joining, by = "gene_name"

    ## Joining, by = c("Day", "gene_name", "number_below_2fold")

    ## Warning: Removed 20 rows containing non-finite values (stat_bin).

![](CRISPRi_analysis_files/figure-markdown_github/alignment_hits_vis-1.png)

This is the same plot as above, execept with log10 scale. Note that since some of the ratios are zero (reads at day 28 were zero), I had to add a "psuedo" ratio min value to make the log10 scale work.

    ## Joining, by = "gene_name"

    ## Joining, by = c("Day", "gene_name", "number_below_2fold")

    ## Warning: Removed 20 rows containing non-finite values (stat_bin).

    ## Warning: Removed 2 rows containing missing values (geom_bar).

![](CRISPRi_analysis_files/figure-markdown_github/alignment_hits_vis_log10-1.png)
