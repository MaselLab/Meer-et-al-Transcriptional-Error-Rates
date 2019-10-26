---
title: "Transcription Error rates R Notebook"
output: html_notebook
---

```{r}
#load data

ecoli_noCU <- read.delim("inal_ecoli_binned_coding_table.txt")
ecoli_CU<-read.delim("CU_binned_8_table.txt")
ecoli_CU<-subset(ecoli_CU, select = -c(T, lifetime))
ecoli_CU<-merge(ecoli_CU, unique(ecoli_noCU[c("gene", "groEL")]), by = "gene")
ecoli<-rbind(ecoli_CU, ecoli_noCU)
```


```{r}
#Testing whether the slope with abundance improves the model
ecoli_noCU$RlogA<-ecoli_noCU$R*log(ecoli_noCU$abundance)
starting_vals = c(10^-5)
summary(no_slope_glm<-glm(E ~R+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
starting_vals = c(10^-5, -10^-6)
summary(single_intercept_glm<-glm(E ~R + RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
anova(no_slope_glm, single_intercept_glm, test = "Chisq")

```


```{r}
#Testing whether each error type should have a separate intercept
starting_vals = c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5)
summary(no_cond_inter_glm<-glm(E ~subtype:R+ RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
anova(single_intercept_glm, no_cond_inter_glm, test = "Chisq")
```


```{r}
#Testing whether each condition should have a separate intercept
starting_vals = c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5)
summary(no_cond_inter_glm<-glm(E ~subtype:R+ RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
starting_vals = c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5)
summary(cond_w_inter_glm<-glm(E ~subtype:R+condition:R + RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
anova(no_cond_inter_glm, cond_w_inter_glm, test = "Chisq")
```

```{r}
#Testing whether each condition should have a separate slope
starting_vals = c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5 , 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, -10^-6, -10^-6, -10^-7, -10^-8)
summary(cond_slope_glm<-glm(E ~subtype:R+condition:R + condition:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
anova(cond_w_inter_glm, cond_slope_glm, test = "Chisq")
```


```{r}
#testing GA
starting_vals = c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5, -10^-6, -10^-6)
ecoli_noCU$isGA<-ifelse(ecoli_noCU$subtype == "GA", "GA", "nonGA")
ga_slope_glm<-glm(E ~subtype:R+condition:R + isGA:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU)
starting_vals = c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5)
pooled_glm<-glm(E ~subtype:R+condition:R + RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU)
anova(pooled_glm, ga_slope_glm, test = "Chisq")

```

```{r}
#getting slopes for GA and non-GA
starting_vals = c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5)
summary(ga_glm<-glm(E ~condition:R + RlogA+0,start = starting_vals, family = poisson(link = identity), data = subset(ecoli_noCU, subtype == 'GA')))
starting_vals = c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, -10^-6, -10^-6)
summary(no_ga_glm<-glm(E ~subtype:R+condition:R + subtype:RlogA+0,start = starting_vals, family = poisson(link = identity), data = subset(ecoli_noCU, subtype != 'GA')))
```


```{r}
#testing whether the model is improved by giving the remaining 10 error types their own slopes.
starting_vals_sep<- c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5 , 10^-5, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6, -10^-6 , -10^-6)
separate_slopes_glm<-glm(E ~subtype:R+condition:R + subtype:RlogA+0,start = starting_vals_sep, family = poisson(link = identity), data = ecoli_noCU)
summary(separate_slopes_glm)
anova(ga_slope_glm, separate_slopes_glm, test = "Chisq")

```


```{r}
#testing Synonymous vs. non-synonymous
starting_vals = c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5, -10^-6, -10^-6)
summary(ga_syn_glm<-glm(E ~subtype:R+condition:R+SNS:R + isGA:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
anova(ga_slope_glm, ga_syn_glm, test = "Chisq")

```

```{r}
#testing groEL
starting_vals = c(10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5,10^-5,10^-5, -10^-6, -10^-6)
summary(ga_groel_glm<-glm(E ~subtype:R+condition:R +groEL:R + isGA:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
anova(ga_slope_glm, ga_groel_glm, test = "Chisq")

```

```{r}
#testing gene length and locus position
max_locus <- aggregate(loci~gene, data = ecoli_noCU, FUN = "max")
min_locus <- aggregate(loci~gene, data = ecoli_noCU, FUN = "min")
min_locus$gene_length<-max_locus$loci - min_locus$loci
min_locus$start_locus<-max_locus$loci
min_locus<-subset(min_locus, select = -loci)

ecoli_noCU<-merge(ecoli_noCU, min_locus, by = 'gene')
ecoli_noCU$loci_position_abs<-ecoli_noCU$loci - ecoli_noCU$start_locus
ecoli_noCU$loci_position_rel<-(ecoli_noCU$loci - ecoli_noCU$start_locus)/ecoli_noCU$gene_length
```

```{r}
starting_vals= c(-10^-6, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 10^-5, 0, 0)
summary(length_glm<-glm(E ~subtype:R+ gene_length:R+RlogA+gene_length:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
summary(rel_pos_glm<-glm(E ~subtype:R+ loci_position_rel:R+RlogA+loci_position_rel:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
summary(abs_pos_glm<-glm(E ~subtype:R+ loci_position_abs:R+RlogA+loci_position_abs:RlogA+0,start = starting_vals, family = poisson(link = identity), data = ecoli_noCU))
```


