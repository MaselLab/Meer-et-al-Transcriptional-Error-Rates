---
title: "Saccharomyces transcription errors Notebook"
output:
  html_notebook: default
  pdf_document: default
---

```{r}
yeast_noCU <- read.delim("~/Downloads/final_sacc_binned_table.txt")
yeast_CU<-read.delim("~/Downloads/CU_sacc_binned_6_table.txt")
yeast_CU<-subset(yeast_CU, select = -T)
yeast<-rbind(yeast_CU, yeast_noCU)
```

```{r}
yeast$RlogA<-yeast$R*log(yeast$abundance)
summary(tot_err_glm<-glm(E~type:R+RlogA+ 0, family = poisson(link = identity), data = yeast))
summary(CU_err_glm<-glm(E~type:R+type:RlogA+ 0, family = poisson(link = identity), data = yeast))
anova(tot_err_glm, CU_err_glm, test = "Chisq")

```




```{r}
#all intercepts pooled vs. separate intercepts 
yeast_noCU$RlogA<-yeast_noCU$R*log(yeast_noCU$abundance)
summary(sacc_single_intercept_glm<-glm(E~ R + RlogA+0, family = poisson(link=identity), data = yeast_noCU))
summary(sacc_sep_intercepts_glm<-glm(E~ subtype:R + RlogA+0, family = poisson(link=identity), data = yeast_noCU))
anova(sacc_sep_intercepts_glm, sacc_single_intercept_glm, test = "Chisq")

```

```{r}
#all slopes pooled vs. GA has it's own slope
yeast_noCU$type_GA<-ifelse(yeast_noCU$subtype=='GA', yeast_noCU$RlogA, 0)
summary(sacc_ga_glm<-glm(E~ subtype:R +RlogA +type_GA +0, family = poisson(link=identity), data = yeast_noCU))
anova(sacc_ga_glm, sacc_sep_intercepts_glm, test = "Chisq")

```




```{r}
#all slopes pooled vs. separate slopes for each type
starting_vals<-c(10^-6, 10^-6, 10^-6, 10^-6, 10^-6, 10^-6, 10^-6, 10^-6, 10^-6, 10^-6, 10^-6, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8, -10^-8)
summary(sacc_sep_slopes_glm<-glm(E~ subtype:R +subtype:RlogA +0, family = poisson(link=identity), start = starting_vals, data = yeast_noCU))
anova(sacc_sep_slopes_glm, sacc_sep_intercepts_glm, test = "Chisq")

```



