# Part of the code used in:
# Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
# 
# From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny
# MIT License


##  R file for reanalysis of Knowles et al. 2016 figure S4a: Microbial cells vs. % prophage-like reads


### Dependencies ###

library(MASS)  # requried for robust regression
# Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with
#  S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0

library(resample)  # required for bootstrap sampling
# Tim Hesterberg (2015). resample: Resampling Functions. R package
#  version 0.4. https://CRAN.R-project.org/package=resample



### Import data ###

# This is downloaded from the SI in Knowles et al. 2016 (see "Extended Data Figure 4: Temperateness of viral communities increases with host density and viral functional composition change. (189 KB)" from http://dx.doi.org/10.1038/nature17193)

	dataSI4a <- read.csv("nature17193-sf4.csv", header=TRUE)


#column 1 is the site location
#column 2 is the island chain the site is located at
#column 3 is the microbial cell density per ml
#column 4 is the % provirus like reads


#Only want to look at microbial cell density and provirus like reads. Store these
	DATA = as.data.frame(cbind(dataSI4a$Microbes.per.ml,dataSI4a$X..prophage.reads))

#Rename headers of the dataset
	COLNAMES = c("Cells","PProphage")
	colnames(DATA) = COLNAMES



### Reanalyse confidence intervals for data ###



#Samples
SAMPLES <- 10000 #want to use 10,000 bootstrap samples

#confidence levels
#95%
CONF_LEVEL_95 <- 0.95
PROBS_95 <- (1 + c(-1, 1) * CONF_LEVEL_95) / 2 




## Standard approaches

#1) Pearson
	
	FPearson <- function(data){
		cor.test(log10(data$Cells),data$PProphage, method="pearson")$estimate
	}

	set.seed(1) # for reproductibility
	PearsonBS = bootstrap(DATA, FPearson, SAMPLES) #perform a 10,000 sample bootstrap
	
	Pearson95 = quantile(PearsonBS$replicates, probs = PROBS_95)

	print("Pearson 95% CI:")
	print(Pearson95)

#2) Kendall

	FKendall <- function(data){
		cor.test(log10(data$Cells),data$PProphage, method="kendall")$estimate
	}

	set.seed(1) # for reproductibility
	KendallBS = bootstrap(DATA, FKendall, SAMPLES) #perform a 10,000 sample bootstrap
		
	Kendall95 = quantile(KendallBS$replicates, probs = PROBS_95)

	print("Kendall 95% CI:")
	print(Kendall95)

#3) Spearman
	
	FSpearman <-function(data){
		cor.test(log10(data$Cells),data$PProphage, method="spearman")$estimate
	}

	set.seed(1) # for reproductibility
	SpearmanBS = bootstrap(DATA, FSpearman, SAMPLES) #perform a 10,000 sample bootstrap
		
	Spearman95 = quantile(SpearmanBS$replicates, probs = PROBS_95)

	print("Spearman estimation 95% CI:")
	print(Spearman95)


## Robust regression techniques

#1) Tukey bisquare

	FTbisquare <- function(data){
		A<-rlm(PProphage ~ log10(Cells), data, method="M", maxit=1000,psi="psi.bisquare")
		return( c(A$coefficients, A$converged) )
	}

	set.seed(1) # for reproductibility
	TbisBS = bootstrap(DATA, FTbisquare, SAMPLES) #perform a 10,000 sample bootstrap
	
	USE = TbisBS$replicates[which(TbisBS$replicates[,3]==1),] #ONLY look at converged values
	STbi = dim(USE)

	Tbis95 = quantile(USE[,2], probs = PROBS_95)

	print("Tukey bisquare 95% CI:")
	print(Tbis95)

#2) Huber

	FHuber <- function(data){
		A<-rlm(PProphage ~ log10(Cells), data, method="M", maxit=1000,psi="psi.huber")
		return( c(A$coefficients, A$converged) )
	}

	set.seed(1) # for reproductibility
	HuberBS = bootstrap(DATA, FHuber, SAMPLES) #perform a 10,000 sample bootstrap
	
	USE = HuberBS$replicates[which(HuberBS$replicates[,3]==1),] #ONLY look at converged values
	SHuber = dim(USE)

	Huber95 = quantile(USE[,2], probs = PROBS_95)

	print("Huber estimation 95% CI:")
	print(Huber95)

#3) Hampel
	
	FHampel <- function(data){
		A <- rlm(PProphage ~ log10(Cells), data, method="M", maxit=1000,psi="psi.hampel")
		return( c(A$coefficients, A$converged) )
	}

	set.seed(1) # for reproductibility
	HampelBS = bootstrap(DATA, FHampel, SAMPLES) #perform a 10,000 sample bootstrap
	
	USE = HampelBS$replicates[which(HampelBS$replicates[,3]==1),] #ONLY look at converged values
	SHampel = dim(USE)

	Hampel95 = quantile(USE[,2], probs = PROBS_95)

	print("Hampel estimation 95% CI:")
	print(Hampel95)


