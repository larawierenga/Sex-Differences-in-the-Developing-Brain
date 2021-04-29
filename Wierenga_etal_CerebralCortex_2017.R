# Wierenga et al., Cerebral Cortex 2017
# "A Key Characteristic of Sex Differences
# in the Developing Brain: Greater Variability in Brain"
# Structure of Boys than Girls 
# RCode by Lara Wierenga 

# ---------------------------------------------------
# Load Libraries
# ---------------------------------------------------
# Most libraries you dont need but these are my standard libraries 
library(lattice)
library(lme4)
library(ggplot2)
library(reshape2)
library(plyr)
library(abind)
library(Hmisc)
library(RColorBrewer)
library(mgcv)
library(gtools)
library(scales)
library(grid)
library(methods)
library(igraph)
library(reshape2)
library(devtools)
library(edgebundleR)
library(data.table)
library(Rmpfr)
library(circlize)
library(qvalue)
library(pgirmess)
library(boot)
library(R.matlab)
require(MatchIt)
library(lawstat)
library(car)
library(randomForest)
library(GGally)
library(plotly)
library(randomForest)
library(quantregForest)
library(data.table)
library(reprtree)
library(qqnorm)
library(ggbeeswarm)
library(rogme)
library(cowplot)
library(foreign)

# ---------------------------------------------------
# Load Data File
# ---------------------------------------------------


xx <-  data.frame(read.table("~/Datasets/PINGdata_ROIs_allScansites2.csv",sep=",",header=TRUE))



# ---------------------------------------------------
#  Variance difference (actually log(var_M/var_F) )
# ---------------------------------------------------

## Create data frame with residuals (controlling for Gender, Age and Device) 
vnames=names(xx)[c(2:11,14)]
resid_DF=NULL
for(vname in vnames)
{
	df=data.frame(y=xx[[vname]],gender=as.factor(xx$Gender),age=xx$Age,dev=as.factor(xx$Device))
	ff=randomForest(y~.,data=df,ntrees=500,proximity=TRUE,importance=TRUE, do.trace=100)
	residuals=ff$y -predict(ff)
	resid_DF=cbind(resid_DF,residuals)
}

resid_DF=xx$Age 
names(resid_DF)=vnames
gender=xx$Gender

## Permute 
B=10000
perm=NULL
n=nrow(xx)
for(b in 1:B)
{
	gender_permute=gender[sample(n,n)]
	aa=aggregate(resid_DF,by=list(gr=gender_permute),var)
	ratio=log(aa[aa$gr=="M",-1]/aa[aa$gr=="F",-1])
	perm=rbind(perm,ratio)

}

## p.values
aa=aggregate(resid_DF,by=list(gr=gender),var)
ratio_obs=log(aa[aa$gr=="M",-1]/aa[aa$gr=="F",-1])
  
p_value=NULL
for(k in 1:length(ratio_obs))
	p_value=c(p_value,mean(perm[,k]>=as.numeric(ratio_obs[k])))


# ----------------------------
# Shift plots Rousselet
# ----------------------------

vnames=names(xx)[c(2:11,14)]

resid_DF=NULL
for(vname in vnames)
{
  df=data.frame(y=xx[[vname]],gender=as.factor(xx$Gender),age=xx$Age,dev=as.factor(xx$Device))
  ff=randomForest(y~.,data=df,ntrees=500,proximity=TRUE,importance=TRUE, do.trace=100)
  residuals=ff$y -predict(ff) #(df$y-ff$predicted)
  resid_DF=cbind(resid_DF,residuals)
}

resid_DF=data.frame(resid_DF)
names(resid_DF)=vnames
resid_DF_all_temp <- cbind(xx$Gender,resid_DF)
names(resid_DF_all_temp)[1] <- "Gender"

font_size <- 8

for(vname in vnames){
  resid_DF_all <- resid_DF_all_temp
  resid_DF_all$y <- as.integer(round(resid_DF_all[[vname]],0))
  resid_DF_all <- resid_DF_all[,c(1,ncol(resid_DF_all))]
  
  scaleFUN <- function(x) sprintf("%.0f", x) # te use round labels on x axis
  
  title_name <- names2[which(vnames==vname)]

  ps <- plot_scat2(resid_DF_all,
                   xlabel = vname,
                   ylabel = "Adjusted volume",
                   alpha = .9,
                   shape = 21,
                   size = 2
                   # colour = ,
                   # fill = "grey90"
  ) # scatterplots
  ps <- ps +  
    scale_colour_manual(values=rev(colour_gender)) +  scale_fill_manual(values=rev(colour_gender))+
    scale_y_continuous(labels = scaleFUN) +
    coord_flip() + 
    ggtitle(vname) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.text.x = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.text.y = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.title.y = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.title.x = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.ticks=element_blank(),
          plot.title=element_text (size =16,color="#888888",hjust=.5),
          legend.position="none")
  # print(ps)
  
  # compute shift function
  sf <- shifthd(data=resid_DF_all, formula = y ~ Gender, nboot = 200)
  sf <- round(sf,0)
  sf_rev <- sf*-1
  # plot shift function
  psf <- plot_sf(sf, plot_theme = 2,symb_fill=(colour_gender))

  # change axis labels
  psf <- psf +
    labs(x = "Quantiles of volume",
         y = "Differences M - F")
  

  psf <- add_sf_lab(psf, sf, y_lab_nudge = .1,link_col = c(colour_gender)) + 
    scale_y_continuous(labels = scaleFUN) +
    ggtitle("")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.text.x = element_text (size = font_size,color="#888888"),#,family="Arial"),
          axis.text.y = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.title.y = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.title.x = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.ticks=element_blank(),
          plot.title=element_text (size =16,color="#888888",hjust=.5),
                # legend.text = element_text (size = 12,color="#888888",family="Helvetica"),
          legend.position="none")
  print(psf)
  
  # scatterplot + deciles + colour coded decile differences
  
  p <- plot_scat2(resid_DF_all,
                  xlabel = vname,
                  ylabel = "Adjusted volume",
                  alpha = .9,
                  shape = 21,
                  size = 2
                  #                   colour = "grey10",
                  #                   fill = "grey90"
  )  # scatterplots
  p <- plot_dec_links(p, sf,
                      dec_size = 1,
                      md_size = 1.5,
                      add_rect = TRUE,
                      rect_alpha = 0.2,
                      rect_col = "grey50",
                      link_col = c(colour_gender),
                      add_lab = TRUE) # superimposed deciles + rectangle
  p <- p + coord_flip() +  scale_colour_manual(values=rev(colour_gender)) +  scale_fill_manual(values=rev(colour_gender)) +
    ggtitle(title_name) + 
    scale_y_continuous(labels = scaleFUN) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.text.x = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.text.y = element_blank(),#element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.title.y = element_blank(),#element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.title.x = element_text (size = font_size,color="#888888"),#,family="Arial"), 
          axis.ticks=element_blank(),
          plot.title=element_text (size =16,color="#888888",hjust=.4),
           legend.position="none")# flip axes

  print(cowplot::plot_grid(p, psf, labels=c("A", "B"),ncol = 1, nrow = 2, rel_heights = c(1, 1, 1), label_size = 10, hjust = -0.5, scale=.95))
  
}



# ------------------------------------
#  Gender differences in correlation 
# ------------------------------------

resid_males=scale(xx[xx$Gender=="M",c(2:11,14)])
resid_females=scale(xx[xx$Gender=="F",c(2:11,14)])

zi=rbind(resid_males,resid_females)
gender=c(rep("M",nrow(resid_males)),rep("F",nrow(resid_females)))

cm=cor(zi[gender=="M",])
cf=cor(zi[gender=="F",])

dd=cm-cf
vv=as.numeric(dd)

### Permutations 
ha=NULL
for(b in 1:10000)
{
	GG=sample(gender,size=length(gender))
	cmb=cor(zi[GG=="M",])
	cfb=cor(zi[GG=="F",])
	hh=as.numeric(cmb-cfb)
	ha=c(ha,hh)
}

mm=matrix(ha,ncol=length(vv),byrow=TRUE)
mm=rbind(mm,vv)

# calculate p-values
pvalue_correlation_stat_all <- NULL
for (i in 1:(11*11)){
  temp=mean(mm[,i]>=vv[i])

  pvalue_correlation_stat_all <- cbind(pvalue_correlation_stat_all,temp)
}

matrix_p_values <- matrix(pvalue_correlation_stat_all,nrow=11,ncol=11)

#### Depending on your hypohtesis you can also estiamte P_value females > males

# For illustration you can plot dd with levelplot
levelplot(dd,
          at=col.seq.all,
          scales=list(tck=0, x=list(rot=90,alternating=2),cex=.5,col="#888888"),
          col.regions=(rgb.palette(120)),
          main=paste(title, sep=""),
          ylab=(list(col="#888888")),
          par.settings=list(axis.line=list(col="transparent")),
          shrink = c(0.9, 0.9)
)

