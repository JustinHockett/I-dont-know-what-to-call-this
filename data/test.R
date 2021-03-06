## NMDS for Verde 2 year dataset with preliminary trib data 


---------------------------------------------------------------------------
  # Clear workspace and close all graphics ----------------------------------

rm(list=ls())
graphics.off()






---------------------------------------------------------------------------
  # Set WD ------------------------------------------------------------------


##This sets it to the data file for the analytical workflows class repo
setwd("C:/GIT/I-dont-know-what-to-call-this/data/")





---------------------------------------------------------------------------
  # Download/retrieve necessary packages ------------------------------------



#install.packages("ade4")
library(ade4)
#install.packages("vegan")
library(vegan) 
#install.packages("ggplot2")
library(ggplot2)
#install.packages("raster")
library(raster)
#install.packages("sp")
library(sp)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("gclus")
library(gclus)
#install.packages("ape")
library(ape)




----------------------------------------------------------------------------------------------------------
  # Bring in species data matrix and prepare for analyses --------------------------------------------------


verdsp_yr2<-read.csv("verde_tribs_spdata_f17_s19.csv")

#head() shows you the columns and first six rows of the imported
head(verdsp_yr2)

#observe the dimensions of the imported matrix#
dim(verdsp_yr2) #126 rows by 139 columns#

##subset the dataframe on which to base the ordination. i.e. subset only the spp data. site names will be assigned again later from second matrix. This step is necessary for ordinations. Data must only be data without ID/groupings.##

#this code subsets columns 2-139
spdat_yr2<-verdsp_yr2[,2:139] 

head(spdat_yr2) # double checks that the correct data was subset

##LOG(X+1) transform matrix. Will run both to try to obtain lowest stress final solution. For this data it fulfills the requirements to be able to do a log(x+1) transformation which is (lowest number must be 0 or 1 (e.g. cannot be 0.3).
##This will provide a lower stress solution in the long run#

sptrns<-log(spdat_yr2+1)





----------------------------------------------------------------------
  #  Bring in the Habitat Matrix which will provide site names. --------

##Note: the order of the data in the habitat matrix must match the order of the data from the abundance data in order to properly correlate site names and habitat variables#
verdhab_yr2<-read.csv("verde_tribs_habdata_f17_s19.csv")

head(verdhab_yr2)

colnames(verdhab_yr2)

dim(verdhab_yr2) ##126 rows by 21 columns





-------------------------------------------------------------------------------
  # Compare raw data to transformed data to see which has a lower stress --------

##ordination by NMDS with RAW Data. Run this if you only want to view the solutions related to the untransformed data. This data has a higher stress solutionn when compared to the transformed data#
NMDS1 <- metaMDS(spdat_yr2, distance = "bray", k = 2)
#head(NMDS1)
NMDS1$stress #stress=0.19

##creates a plot of distance as a function of dissimilarity and shows the fit of the line	
stressplot(NMDS1)				

#Function from Vegan --> metaMDS
#metaMDS(matrix to run NMS, distance metric (here we use Bray-Curtis/Sorensen), k --> seed to start

#ordination by NMDS with log(x+1) transformed Data				
NMDS2 <- metaMDS(sptrns, distance = "bray", k = 2)

# R will give you stats on the runs

head(NMDS2) 

##creates a plot of distance as a function of dissimilarity and shows the fit of the line			
NMDS2$stress #stress=0.1586

stressplot(NMDS2) #creates a plot of distance as a function of dissimilarity and shows the fit of the line	

##Decide to use NMDS1 or NMDS2. Transformed data has lower column and row variance as well as 
##results in a lower final stress









---------------------------------------------------------------------------			
  # Environmental correlations overlay (will be inserted on the ordi --------


head(verdhab_yr2)

##extract only the variables we want to measure within the habitat matrix
envhab_yr2<-verdhab_yr2[,c(7:21)] 

joint<-envfit(NMDS1,envhab_yr2,permutations=999,strata=NULL,choices=c(1,2),scaling="sites")

#joint<-envfit(NMDS2,verdhab_yr2,permutations=999,strata=NULL,choices=c(1,2),scaling="species")
joint

colnames(joint)

class(joint)

scores(joint,"vectors")

envpts<-scores(joint,"vectors")

class(joint)		

ordiArrowMul(joint)	


envscores<-as.data.frame(scores(joint,display="vectors"))
#this provides the x and y points related to the arrows that will represent correlations of the species/samples	
envscores 





---------------------------------------------------------------------------
  # Data visualization with NMDS --------------------------------------------

##Tells R how to group the shapes, we want them grouped by microhabitat type and by year, which is how we've ordered the data matrix
###micro.season organizes the data into microhabitats by season and by year for the mainstem (fall riffle yr1, fall pool yr1, fall run yr1, spring riffle yr1, spring pool yr1, spring run yr1, fall riffle yr2...)
###For the tribs it organizes by trib and by microhabitat, since only one season of data is presented (Tangle riffle, tangle pool, west clear riffle, west clear pool, fossil riffle, fossil pool)
aspect<- factor(verdhab_yr2$micro.season)

##Tells R how we want to color each point, which is by microhabitat type
color<- factor(verdhab_yr2$micro.season)

##Tells R what specific colors we want. We want to color code by microhabitat
###In this case, it is aphabetical by habitat type, pool is black, riffle is magenta, and run is blue
###For the added tribs, tangle is red, clear is purple, and fossil is green
co<-c("black", "magenta", "blue","black", "magenta", "blue","black", "magenta", "blue","black", "magenta", "blue",
      "red", "red", "brown","brown","green", "green")

##shape of the individual points (samples) in order as they appear on the matrices, 
###specifically the habitat matrix but the order should match the species matrix
###We want closed shapes for fall and closed shapes for spring, 0,1,2 are the open versions of 15,16,17
####Riffle fall is 0, pool fall is 1, run fall is 2.
####For the tribs, 3 is riffle and 4 is pool. 
shape<-c(0,1,2,15,16,17,0,1,2,15,16,17,3,4,3,4,3,4)

##This plots our points using the untransformed data
plot(NMDS1$points, col=co[color],asp=1,pch = shape[aspect], cex=1.2,  xlab = "NMDS1", ylab = "NMDS2")


#Add legend and additional text


txt <- c("Riffle Fall", "Riffle Spring", "Pool Fall", "Pool Spring", "Run Fall", "Run Spring", "Tangle Riffle", "Tangle Pool",
         "Clear Riffle", "Clear Pool", "Fossil Riffle", "Fossil Pool")


legend('bottomright', txt , pch=c(0,15,1,16,2,17,3,4,3,4,3,4),
       col=c("black","black", "magenta","magenta" ,"blue","blue","red", "red", "brown","brown", "green", "green") ,cex=0.9, bty = "y")


text(-2,-1.5,pos=1,"Stress=.1902",cex=.9)

