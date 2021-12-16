#### NMDS for Verde Year 1 ####


# Clear workspace and close all graphics ----------------------------------

rm(list=ls())
graphics.off()


# Set WD ------------------------------------------------------------------


##This sets it to the data file for the analytical workflows class repo
setwd("C:/GIT/I-dont-know-what-to-call-this/data/")



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
#install.packages("cluster")
library("cluster")
#install.packages("gclus")
library(gclus)
#install.packages("ape")
library(ape)



# Bring in data matrices for species and habitat data and prepare them for analyses --------------------------------------------------

##Species matrix
verdsp<-read.csv("verde_yr1.csv")

#head() shows you the columns and first six rows of the imported
head(verdsp)

#observe the dimensions of the imported matrix#
dim(verdsp) 

##subset the dataframe on which to base the ordination. i.e. subset only the spp data. 
##site names will be assigned again later from second matrix. 
##This step is necessary for ordinations. Data must only be data without ID/groupings.

#this code subsets columns 2-188, this is the number of species
spdat<-verdsp[2:118] 

head(spdat) # double checks that the correct data was subset

##LOG(X+1) transform matrix. Will run both to try to obtain lowest stress final solution. 
##For this data it fulfills the requirements to be able to do a log(x+1) transformation which is 
##(lowest number must be 0 or 1 (e.g. cannot be 0.3). This will provide a lower stress solution in the long run#

sptrns<-log(spdat+1)


## Bring in the Habitat Matrix which will provide site names. columns that contains the descriptive/environmental data ##

#Note: the order of the data in the habitat matrix must match the order of the data 
##from the abundance data in order to properly correlate site names and habitat variables#

verdhab<-read.csv("verde_hab_yr1.csv")

head(verdhab)

colnames(verdhab)

dim(verdhab) 


# Compare raw data to transformed data to see which has a lower stress --------

##ordination by NMDS with RAW Data. Run this if you only want to view the solutions related to the untransformed data. 
##This data has a higher stress solutionn when compared to the transformed data#
NMDS1 <- metaMDS(spdat, distance = "bray", k = 2)
#head(NMDS1)
NMDS1$stress

#ordination by NMDS with log(x+1) transformed Data

#Function from Vegan --> metaMDS
#metaMDS(matrix to run NMS, distance metric (here we use Bray-Curtis/Sorensen), k --> seed to start
#Check the transformed data, NMDS2

NMDS2 <- metaMDS(sptrns, distance = "bray", k = 2)

# R will give you stats on the runs

head(NMDS2) 

NMDS2$stress #provides the lowest stress solution

stressplot(NMDS2) #creates a plot of distance as a function of dissimilarity and shows the fit of the line	

##Decide to use NMDS1 or NMDS2. Transformed data has lower column and row variance as well as 
##results in a lower final stress##	



# Environmental correlations overlay (will be inserted on the ordi --------


###Environmental Correlations Overlay##########################	

head(verdhab)

##extract only the variables we want to measure within the habitat matrix
envhab<-verdhab[,c(7:21)] 

##stores data for use in the ordination, not entirely sure how it works
##This is using untransformed data (NMDS1)
joint<-envfit(NMDS1,envhab,permutations=999,strata=NULL,choices=c(1,2),scaling="sites")

##This is alternative code using transformed data, the raw habitat data, and scaling=species
##Not entirely sure what a lot of these unputs do
#joint<-envfit(NMDS2,verdhab,permutations=999,strata=NULL,choices=c(1,2),scaling="species")

joint
colnames(joint)
class(joint)
scores(joint,"vectors")

##Creates a new variable for later use
envpts<-scores(joint,"vectors")

class(joint)		
ordiArrowMul(joint)	

##this provides the x and y points related to the arrows that will represent correlations of the species/samples
envscores<-as.data.frame(scores(joint,display="vectors"))

envscores 

##Alternative code using env variables with R^2 > 0.5
envhab2<-verdhab[,c(7,8,11)]
joint2<-envfit(NMDS1,envhab2,permutations=999,strata=NULL,choices=c(1,2),scaling="sites")




# Create an ordination using NMDS1 (there is also code for NMDS2 for comparisons) -------------------


#plot(NMDS1)
##Tells R how to group the shapes, we want them grouped by microhabitat type and by year, which is how we've ordered the data matrix
aspect<- factor(verdhab$order)

##Tells R how we want to color each point, which is by microhabitat type
color<- factor(verdhab$Hab)

##Tells R what specific colors we want. We want to color code by microhabitat
##In this case, it is aphabetical by habitat type, pool is black, riffle is magenta, and run is blue
co<-c("black", "magenta", "blue")

##shape of the individual points (samples) in order as they appear on the matrices, 
##specifically the habitat matrix but the order should match the species matrix
##We want closed shapes for fall and open shapes for spring, 0,1,2 are the closed versions of 15,16,17
shape<-c(0,1,2,15,16,17)

##This plots our points using the untransformed data
plot(NMDS1$points, col=co[color],asp=1,pch = shape[aspect], cex=1.2,  xlab = "NMDS1", ylab = "NMDS2")

##This plots our points using the transformed data
#plot(NMDS2$points, col=co[color],asp=1,pch = shape[aspect], cex=1.2,  xlab = "NMDS1", ylab = "NMDS2")




# This links all points to centroid (old code, would have to modify)  --------------------------------------

#Connect the points that belong to the same treatment with ordispider and create convex hulls around sites habitat (riffle, pool, run)
#ordispider(NMDS2, groups = verdhab$Hab,  label = TRUE,cex=.8,col=c(17,76,54))



# This creates a polygon around each grouping of interest --------

##We are ordering by microhabitat (Micro column of hab matrix), not sure what display does, the border order is based on microhabitat, so riffles, then pools,
##then runs (based on Micro column of hab matrix)
#ordihull(NMDS1, verdhab$Micro, display= "sites", draw= c("polygon"), 
#         col=NULL, border=c("magenta", "magenta", "magenta", "black", "black", "black", "blue", "blue", "blue" ) ,
#        lty= c(1), lwd=2.5) 

##Old code
#ordihull(NMDS2,groups=verdhab$order,label=FALSE,cex=.6,col=c(17,76,54, 2, 30, 110, 92, 128, 150),draw="polygon",alpha=.05,lty=3)



# This adds the environmental overlay  ------------------------------------

##Add joint plot overlay of env correlations. We created this above#



##For env overlay
plot(joint2, choices = c(1,2), at = c(0,0),axis = FALSE, p.max = 0.05, col ="gray40", add = TRUE,cex=.8)

##For transformed data
#plot(joint2, choices = c(1,2), at = c(0,0),axis = FALSE, p.max = 0.05, col ="gray40", add = TRUE,cex=.5)

##old code
#plot(joint,cex=.5,col="black",asp=1,p=0.05,family="serif")


# Rotate ordination -------------------------------------------------------


#vectorfit(joint)
#Rotate ordination so that axis lines up with env variable
#MDSrotate(...)




# Add a legend to the ordination --------------------------------------------------------------


##Add legend and additional text
##Make sure that colors and shapes are consistent with what's on the ordination
txt <- c("Riffle Fall", "Riffle Spring", "Pool Fall", "Pool Spring", "Run Fall", "Run Spring")
legend('bottomright', txt , pch=c(0,15,1,16,2,17),col=c("magenta","magenta","black","black","blue","blue"),
       cex=1, bty = "y")

text(-1.2,-.9,pos=1,"Stress=0.1495",cex=.9)



# This tests for differences between groups -------------------------------



#Bootstrapping and testing for differences between the groups
#fit <- adonis(sptrns ~ Hab, data=verdhab, permutations=999, method="bray")
#fit
#fit1 <- adonis(spdat ~ Hab, data=verdhab, permutations=999, method="bray")
#fit1



# This checks our assumptions ---------------------------------------------


#Check assumption of homogeneity of multivariate dispersion
#distances_data <- vegdist(sptrns)
#anova(betadisper(distances_data, verdhab$Hab))		

#distances_data1 <- vegdist(spdat)
#anova(betadisper(distances_data1, verdhab$Hab))




# This saves our plot as a high res jpeg file -----------------------------

#Plot ordination so that points are colored and shaped according to the groups of interest


#aspect <- factor(verdhab$hs)
#color <- factor(verdhab$hs)
#co<-c("darkorchid1", "darkmagenta", "deepskyblue", "deepskyblue4", 
#"black", "black", "blue", "blue",
#"magenta", "magenta", "chocolate4", "chocolate3")
#shape<-c(0,15,2,17,1,16)
#shape<-c(18,18,18,18,0,15,2,17,1,16,18,18)

#This will initiate the sequence of saving the plot as a jpeg, make sure to run all necessary code before dev.off
#dev.set()

#jpeg("tribplot1.jpeg", width = 6, height = 6, units = 'in', res = 300)
#jpeg(filename="tribplot3.jpeg",res=600,height=8,width=12,units="in")

#plot(NMDS2$points, col=co[color],asp=1,pch = shape[aspect], cex=1.2,  xlab = "NMDS1", ylab = "NMDS2")

#This will save the image
#dev.off()
