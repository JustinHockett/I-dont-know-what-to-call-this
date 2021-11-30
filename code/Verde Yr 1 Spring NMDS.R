#### NMDS for Verde Year 1 SPRING ONLY (2018)####
##CODE NEEDS TO BE CHECKED FOR ERRORS##
##Code copied and variables changed from Year 1 NMDS

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
#install.packages("gclus")
library(gclus)
#install.packages("ape")
library(ape)



# Bring in data matrices for species and habitat data and prepare them for analyses --------------------------------------------------

##Species matrix
verdsp_sp18<-read.csv("verde_yr1_spring.csv")

#head() shows you the columns and first six rows of the imported
head(verdsp_sp18)

#observe the dimensions of the imported matrix#
dim(verdsp_sp18) 

##subset the dataframe on which to base the ordination. i.e. subset only the spp data. 
##site names will be assigned again later from second matrix. 
##This step is necessary for ordinations. Data must only be data without ID/groupings.

##this code subsets columns 2-118, this is the number of species
spdat_sp18<-verdsp_sp18[2:118] 

head(spdat_sp18) # double checks that the correct data was subset

##LOG(X+1) transform matrix. Will run both to try to obtain lowest stress final solution. 
##For this data it fulfills the requirements to be able to do a log(x+1) transformation which is 
##(lowest number must be 0 or 1 (e.g. cannot be 0.3). This will provide a lower stress solution in the long run#

sptrns<-log(spdat_sp18+1)


## Bring in the Habitat Matrix which will provide site names. columns that contains the descriptive/environmental data ##

#Note: the order of the data in the habitat matrix must match the order of the data 
##from the abundance data in order to properly correlate site names and habitat variables#

verdhab_sp18<-read.csv("verde_hab_yr1_spring.csv")

head(verdhab_sp18)

colnames(verdhab_sp18)

dim(verdhab_sp18) 


# Compare raw data to transformed data to see which has a lower stress --------

##ordination by NMDS with RAW Data. Run this if you only want to view the solutions related to the untransformed data. 
##This data has a higher stress solutionn when compared to the transformed data#
NMDS1 <- metaMDS(spdat_sp18, distance = "bray", k = 2)
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
##results in a lower final stress



# Environmental correlations overlay (will be inserted on the ordi --------



head(verdhab_sp18)

##extract only the variables we want to measure within the habitat matrix
envhab<-verdhab_sp18[,c(8:21)] 

##stores data for use in the ordination, not entirely sure how it works
##This is using untransformed data (NMDS1)
joint<-envfit(NMDS1,envhab,permutations=999,strata=NULL,choices=c(1,2),scaling="sites")

##This is alternative code using transformed data, the raw habitat data, and scaling=species
##Not entirely sure what a lot of these unputs do
#joint<-envfit(NMDS2,verdhab_sp18,permutations=999,strata=NULL,choices=c(1,2),scaling="species")


##check significance of environmental variables
joint 


colnames(joint)

class(joint)

scores(joint,"vectors")

##Creates a new variable for later use
envpts<-scores(joint,"vectors")

##Not sure what this tells us	
ordiArrowMul(joint)	

##this provides the x and y points related to the arrows that will represent correlations of the species/samples
envscores<-as.data.frame(scores(joint,display="vectors"))

envscores 

##Alternative code using only SIGNIFICANT env variables (less hectic on the ordination)
envhab2_spring18<-verdhab_sp18[,c(8,11,17,19)]
joint2_spring18<-envfit(NMDS1,envhab2_spring18,permutations=999,strata=NULL,choices=c(1,2),scaling="sites")





# Create an ordination using NMDS1 (there is also code for NMDS2 for comparisons) -------------------


#plot(NMDS1)
##Tells R how to group the shapes, we want them grouped by microhabitat type and by year, which is how we've ordered the data matrix
aspect<- factor(verdhab_sp18$order)

##Tells R how we want to color each point, which is by microhabitat type
color<- factor(verdhab_sp18$Hab)

##Tells R what specific colors we want. We want to color code by microhabitat
##In this case, it is aphabetical by habitat type, pool is black, riffle is magenta, and run is blue
co<-c("black", "magenta", "blue")

##shape of the individual points (samples) in order as they appear on the matrices, 
##specifically the habitat matrix but the order should match the species matrix
##We want closed shapes for fall and open shapes for spring, 0,1,2 are the closed versions of 15,16,17
shape<-c(15,16,17)

##This plots our points using the untransformed data
plot(NMDS1$points, col=co[color],asp=1,pch = shape[aspect], cex=1.2,  xlab = "NMDS1", ylab = "NMDS2")

##This plots our points using the transformed data
#plot(NMDS2$points, col=co[color],asp=1,pch = shape[aspect], cex=1.2,  xlab = "NMDS1", ylab = "NMDS2")




# This links all points to centroid (old code, would have to modify)  --------------------------------------

#Connect the points that belong to the same treatment with ordispider and create convex hulls around sites habitat (riffle, pool, run)
#ordispider(NMDS2, groups = verdhab_sp18$Hab,  label = TRUE,cex=.8,col=c(17,76,54))



# This creates a polygon around each grouping of interest --------

##We are ordering by microhabitat (Micro column of hab matrix), not sure what display does, the border order is based on microhabitat, so riffles, then pools,
##then runs (based on Micro column of hab matrix)
ordihull(NMDS1, verdhab_sp18$Micro, display= "sites", draw= c("polygon"), 
      col=NULL, border=c("magenta", "magenta", "magenta", "black", "black", "black", "blue", "blue", "blue" ) ,
       lty= c(1), lwd=2.5) 

##Old code
#ordihull(NMDS2,groups=verdhab_sp18$order,label=FALSE,cex=.6,col=c(17,76,54, 2, 30, 110, 92, 128, 150),draw="polygon",alpha=.05,lty=3)



# This adds the environmental overlay  ------------------------------------

##Add joint plot overlay of env correlations. We created this above#



##For untransformed data
###cex changes size of text
#plot(joint, choices = c(1,2), at = c(0,0),axis = FALSE, p.max = 0.05, col ="gray40", add = TRUE,cex=.5)

##Significant variables only
plot(joint2_spring18, choices = c(1,2), at = c(0,0),axis = FALSE, p.max = 0.05, col ="gray40", add = TRUE, cex=.8)

##old code
#plot(joint,cex=.5,col="black",asp=1,p=0.05,family="serif")


# Rotate ordination -------------------------------------------------------


#vectorfit(joint)
#Rotate ordination so that axis lines up with env variable
#MDSrotate(...)




# Add a legend to the ordination --------------------------------------------------------------


##Add legend and additional text
##Make sure that colors and shapes are consistent with what's on the ordination
txt <- c("Riffle", "Pool", "Run")
legend('bottomright', txt , pch=c(15,16,17),col=c("magenta","black","blue"),
       cex=1, bty = "y")
text(-1.2,-.7,pos=1,"Stress=15.86",cex=.9)



# This tests for differences between groups -------------------------------



#Bootstrapping and testing for differences between the groups
#fit <- adonis(sptrns ~ Hab, data=verdhab_sp18, permutations=999, method="bray")
#fit
#fit1 <- adonis(spdat_sp18 ~ Hab, data=verdhab_sp18, permutations=999, method="bray")
#fit1



# This checks our assumptions ---------------------------------------------


#Check assumption of homogeneity of multivariate dispersion
#distances_data <- vegdist(sptrns)
#anova(betadisper(distances_data, verdhab_sp18$Hab))		

#distances_data1 <- vegdist(spdat_sp18)
#anova(betadisper(distances_data1, verdhab_sp18$Hab))




# This saves our plot as a high res jpeg file -----------------------------

#Plot ordination so that points are colored and shaped according to the groups of interest


#aspect <- factor(verdhab_sp18$hs)
#color <- factor(verdhab_sp18$hs)
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
