rm(list=ls())
source("ClustIndices.R")         #Required for external and internal evaluation of clustering results
library(clue)                    #required for evaluation of clustering accuracy






DataSet="3Sources" 	#Data set name (if any)
n=169			#Number of samples in the data set
K=6				#Number of clusters in the data set

dataDir="DataSets/3Sources/Similarity Matrices"    #Path to data folder

#Reading ground truth class labels if present, for evaluation, not required for proposed model
files=read.table(file.path(dataDir,"Labels"),header=TRUE)
n=dim(files)[1]
true=files[,2]
K=length(unique(true))



#########################################         In Similarity Space   #########################################
#3Sources data set has pairwise Similarity matrix based representation of views
Data<-list()
Data[[1]] <- as.matrix(read.table(file.path(dataDir,"bbc")))       #Load View1 similarity matrix
Data[[2]] <- as.matrix(read.table(file.path(dataDir,"guardian")))  #Load View2 similarity matrix
Data[[3]] <- as.matrix(read.table(file.path(dataDir,"reuters")))   #Load View3 similarity matrix


# #########################################         For Feature Space baed representation  #########################################
# Data<-list()
# Data[[1]] <- as.matrix(read.table("Data/bbc"))
# Data[[2]] <- as.matrix(read.table("Data/guardian"))
# Data[[3]] <- as.matrix(read.table("Data/reuters"))





M=length(Data)
mod=c("bbc","guardian","reuters")
nkmeans=10
params=NULL;

#############################################################               Grassmannian-SPD-Optimization                         ##############################################################################

rank=45                          #Specify rank for a fixed user defined Rank
# rank =NULL                        #If rank not specified then algorithm runs for rank K to 50
                                  #In that case, perform optimal rank tuning as discussed in paper from results obtained in Results folder

Algo="GeARS"
source("Grassmann-SPD-Optimize.R")
ManifoldGrassmannMinimize(Data,K=K,true=true,modname=mod,DataSet=DataSet,simFromFile=1, rank=rank)
#For feature space based representation of views set simFromFile=0 in ManifoldGrassmannMinimize function
#This function also takes the ground truth cluster labels in the variable 'true'
#Ground truth annotations are used only for external evaluation of clustering results not in the proposed model during clustering which is unsupervised.

#Open files in Results folder using a spreadsheet software like MS-Excel or Libre-Calc using space ' ' as delimeter
#These spreadsheets contains proposed algorithm results for different values of rank r for rank tunning




if(!is.null(rank))    #In case a particular  rank is specified then joint subspace corresponding to that rank is written to file and k-means on joint subspace is performed
{
    #Perform K-means clustering on joint subspace
    UjointSub=as.matrix(read.table("UjointStar.txt",sep=" ",header=FALSE))
    cat("\n First few rows of Ujoint* subspace:\n")
    print(UjointSub[1:5,1:K])
    cat("\n Subspace Dimension: ",dim(UjointSub)[1]," rows",dim(UjointSub)[2]," columns")
    cat("\nClustering on First k columns")
    UjointSubK=UjointSub[,1:K]
    km=kmeans(UjointSubK,K)$cluster
    df=data.frame(cbind(km))
    write.table(df,quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste0(DataSet,"-ClusterAssignment.txt"))
    cat("\n\nFinal cluster assignments written to file:",paste0(DataSet,"-ClusterAssignment.txt\n\n"))
}
