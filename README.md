# ReGraph
ReGraph algorithm for Multi-View Clustering based on Grassmannian and Symmetric Positive Definite Manifold Optimization


The GrassGO algorithm permforms integrative clustering on high-dimensional multimodal data sets. A multimodal data set consists of ``M`` modalities X<sub>1</sub>, ..., X<sub>m</sub>, ..., X<sub>M</sub>. Each modality X<sub>m</sub> represents the observations for same set of ``n`` samples from the ``m``-th data source.

Inorder to execute the R code for the 3Sources data set,  within the ``R`` environment execute:
>source("3Sources-demo.R")

The 3Sources is a multi-view multi-source news article clustering data set. It consists of 3 views, that is, news articles from three different news sources, namely, BBC News, The Guardian, and Reuters. The objective here is to cluster the news atricle considering information from multiple news sources.


Optimal low-rank joint subspace corresponding is written to file : ``UjointStar.txt``   
UjointStar.txt contains a ``(n x r)`` matrix.   
Here ``n`` is the number of samples in the data set and ``r`` is the optimal/required rank of the joint subspace.   

``k``-means clustering can be performed on the rows of UJointStar matrix to get the clusters in the data set.   
The cluster assignments are written to the file ``3Sources-ClusterAssignment.txt`` for the 3Sources data set.  

The file ``Grassmann-SPD-Optimize.R`` contains the ``R`` implementation of the MiMIC algorithm as a function `ManifoldJointMinimize`. 
Details of the fuctions is as follows:

Function Name: `ManifoldGrassmannMinimize`

###### #Usage 
`ManifoldGrassmannMinimize(Data,K=K,true=true,modname=mod,DataSet=DataSet,simFromFile=1, rank=rank)
`


Arguments
``Data``:  A list object containing ``M`` data matrices representing ``M`` different omic data types measured in a set of ``n`` samples.    
For each matrix, the rows represent samples, and the columns represent genomic features.
The matrices in the list can have variable numbe of columns(features), but all must have the same number of *n* rows(samples).

``K``: The number of clusters in the data set.

``rank``: The rank of the individual and joint Laplacian. 
Default value: ``NULL``.
if ``rank=NULL``, the algorithm varies the rank between ``K`` to 50 and selects the optimal rank of the subspace.

``simFromFile``: Boolean value in terms of ``1`` or ``0``   
if `0`, algorithm considers the matrices in the `Data` list as feature-based representation, i.e., as `(n x d_m)` data matrices, and computes the graph Laplacian from data matrices.   
if ``1``, algorithm considers the matrices in the `Data` list as graph-based representation, i.e., as `(n x n)` similarity matrices, and computes the graph Laplacian from similarity matrices.
Default value: ``FALSE``.

`mod`: Array containing names of modalities
Default value: `NULL`.
if `mod=NULL`, the algorithm names the modalities as  `1,2,...M`.




# Example call:

```r
DataSet="3Sources" 	 #Data set name (if any)
n=169			          #Number of samples in the data set
K=6				          #Number of clusters in the data set

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


M=length(Data)
mod=c("bbc","guardian","reuters")
nkmeans=10
params=NULL;   #Any additional parameters if required

#############################################################               Grassmannian-SPD-Optimization                         ##############################################################################

# rank=45                        #Specify rank for a fixed user defined Rank
rank =NULL                       #If rank not specified then algorithm runs for rank K to 50
                                 #In that case, perform optimal rank tuning as discussed in paper from results obtained in Results folder

Algo="ReGraph"
source("Grassmann-SPD-Optimize.R")
ManifoldGrassmannMinimize(Data,K=K,true=true,modname=mod,DataSet=DataSet,simFromFile=1, rank=rank)
```
