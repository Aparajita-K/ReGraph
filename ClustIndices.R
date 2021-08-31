##************************************************************************** Detailed Survival Analysis of Subtypes ********************************************
detailedSurv<-function(clust,survInfo,foundsample,Algo)
{
		K=length(unique(clust))
        sink(file=paste0("Results/Survival/",Algo,"DetailedSurvival"),append=FALSE)
        group=clust[foundsample]
        fit <- survfit(Surv(times, patient.vital_status) ~ group, data = survInfo)
        sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo)
        pval0 = 1 - pchisq(sd$chisq, length(sd$n) - 1)
        options(width=1000)
        cat("Log-Rank p-value of survival difference : ",pval0,"\n")
	    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo)
	    pvalLR = 1 - pchisq(sd$chisq, length(sd$n) - 1)
	    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,rho=1)
	    pvalGW = 1 - pchisq(sd$chisq, length(sd$n) - 1)
        cat("\nLog Rank=",pvalLR)
        cat("\nWilcoxon=",pvalGW)
        cat("\n")
        print(summary(fit,scale=365,time=c(2,5,7)))
        print(fit,scale=365)
        if(K<=4)
        {
            gr12=which(group==1 | group==2)
		    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,subset=gr12)
		    pval = 1 - pchisq(sd$chisq, length(sd$n) - 1)
		    cat("Subtype1 and Subtype2 p-value : ",pval,"\n")
        }
        if(K<=4 && K>2 && K<=3)
		{
			gr13=which(group==1 | group==3)
		    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,subset=gr13)
		    pval = 1 - pchisq(sd$chisq, length(sd$n) - 1)
		    cat("Subtype1 and Subtype3 p-value : ",pval,"\n")
		    gr23=which(group==2 | group==3)
		    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,subset=gr23)
		    pval = 1 - pchisq(sd$chisq, length(sd$n) - 1)
		    cat("Subtype2 and Subtype3 p-value : ",pval,"\n")
		}
		if(K>3 && K<=4)
		{
			gr14=which(group==1 | group==4)
		    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,subset=gr14)
		    pval = 1 - pchisq(sd$chisq, length(sd$n) - 1)
		    cat("Subtype1 and Subtype4 p-value : ",pval,"\n")
		    gr24=which(group==2 | group==4)
		    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,subset=gr24)
		    pval = 1 - pchisq(sd$chisq, length(sd$n) - 1)
		    cat("Subtype2 and Subtype4 p-value : ",pval,"\n")
		    gr34=which(group==3 | group==4)
		    sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,subset=gr34)
		    pval = 1 - pchisq(sd$chisq, length(sd$n) - 1)
		    cat("Subtype3 and Subtype4 p-value : ",pval,"\n")
		}
        sink()
}

##******************************************************************** Result printing functions *************************************************************************

printHeader<-function(file,surv=F,internal=T)
{
    cat("\nAlgorithm","Clusters","Dimension","Modalities","Accuracy","Fmeasure","ARI","Purity","NMI","Rand","Jaccard","Dice",file=file,append=TRUE)
    if(surv==T)
        cat(" Log-Rank","Wilcoxon",file=file,append=TRUE)
    if(internal==T)
	    cat(" Silhouette->MAX","CH->MAX","DB->MIN","Dunn->MAX","C->MIN","XB->MIN","PBM->MAX","WithinSS->MIN","BetweenSS->MAX",file=file,append=TRUE)
	cat(" TIME","Parameters",file=file,append=TRUE)
}
printResult<-function(Dmat,clust,true,M,Algo,survInfo=NULL,foundsample=NULL,params=NULL,time=0,file,plot=FALSE,parDir=NULL,writeData=F)
{
        colours=c("deeppink4","darkturquoise","darkslateblue","orange")
        lines=c(2,4,1,3)
        K=max(clust)
        ext=confusion(true,clust)
	    npclust=rep(0,K)
        for(k in 1:K)
            npclust[k]=length(which(clust==k))
        pvalLR=pvalGW=-1
        if(!is.null(survInfo))
        {
            group=clust[foundsample]
            sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo)
	        pvalLR = 1 - pchisq(sd$chisq, length(sd$n) - 1)
	        sd=survdiff(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo,rho=1)
	        pvalGW = 1 - pchisq(sd$chisq, length(sd$n) - 1)
	    }
        # params=paste(params),"#no_points_in_Clusters=",paste(npclust,collapse=","))
        ic=list()
        if(!is.null(Dmat))
        {
        	if(writeData==TRUE)
        	{
		    	dir.create(file.path(parDir,"LowRankOutput"), showWarnings = F)
		        dfname=paste0(parDir,"/LowRankOutput/",Algo)
				# df=data.frame(Dmat,clust)
				df=data.frame(Dmat,true)     #####  Write data with ground truth labels in last column for scatter plot
				write.table(df,col.names=FALSE,row.names=FALSE,quote=FALSE,file=dfname)
            }
            ic=InternalIndices(Dmat,clust,indices=c("Silhouette","Calinski-Harabasz","Davies-Bouldin","Dunn","C-Index","Xie-Beni","PBM"))
            if(is.null(survInfo))
                cat(paste0("\n",Algo),K,dim(Dmat)[2],M,ext$ACC,ext$Fmeasure,ext$ARI,ext$Purity,
                    ext$NMI,ext$Rand,ext$Jac,ext$Dice,
                    ic$Sil,ic$CH,ic$DB,ic$Dunn,ic$CI,ic$XB,ic$PBM,ic$WSS,ic$BSS,time,params,file=file,append=TRUE)
            else
                cat(paste0("\n",Algo),K,dim(Dmat)[2],M,ext$ACC,ext$Fmeasure,ext$ARI,ext$Purity,
                    ext$NMI,ext$Rand,ext$Jac,ext$Dice,pvalLR,pvalGW,
                    ic$Sil,ic$CH,ic$DB,ic$Dunn,ic$CI,ic$XB,ic$PBM,ic$WSS,ic$BSS,time,params,file=file,append=TRUE)
        }
        else if(!is.null(survInfo))
        {
            cat(paste0("\n",Algo),K,"-",M,ext$ACC,ext$Fmeasure,ext$ARI,ext$Purity,ext$NMI,
                ext$Rand,ext$Jac,ext$Dice,pvalLR,pvalGW,time,params,file=file,append=TRUE)
        }
		else
        {
            cat(paste0("\n",Algo),K,"-",M,ext$ACC,ext$Fmeasure,ext$ARI,ext$Purity,ext$NMI,
                ext$Rand,ext$Jac,ext$Dice,time,params,file=file,append=TRUE)
        }
        if(plot==TRUE)
        {
            if(max(clust)>1)
	        {
		        subtLabel=NULL
		        for(s in 1:max(clust))
			        subtLabel[s]=paste0("Subtype",s)
			    palette=colours[1:K]
	            linetype=lines[1:K]
			    group=clust[foundsample]
		       	fit <- survfit(Surv(times, patient.vital_status) ~ as.factor(group), data = survInfo)
	            ggsurvplot(fit, title = paste(DataSet,"-",Algo), pval=FALSE, ylab="Survival Probability", xlab="Time (Days)", conf.int = FALSE, font.main = c(28, "bold"),
	                font.x = c(24, "bold"),font.y = c(24, "bold"),font.tickslab = c(20, "bold"), font.legend = c(15, "bold"), ggtheme = theme_bw(),
	                data=survInfo,palette=palette,linetype=linetype,legend.labs= subtLabel,legend=c(0.86,0.86),
	                legend.title ="",pval.coord=c(2500,0.95), pval.size=8)
				dir.create(file.path(parDir,"Survival"), showWarnings = F)
				fname=paste0(parDir,"/Survival/",Algo,".eps")
		        ggsave(fname,device="eps",scale=1,color=TRUE)
	        }
        }
        return(list(ACC=ext$ACC,FM=ext$Fmeasure,ARI=ext$ARI,Purity=ext$Purity,NMI=ext$NMI,Rand=ext$Rand,Jac=ext$Jac,Dice=ext$Dice,
                    pvalLR=pvalLR,pvalGW=pvalGW,Sil=ic$Sil,CH=ic$CH,DB=ic$DB,Dunn=ic$Dunn,CI=ic$CI,XB=ic$XB,PBM=ic$PBM,
                    WSS=ic$WSS,BSS=ic$BSS,Time=time))
}

kmeansMy<-function(Dmat,K,iter.max=100,nstart=30,algorithm="Lloyd")
{
    Cent=as.matrix(read.table("Centers"))
    n=dim(Dmat)[1]
    nstart=dim(Cent)[1]
    centMat=Dmat[Cent[1,],]
    kmbest=kmeans(Dmat,centMat,iter.max=100)
    for(i in 2:nstart)
    {
        centMat=Dmat[Cent[i,],]
        km=kmeans(Dmat,centMat,iter.max=100)
        if(km$tot.withinss<kmbest$tot.withinss)
            kmbest=km
    }
    return(kmbest)
}
BayesInf <- function(fit)
{
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    return(D + 2*m*k)
}
#Internal Indices
sumsq<-function(v1,v2)
{
	eqd=sum((v1-v2)^2)
	return(eqd)
}
proximity<-function(Dmat)
{
	n=dim(Dmat)[1]
	proxim=matrix(0,n,n)
	for(i in 1:n)
	{
		for(j in 1:(i-1))
		{
			proxim[i,j]=sqrt(sumsq(Dmat[i,],Dmat[j,]))
			proxim[j,i]=proxim[i,j]
		}
	}
	return(proxim)
}
InternalIndices<-function(Dmat,clust,indices=c("Silhouette","Calinski-Harabasz","Davies-Bouldin","Dunn","C-Index","Xie-Beni","PBM"))
{

    proxim=proximity(Dmat)
    K=max(clust)
	cl=seq(1:K)
	n=dim(Dmat)[1]
	d=dim(Dmat)[2]
	Xcent=scale(Dmat,center=TRUE,scale=FALSE)
	TSS=sum(Xcent^2)
	sil=dunn=ch=ci=db=xb=pbm=-Inf
	centers=matrix(0,K,d)
	WSS=rep(0,K)
	Ni=rep(0,K)
	for(i in 1:K)
	{
	    Ni[i]=length(which(clust==i))
	    centers[i,]=colMeans(Dmat[which(clust==i),,drop=F])
	    WSS[i]=sum((Dmat[which(clust==i),]-(rep(1,Ni[i])%o%centers[i,]))^2)
	}
	BSS=TSS-sum(WSS)
	if("Silhouette" %in% indices)
	{
	    	s<-vector(mode="numeric",length=n)
	        bvec<-vector(mode="numeric",length=K-1)
	        for(i in 1:n)
	        {
		        own=clust[i]
		        temp=proxim[i,]
		        restown=setdiff(which(clust==own),i)
		        if(length(restown)>0)
		        {
			        owndist=temp[restown]
			        ai=mean(owndist)
		        }
		        else
			        ai=0
		        oc=setdiff(cl,own)
		        for(j in 1:length(oc))
		        {
			        oth=oc[j]
			        otherclust=which(clust==oth)
			        othdist=temp[otherclust]
			        bvec[j]=mean(othdist)
		        }
		        bi=min(bvec)
		        s[i]=(bi-ai)/max(ai,bi)
	        }
	        clustsil=rep(0,K)
	        for(i in 1:K)
	        {
		        temp=which(clust==i)
		        clsil=s[temp]
		        clustsil[i]=mean(clsil)
	        }
            sil=mean(s)
	}
	if("Dunn" %in% indices)
	{
	    maxIntra=-Inf
	    for(i in 1:K)
	    {
	        clusti=which(clust==i)
	        maxi=max(proxim[clusti,clusti])
	        if(maxi>maxIntra)
	            maxIntra=maxi
	    }
	    clustPairs=as.matrix(combn(K,2))
	    minInter=Inf
	    for(id in 1:dim(clustPairs)[2])
	    {
	        i=clustPairs[1,id]
	        j=clustPairs[2,id]
	        clusti=which(clust==i)
	        clustj=which(clust==j)
	        minij=min(proxim[clusti,clustj])
	        if(minij<minInter)
	            minInter=minij
	    }
	    dunn=minInter/maxIntra
	}
	if("Calinski-Harabasz" %in% indices)
	{
	    Xmean=colMeans(Dmat)
	    num=0
	    for(i in 1:K)
	        num=num+Ni[i]*sumsq(centers[i,],Xmean)
	    TWSS=sum(WSS)
	    ch=((n-K)/(K-1))*(num/TWSS)
	}
	if("C-Index" %in% indices)
	{
	    Nw=0
	    S=0
	    for(i in 1:K)
	    {
	        clusti=which(clust==i)
	        Subi=proxim[clusti,clusti]
	        S=S+sum(Subi[lower.tri(Subi)])
	        Nw=Nw+Ni[i]*(Ni[i]-1)/2
	    }
	    allPdist=proxim[lower.tri(proxim)]
	    Smin=sum(sort(allPdist)[1:Nw])
        Smax=sum(sort(allPdist,decreasing=T)[1:Nw])
        ci=(S-Smin)/(Smax-Smin)
	}
	if("Davies-Bouldin" %in% indices)
	{
		Si=rep(0,K)
		for(i in 1:K)
			Si[i]=mean(sqrt(rowSums((Dmat[which(clust==i),,drop=F]-(rep(1,Ni[i])%o%centers[i,]))^2)))
		ProxCent=proximity(centers)
		sum=0
		for(i in 1:K)
		{
		    a1=Si[-i]+Si[i]
		    a2=ProxCent[i,-i]
		    a3=a1/a2
			sum=sum+max(a3)
		}
		db=sum/K
	}
	if("Xie-Beni" %in% indices)
	{
		ProxCent=proximity(centers)^2
		minCentDist=min(ProxCent[lower.tri(ProxCent)])
		xb=sum(WSS)/(n*minCentDist)

		#According To clusterCrit package definition*************************
		clustPairs=combn(K,2)
	    minInter=Inf
	    for(id in 1:dim(clustPairs)[2])
	    {
	        i=clustPairs[1,id]
	        j=clustPairs[2,id]
	        clusti=which(clust==i)
	        clustj=which(clust==j)
	        minij=min(proxim[clusti,clustj])
	        if(minij<minInter)
	            minInter=minij
	    }
	    xb=sum(WSS)/(n*minInter*minInter)
	}
	if("PBM" %in% indices)
	{
		ProxCent=proximity(centers)
		maxCentDist=max(ProxCent[lower.tri(ProxCent)])
		Si=rep(0,K)
		for(i in 1:K)
			Si[i]=sum(sqrt(rowSums((Dmat[which(clust==i),,drop=F]-(rep(1,Ni[i])%o%centers[i,]))^2)))
		Ew=sum(Si)
		Et=sum(sqrt(rowSums(scale(Dmat,center=T,scale=F)^2)))
		pbm=((1/K)*(Et/Ew)*maxCentDist)^2
	}
	return(list(Sil=sil,Dunn=dunn,CH=ch,CI=ci,DB=db,XB=xb,PBM=pbm,WSS=sum(WSS),BSS=BSS))
}
silhouette<-function(Dmat,kmclust)
{
    proxim=proximity(Dmat)
    K=max(kmclust)
	cl=seq(1:K)
	n=dim(Dmat)[1]
	s<-vector(mode="numeric",length=n)
	bvec<-vector(mode="numeric",length=K-1)
	for(i in 1:n)
	{
		own=kmclust[i]
		temp=proxim[i,]
		restown=setdiff(which(kmclust==own),i)
		if(length(restown)>0)
		{
			owndist=temp[restown]
			ai=mean(owndist)
		}
		else
			ai=0
		oc=setdiff(cl,own)
		for(j in 1:length(oc))
		{
			oth=oc[j]
			otherclust=which(kmclust==oth)
			othdist=temp[otherclust]
			bvec[j]=mean(othdist)
		}
		bi=min(bvec)
		s[i]=(bi-ai)/max(ai,bi)
	}
	clustsil=rep(0,K)
	for(i in 1:K)
	{
		temp=which(kmclust==i)
		clsil=s[temp]
		clustsil[i]=mean(clsil)
	}
    sil=mean(s)
	return(sil)
}
CH<-function(Dmat,clust)
{
	    K=max(clust)
		n=dim(Dmat)[1]
		d=dim(Dmat)[2]
		centers=matrix(0,K,d)
		WSS=rep(0,K)
		Ni=rep(0,K)
		for(i in 1:K)
		{
			Ni[i]=length(which(clust==i))
			centers[i,]=colMeans(Dmat[which(clust==i),,drop=F])
			WSS[i]=sum((Dmat[which(clust==i),]-(rep(1,Ni[i])%o%centers[i,]))^2)
		}
	    Xmean=colMeans(Dmat)
	    num=0
	    for(i in 1:K)
	        num=num+Ni[i]*sumsq(centers[i,],Xmean)
	    TWSS=sum(WSS)
	    Xcent=scale(Dmat,center=TRUE,scale=FALSE)
		TSS=sum(Xcent^2)
	    ch=((n-K)/(K-1))*(num/TWSS)
	    ch=(num/TWSS)
	    return(ch)
}
Dunn<-function(Dmat,clust)
{
    proxim=proximity(Dmat)
    K=max(clust)
	cl=seq(1:K)
	n=dim(Dmat)[1]
    maxIntra=-Inf
    for(i in 1:K)
    {
        clusti=which(clust==i)
        maxi=max(proxim[clusti,clusti])
        if(maxi>maxIntra)
            maxIntra=maxi
    }
    clustPairs=as.matrix(combn(K,2))
    minInter=Inf
    for(id in 1:dim(clustPairs)[2])
    {
        i=clustPairs[1,id]
        j=clustPairs[2,id]
        clusti=which(clust==i)
        clustj=which(clust==j)
        minij=min(proxim[clusti,clustj])
        if(minij<minInter)
            minInter=minij
    }
    dunn=minInter/maxIntra
	return(dunn)
}


#External Indices
adjustedRand<-function(true,clust)
{
	n=length(true)
	Kt=max(true)
	Kc=max(clust)
	ContgCT=matrix(0,Kc,Kt)
	for(i in 1:Kc)
	{
		for(j in 1:Kt)
		{
			ci=which(clust==i)
			tj=which(true==j)
			Iij=intersect(ci,tj)
			ContgCT[i,j]=length(Iij)
		}
	}

	t1=0
	for(i in 1:Kc)
	{
		ci=which(clust==i)
		lci=length(ci)
		t1=t1+choose(lci,2)
	}
	t2=0
	for(j in 1:Kt)
	{
		tj=which(true==j)
		ltj=length(tj)
		t2=t2+choose(ltj,2)
	}
	t3=2*t1*t2/(n*(n-1))
	adjR=0

	sm=0
	for(i in 1:Kc)
	{
		for(j in 1:Kt)
		{
			sm=sm+choose(ContgCT[i,j],2)
		}
	}
	num=sm-t3
	den=1/2*(t1+t2)-t3
	adjR=num/den
	return(adjR)
}

randind<-function(true, clust)
{
	n=length(true)
	n11=0
	n10=0
	n01=0
	n00=0
	for(i in 1:(n-1))
	{
		for(j in (i+1):n)
		{
			if((true[i]==true[j])&(clust[i]==clust[j]))
				n11=n11+1
			if((true[i]==true[j])&(clust[i]!=clust[j]))
				n10=n10+1
			if((true[i]!=true[j])&(clust[i]==clust[j]))
				n01=n01+1
			if((true[i]!=true[j])&(clust[i]!=clust[j]))
				n00=n00+1
		}
	}
	Rand=(n11+n00)/(n11+n10+n01+n00)
	ARI=(2*(n00*n11-n01*n10))/((n00+n01)*(n01+n11)+(n00+n10)*(n10+n11))
	return(Rand)
}

confusion<-function(true,clust)
{
	n=length(true)
	Kt=max(true)
	Kc=max(clust)
	n11=0
	n10=0
	n01=0
	n00=0
	for(i in 1:(n-1))
	{
		for(j in (i+1):n)
		{
			if((true[i]==true[j])&(clust[i]==clust[j]))
				n11=n11+1
			if((true[i]==true[j])&(clust[i]!=clust[j]))
				n10=n10+1
			if((true[i]!=true[j])&(clust[i]==clust[j]))
				n01=n01+1
			if((true[i]!=true[j])&(clust[i]!=clust[j]))
				n00=n00+1
		}
	}
	Rand=(n11+n00)/(n11+n10+n01+n00)
	ARI=(2*(n00*n11-n01*n10))/((n00+n01)*(n01+n11)+(n00+n10)*(n10+n11))
	Dice=(2*n11)/(2*n11+n01+n10)
	Jac=n11/(n11+n10+n01)


	ContgCT=matrix(0,Kc,Kt)
	for(i in 1:Kc)
	{
		for(j in 1:Kt)
		{
			ci=which(clust==i)
			tj=which(true==j)
			Iij=intersect(ci,tj)
			ContgCT[i,j]=length(Iij)
		}
	}

	t1=0
	for(i in 1:Kc)
	{
		ci=which(clust==i)
		lci=length(ci)
		t1=t1+choose(lci,2)
	}
	t2=0
	for(j in 1:Kt)
	{
		tj=which(true==j)
		ltj=length(tj)
		t2=t2+choose(ltj,2)
	}
	t3=2*t1*t2/(n*(n-1))
	adjR=0

	sm=0
	for(i in 1:Kc)
	{
		for(j in 1:Kt)
		{
			sm=sm+choose(ContgCT[i,j],2)
		}
	}
	num=sm-t3
	den=1/2*(t1+t2)-t3
	adjR=num/den

	Purity=0
	for(i in 1:Kc)
	{
		maxtj=max(ContgCT[i,])
		Purity=Purity+maxtj
	}
	Purity=Purity/n

    ACC=0
	count=0
	assign=minWeightBipartiteMatching(true,clust)
	for(i in 1:n)
	{
		yi=true[i]
		if(clust[i]==assign[yi])
			count=count+1
	}
	ACC=count/n


	Fmeasure=0
	Fij=matrix(0,Kc,Kt)
	for(i in 1:Kc)
	{
		for(j in 1:Kt)
		{
			Fij[i,j]=(2*ContgCT[i,j])/(sum(ContgCT[i,])+sum(ContgCT[,j]))
		}
	}
	for(j in 1:Kt)
	{
		nj=sum(ContgCT[,j])
		maxFi=max(Fij[,j])
		Fmeasure=Fmeasure+nj*maxFi
	}
	Fmeasure=Fmeasure/n

	InfCT=0
	HC=0
	HT=0
	for(i in 1:Kc)
	{
		for(j in 1:Kt)
		{
		  wl=(n*ContgCT[i,j])/(sum(ContgCT[i,])*sum(ContgCT[,j]))
		  if(wl!=0)
				InfCT=InfCT+(ContgCT[i,j]/n)*log(wl)
		}
	}
	if(InfCT<0)
		InfCT=abs(InfCT)
	for(i in 1:Kc)
	{
		ci=sum(ContgCT[i,])/n
		if(ci!=0)
			HC=HC+ci*log(ci)
	}
	HC=-HC
	for(j in 1:Kt)
	{
		tj=sum(ContgCT[,j])/n
		if(tj!=0)
			HT=HT+tj*log(tj)
	}
	HT=-HT
	NMI=2*InfCT/(HC+HT)
    NMIsqrt=InfCT/(sqrt(HC*HT))
	return(list(ACC=ACC,Fmeasure=Fmeasure,ARI=ARI,Purity=Purity,NMI=NMI, adjR=adjR,Rand=Rand,Dice=Dice,Jac=Jac))
}

# labels from cluster A will be matched on the labels from cluster B
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
    require(clue)
    idsA <- unique(clusteringA)  # distinct cluster ids in a
    idsB <- unique(clusteringB)  # distinct cluster ids in b
    nA <- length(clusteringA)  # number of instances in a
    nB <- length(clusteringB)  # number of instances in b
    if (length(idsA) != length(idsB) || nA != nB) {
        stop("number of cluster or number of instances do not match")
    }

    nC <- length(idsA)
    tupel <- c(1:nA)

    # computing the distance matrix
    assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
    for (i in 1:nC) {
        tupelClusterI <- tupel[clusteringA == i]
        solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
            nA_I <- length(tupelA_I)  # number of elements in cluster I
            tupelB_I <- tupel[clusterIDsB == i]
            nB_I <- length(tupelB_I)
            nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
            return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
        }, clusteringB, tupelClusterI)
        assignmentMatrix[i, ] <- solRowI
    }

    # optimization
    result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
    attr(result, "assignmentMatrix") <- assignmentMatrix
    return(result)
}










# POD<-function(Dmat,clust,dis=0)
# {
# 	K=max(clust)
# 	n=dim(Dmat)[1]
# 	d=dim(Dmat)[2]
# 	clsort=sort(clust,index.return=TRUE)
# 	sort.assign=clsort$x
# 	sort.index=clsort$ix
# 	Dmatorder=Dmat[sort.index,]
# 	B=Dmatorder%*%t(Dmatorder)
# 	Bstand=matrix(0,n,n)
# 	for(i in 1:n)
# 	{
# 		for(j in 1:n)
# 		{
# 			Bstand[i,j]=B[i,j]/sqrt(B[i,i]*B[j,j])
# 			if(Bstand[i,j]<0)
# 				Bstand[i,j]=0
# 		}
# 	}
# 	Z=model.matrix(~0+as.factor(sort.assign))
# 	Bdiag=Z%*%t(Z)
#     Bdiag=Bdiag[sort.index,sort.index]
#     pod=sum(abs(Bstand-Bdiag))/(n*n)
#
# #   s=25
# #   t=160
# #   #flip matrix for plot orientation
# #   f.a=t(as.matrix(rev(as.data.frame(t(Bstand)))))
# #   image(1:(ncol(f.a)+1), 1:(nrow(f.a)+1),t(f.a), axes=FALSE, col=gray(25:0/25),ylab="Samples",xlab="Samples",cex.lab=1.3,cex.axis=1.2,font.lab=2)
# #   if(!is.null(clust)){
# #     axis(side=1, at=c(seq(1,t,s),n), c(1,seq(s,t,s),n),las=2,cex.axis=0.8,cex.lab=2,font=2)
# #     axis(side=2, at=c(seq(1,t,s),n), c(1,seq(s,t,s),n),las=2,cex.axis=0.8,cex.lab=2,font=2)
# #     }else{
# #     axis(side=1, at=1:n, 1:n)
# #     axis(side=2, at=1:n, 1:n)
# #     }
# #   mtext(side=3,line=1,paste(dis,"k=",max(clust),"POD=",round(pod, digits=4),sep=" "),font=2, cex=2)
# #	mtext(side=3,line=1,paste(dis,"POD=",round(pod, digits=4),sep="   "),font=2, cex=2)
# #    box()
#    return(pod)
# }


# #################################################################Clustering Eigenvectors without k-means #####################################################
# discretisationEigenVectorData <- function(eigenVector) {
#   Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
#   maxi <- function(x) {
#     i = which(x == max(x))
#     return(i[1])
#   }
#   j = apply(eigenVector,1,maxi)
#   Y[cbind(1:nrow(eigenVector),j)] = 1
#   return(Y)
# }
#
# discretisationLabels <- function(eigenVectors,ntimes=30) {
# #  normalize <- function(x) x / sqrt(sum(x^2))
# #  eigenVectors = t(apply(eigenVectors,1,normalize))
#   n = nrow(eigenVectors)
#   k = ncol(eigenVectors)
#   seed=sample(n,ntimes)
#   pvemax=0
#   for(t in 1:ntimes)
#   {
#       R = matrix(0,k,k)
#       R[,1] = t(eigenVectors[seed[t],])
#       mini <- function(x) {
#         i = which(x == min(x))
#         return(i[1])
#       }
#       c = matrix(0,n,1)
#       for (j in 2:k) {
#         c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
#         i = mini(c)
#         R[,j] = t(eigenVectors[i,])
#       }
#         eigenDiscrete = discretisationEigenVectorData(eigenVectors %*% R)
#
#
#
#         cluster = apply(eigenDiscrete,1,which.max)
#         K=max(cluster)
#         totss=sum(scale(eigenVectors,center=T,scale=F)^2)
#         withinss=0
#         for(i in 1:K)
#         {
#             clusti=which(cluster==i)
#             withinss=withinss+sum(scale(eigenVectors[clusti,],center=T,scale=F)^2)
#         }
#         betweenss=totss-withinss
#         pve=betweenss/totss
#         if(pve>pvemax)
#         {
#             clustBest=cluster
#             betweenssBest=betweenss
#             totssBest=totss
#             pvemax=pve
#         }
#     }
#   return(list(cluster=clustBest,betweenss=betweenssBest,totss=totssBest))
# }
