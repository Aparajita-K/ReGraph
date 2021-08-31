library(pracma)
LaplacianFromFile<-function(W)
{
     n=dim(W)[1]
     Dg=rowSums(W)
     DH=diag(Dg^(-0.5))
     AN=DH%*%W%*%DH
     I=diag(n)
     L=I+AN
     return(list(L=L,W=W,Dg=Dg))
}

GraphLaplacian<-function(Dmat,rbfsg=1000,mod=1)
{
     show=5
     n=dim(Dmat)[1]
     W=matrix(0,n,n)
     for(i in 1:n)
     {
        for(j in 1:i)
        {
             W[i,j]=W[j,i]=exp(-sum((Dmat[i,]-Dmat[j,])^2)/(2*rbfsg*rbfsg))
        }
     }
     Dg=rowSums(W)
     DH=diag(Dg^(-0.5))
     Wn=DH%*%W%*%DH
     I=diag(n)
     L=I+Wn
     return(list(L=L,W=W))
}

computeFobj<-function(LJoint,LList,UList,UJoint)
{
    M=length(LList)
    fclust=-sum(diag(t(UJoint)%*%(LJoint)%*%UJoint))
    fdag=0
    for(m in 1:M)
        fdag=fdag-sum(diag(t(UList[[m]])%*%((LList[[m]]))%*%UList[[m]]))-sum(diag(UJoint%*%t(UJoint)%*%UList[[m]]%*%t(UList[[m]])))
    fpair=0
    for(i in 1:M)
    {
        for(j in 1:M)
        {
            if(i!=j)
                fpair=fpair-sum(diag(UList[[i]]%*%t(UList[[i]])%*%UList[[j]]%*%t(UList[[j]])))
        }
    }
    f=(fclust+(fdag/M)+fpair/(M*M-1))/2
    return(f)
}

maximize_GrassmannManifold<-function(L,Ut,alpha=0.05)
{
	n=nrow(Ut)
	c=ncol(Ut)
	Qt=L%*%Ut
	Zt=(diag(n)-Ut%*%t(Ut))%*%Qt
	Zt1=Ut+(alpha*Zt)
	sv=svd(Zt1)
	Ut1R=sv$u[,1:c]%*%t(sv$v[,1:c])
	Ut1=Ut1R
	return(list(Ut1=Ut1,gradNorm=Qt))
}

symm<-function(B)
{
    Bsymm=(B+t(B))/2
    return(Bsymm)
}

maximize_SPDManifold<-function(Qt,L,alpha=0.05)
{
    symmQt=symm(Qt)
    Zt=L%*%symmQt%*%L
    Zt1=L+alpha*Zt
    Lt1=Zt1
	return(Lt1)
}

alpha_update<-function(UJoint, LList, kappa=2)
{
    nv=length(LList)
    gm=rep(0,nv)
    pow=1/(1-kappa)
    for(i in 1:nv)
        gm[i]=(0.5*sum(diag(t(UJoint)%*%LList[[i]]%*%UJoint)))^pow
    alpha_Temp=gm/sum(gm)
    return(alpha_Temp)
}

alpha_update_gradient<-function(alpha_t, UJoint, LList, kappa=2, step=0.001)
{
    nv=length(LList)
    alpha_t1=rep(0,nv)
    for(i in 1:nv)
    {
        nGrad=kappa*(alpha_t[i]^(kappa-1))*(0.5*sum(diag(t(UJoint)%*%LList[[i]]%*%UJoint)))
        alpha_t1[i]=alpha_t[i]+step*nGrad
    }
    alpha_t1=alpha_t1/sum(alpha_t1)
    return(alpha_t1)
}

ManifoldGrassmannMinimize<-function(Data,K,true,modname="View1",simFromFile=0,DataSet="3Sources",rank=NULL)
{
    set.seed(NULL)
    options(digits=7)
    clearFolder=1
    writeSimilarity=FALSE
	maxiter=100
  	Beta=0.5
  	dfEps=0.01
  	sigma_decrease=1e-05
    dampF=1
    if(is.null(rank))
    {
        startRk=K
        maxRk=50
    }
    else
        startRk=maxRk=rank
    step=1
    nkmeans=10
    kappa=2
    msg=NULL
    surv=FALSE; survInfo=NULL; foundsample=NULL;  #Survival analysis parameters for multi-omics cancer subtype analysis
	parDir=getwd()
	subDir=paste0(DataSet,"-Results")
	parDir=file.path(parDir, subDir)
	dir.create(parDir, showWarnings = F)
	M=length(Data)
	app=FALSE
	nstart=100
	rkSeq=seq(startRk,maxRk,step)
    if(clearFolder==1)
        file.remove(file.path(parDir,list.files(parDir)))
    fileAvg=file.path(parDir,paste0(msg,"Rank=",startRk,"-",maxRk," Diff-f=",dfEps," dampF=",dampF," Rank-Wise ","Average","."))
    cat(file=fileAvg,append=F)
    cat("\n",msg,file=fileAvg,append=T)
    cat("\nStep_Length= 0.01 for Gr, 0.001 for SPD","StepReduce_Beta=",Beta,"kappa=",kappa,"eps=",dfEps,"Ca parameter=",sigma_decrease,"dampF=",dampF,file=fileAvg,append=T)
    printHeader(file=fileAvg,surv=surv,internal=T)
    fileStd=file.path(parDir,paste0(msg,"Rank=",startRk,"-",maxRk," Diff-f=",dfEps," dampF=",dampF," Rank-Wise ","Std","."))
    cat(file=fileStd,append=F)
    cat("\n",msg,file=fileStd,append=T)
    cat("\nStep_Length= 0.01 for Gr, 0.001 for SPD","StepReduce_Beta=",Beta,"kappa=",kappa,"eps=",dfEps,"Ca parameter=",sigma_decrease,"dampF=",dampF,file=fileStd,append=T)
    printHeader(file=fileStd,surv=surv,internal=T)
	fileL=file.path(parDir,paste0(msg,"Rank=",startRk,"-",maxRk," Diff-f=",dfEps," dampF=",dampF))
	cat(file=fileL,append=F)
    cat("\n",msg,file=fileL,append=T)
    cat("\nStep_Length= 0.01 for Gr, 0.001 for SPD","StepReduce_Beta=",Beta,"kappa=",kappa,"eps=",dfEps,"Ca parameter=",sigma_decrease,"dampF=",dampF,file=fileL,append=T)
	show=5
	ShiLm=list()
	Datas=list()
	Rel=vector(mode="numeric",length=M)
	Alpha=vector(mode="numeric",length=M)
	order=vector(mode="numeric",length=M)
	Dist=vector(mode="numeric",length=M)
	sigma=vector(mode="numeric",length=M)
	n=dim(Data[[1]])[1]
	for(rkFrac in rkSeq)
	{
		cat("\n",DataSet,"At Rank ",rkFrac)
		nsc=rkFrac
        dK=min(K,nsc)
	    for(m in 1:M)
	    {
		    if(simFromFile==1)
		       Lp=LaplacianFromFile(Data[[m]])
		    else
		    {
		        Datas[[m]]=scale(Data[[m]],center=T,scale=F)
		        sgFrac=0.5
		        Dist[m]=max(as.numeric(dist(Datas[[m]])))
		        sigma[m]=sgFrac*Dist[m]
	            Lp=GraphLaplacian(Datas[[m]],rbfsg=sigma[m],mod=m)
	        }
            if(writeSimilarity==TRUE)
            {
                Wfile=file.path(parDir,paste0("Similarity Matrices/",modname[m]))
                write.table(Lp$W,file=Wfile,row.names=FALSE,col.names=FALSE,quote=FALSE)
            }
	        Lm=Lp$L
            Wm=Lp$W
	        evi=eigen(Lm)
	        evi$values[ abs(evi$values)<1e-10 ] <- 0
	        evi$values<-round(evi$values,digits=10)
	        ln=length(evi$values)
	        evind=which(evi$values!=2)[1:nsc]
	        Um=evi$vectors[,evind,drop=FALSE]
	        Dm=evi$values[evind]
            ind=seq(1,K)
            cat("\nEigenvalues",m,evi$values[1:K],sum(evi$values[ind]))
            Rel[m]=sum(evi$values[ind])
		    Lmr=Um%*%diag(Dm)%*%t(Um)
		    ShiLm[[m]]=list(L=Lm,U=Um[,1:K,drop=FALSE],Lr=Lmr)
	    }
	    tstart=proc.time()[3]
	    order=sort(Rel,index.return=TRUE,decreasing=TRUE)$ix
        alphaW=rep(0,M)
    	for(m in 1:M)
        	alphaW[order[m]]=Rel[order[m]]*(1/(dampF^m))
    	Alpha=alphaW/sum(alphaW)
        # Alpha=rep(1,M)/M
        cat("\nRelevance=",Rel)
        cat("\norder=",order)
        cat("\nAlpha^k=",Alpha^kappa)
    	LJoint=matrix(0,n,n)
        for(m in 1:M)
            LJoint=LJoint+(Alpha[m]^kappa)*ShiLm[[m]]$Lr
######################### Gradient Descent Manifold Optimization #################################################
        eigLJoint=eigen(LJoint)
        UJoint=eigLJoint$vectors[,1:nsc,drop=FALSE]                   #To change to 1:K in revision
		UList=list()
		UListTemp=list()
        LList=list()
        LListTemp=list()
		fprev=0
		for(i in 1:M)
        {
			UList[[i]]=ShiLm[[i]]$U
            LList[[i]]=ShiLm[[i]]$Lr
        }
		JointKernel=UJoint%*%t(UJoint)
	    SumKernel=matrix(0,n,n)
	    for(i in 1:M)
	        SumKernel=SumKernel+UList[[i]]%*%t(UList[[i]])
	    cat("\n\n",file=fileL,append=TRUE)
	    cat("\nStep_Length= 0.01 for Gr, 0.001 for SPD","StepReduce_Factor=",Beta,"eps=",dfEps,file=fileL,append=T)
	    printHeader(file=fileL,surv=surv)
		Algo=paste0("UJoint-0-Rank-",nsc,msg)
	    km=kmeans(UJoint[,1:K],K,iter.max=100,nstart=nstart)$cluster
        # km=kmeans_centroid(UJoint[,1:K],K=K,true=true)
        time=(proc.time()[3]-tstart)
		f=fprev=computeFobj(LJoint,LList,UList,UJoint)
    	params=paste("f=",f)
        cat("\niteration 0 : f=",f)
        res=printResult(Dmat=UJoint[,1:dK],clust=km,true=true,M=length(Data),Algo=Algo,survInfo=survInfo,foundsample=foundsample,params=params,file=fileL,parDir=parDir,plot=F,writeData=F,time=time)
#### Manifold Optimization Initialization Ends Here
        prevdf=Inf
	    t=0
        stepGr=0.01
        stepSPD=0.001
	    while(stepGr>1e-03)
	    {
            cat("\niteration",t+1,":")

	      	#Joint Subspace Optimization
            JointOpt=maximize_GrassmannManifold(LJoint+SumKernel/M,UJoint,alpha=stepGr)
            UJointTemp=JointOpt$Ut1
			JointKernel=UJointTemp%*%t(UJointTemp)

            #Individual Subspace Optimization
			for(i in 1:M)
			{
                # UListTemp[[i]]=maximize_GrassmannManifold((LList[[i]]+JointKernel),UList[[i]],alpha=stepGr)$Ut1

                UjNoti=matrix(0,n,n)
                for(j in 1:M)
                {
                    if(j!=i)
                        UjNoti=UjNoti+(UList[[j]]%*%t(UList[[j]]))
                }
                UListTemp[[i]]=maximize_GrassmannManifold((LList[[i]]+JointKernel+UjNoti/(M-1)),UList[[i]],alpha=stepGr)$Ut1
			}

            #Individual Graph Laplacian Optimization
            for(m in 1:M)
            {
                 UjUjt=(Alpha[m]^kappa)*JointKernel
                 UUt=UListTemp[[m]]%*%t(UListTemp[[m]])
                 LListTemp[[m]]=maximize_SPDManifold(UjUjt+UUt,LList[[m]],alpha=stepSPD)
            }

            #Graph Weight Update
            AlphaTemp=alpha_update(UJointTemp, LListTemp, kappa=kappa)
            # AlphaTemp=alpha_update_gradient(Alpha, UJointTemp, LListTemp, kappa=kappa, step=0.001)

            #Joint Graph Laplacian Update
            LJointTemp=matrix(0,n,n)
            for(m in 1:M)
                 LJointTemp=LJointTemp+(AlphaTemp[m]^kappa)*LListTemp[[m]]
            LJointTemp=diag(n)+LJointTemp


            #Updated Objective Computation
            f=computeFobj(LJointTemp,LListTemp,UListTemp,UJointTemp)
			df=fprev-f; Ca=0.0
			Ca=fprev-f-sigma_decrease*stepGr*sum(diag(t(JointOpt$gradNorm)%*%JointOpt$gradNorm))
            cat(" f=",f)
            # readline(prompt="Press [enter] to continue")
			if(df>=dfEps && Ca>=0 && t<maxiter)
			{
				t=t+1
				UJoint=UJointTemp
                Alpha=AlphaTemp
                LJoint=matrix(0,n,n)
				for(m in 1:M)
                {
					UList[[m]]=UListTemp[[m]]
                    LList[[m]]=LListTemp[[m]]
                    LJoint=LJoint+(Alpha[m]^kappa)*LList[[m]]
                }

			  	Algo=paste0("UJoint-",t,"-Rank-",nsc,msg)
                params=paste("f=",f,"Diff-f=",df,"Ca=",Ca)
				km=kmeans(UJoint[,1:K],K,iter.max=100,nstart=nstart)$cluster
				time=(proc.time()[3]-tstart)
				res=printResult(Dmat=UJoint[,1:K],clust=km,true=true,M=length(Data),Algo=Algo,survInfo=survInfo,foundsample=foundsample,params=params,file=fileL,plot=F,parDir=parDir,writeData=F,time=time)
				SumKernel=matrix(0,n,n)
				for(i in 1:M)
					SumKernel=SumKernel+UList[[i]]%*%t(UList[[i]])
				fprev=f
                prevdf=df
			}
			else
			{
				stepGr=Beta*stepGr; stepSPD=Beta*stepSPD
				cat("\nf=",f,"Diff-f=",df,"Ca=",Ca,"Objective not improving, reducing step=",stepGr,file=fileL,append=T)
			}
		}  #End While

        UJointStar=UJoint[,1:dK]
        cat("\n\nk-Means on UJoint*",file=fileL,append=TRUE)
        printHeader(file=fileL,surv=surv,internal=T)
        resList<-list()
        for(test in 1:nkmeans)
        {
            clust=kmeans(UJointStar,K,iter.max=100,nstart=30)$cluster
            res=printResult(Dmat=UJointStar,clust=clust,true=true,M=M,Algo=Algo,survInfo=survInfo,foundsample=foundsample,time=time,file=fileL,parDir=parDir,writeData=F,plot=F,params=params)
            resList[[test]]=res
        }
        if(surv==T)
        {
            rm=apply(as.matrix(sapply(resList, unlist)),1,mean)
            Algo2=paste0(Algo,"-Mean")
            cat(paste0("\n\n",Algo2),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],rm["pvalLR"],rm["pvalGW"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileL,append=TRUE)
            cat(paste0("\n",Algo),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],rm["pvalLR"],rm["pvalGW"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileAvg,append=TRUE)



            rm=apply(as.matrix(sapply(resList, unlist)),1,sd)
            Algo2=paste0(Algo,"-Std")
            cat(paste0("\n",Algo2),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],rm["pvalLR"],rm["pvalGW"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileL,append=TRUE)
            cat(paste0("\n",Algo),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileStd,append=TRUE)
        }
        else
        {
            rm=apply(as.matrix(sapply(resList, unlist)),1,mean)
            Algo2=paste0(Algo,"-Mean")
            cat(paste0("\n\n",Algo2),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileL,append=TRUE)
            cat(paste0("\n",Algo),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileAvg,append=TRUE)



            rm=apply(as.matrix(sapply(resList, unlist)),1,sd)
            Algo2=paste0(Algo,"-Std")
            cat(paste0("\n",Algo2),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileL,append=TRUE)
            cat(paste0("\n",Algo),K,dim(UJointStar)[2],M,rm["ACC"],rm["FM"],rm["ARI"],rm["Purity"],rm["NMI"],rm["Rand"],rm["Jac"],rm["Dice"],
                    rm["Sil"],rm["CH"],rm["DB"],rm["Dunn"],rm["CI"],rm["XB"],rm["PBM"],rm["WSS"],rm["BSS"],rm["RWB"],rm["Time.elapsed"],params,file=fileStd,append=TRUE)
        }

    } #Rank Loop End
    if(!is.null(rank))
    {
        cat("\n\nFor a user provided fixed rank=",rank)
        cat("\nOptimal Subspace Ujoint* corresponding to Joint view written to file: UjointStar.txt")
        write.table(UJointStar,row.names=FALSE,col.names=FALSE,quote=FALSE,file="UjointStar.txt")
    }
    warnings()
}
