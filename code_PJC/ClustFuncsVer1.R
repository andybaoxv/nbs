# Functions defined for the COPD subtype clustering analysis.

# define class for the results object of a clustering function.
# The class contains:
#     -dataset: IDs + clustering assignment
#     -model: the model learned from the respective clustering function
#     -fit metrics: internal measures of cluster quality

require(cluster)
require(randomForest)
setOldClass("silhouette")
#Define a simple class to store data about pruning
pruned<-setClass("pruned",representation(threshold="numeric",percPruned="numeric",drops="vector"))
#Define a class to store k-means cluster results
setOldClass("randomForest")  #allow storage of random forest object in new class
setOldClass("randomForest.formula")
kmclus<-setClass("kmclus",representation(assign="numeric",k="numeric",model="matrix",mbmodel="vector",qual="matrix",bic="vector",avgSil="numeric",clusAvgSil="vector",prune="pruned",cdata="data.frame",cpred="numeric",cpredobj="randomForest.formula",cluschar="data.frame",phenoAssoc="list",phenoAssocAdj="list",SNPAssoc="list",SNPAssocAdj="list",valassign="numeric",vpredobj="randomForest.formula",vqual="matrix",vdata="data.frame",valcluschar="data.frame",valphenoAssoc="list",valphenoAssocAdj="list",valSNPAssoc="list",valSNPAssocAdj="list"))


#This function generates a takes a dataframe of cluster results with sid and merges it to a phenotype dataframe to give mean cluster results
# INPUTS
#  clusres - vector of cluster assignments with sid
#  bdata - data with "characteristic" variables with sid
#  clusvar - optional variable to specify name tat you want the clusters to be called
#
# OUTPUTS
#  data frame with mean cluster characteristics
#

#clusvar is not name of the cluster, it is what you want to cluster to be named
cluscharsd<-function(clusres,bdata,idname="sid",clusvar="cluster",N="Y",SD="Y"){
  temp2<-merge(clusres,bdata,by=idname)
  x<-describeBy(temp2[,3:length(names(temp2))],temp2[,clusvar],skew=F,ranges=F)
  y<-cbind(x[[1]]$n[1],x[[1]]$mean,x[[1]]$sd)
   for (a in 2:length(x)){ 
     y=cbind(y,x[[a]]$n[1],x[[a]]$mean,x[[a]]$sd)
   }
   y=as.data.frame(y)
  rownames(y)=rownames(x[[1]])
  for (a in 1:(length(names(y))/3)){
    names(y)[3*a-2]=paste("C",a,":N",sep="")
    names(y)[3*a-1]=paste("C",a,":Mean",sep="")
    names(y)[3*a]=paste("C",a,":SD",sep="")
  }
  if(N=="N") y<-y[,-(grep("N",names(y)))]
  if(SD=="N") y<-y[,-(grep("SD",names(y)))]
  y<-round(y,2)
  y[,grep("N",names(y))]<-round(y[,grep("N",names(y))],0)
  y  
}


#This function computes NMI between two single clustering results
micompclust<-function(clus1,clus2){
    require(entropy)
      tab<-table(clus1,clus2,useNA="no")
      mi.plugin(tab)/sqrt(entropy(table(clus1,useNA="no"))*entropy(table(clus2,useNA="no")))
  }

#This function takes a dataframe of clustering results (each column represents a clustering solution)
#and return a matrix of NMI values
#
#  DEPENDS: micompclust
#

nmimat<-function(clustdf){
  resmat<-matrix(1:(ncol(clustdf)*ncol(clustdf)),ncol=ncol(clustdf))
  for(cl in 1:ncol(clustdf)){
    for(rw in cl:ncol(clustdf)){      
      resmat[rw,cl]<-micompclust(clustdf[,cl],clustdf[,rw])
    }
  }
  rownames(resmat)<-names(clustdf)
  colnames(resmat)<-names(clustdf)
  resmat
}


makeassoc<-function(ccsall,ccssome,cdf,gdata,clusvar,SD="N",k,lambda,printtab=FALSE){   #ccsall = dataframe with all clinical data + cluster results, ccssome = dataframe with clinical characteristics for printout, cdf = dataframe with cluster results and sid only, gdata = dataframe with additively coded genotypes and sids, clusvar = name of cluster variable (vector with all cluster assignments), k = number of clusters, lambda = lambda for cc algorithm
  if(min(table(cdf[,clusvar]))<10 | length(table(cdf[,clusvar]))<k ) {   # This chunk prevents the ANOVA from failing when one cluster size is too small
    aovout<-vector()
    for(a in 1:5){
      aovout[a]<-1
    }
  } else {
    char<-cluscharsd(ccsall[,c("sid",clusvar)],bdata=ccssome,clusvar=clusvar,SD=SD)
    ref<-names(char)[2*which.max(char["FEV1pp_utah",seq(2,ncol(char),2)])]   # get the number of the cluster with best FEV1, this will be reference
    ref<-strsplit(ref,":")[[1]][1]
    ref<-strsplit(ref,"C")[[1]][2]
    if(printtab) print(xtable(char,caption=paste("Results for 4 bins, k=",k,", lambda=",lambda,sep=""),digits=c(0,(rep(c(0,2),ncol(char)/2)))),floating.environment="sidewaystable")
    gd<-merge(ccssome,gdata,by="sid")
    gd<-merge(gd,cdf[,c("sid",clusvar)])
    x<-snpglm(gd,snps=c(names(gdata)[-1]),ref=ref,cname=clusvar,a=k,rowname=TRUE,adjust=NULL)
    if(printtab) t<-snpORPrint(snpres=x,k=k,ref=ref,tables=NULL,file=NULL,print=TRUE)
    aovout<-vector()
    for(a in names(gdata)[-1]){
      aovout[a]<-summary(aov(gd[,a]~as.factor(gd[,clusvar])))[[1]][[5]][1]
    }
  }
  aovout
}

#function renames cluster variable so that clusters are ordered by FEV1, with cluster a having highest FEV1 % predicted

renameClusFEV<-function(mdata,clusvar="cluster"){
  require(reshape)
  x=describeBy(mdata[,"FEV1pp_utah"],mdata[,clusvar])
  vec=vector()
  #get the means for each cluster
  for(a in 1:length(x)){
    vec[a]=x[[a]]$mean
  }
  #associate the name of each cluster with its respective mean FEV1
  names(vec)<-seq(range(mdata[,clusvar])[1],range(mdata[,clusvar])[2])
  #sort, so that the highest FEV1 cluster in first, etc
  y=sort(vec,decreasing=TRUE)
  #save the original cluster values for comparison
  mdata$origclus<-mdata[,clusvar]
  #rename clusters so that cluster with highest mean FEV is cluster 1, 2nd highest is cluster 2, etc
  index=1
  mdata[,clusvar]=NA
  for(a in names(y)){
    mdata[which(mdata$origclus==a),clusvar]<-index
    index=index+1
  }
  mdata
}


############################## OLD ##############

#wrapper for GEE Assoc that allows for testing multiple variables and
#outputting results objects to a list
#to do ordered logistic regression on categorical variables, these need to be factors in the mdata dataset

glmAssocMult<-function(clusres,bdata,idname="sid",outcomes,adjustvar=NULL,adjfactors=NULL,adjfactlevels=NULL,clusvar="cluster",gee=FALSE,ord=NULL){
    mdata<-merge(clusres,bdata,by=idname)
    res<-list()
      #make ordered categorical factors if specified by the ord argument, this results in ordinal logistic regression being done
      if(!is.null(ord)){
            if(!is.vector(ord)) cat("\nord must be a vector. Error in glmAssocMult.\n")
                for(o in ord){
                        mdata[,o]<-ordered(mdata[,o])
                      }
          }
      for (a in outcomes){
            if(is.factor(mdata[,a])) {
                    res[[a]]<-clmAssoc(mdata,a,adjustvar=adjustvar,adjfactors=adjfactors,clusvar=clusvar)
                  } else {
                          res[[a]]<-glmAssoc(mdata,a,adjustvar=adjustvar,adjfactors=adjfactors,adjfactlevels=adjfactlevels,clusvar=clusvar,gee=gee)
                        }
          }
      res
  }

#function takes in a data frame that has cluster results and other variables
#(such as created with mergetoorig), performs GEE to relate the cluster assignments
#to the variables specified in assocvar.
# INPUTS
#    mdata = data frame that has cluster results and other variables(such as made by mergetoorig)
#    outcome = variables name to test for association with cluster as the dependent variable
#    adjustvar = vector of variables names for adjustment variables for GEE
#    adjfactors = vector of variable names that should be analyzed as factors
#    adjfactlevels = vector of vectors, in same order as adjfactors, has the order of levels of the factor - first level in each vector will be the reference
#    clusvar = name of the clustering variable, default is "cluster"
# OUTPUT
#    list of geeglm object for each element on assocvar
# ASSUMES
#    all elements of outcome are appropriate for either the binomial or gaussian family

glmAssoc<-function(mdata,outcome,adjustvar=NULL,adjfactors=NULL,adjfactlevels=NULL,clusvar="cluster",cluslevels=NULL,gee=FALSE){
    #code cluster variable as a factor
    if(!is.null(cluslevels)) {
          mdata[,clusvar]<-factor(mdata[,clusvar],levels=cluslevels)
        } else {
              mdata[,clusvar]<-as.factor(mdata[,clusvar])
            }
      #make adjustment variables factors as necessary
      if(!is.null(adjfactors)){
            if(length(adjfactors>0)){
                    mdata<-adjfactor(mdata,factors=adjfactors,level=adjfactlevels)     #depednds on helper function written elsewhere here
                  }
          }
      #order dataset by clustering variable (necessary for GEE)
      mdata<-mdata[order(mdata[,clusvar]),]
      #generate formula object
      if(length(adjustvar)>0){
            begin<-paste(outcome,clusvar,sep="~")
                avars<-paste(adjustvar,collapse="+")
                model<-formula(paste(begin,avars,sep="+"))
          }
      else model<-formula(paste(outcome,clusvar,sep="~"))
      #change family to binomial if necessary, NOTE - this won;t be approriate for proportion variables
      orange<-range(mdata[,outcome],na.rm=T)
      if(orange[2]-orange[1]==1) {   #use binomial family
            #run glm with option for gee
            if(gee){
                    geeglm(model,data=mdata,id=mdata[,clusvar],family=binomial)
                  } else {
                          glm(model,data=mdata,family=binomial)
                        }
          } else {    #use gaussian family
                if(gee){
                        geeglm(model,data=mdata,id=mdata[,clusvar],family=gaussian)
                      } else {
                              glm(model,data=mdata,family=gaussian)
                            }
              }
  }




#snpOR helper function - takes variables and makes them factors
#
#INPUT
#
#OUTPUT
#

 adjfactor<-function(mdata,factors,level=NULL){
      counter=1
         for (a in factors){
                if(!is.null(level)) {
                         if (!is.null(level[a])) {
                                    mdata[,a]<-factor(mdata[,a],levels=level[[a]])
                                  } else {
                                             mdata[,a]<-as.factor(mdata[,a])
                                           }
                       } else {
                                mdata[,a]<-as.factor(mdata[,a])
                              }
              }
         mdata
    }


#snpOR helper function
#
#INPUT
#
#OUTPUT

snpglm<-function(mdata,snps,ref,cname,a,rowname=TRUE,adjust=NULL){
  reslist<-list()
  for(s in snps){
    counter=1
                                        #Make matrix object to hold results
    res<-matrix(1:(3*(a-1)),nrow=(a-1))
    for(c in range(mdata[,cname])[1]:range(mdata[,cname])[2]){
      if(c != ref){
        temp<-mdata[which(mdata[,cname]==ref | mdata[,cname]==c),]
        temp[,cname]<-as.factor(temp[,cname])   #this makes cluster variables work in logistic regression, reference will be the smaller number
        if(is.null(adjust)){         #two different formulas depending on whether adjustment variable needed
          f<-formula(paste(cname,s,sep="~"))
        } else {
          lhs<-paste(s,paste(adjust,collapse="+"),sep="+")
          f<-formula(paste(cname,lhs,sep="~"))
        }
        glmres<-glm(f,data=temp,family=binomial)
                                        #fill results row
        res[counter,]<-summary(glmres)$coefficients[grep(s,dimnames(summary(glmres)$coefficients)[[1]]),c(1,2,4)]
        counter=counter+1
      }
    }
                                        #exponentiate beta coefficient
    res[,1]<-exp(res[,1])
    res<-as.data.frame(res)
    names(res)<-c("OR","SE","P")
    clusnames<-seq(range(mdata[,cname])[1]:range(mdata[,cname])[2])
    clus<-paste("C",clusnames[which(!clusnames %in% ref)],sep=":")
    if(rowname) rownames(res)<-clus        #this line makes the code work for GOLD, otherwise it presumes cluster names start at 1
    reslist[[s]]<-res
  }
  reslist
}

#as a side effect, it turns the list into a dataframe
#INPUT
#
#OUTPUT
#

snpORPrint<-function(snpres,k,ref,tables=NULL,file=NULL,print=TRUE){
  require(xtable)
  if(length(snpres)>1){
    namevec<-names(snpres)
    res=snpres[[1]]
    rownames(res)<-paste(rep(namevec[1],each=(k-1)),rownames(res),sep=":")
    for(a in 2:length(snpres)){
      rownames(snpres[[a]])<-paste(rep(namevec[a],each=(k-1)),rownames(snpres[[a]]),sep=":")
      res<-rbind(res,snpres[[a]])
    }
  }
  if(print){
                                        #print results
    print(xtable(res,caption=paste("SNP Assocation with Cluster Membership, Reference = Cluster ",ref,sep=""),display=c("s","f","f","g")),hline.after=seq((k-1),nrow(res),by=(k-1)))
  }
  if(!is.null(tables) && !is.null(file)){      #tables and file must have values to write file
    write.table(res,file=paste(tables,file,sep=""),quote=F,sep="\t")
  }
  res
}


#function prints a transposed version of SNP results table
#see above for details
#this is the one used for the paper results

snpORPrint2<-function(snpres,k,ref,tables=NULL,file=NULL,print=TRUE,sideways=FALSE){
  require(xtable)
  if(length(snpres)>1){
            #Make new data container
    resmat<-matrix(1:(2*length(snpres)*(k-1)),ncol=(2*(k-1)))
    dimnames(resmat)[[1]]<-names(snpres)
    res=snpres[[1]]
    dimnames(resmat)[[2]]<-paste(rep(rownames(res),each=2),c("beta/OR","pval"),sep=":")
    for(a in names(snpres)){
      res<-snpres[[a]]
      for(b in 1:(k-1)){
        resmat[a,(b*2-1)]<-res[b,1]
        resmat[a,(b*2)]<-res[b,3]
      }
    }
  }
  if(print){
                                        #print results
    if(sideways){
      print(xtable(resmat,caption=paste("SNP Assocation with Cluster Membership, Reference = Cluster ",ref,sep=""),display=c("s",rep(c("f","g"),(k-1)))),floating.environment="sidewaystable",size="\\tiny")
    } else {
      print(xtable(resmat,caption=paste("SNP Assocation with Cluster Membership, Reference = Cluster ",ref,sep=""),display=c("s",rep(c("f","g"),(k-1)))))
    }
  }
  if(!is.null(tables) && !is.null(file)){      #tables and file must have values to write file
    write.table(resmat,file=paste(tables,file,sep=""),quote=F,sep="\t")
  }
  
  resmat
}
