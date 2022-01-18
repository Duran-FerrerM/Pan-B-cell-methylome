
########## R script that exemplifies how the predictors 1, 2, 3 and 5 of the classifier algorithm for B cell tumor subtypes are built
########## The methodology is fully described in Martí Duran-Ferrer et al. The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome

##### PACKAGES REQUIRED

library(parallel)
library(e1071)
library(BioQC)

##### FUNCTIONS REQUIRED TO TRAIN THE SVM PREDICTOR

# function that computes the p-value of the one.vs.all differential methylation (based on wilcox.test)

diff.meth=function(betas,g,cores){
  
  # betas = methylation values (row=samples, col=cpgs/variables)
  # g = class of each sample (factor)
  # cores = The number of cores to use (see mc.cores argument in mclapply)
  
  if(class(g)!="factor") stop("g must be a factor")
  if(class(betas)!="matrix") stop("betas must be a matrix")
  if(is.null(rownames(betas))|is.null(colnames(betas))) stop("betas matrix must have row and column names")
  
  output=list()

  for (k in 1:nlevels(g)){
      c1=levels(g)[k]
      splits=split(1:ncol(betas),sort(rep(1:cores,length.out=ncol(betas))))
      
      p=do.call("c",
          mclapply(splits,function(x,betas,g,c1){
            wmwTest(betas[,x], g %in% c1,valType="p.two.sided")
          },betas=betas,g=g,c1=c1,mc.cores=cores)
      )
      
      names(p)=colnames(betas)
      output[[k]]=p
  }
  
  names(output)=levels(g)
  return(output)
}

# function that selects the cpgs with lower p-value

sel.top=function(score,M){
  
  score=sort(score)
  score.topM=score[1:M]
  
  score.draw=score[score==max(score.topM)]
  score.draw=score.draw[!names(score.draw) %in% names(score.topM)]
  
  return(names(c(score.topM,score.draw)))
  
}

# computation of the signatures of each class (training data)

calc.signatures=function(betas,g,num.cpgs,pval.diff){
  
  if(is.null(names(num.cpgs))) stop("num.cpgs must be a named vector with one value for each factor g class")
  if(!all(levels(g) %in% names(num.cpgs))) stop("num.cpgs must have a value for each factor g class")
  
  signatures=matrix(nrow=nrow(betas),ncol=length(num.cpgs))
  colnames(signatures)=levels(g)
  rownames(signatures)=rownames(betas)
  
  sign=list()
  
  for (i in names(num.cpgs)){
    
    if (num.cpgs[i]>0){
      
      # select the cpgs with lower p-values
      sel=sel.top(pval.diff[[i]],M=num.cpgs[i])
      
      # use methylation difference as tie-breaker of the p-values
      if (length(sel)>num.cpgs[i]){
        which.draw=names(which(pval.diff[[i]][sel]==max(pval.diff[[i]][sel])))
        n.to.remove=length(sel)-num.cpgs[i]
        diffmeans=abs(apply(simplify2array(by(betas[,which.draw],g %in% i,colMeans)),1,diff))
        remove=names(sort(diffmeans)[1:n.to.remove])
        sel=sel[!sel %in% remove]
      }
      
      # compute a signature for each class
      aux=betas[,sel,drop=F]
      if (ncol(aux)>1){
        sign[[i]]=ifelse(apply(simplify2array(by(aux,g %in% i,colMeans)),1,diff)>0,0,1)
        signatures[,i]=rowMeans(sweep(aux,MARGIN=2,STATS=sign[[i]],FUN=function(x,y) abs(y-x)))
      } else {
        sign[[i]]=ifelse(diff(tapply(aux[,1],g %in% i,mean))>0,0,1)
        names(sign[[i]])=colnames(aux)
        signatures[,i]=abs(sign[[i]]-aux)
      }
    }
  }
  
  signatures=signatures[,names(num.cpgs)[num.cpgs>0],drop=F]
  
  # retain if the cpg is hypo or hypermethylated in the corresponding class
  attributes(signatures)$sign=sign
  
  return(signatures)
}

# computation of the signatures of each class (new data)

calc.signatures.newdata=function(newdata,train){
  
  if(class(newdata)!="matrix") stop("newdata must be a matrix")
  if(is.null(rownames(newdata))|is.null(colnames(newdata))) stop("newdata matrix must have row and column names")
  
  signatures=matrix(nrow=nrow(newdata),ncol=ncol(train))
  colnames(signatures)=colnames(train)
  rownames(signatures)=rownames(newdata)
  
  sign=attributes(train)$sign
  
  for (i in colnames(train)){
    sel=names(sign[[i]])
    aux=newdata[,sel,drop=F]
    
    if (ncol(aux)>1){
      signatures[,i]=rowMeans(sweep(aux,MARGIN=2,STATS=sign[[i]],FUN=function(x,y) abs(y-x)))
    } else {
      signatures[,i]=abs(sign[[i]]-aux)
    }
  }
  
  return(signatures)
}


##### END FUNCTIONS REQUIRED TO TRAIN THE SVM PREDICTOR


##### TRAIN OF THE SVM PREDICTOR (with example data)

# example data

betas=matrix(runif(30000),ncol=1000) # matrix with the methylation values
colnames(betas)=paste0("cpg",1:ncol(betas))
rownames(betas)=paste0("S",1:nrow(betas))

g=factor(rep(letters[1:3],each=10)) # factor indicating the class of each sample

# parameters required to build the svm predictor (set to an arbitrary value in this example) 
# WARNING: A proper methodology should be used to set or optimize these parameters (for example, using cross-validation). 
# Our crossvalidated strategy will be uploaded soon.

num.cpgs=c(a=2,b=5,c=1) # number of cpgs specific of each class included in the predictor
cost=1 # cost parameter of the svm

# p-value of the differential methylation comparing each class vs all other classes merged

cores=2 # cpu cores used
pval.diff=diff.meth(betas,g=g,cores=cores)

# computation of the signature score for each class

train=calc.signatures(betas,g,num.cpgs,pval.diff)

# train the svm predictor

weights=length(g)/(nlevels(g)*table(g))
predictor=svm(train,g,kernel="linear",cost=cost,probability=T,class.weights=weights)

# apply the svm predictor to new data

newdata=matrix(runif(30000),ncol=1000) # matrix with the methylation values
colnames(newdata)=paste0("cpg",1:ncol(newdata))
rownames(newdata)=paste0("T",1:nrow(newdata))

newdata.signatures=calc.signatures.newdata(newdata,train)

predictions=predict(predictor,newdata=newdata.signatures,probability=T)
predictions


