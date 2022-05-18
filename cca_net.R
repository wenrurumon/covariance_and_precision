rm(list=ls())
library(impute)
library(dplyr)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(igraph)
library(qgraph)
#Functions
cca <- function(A,B){
  n = nrow(A);
  p = mrank(A);
  q = mrank(B);
  if (p <= q){
    X = A;
    Y = B;
  }else{
    X = B;
    Y = A;
  }
  R = p_ginv_sq(cov(X),0.5) %*% cov(X,Y) %*% p_ginv_sq(cov(Y),1) %*% cov(Y,X) %*% p_ginv_sq(cov(X),0.5);
  k = mrank_sq(R);
  d = Re(eigen(R)$values);
  rho = d[1:k]^(0.5);
  rho[rho >= 0.9999]=0.9;
  chisq_p = CCA_chisq_test(rho,n,p,q);
  return(c("chisq_p"=chisq_p,"df"=p*q));
}
mrank <- function(X){
  X.svd = svd(X);
  X.rank = sum(X.svd$d>1e-6);
  return(X.rank);
}
mrank_sq <- function(X){
  X.eigen = eigen(X);
  X.rank = sum(Re(X.eigen$values)>1e-6);
  return(X.rank);
}
CCA_chisq_test <- function(rho,n,p,q){
  tstat = -1*n*sum(log(1-rho^2));
  p_value = pchisq(tstat,(p*q),lower.tail=FALSE);
  return(p_value);          
}
p_ginv_sq <- function(X,p){
  X.eigen = eigen(X);
  X.rank = sum(X.eigen$values>1e-8);
  X.value = X.eigen$values[1:X.rank]^(-1*p);
  if (length(X.value)==1){
    D = as.matrix(X.value);
  }else{
    D = diag(X.value);
  }
  rlt = X.eigen$vectors[,1:X.rank] %*% D %*% t(X.eigen$vectors[,1:X.rank]);
  return(rlt);
}
qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
#process data
setwd('/Users/wenrurumon/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/da60239cab44a6cfe3abf90382713324/Message/MessageTemp/e2276c0e4ac5589adf520cf6316eef4f/File')
pheno_data <- read.table("pheno_clean_data.csv",sep = ",",stringsAsFactors = FALSE,header = TRUE,row.names=1)
pheno_class <- read.table("Pheno_class.csv",sep = ",",stringsAsFactors = FALSE,header = TRUE)
Y <- pheno_data %>% as.matrix
Y <- apply(Y,2,function(x){
  x[is.na(x)] <- rnorm(sum(is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T))
  scale(x)
})
#analysis by stage
rlt <- do.call(rbind,lapply(unique(substr(rownames(pheno_data),1,1)),function(stagei){
  sel <- (substr(rownames(pheno_data),1,1))==stagei
  Y <- Y[sel,,drop=F]
  Y.class <- lapply(unique(pheno_class$Class),function(classi){
    scale(Y[,colnames(Y)%in%(pheno_class %>% filter(Class==classi))$paper_ID,drop=F])
  })
  stage <- outer(substr(rownames(pheno_data),1,1),unique(substr(rownames(pheno_data),1,1)),'==')+0
  Y.class2 <- lapply(Y.class,function(y){
    scale(y-predict(lm(y~Y.class[[1]])))
  })[-2]
  rlt <- sapply(Y.class2,function(i){
    sapply(Y.class2,function(j){
      cca(i,j)[1]
    })
  })
  dimnames(rlt) <- list(unique(pheno_class$Class)[-1],unique(pheno_class$Class)[-1])
  data.frame(stage=stagei,melt(rlt))
}))
rlt$pvalue <- p.adjust(rlt$value,method='bonferroni')
ggplot(data=rlt) + geom_tile(aes(x=Var1,y=Var2,fill=(pvalue<0.05)),colour='black') + facet_grid(.~stage) +
  theme(axis.text.x = element_text(angle=90))
#Build Graph
g <- lapply(unique(rlt$stage),function(stagei){
  graph_from_adjacency_matrix(((filter(rlt,stage==stagei) %>% acast(Var1~Var2,value.var='pvalue'))<0.05)+0,mode='undirected')
})
#Graph Attribute
sapply(g,function(gi){
  c(unlist(smallworldIndex(gi)),roburst=mean(centr_degree(gi)$res^2)/mean(centr_degree(gi)$res))
})
