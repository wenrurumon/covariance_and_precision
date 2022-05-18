
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
  x[is.na(x)] <- rnorm(sum(is.na(x)),mean(x,na.rm=T),sd(x,na.rm=T)/2)
  scale(x)
})
Y <- scale(Y-predict(lm(Y~cbind(Y[,1:4],outer(substr(rownames(pheno_data),1,1),unique(substr(rownames(pheno_data),1,1)),'=='))-1)))
Y <- data.frame(stage=substr(rownames(pheno_data),1,1),
           id=substr(rownames(pheno_data),2,nchar(rownames(pheno_data))),
           Y) %>% melt(id=1:2)
Y <- merge(Y,pheno_class %>% select(variable=1,class=2)) %>% filter(class!='Demographic data')
system.time(Y2 <- merge(Y,Y,by=c('id','stage')))
Y2.cov <- Y2 %>% group_by(stage,class.x,class.y,variable.x,variable.y) %>% summarise(cov=mean(value.x*value.y)) %>%
  mutate(cov=ifelse(variable.x==variable.y,1,cov))
#by stage
for(i in c('A','B','C')){
  temp <- filter(Y2.cov,stage==i) %>% acast(paste(class.x,variable.x)~paste(class.y,variable.y),value.var='cov')
  pheatmap((abs(solve(temp))>=quantile(abs(solve(temp)),0.75))+0,fontsize = 5,main=i,cluster_rows = F,cluster_cols = F)
}
sapply(c('A','B','C'),function(i){
  temp <- filter(Y2.cov,stage==i) %>% acast(paste(class.x,variable.x)~paste(class.y,variable.y),value.var='cov')
  g <- graph_from_adjacency_matrix((abs(solve(temp))>=quantile(abs(solve(temp)),0.75))+0)
  c(unlist(smallworldIndex(g)),roburst=mean(centr_degree(g)[[1]]^2)/mean(centr_degree(g)[[1]]))
})
#gap
system.time(Y2 <- merge(Y,Y,by=c('id','variable','class')) %>% mutate(gap=value.y-value.x) %>%
  select(id,stage.x,stage.y,class,variable,gap))
system.time(Y2.cov <- merge(Y2,Y2,by=c('id','stage.x','stage.y')) %>% 
              filter(stage.x!=stage.y) %>%
              group_by(stage.x,stage.y,class.x,class.y,variable.x,variable.y) %>% 
              summarise(r=cor(gap.x,gap.y),p=cor.test(gap.x,gap.y)$p.value))
Y2.cov$padjust <- p.adjust(Y2.cov$p,method='fdr')#method='bonferroni')
ggplot() + 
  geom_tile(data=Y2.cov,aes(x=paste(class.x,variable.x),y=paste(class.y,variable.y),
                            fill=ifelse(padjust<0.05,r,0)),
            colour='grey') + 
  facet_grid(stage.x~stage.y) + 
  scale_fill_gradientn(colors=c('blue','white','red'),limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle=90),text=element_text(size=5))
