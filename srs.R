biocLite<-function(pkgs,groupName='lite', ...){
  if (missing(pkgs))
		biocinstall(groupName=groupName, ...)
	else
		biocinstall(pkgs=pkgs,groupName=groupName, ...)
}

#######################################################################
###FUNCION NORMALIZACION
#######################################################################
normalization<-function(file, methods,baseline){
if (!("limma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("limma")
    }
library('limma')
a<-read.table(file,sep='\t',header=T)
m<-colnames(a)
r<-NULL
for(i in 1:length(m)){
	if (baseline==m[i]){
		r<-i
		}
	}
r<-(r-1)	
b<-read.table(file,sep='\t')
b<-as.matrix(b)
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
y<-(a[2:length(a)])
y[is.na(y)]<-0
y[is.nan(y)]<-0
q<-cbind((b[,1][2:length(b[,1])]),y)
colnames(q)<-c('GeneName',colnames(y))
colnames(b)<-NULL
x<-NULL
w<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
			w<-cbind(w,log(as.double(as.vector(b[,j][2:length(b[,j])]))))
		}
	}
}

colnames(x)<-c
e<- b[2:length(b[,1]),]
colnames(e)<-b[1,]
if (methods=='scale'||methods=='quantile'){
	x<-normalizeBetweenArrays(x, method=methods)
	w<-normalizeBetweenArrays(w, method=methods)
	}
if(methods=='z_score'){
	stdv<-apply(x,1,sd,na.rm=TRUE)
	average<-apply(x,1,mean,na.rm=TRUE)
	for (i in 1:nrow(x)){
		for (j in 1:ncol(x)){
			if (is.na(x[i,j])==TRUE){
				x[i,j]<-NA
				}
			if (is.na(x[i,j])==FALSE && is.na(stdv[i])==FALSE){
				x[i,j]<-(x[i,j]-average[i])/stdv[i]
				}
			if (is.na(x[i,j])==FALSE && is.na(stdv[i])==TRUE && is.na(average[i])==FALSE){
				x[i,j]<-NaN
				}
			}
		}
	}
if(methods=='fold_change'){
	for (i in 1:nrow(x)){
		index<-x[i,r]
		for (j in 1:ncol(x)){
			if (is.na(x[i,j])==TRUE){
				x[i,j]<-NA
				}
			if (is.na(x[i,j])==FALSE){
				x[i,j]<-(x[i,j]/index)
				}
			}
		}
	}

x[is.nan(x)==TRUE]<-0
w[is.nan(w)==TRUE]<-0
x[is.na(x)==TRUE]<-0
w[is.na(w)==TRUE]<-0
data<-x
datas<-w
normdata<-cbind(e[,1],datas)
affy_normdata<-cbind(e[,1],data)
two_channel<-cbind(e[,1],data)
colnames(data)<-paste(colnames(data),'.normalized',sep="")
colnames(datas)<-paste(colnames(data))
colnames(normdata)<-c(colnames(e)[1],colnames(data))
colnames(affy_normdata)<-c(colnames(e)[1],colnames(data))
data<-cbind(e,data)
datas<-cbind(q,datas)
write.table(data,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(datas,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
write.table(datas,'c:\\PG\\processes\\data\\illumina_result.txt',sep='\t',row.names=FALSE)
write.table(data,'c:\\PG\\processes\\data\\affymetrix_result.txt',sep='\t',row.names=FALSE)
write.table(normdata,'c:\\PG\\processes\\data\\normalized_result.txt',sep='\t',row.names=FALSE)
write.table(affy_normdata,'c:\\PG\\processes\\data\\affymetrix_normalized_result.txt',sep='\t',row.names=FALSE)
write.table(affy_normdata,'c:\\PG\\processes\\data\\twochannel_normalized_result.txt',sep='\t',row.names=FALSE)
}

#######################################################################
###FUNCION NORMALIZACIÓN VALORES A
#######################################################################
anormalization<-function(file, methods,baseline){

if (!("limma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("limma")
    }

library('limma')
a<-read.table(file,sep='\t',header=T)
m<-colnames(a)
r<-NULL
for(i in 1:length(m)){
	if (baseline==m[i]){
		r<-i
		}
	}
r<-(r-1)	
b<-read.table(file,sep='\t')
b<-as.matrix(b)
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
y<-(a[2:length(a)])
y[is.na(y)]<-0
y[is.nan(y)]<-0
q<-cbind((b[,1][2:length(b[,1])]),y)
colnames(q)<-c('GeneName',colnames(y))
colnames(b)<-NULL
x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}

colnames(x)<-c
e<- b[2:length(b[,1]),]
colnames(e)<-b[1,]
if (methods=='scale'||methods=='quantile'){
	x<-normalizeBetweenArrays(x, method=methods)
	}
if(methods=='z_score'){
	stdv<-apply(x,1,sd,na.rm=TRUE)
	average<-apply(x,1,mean,na.rm=TRUE)
	for (i in 1:nrow(x)){
		for (j in 1:ncol(x)){
			if (is.na(x[i,j])==TRUE){
				x[i,j]<-NA
				}
			if (is.na(x[i,j])==FALSE && is.na(stdv[i])==FALSE){
				x[i,j]<-(x[i,j]-average[i])/stdv[i]
				}
			if (is.na(x[i,j])==FALSE && is.na(stdv[i])==TRUE && is.na(average[i])==FALSE){
				x[i,j]<-NaN
				}
			}
		}
	}
if(methods=='fold_change'){
	for (i in 1:nrow(x)){
		index<-x[i,r]
		for (j in 1:ncol(x)){
			if (is.na(x[i,j])==TRUE){
				x[i,j]<-NA
				}
			if (is.na(x[i,j])==FALSE){
				x[i,j]<-(x[i,j]/index)
				}
			}
		}
	}

x[is.nan(x)==TRUE]<-0
x[is.na(x)==TRUE]<-0
data<-x
colnames(data)<-paste(colnames(data),'.normalized',sep="")
write.table(data,'c:\\PG\\processes\\data\\normalized_a_values.txt',sep='\t',row.names=FALSE)
}


#######################################################################
##FUNCION NORMALIZACION POR ENDOGENOS CONTROL
#######################################################################
ECnormalization<-function(file,gene_list){
a<-read.table(file,sep='\t',header=T)
b<-read.table(file,sep='\t')
f<-a
b<-as.matrix(b)
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
x<-cbind(a[1],x)
d<-c
colnames(x)<-c('gene_id',d)
a<-x
gene_names<-read.table(gene_list,sep='\t',header=FALSE)
gene_names<-as.vector(gene_names$V1)
genes<-NULL
for (i in 1:length(a[,1])){
	for (j in 1:length(gene_names)){
		if(gene_names[j]==as.vector(a[i,1])){
			genes<-rbind(genes,a[i,])
			}
		}
	}
computed_gene<-NULL
for (i in 2:length(genes)){
	computed_gene<-cbind(computed_gene,mean(genes[i]))
	}
colnames(computed_gene)<-colnames(genes)[2:length(colnames(genes))]
normalized_data<-a[2:length(a)]-computed_gene
colnames(normalized_data)<-paste(colnames(a[2:length(a)]),'.normalized',sep="")
data<-cbind(f,normalized_data)
datas<-cbind(a[1],normalized_data)
write.table(data,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(datas,'c:\\PG\\processes\\data\\normalized_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION IMPUTACION
#######################################################################
imputation<-function(file, methods, value){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)

colnames(b)<-NULL

x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
colnames(x)<-c
e<- b[2:length(b[,1]),]
colnames(e)<-b[1,]

for(i in 1:length(x[1,])){
	for(j in 1:length(x[,1])){
		if ((x[j,i])==0 && methods=='Row average'){
			x[j,i]<-mean(x[i,],na.rm=TRUE)
			}
		if ((x[j,i])==0 && methods=='Column average'){
			x[j,i]<-mean(x[,i],na.rm=TRUE)
			}
		if ((x[j,i])==0 && methods=='Constant'){
			x[j,i]<-(value)
			}	
		}
	}
data<-x
z<-length(b[,1])
colnames(data)<-paste(colnames(data),'.imputed',sep="")
data<-cbind(e,data)
write.table(data,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(data,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
}

#FUNCION MATRIX TRANSPOSITION
transpose<-function(file){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]

x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
colnames(x)<-c

data<-x
z<-length(b[,1])
y<-cbind(as.matrix(b[1][2:z,]))
data<-cbind(y,data)
data<-t(data)
rownames(data)[1]<-'GeneID'
write.table(data,'c:\\PG\\processes\\data\\result.txt',sep='\t',col.names=FALSE)

}

#######################################################################
###FUNCION K-MEANS
#######################################################################
kmeanscluster<-function(file, groups){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)

colnames(b)<-NULL

x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
colnames(x)<-c
fits<-kmeans(x,groups)
kcluster<-NULL
for (i in 2:length(b[,1])){
	kcluster<-rbind(kcluster,b[i,])
	}
clusters<- cbind(kcluster, fits$cluster)
colnames(clusters)<-c(b[1,],'group')
write.table(clusters,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION ANALISIS COMPONENTES PRINCIPALES
#######################################################################
pcanalysis<-function(file){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)

colnames(b)<-NULL

x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
colnames(x)<-c
fits<-princomp(x,cor=TRUE)
pca<-NULL
for (i in 2:length(b[,1])){
	pca<-rbind(pca,b[i,])
	}
pca_names<-paste(c,'.score',sep="")
pca_result<- cbind(pca, fits$scores)
colnames(pca_result)<-c(b[1,],pca_names)
write.table(pca_result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION SUMMARIZATION
#######################################################################
summarization<-function(file, methods){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)

colnames(b)<-NULL

x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
colnames(x)<-c
e<- b[2:length(b[,1]),]
colnames(e)<-b[1,]
data<-NULL

if (methods=='Average'){
	for (i in 1:length(x[,1])){
		data<-append(data,mean(as.numeric((x[i,][1:length(x[i,])]))))
		}
	}
if (methods=='Median'){
	for (i in 1:length(x[,1])){
		data<-append(data,median(as.numeric((x[i,][1:length(x[i,])]))))
		}
	}
if (methods=='Standard deviation'){
	for (i in 1:length(x[,1])){
		data<-append(data,sd(as.numeric((x[i,][1:length(x[i,])]))))
		}
	}
if (methods=='Variance'){
	for (i in 1:length(x[,1])){
		data<-append(data,var(as.numeric((x[i,][1:length(x[i,])]))))
		}
	}
data<-as.matrix(data)
z<-length(b[,1])
if (methods=='Average'){
	colnames(data)<-paste(colnames(data),'Average',sep="")
	}
if (methods=='Median'){
	colnames(data)<-paste(colnames(data),'Median',sep="")
	}
if (methods=='Standard deviation'){
	colnames(data)<-paste(colnames(data),'Stdev',sep="")
	}
if (methods=='Variance'){
	colnames(data)<-paste(colnames(data),'Variance',sep="")
	}
data<-cbind(e,data)
write.table(data,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}




#######################################################################
###FUNCION LIMMA
#######################################################################
limma<-function(technology,file,fdr_method){
if (!("limma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("limma")
    }

library('limma')
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
w<-read.table('c:\\PG\\processes\\data\\design.txt',sep='\t')
design<-c(w)
if(technology=='Affymetrix'|| technology=='All-Affymetrix'|| technology=='Custom'|| technology=='Illumina' || technology=='qPCR'){
	design<-as.numeric(design$V1)
	design<-cbind(design,design[length(design):1])
	design <- model.matrix(~ -1+factor(c(w$V1)))
	colnames(design)<-c('control','target')
}
if(technology=='Genepix' || technology=='Agilent' || technology=='Imagene' || technology=='Scanarray' || technology=='Quantarray' || technology=='SMD'){
	design<-as.numeric(design$V1)
	}
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)
colnames(b)<-NULL


x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
means<-rowMeans(x)
colnames(x)<-c
if (technology=='Affymetrix'||technology=='All-Affymetrix'){
	e<- cbind(b[2:length(b[,1]),],means)
	colnames(e)<-c(b[1,],"AveExpr")
	}
if (technology!='Affymetrix' && technology!='All-Affymetrix'){
	e<- cbind(b[2:length(b[,1]),])
	colnames(e)<-c(b[1,])
	}
z<-length(b[,1])
fit<-lmFit(x,design)
fit$genes<-b[,1][2:length(b[,1])]
if (technology=='Affymetrix'||technology=='All-Affymetrix'||technology=='Custom'||technology=='Illumina'||technology=='qPCR'){
	contrast.matrix <- makeContrasts(target-control, levels=design)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	}
if(technology=='Genepix'||technology=='Agilent' || technology=='Imagene' || technology=='Scanarray'|| technology=='Quantarray' || technology=='SMD'){
	fit2 <- fit
	}
fit2 <- eBayes(fit2)
y<-topTable(fit2,adjust=fdr_method,number=length(e[,1]),coef=1)
data<-y
result<-cbind(e,data)
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION SIGNIFICANT ANALYSIS FOR MICROARRAYS
#######################################################################
sam<-function(file, technology, method, statistics){

if (!("samr" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("samr")
    }
if (!("multtest" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("multtest")
    }
library("multtest")
library("samr")
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
w<-read.table('c:\\PG\\processes\\data\\design.txt',sep='\t')
design<-c(w)
if(technology=='Affymetrix'|| technology=='Custom'|| technology=='Illumina'||technology=='Genepix'||technology=='Agilent'){
	design<-as.numeric(design$V1)
	design<-cbind(design,design[length(design):1])
	design <- model.matrix(~ -1+factor(c(w$V1)))
	colnames(design)<-c('control','target')
}

c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)
colnames(b)<-NULL

for (i in 1:length(design)){
	if (design[i]==0){
		design[i]<-2
		}
	}

x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
			}
		}
	}

data=list(x=x,y=design[,1], geneid=b[,1], logged2=TRUE)
samr.obj<-samr(data,  resp.type="Two class unpaired", testStatistic=statistics, nperms=25)
p.val<-samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
vector<-as.vector(p.val)
adjusted<-mt.rawp2adjp(vector,proc=method)
adj.p.val<-adjusted$adjp[,2]
sam_result<-cbind(b[,1][2:length(b[,1])],samr.obj$foldchange,p.val,adj.p.val)
result<-cbind(b[2:length(b[,1]),],sam_result)
colnames(result)<-c(b[1,],'ID','logFC','P.Value','adj.P.Val')
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}




#######################################################################
###FUNCION ANOVA
#######################################################################
anovas<-function(file){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
w<-read.table('c:\\PG\\processes\\data\\design.txt',sep='\t')
design<-c(w)
design<-as.numeric(design$V1)
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)

colnames(b)<-NULL
x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
colnames(x)<-c
e<- b[2:length(b[,1]),]
colnames(e)<-b[1,]
z<-length(b[,1])

q<-c(x[1:length(x[,1]),])
n<-rep(length(x[,1]),length(x[1,]))
group<-rep(design,n)
data<-data.frame(y=q,group=factor(group))
fit<-lm(y~group,data)
data<-anova(fit)
data<-as.data.frame(data)
write.table(data,'c:\\PG\\processes\\data\\result_anova.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION CUSTOM
#######################################################################
custom<-function(file,log_scale,extension){
if (extension=='txt'){
	a<-read.table(file,sep='\t',header=T,check.names=TRUE)
	}
if (extension=='vsc'){
	a<-read.csv(file,sep=',',header=T,check.names=TRUE)
	}
if (log_scale!='TRUE'){
	x<-log(a[2:length(a)])
	x[is.na(x)]<-0
	x[is.nan(x)]<-0
	x[is.infinite(x)]<-0
	}
if (log_scale=='TRUE'){
	x<-a[2:length(a)]
	x[is.na(x)]<-0
	x[is.nan(x)]<-0
	x[is.infinite(x)]<-0
	}
for(i in 1:length(x[1,])){
	for(j in 1:length(x[,1])){
		if (is.na(x[j,i])==TRUE || is.nan(x[j,i])==TRUE){
			x[j,i]<-0
			}
		if ((x[j,i]<=0)==TRUE && log_scale!='TRUE'){
			x[j,i]<-1
			}
		}
	}
result<-x
result<-cbind(a[1],result)
write.table(result,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
if (extension=='vsc'){
	write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
	}
}



#######################################################################
###FUNCION PROJECT
#######################################################################
projectcustom<-function(file,log_scale,extension,technology){
if (extension=='txt'){
	a<-read.table(file,sep='\t',header=T,check.names=TRUE)
	}
if (extension=='vsc'){
	a<-read.csv(file,sep=',',header=T,check.names=TRUE)
	}
if (log_scale!='TRUE'){
	x<-log(a[2:length(a)])
	x[is.na(x)]<-0
	x[is.nan(x)]<-0
	x[is.infinite(x)]<-0
	}
if (log_scale=='TRUE'){
	x<-a[2:length(a)]
	x[is.na(x)]<-0
	x[is.nan(x)]<-0
	x[is.infinite(x)]<-0
	}
for(i in 1:length(x[1,])){
	for(j in 1:length(x[,1])){
		if (is.na(x[j,i])==TRUE || is.nan(x[j,i])==TRUE){
			x[j,i]<-0
			}
		}
	}
n=0
m=0
for (i in 1:length(x[1,])){
	if (colnames(x)[i]=='ID'){
		n=i
	}
}
if (technology=='Affymetrix'||technology=='Agilent'|| technology=='Genepix'||technology=='Quantarray'||technology=='Scanarray'||technology=='Imagene'||technology=='SMD'){
	x<-2^x
	}
if (technology=='Custom'){
	x<-(2.73)^x
	}
if (technology=='Illumina'){
	x<-log(x)
	}
result<-x
result<-cbind(a[1],result)
if (n>0){
	for (i in 1:length(x[1,])){
		if (colnames(x)[i]=='AveExpr'){
			m=i
			}
		}
	if (m>0){
		if (technology!='Illumina'){
			results<-result[1:n-1]
			}
		if (technology=='Illumina'){
			results<-result[1:n]
			}
	}
	if (m==0){
		results<-result[1:n]
	}
}

if (n==0){
	n=length(x[1,])
	results<-result[1:n+1]
	results<-cbind(a[1],results)
}
write.table(result,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
write.table(results,'c:\\PG\\processes\\data\\project_data.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION CONTENIDO PROYECTO
#######################################################################
project<-function(file,file_name){
a<-read.table(file,sep='\t',header=T)
if (file_name!=''){
	hyb_name<-paste(file_name,'_hyb.txt',sep='')
	hyb_file<-paste('c:\\PG\\processes\\data\\projects\\',hyb_name,sep='')
	hyb_content<-read.table(hyb_file,sep='\t',header=T)
	write.table(hyb_content,'c:\\PG\\processes\\data\\hyb_result.txt',sep='\t',row.names=FALSE)
	}
b<-grep('normalized',colnames(a))
q<-NULL
p<-NULL
p<-a[1]
q<-a[1]
#RAW DATA
if (length(b)>0){
	for (i in 3:b[1]-1){
		p<-cbind(p,a[i])
		}
#NORMALIZED
	for (i in 1:length(b)){
		q<-cbind(q,a[b[i]])
		}
	}
if (length(b)==0){
	b<-grep('AveExpr',colnames(a))
	#RAW DATA
	for (i in 3:b[1]-1){
		p<-cbind(p,a[i])
		}
	if (length(b)==0){
		b<-grep('ID',colnames(a))
		#RAW DATA
		for (i in 3:b[1]-1){
			p<-cbind(p,a[i])
			}
		}
	}

write.table(q,'c:\\PG\\processes\\data\\affymetrix_normalized_result.txt',sep='\t',row.names=FALSE)
write.table(p,'c:\\PG\\processes\\data\\box_result.txt',sep='\t',row.names=FALSE)
write.table(q,'c:\\PG\\processes\\data\\illumina_result.txt',sep='\t',row.names=FALSE)
write.table(q,'c:\\PG\\processes\\data\\normalized_result.txt',sep='\t',row.names=FALSE)
write.table(p,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(q,'c:\\PG\\processes\\data\\twochannel_normalized_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION ILLUMINA
#######################################################################
illumina<-function(file,skip_lines){
if (!("beadarray" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("beadarray")
    }

library("beadarray")
a<-readBeadSummaryData(dataFile=file,skip=skip_lines,columns = list(exprs = "AVG_Signal"))
b<-exprs(a)
d<-log(exprs(a))
c<-cbind(rownames(exprs(a)),b)
e<-cbind(rownames(exprs(a)),d)
colnames(c)<-c('GeneID',colnames(c)[2:length(c[1,])])
colnames(e)<-c('GeneID',colnames(e)[2:length(e[1,])])
result<-c
results<-e
write.table(results,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(results,'c:\\PG\\processes\\data\\illumina_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION ANALISIS AFFYMETRIX
#######################################################################
affymetrix<-function(file,analysis){
if (!("affy" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite<-function(pkgs,groupName='lite', ...){
			if (missing(pkgs))
				biocinstall(groupName=groupName, ...)
			else
				biocinstall(pkgs=pkgs,groupName=groupName, ...)
			}	
		  biocLite("affy")
    }

library("affy")
design_data<-read.table(file,sep="\t",header=FALSE)
samples<-as.character(design_data[,1])
files<-as.character(design_data[,2])
#Selección del método para tratar los chips (el resultado en guardado en la variable exp_data)
if (analysis=="RMA"){
	exp_data<-ReadAffy(sampleNames=c(samples),filenames=c(files))
	eset<-rma(exp_data)
	gene_id<-geneNames(exp_data)
	}
if (analysis=="MAS5"){
	exp_data<-ReadAffy(sampleNames=c(samples),filenames=c(files))
	eset<-mas5(exp_data)
	gene_id<-geneNames(exp_data)
	}

result<-as.data.frame(exprs(eset))
if (analysis=="MAS5"){
	result=log(result)
	}
result<-cbind(gene_id,result)
colnames(result)<-c('GeneID',samples)
controles<-grep("AFFX.+",geneNames(exp_data))
control_id<-geneNames(exp_data)[controles]
a<-controles[1]
b<-controles[length(controles)]
d<-pm(exp_data)[a:b,]
e<-cbind(control_id,d)

controles_Hyb<-grep("AFFX-Bio.+",geneNames(exp_data))
hyb_id<-geneNames(exp_data)[controles_Hyb]
f<-controles_Hyb[1]
g<-controles_Hyb[length(controles_Hyb)]
h<-pm(exp_data)[f:g,]
i<-cbind(hyb_id,h)

write.table(result,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\affymetrix_result.txt',sep='\t',row.names=FALSE)
write.table(e,'c:\\PG\\processes\\data\\control_result.txt',sep='\t',row.names=FALSE)
write.table(i,'c:\\PG\\processes\\data\\hyb_result.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION CONTROL AFFYMETRIX AL USAR CDF
#######################################################################

affycontrol<-function(){
exp_data<-read.table('c:\\PG\\processes\\data\\result.txt',sep='\t',header=TRUE)
controles<-grep("AFFX.+",exp_data[,1])
control_id<-exp_data[,1][controles]
a<-controles[1]
b<-controles[length(controles)]
d<-exp_data[a:b,]
write.table(d,'c:\\PG\\processes\\data\\hyb_result.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION TECNOLOGIAS DE DOBLE CANAL
#######################################################################
twochannel<-function(file, technology, bg_method,normalization_method,complement){
if (!("limma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("limma")
    }
library("limma")
targets <- readTargets(file)
if(technology=='scanarray'){
	RG <- read.maimages(targets$FileName, source=technology,annotation=c("ID","Name"))
	}
if(technology!='scanarray'){
	RG <- read.maimages(targets$FileName, source=technology)
	}
RGb <- backgroundCorrect(RG, method=bg_method)
MA <- normalizeWithinArrays(RGb,method=normalization_method)

if(technology=='genepix'){
	result<-cbind(MA$genes$ID,MA$M)
	result_2<-cbind(MA$genes$ID,MA$A)
	colnames(result)<-c('GeneID',targets$Samples)
	colnames(result_2)<-c('GeneID',targets$Samples)
	}
if(technology=='agilent'){
	if (complement!='CGH'){
		result<-cbind(MA$genes$ProbeName,MA$M)
		result_2<-cbind(MA$genes$ProbeName,MA$A)
		colnames(result)<-c('GeneID',targets$Samples)
		colnames(result_2)<-c('GeneID',targets$Samples)
		}
	if (complement=='CGH'){
		result<-cbind(MA$genes$SystematicName,MA$M)
		result_2<-cbind(MA$genes$SystematicName,MA$A)
		colnames(result)<-c('GeneID',targets$Samples)
		colnames(result_2)<-c('GeneID',targets$Samples)
		}
	}
if(technology=='quantarray'){
	result<-cbind(MA$genes$Name,MA$M)
	result_2<-cbind(MA$genes$Name,MA$A)
	colnames(result)<-c('GeneID',targets$Samples)
	colnames(result_2)<-c('GeneID',targets$Samples)
	}
if(technology=='scanarray'){
	result<-cbind(MA$genes$Name, MA$M)
	result_2<-cbind(MA$genes$Name,MA$A)
	colnames(result)<-c('GeneID',targets$Samples)
	colnames(result_2)<-c('GeneID',targets$Samples)
	}
if (technology=='smd'){
	result<-cbind(MA$genes$"Gene Symbol",MA$M)
	result_2<-cbind(MA$genes$"Gene Symbol",MA$A)
	colnames(result)<-c('GeneID',targets$Samples)
	colnames(result_2)<-c('GeneID',targets$Samples)
	}
if (technology=='imagene'){
	result<-cbind(as.vector(MA$genes$"Gene ID"),MA$M)
	result_2<-cbind(as.vector(MA$genes$"Gene ID"),MA$A)
	colnames(result)<-c('GeneID',unique(targets$Samples))
	colnames(result_2)<-c('GeneID',unique(targets$Samples))
	}

write.table(result,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(result_2,'c:\\PG\\processes\\data\\box_result.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION ONE CHANNEL
#######################################################################
singlechannel<-function(file,technology,bg_method){

if (!("limma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("limma")
        }
library("limma")

onechannel<-function (files = NULL, source = "generic", path = NULL, ext = NULL, 
    names = NULL, columns = NULL, other.columns = NULL, annotation = NULL, 
    green.only = FALSE, wt.fun = NULL, verbose = TRUE, sep = "\t", 
    quote = NULL, ...) 
{
    if (is.null(files)) {
        if (is.null(ext)) 
            stop("Must specify input files")
        else {
            extregex <- paste("\\.", ext, "$", sep = "")
            files <- dir(path = ifelse(is.null(path), ".", path), 
                pattern = extregex)
            files <- sub(extregex, "", files)
        }
    }
    source <- match.arg(source, c("generic", "agilent", "agilent.median", 
        "arrayvision", "arrayvision.ARM", "arrayvision.MTM", 
        "bluefuse", "genepix", "genepix.mean", "genepix.median", 
        "genepix.custom", "imagene", "imagene9", "quantarray", 
        "scanarrayexpress", "smd.old", "smd", "spot", "spot.close.open"))
    source2 <- strsplit(source, split = ".", fixed = TRUE)[[1]][1]
    if (is.null(quote)) 
        if (source == "agilent") 
            quote <- ""
        else quote <- "\""
    if (source2 == "imagene") 
        return(read.imagene(files = files, path = path, ext = ext, 
            names = names, columns = columns, other.columns = other.columns, 
            wt.fun = wt.fun, verbose = verbose, sep = sep, quote = quote, 
            ...))
    if (is.data.frame(files)) {
        targets <- files
        files <- files$FileName
        if (is.null(files)) 
            stop("targets frame doesn't contain FileName column")
        if (is.null(names)) 
            names <- targets$Label
    }
    else {
        targets <- NULL
    }
    slides <- as.vector(as.character(files))
    if (!is.null(ext)) 
        slides <- paste(slides, ext, sep = ".")
    nslides <- length(slides)
    if (is.null(names)) 
        names <- removeExt(files)
    if (is.null(columns)) {
        if (source2 == "generic") 
            stop("must specify columns for generic input")
        columns <- switch(source, agilent = list(G = "gMeanSignal", 
            Gb = "gBGMedianSignal", R = "rMeanSignal", Rb = "rBGMedianSignal"), 
            agilent.median = list(G = "gMedianSignal", Gb = "gBGMedianSignal", 
                R = "rMedianSignal", Rb = "rBGMedianSignal"), 
            arrayvision = , arrayvision.ARM = list(G = "ARM Dens - Levels", 
                Gb = "Bkgd", R = "ARM Dens - Levels", Rb = "Bkgd"), 
            arrayvision.MTM = list(G = "MTM Dens - Levels", Gb = "Bkgd", 
                R = "MTM Dens - Levels", Rb = "Bkgd"), bluefuse = list(G = "AMPCH1", 
                R = "AMPCH2"), genepix = , genepix.mean = list(R = "F635 Mean", 
                G = "F532 Mean", Rb = "B635 Median", Gb = "B532 Median"), 
            genepix.median = list(R = "F635 Median", G = "F532 Median", 
                Rb = "B635 Median", Gb = "B532 Median"), genepix.custom = list(R = "F635 Mean", 
                G = "F532 Mean", Rb = "B635", Gb = "B532"), quantarray = list(R = "ch2 Intensity", 
                G = "ch1 Intensity", Rb = "ch2 Background", Gb = "ch1 Background"), 
            imagene9 = list(R = "Signal Mean 2", G = "Signal Mean 1", 
                Rb = "Background Median 2", Gb = "Background Median 1"), 
            scanarrayexpress = list(G = "Ch1 Mean", Gb = "Ch1 B Median", 
                R = "Ch2 Mean", Rb = "Ch2 B Median"), smd.old = list(G = "CH1I_MEAN", 
                Gb = "CH1B_MEDIAN", R = "CH2I_MEAN", Rb = "CH2B_MEDIAN"), 
            smd = list(G = "Ch1 Intensity (Mean)", Gb = "Ch1 Background (Median)", 
                R = "Ch2 Intensity (Mean)", Rb = "Ch2 Background (Median)"), 
            spot = list(R = "Rmean", G = "Gmean", Rb = "morphR", 
                Gb = "morphG"), spot.close.open = list(R = "Rmean", 
                G = "Gmean", Rb = "morphR.close.open", Gb = "morphG.close.open"), 
            NULL)
        if (green.only) {
            columns$R <- columns$Rb <- NULL
            nRG <- 1
            E <- FALSE
        }
        else {
            nRG <- 2
            E <- FALSE
        }
        cnames <- names(columns)
    }
    else {
        columns <- as.list(columns)
        cnames <- names(columns)
        if (is.null(cnames)) {
            if (length(columns) == 1) {
                names(columns) <- "E"
                E <- TRUE
                nRG <- 0
            }
            else {
                stop("columns needs to be a named list")
            }
        }
        else {
            names(columns)[cnames == "Gf"] <- "G"
            names(columns)[cnames == "Rf"] <- "R"
            cnames <- names(columns)
            nRG <- sum(c("R", "G") %in% cnames)
            E <- ("E" %in% cnames)
            if (E && nRG > 0) 
                stop("columns can be R,G for two color data, or E for single channel, but not both")
            if (!E && nRG == 0) 
                stop("columns must specify foreground G or R or E")
            if (!all(cnames %in% c("G", "R", "Gb", "Rb", "E", 
                "Eb"))) 
                warning("non-standard columns specified")
        }
    }
    if (is.null(annotation)) 
        annotation <- switch(source, agilent = , agilent.median = c("Row", 
            "Col", "Start", "Sequence", "SwissProt", "GenBank", 
            "Primate", "GenPept", "ProbeUID", "ControlType", 
            "ProbeName", "GeneName", "SystematicName", "Description"), 
            arrayvision = , arrayvision.ARM = , arrayvision.MTM = c("Spot labels", 
                "ID"), bluefuse = c("ROW", "COL", "SUBGRIDROW", 
                "SUBGRIDCOL", "BLOCK", "NAME", "ID"), genepix = , 
            genepix.median = , genepix.custom = c("Block", "Row", 
                "Column", "ID", "Name"), imagene9 = c("Meta Row", 
                "Meta Column", "Row", "Column", "Gene ID"), quantarray = c("Array Row", 
                "Array Column", "Row", "Column", "Name"), scanarrayexpress = c("Array Row", 
                "Array Column", "Spot Row", "Spot Column"), smd = c("Spot", 
                "Clone ID", "Gene Symbol", "Gene Name", "Cluster ID", 
                "Accession", "Preferred name", "Locuslink ID", 
                "Name", "Sequence Type", "X Grid Coordinate (within sector)", 
                "Y Grid Coordinate (within sector)", "Sector", 
                "Failed", "Plate Number", "Plate Row", "Plate Column", 
                "Clone Source", "Is Verified", "Is Contaminated", 
                "Luid"), smd.old = c("SPOT", "NAME", "Clone ID", 
                "Gene Symbol", "Gene Name", "Cluster ID", "Accession", 
                "Preferred name", "SUID"), NULL)
    fullname <- slides[1]
    if (!is.null(path)) 
        fullname <- file.path(path, fullname)
    required.col <- unique(c(annotation, unlist(columns), other.columns))
    text.to.search <- if (is.null(wt.fun)) 
        ""
    else deparse(wt.fun)
    switch(source2, quantarray = {
        firstfield <- scan(fullname, what = "", sep = "\t", flush = TRUE, 
            quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE, 
            allowEscapes = FALSE)
        skip <- grep("Begin Data", firstfield)
        if (length(skip) == 0) stop("Cannot find \"Begin Data\" in image output file")
        nspots <- grep("End Data", firstfield) - skip - 2
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, nrows = nspots, flush = TRUE, ...)
    }, arrayvision = {
        skip <- 1
        cn <- scan(fullname, what = "", sep = sep, quote = quote, 
            skip = 1, nlines = 1, quiet = TRUE, allowEscape = FALSE)
        fg <- grep(" Dens - ", cn)
        if (length(fg) != 2) stop(paste("Cannot find foreground columns in", 
            fullname))
        bg <- grep("^Bkgd$", cn)
        if (length(bg) != 2) stop(paste("Cannot find background columns in", 
            fullname))
        columns <- list(R = fg[1], Rb = bg[1], G = fg[2], Gb = bg[2])
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, flush = TRUE, ...)
        fg <- grep(" Dens - ", names(obj))
        bg <- grep("^Bkgd$", names(obj))
        columns <- list(R = fg[1], Rb = bg[1], G = fg[2], Gb = bg[2])
        nspots <- nrow(obj)
    }, bluefuse = {
        skip <- readGenericHeader(fullname, columns = c(columns$G, 
            columns$R))$NHeaderRecords
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, flush = TRUE, ...)
        nspots <- nrow(obj)
    }, genepix = {
        h <- readGPRHeader(fullname)
        if (verbose && source == "genepix.custom") cat("Custom background:", 
            h$Background, "\n")
        skip <- h$NHeaderRecords
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, flush = TRUE, ...)
        nspots <- nrow(obj)
    }, imagene9 = {
        h <- readImaGeneHeader(fullname)
        skip <- h$NHeaderRecords
        FD <- h$"Field Dimensions"
        if (is.null(FD)) stop("Can't find Field Dimensions in ImaGene header")
        nspots <- sum(apply(FD, 1, prod))
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, flush = TRUE, nrows = nspots, ...)
    }, smd = {
        skip <- readSMDHeader(fullname)$NHeaderRecords
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, flush = TRUE, ...)
        nspots <- nrow(obj)
    }, {
        skip <- readGenericHeader(fullname, columns = columns, 
            sep = sep)$NHeaderRecords
        obj <- read.columns(fullname, required.col, text.to.search, 
            skip = skip, sep = sep, quote = quote, stringsAsFactors = FALSE, 
            fill = TRUE, flush = TRUE, ...)
        nspots <- nrow(obj)
    })
    Y <- matrix(NA, nspots, nslides)
    colnames(Y) <- names
    RG <- columns
    for (a in cnames) RG[[a]] <- Y
    if (!is.null(wt.fun)) 
        RG$weights <- Y
    if (is.data.frame(targets)) {
        RG$targets <- targets
    }
    else {
        RG$targets <- data.frame(FileName = files, row.names = names, 
            stringsAsFactors = FALSE)
    }
    if (!is.null(annotation)) {
        j <- match(annotation, colnames(obj), 0)
        if (any(j > 0)) 
            RG$genes <- data.frame(obj[, j, drop = FALSE], check.names = FALSE)
    }
    RG$source <- source
    if (source2 == "agilent") {
        if (!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
            nr <- length(unique(RG$genes$Row))
            nc <- length(unique(RG$genes$Col))
            if (nspots == nr * nc) 
                RG$printer <- list(ngrid.r = 1, ngrid.c = 1, 
                  nspot.r = nr, nspot.c = nc)
        }
    }
    if (source2 == "genepix") {
        if (!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && 
            !is.null(RG$genes$Column)) {
            RG$printer <- getLayout(RG$genes, guessdups = FALSE)
            nblocks <- RG$printer$ngrid.r * RG$printer$ngrid.c
            if (!is.na(nblocks) && (nblocks > 1) && !is.null(obj$X)) {
                blocksize <- RG$printer$nspot.r * RG$printer$nspot.c
                i <- (1:(nblocks - 1)) * blocksize
                ngrid.r <- sum(obj$X[i] > obj$X[i + 1]) + 1
                if (!is.na(ngrid.r) && nblocks%%ngrid.r == 0) {
                  RG$printer$ngrid.r <- ngrid.r
                  RG$printer$ngrid.c <- nblocks/ngrid.r
                }
                else {
                  warning("Can't determine number of grid rows")
                  RG$printer$ngrid.r <- RG$printer$ngrid.c <- NA
                }
            }
        }
    }
    if (source2 == "imagene9") {
        printer <- list(ngrid.r = FD[1, "Metarows"], ngrid.c = FD[1, 
            "Metacols"], nspot.r = FD[1, "Rows"], nspot.c = FD[1, 
            "Cols"])
        if (nrow(FD) == 1) {
            RG$printer <- printer
        }
        else {
            printer$ngrid.r <- sum(FD[, "Metarows"])
            if (all(printer$ngrid.c == FD[, "Metacols"]) && all(printer$nspot.r == 
                FD[, "Rows"]) && all(printer$nspot.c == FD[, 
                "Cols"])) 
                RG$printer <- printer
        }
    }
    if (!is.null(other.columns)) {
        other.columns <- as.character(other.columns)
        j <- match(other.columns, colnames(obj), 0)
        if (any(j > 0)) {
            other.columns <- colnames(obj)[j]
            RG$other <- list()
            for (j in other.columns) RG$other[[j]] <- Y
        }
        else {
            other.columns <- NULL
        }
    }
    for (i in 1:nslides) {
        if (i > 1) {
            fullname <- slides[i]
            if (!is.null(path)) 
                fullname <- file.path(path, fullname)
            switch(source2, quantarray = {
                firstfield <- scan(fullname, what = "", sep = "\t", 
                  flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE, 
                  multi.line = FALSE, allowEscapes = FALSE)
                skip <- grep("Begin Data", firstfield)
            }, arrayvision = {
                skip <- 1
            }, genepix = {
                skip <- readGPRHeader(fullname)$NHeaderRecords
            }, smd = {
                skip <- readSMDHeader(fullname)$NHeaderRecords
            }, {
                skip <- readGenericHeader(fullname, columns = columns)$NHeaderRecords
            })
            if (verbose && source == "genepix.custom") 
                cat("Custom background:", h$Background, "\n")
            obj <- read.columns(fullname, required.col, text.to.search, 
                skip = skip, sep = sep, stringsAsFactors = FALSE, 
                quote = quote, fill = TRUE, nrows = nspots, flush = TRUE, 
                ...)
        }
        for (a in cnames) RG[[a]][, i] <- obj[, columns[[a]]]
        if (!is.null(wt.fun)) 
            RG$weights[, i] <- wt.fun(obj)
        if (!is.null(other.columns)) 
            for (j in other.columns) {
                RG$other[[j]][, i] <- obj[, j]
            }
        if (verbose) 
            cat(paste("Read", fullname, "\n"))
    }
    if (nRG == 1) {
        n <- names(RG)
        n[n == "G"] <- "E"
        n[n == "Gb"] <- "Eb"
        n[n == "R"] <- "E"
        n[n == "Rb"] <- "Eb"
        names(RG) <- n
    }
	RG
}
targets <- readTargets(file)
RG<-onechannel(targets$FileName,source=technology,green.only=TRUE)
RGb<-backgroundCorrect(RG,method=bg_method)
result<-cbind(RGb$genes$ProbeName,log(RGb$E))
colnames(result)<-c('GeneID',targets$Samples)
write.table(result,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\box_result.txt',sep='\t',row.names=FALSE)
}
#######################################################################
###FUNCION ENRIQUECIMIENTO BIOLOGICO
#######################################################################
enrichment<-function(file,chip_name,annotation,threshold_element,threshold,compare){
if (!("topGO" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("topGO")
    }
input<-read.table(file,sep='\t',header=TRUE)
f<-nchar(chip_name,type='char')
for (i in 1:f){
	b<-sub("affy",'',chip_name,ignore.case=T) 
	chip_name<-b
	b<-sub("_", '', chip_name, ignore.case=T)
	chip_name<-b
	}
for (i in 1:length(colnames(input))){
	if(threshold_element==colnames(input)[i]){
	column<-i
	}
}
new_matrix<-NULL
if (compare=='equal'){
	for (i in 1:length(input[,column])){
		if (as.numeric(as.vector(input[,column][i]))== as.numeric(threshold)){
			new_matrix<-rbind(new_matrix,input[i,])
		}
	}
}
if (compare=='lower'){
	for (i in 1:length(input[,column])){
		if (as.numeric(as.vector(input[,column][i]))<= as.numeric(threshold)){
			new_matrix<-rbind(new_matrix,input[i,])
		}
	}
}
if (compare=='bigger'){
	for (i in 1:length(input[,column])){
		if (as.numeric(as.vector(input[,column][i]))>= as.numeric(threshold)){
			new_matrix<-rbind(new_matrix,input[i,])
		}
	}
}

n<-NULL
for (i in 1:length(colnames(new_matrix))){
	if (colnames(new_matrix)[i]=="ID"){
		n<-i
		}
	}
affyLib<-chip_name
affyLib.db <- paste(chip_name, "db", sep = ".")

# Check and upload the packages required

if (!(affyLib.db %in% installed.packages())) { 
	source("http://www.bioconductor.org/getBioC.R")
        biocLite(affyLib.db)
	}

load <- library(package = affyLib.db, character.only = TRUE, logical.return = TRUE)

if (!load) stop("Biological Enrichment for ", affyLib ," genome is not supported.")

affyLibGO<-paste(affyLib,"GO",sep="")
library(package = affyLib.db, character.only = TRUE)
library('topGO')
allgenes<-eval(parse(text=paste("Lkeys(",affyLibGO,")",sep="")))

interesting <- as.character(new_matrix[,n])
interesting<-append(interesting,'NA')

geneList <- factor(as.integer(interesting %in% allgenes))

names(geneList)<-new_matrix[,n]
sampleGOdata <- new("topGOdata", ontology = annotation, allGenes = geneList,annot = annFUN.db, affyLib = affyLib.db)
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,topNodes=20)
write.table(allRes,'c:\\PG\\processes\\data\\enrichment_result.txt',sep='\t',row.names=FALSE)

}



#######################################################################
###FUNCION OUTER ENRICHMENT
#######################################################################
enrichmentouter<-function(file,chip_name,annotation,threshold_element,threshold,compare){
input<-read.table(file,sep='\t',header=FALSE)
input<-as.vector(input$V1)
f<-nchar(chip_name,type='char')
for (i in 1:f){
	b<-sub("affy",'',chip_name,ignore.case=T) 
	chip_name<-b
	b<-sub("_", '', chip_name, ignore.case=T)
	chip_name<-b
	}
affyLib<-chip_name
affyLib.db <- paste(chip_name, "db", sep = ".")

# Check and upload the packages required

if (!(affyLib.db %in% installed.packages())) { 
	source("http://www.bioconductor.org/getBioC.R")
        biocLite(affyLib.db)
	}
if (!(topGO %in% installed.packages())) { 
	source("http://www.bioconductor.org/getBioC.R")
        biocLite(topGO)
	}
load <- library(package = affyLib.db, character.only = TRUE, logical.return = TRUE)

if (!load) stop("Biological Enrichment for ", affyLib ," genome is not supported.")

affyLibGO<-paste(affyLib,"GO",sep="")
library(package = affyLib.db, character.only = TRUE)
library('topGO')
allgenes<-eval(parse(text=paste("Lkeys(",affyLibGO,")",sep="")))

interesting <- input
interesting<-append(interesting,'NA')

geneList <- factor(as.integer(interesting %in% allgenes))

names(geneList)<-input
sampleGOdata <- new("topGOdata", ontology = annotation, allGenes = geneList,annot = annFUN.db, affyLib = affyLib.db)
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,topNodes=20)
write.table(allRes,'c:\\PG\\processes\\data\\enrichment_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION CORRELACION
#######################################################################
correlation<-function(file){
a<-read.table(file,sep='\t',header=T,check.names=FALSE)
y<-colnames(a)
y<-as.matrix(y[2:length(y)])
a<-a[2:length(a)]
d<-colnames(a)
x<-NULL
for (i in 1:length(y)){
	for (j in 1:length(d)){
		if (y[i]==d[j]){
			x<-cbind(x,as.double(as.vector(a[,j][2:length(a[,j])])))
		}
	}
}
x[is.infinite(x)]<-0
b<-cor(x,use="pairwise.complete.obs")
c<-cbind(y,b)
colnames(c)<-c('SampleName',y)
write.table(c,'c:\\PG\\processes\\data\\cor_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION BOXPLOT
#######################################################################
box<-function(file,method,technology){
a<-read.table(file,sep='\t',header=T,check.names=FALSE)
y<-colnames(a)
if (technology!='two_channel'){
	y<-as.matrix(y[2:length(a)])
	}
if (technology=='two_channel'){
	y<-as.matrix(y[1:length(a)])
	}
a<-a[2:length(a)]
d<-colnames(a)
x<-NULL
for (i in 1:length(y)){
	for (j in 1:length(d)){
		if (y[i]==d[j]){
			x<-cbind(x,as.double(as.vector(a[,j][2:length(a[,j])])))
		}
	}
}
x[is.infinite(x)]<-0
b<-mean(as.data.frame(x),na.rm=TRUE)
z<-max(as.data.frame(b))
if (method=="Standard deviation"){
	c<-sd(x,na.rm=TRUE)
	}
if (method=="Variance"){
	c<-var(x,na.rm=TRUE)
	c<-c[1,]
	}
d<-cbind(y,b,c)
colnames(d)<-c('SampleName','Mean','StandardDeviation')
write.table(d,'c:\\PG\\processes\\data\\boxplot.txt',sep='\t',row.names=FALSE)
write.table(z,'c:\\PG\\processes\\data\\boxplot_max.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION ALL AFFYMETRIX
#######################################################################
all.affymetrix<-function(file,dir, analysis){
if (!("affy" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("affy")
    }

if (!("gcrma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("gcrma")
    }
library('affy')
library('gcrma')
design_data<-read.table(file,sep="\t",header=FALSE)
samples<-as.character(design_data[,1])
files<-paste(samples,'CEL',sep='.')
setwd(dir)
if (analysis=="RMA"){
	exp_data<-justGCRMA(sampleNames=c(samples),filenames=c(files))
	gene_id<-featureNames(exp_data)
	}
result<-as.data.frame(exprs(exp_data))
result<-cbind(gene_id,result)
controles<-grep("AFFX.+",featureNames(exp_data))
control_id<-featureNames(exp_data)[controles]
a<-controles[1]
b<-controles[length(controles)]
d<-(exp_data)[a:b,]
e<-cbind(control_id,2^(exprs(d)))

controles_Hyb<-grep("AFFX-Bio.+",featureNames(exp_data))
hyb_id<-featureNames(exp_data)[controles_Hyb]
f<-controles_Hyb[1]
g<-controles_Hyb[length(controles_Hyb)]
h<-(exp_data)[f:g,]
i<-cbind(hyb_id,2^(exprs(h)))

write.table(result,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(result,'c:\\PG\\processes\\data\\affymetrix_result.txt',sep='\t',row.names=FALSE)
write.table(e,'c:\\PG\\processes\\data\\control_result.txt',sep='\t',row.names=FALSE)
write.table(i,'c:\\PG\\processes\\data\\hyb_result.txt',sep='\t',row.names=FALSE)

}


#######################################################################
###FUNCION HIERARCHICAL CLUSTERING
#######################################################################
cluster<-function(file,cluster_method,distance_method,columns,filter_element,threshold,high,medium,low){
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
c<-as.matrix(c)
d<-b[1,]
b<-as.matrix(b)
n<-NULL
m<-NULL
threshold<-as.numeric(threshold)
for (i in 1:length(d)){
	if (d[i]==filter_element){
		n<-i
		}
	}
for (i in 1:length(d)){
	if (d[i]=="ID"){
		m<-i
		}
	}
x<-NULL
for (i in 1:length(c)){
	for (j in 1:length(d)){
		if (c[i]==d[j]){
			x<-cbind(x,as.double(as.vector(b[,j][2:length(b[,j])])))
		}
	}
}
x[is.na(x)==TRUE]<-0
gene_name<-NULL
for (i in 1:length(b[,1])){
	if (b[i,n]<=threshold){
		gene_name<-rbind(gene_name,b[i,][m])
		}
	}
y<-NULL
a<-length(b[,1])
p=length(gene_name)
gene_list<-NULL
if (p>1000){
	p=1000
	}
for (i in 1:(a-1)){
	for (j in 1:p){
		if (gene_name[j]==as.vector(b[i,1])){
			if (length(which(gene_list==gene_name[j]))!=1){
				gene_list<-append(gene_list,gene_name[j])
				y<-rbind(y,x[i,])
				}
			}
		}
	}
if (distance_method=='pearson'){
	h<-as.dendrogram(hclust(as.dist(1-cor(t(y),method="pearson")),method=cluster_method))
	if (columns=='TRUE'){
		hc<-as.dendrogram(hclust(as.dist(1-cor((y),method="pearson")),method=cluster_method))
		}
	}
if (distance_method!='pearson'){
	h<-as.dendrogram(hclust(dist (y,is.na(y),method=distance_method),method=cluster_method))
	if (columns=='TRUE'){
		hc<-as.dendrogram(hclust(dist(t(y),is.na(y),method=distance_method),method=cluster_method))
		}
	}
try((rownames(y)<-gene_list),silent = TRUE)

colnames(y)<-c
if (columns!='TRUE'){
	heatmap(y,Rowv=h,Colv=NA,na.rm=TRUE,col=colorRampPalette(c(high,medium,low))(256),keep.dendro=TRUE,cexCol=0.5,cexRow=0.7)
	}
if (columns=='TRUE'){
	heatmap(y,Rowv=h,Colv=hc,na.rm=TRUE,col=colorRampPalette(c(high,medium,low))(256),cexCol=0.5,cexRow=0.7)	
	}
}



#######################################################################
###FUNCION DATABASE ANNOTATION
#######################################################################

databaseannotation<-function(file,id_column,organism,filter,attribute){
if (!("biomaRt" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("biomaRt")
    }
library('biomaRt')
a<-read.table(file,sep='\t',header=T)
b<-colnames(a)
data<-NULL
for (i in 1:length(b)){
	if (b[i]==id_column){
		data<-a[,i]
		}
	}
data<-as.vector(data)
if (organism!="athaliana_eg_gene"){
	ensembl = useMart("ensembl")
	ensembl = useDataset(organism,mart=ensembl)
	result<-getBM(attributes=c(filter,attribute),filters=filter,values=c(data), mart=ensembl)
	}
if (organism=="athaliana_eg_gene"){
	ensembl = useMart("ENSEMBL_MART_PLANT")
	ensembl = useDataset(organism,mart=ensembl)
	result<-getBM(attributes=c(filter,attribute),filters=filter,values=c(data), mart=ensembl)
	}
final_data<-NULL
for (i in 1:(length(result[,1])-1)){
	if (result[i,1]!=result[i+1,1]){
		final_data<-rbind(final_data,result[i,])
		}
	}
write.table(final_data,'c:\\PG\\processes\\data\\annotation.txt',sep='\t',row.names=FALSE)
}




#######################################################################
###FUNCION T_TEST
#######################################################################
ttest<-function(file){
if (!("fdrtool" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("fdrtool")
	}
library('fdrtool')
b<-read.table(file,sep='\t')
c<-read.table('c:\\PG\\processes\\data\\samples.txt',sep='\t')
w<-read.table('c:\\PG\\processes\\data\\design.txt',sep='\t')
a<-length(w[,1])
b<-as.matrix(b)
d<-b[2:length(b[,1]),]
e<-d[,2:length(b[1,])]
control<-NULL
target<-NULL
for (i in 1:a){
	for (j in 1:length(b[1,])){
		if (c[i,]==b[1,][j]){
			if (w[i,1]==0){
				control<-cbind(control,e[,i])
				}
			if (w[i,1]==1){
				target<-cbind(target,e[,i])
				}
			}
		}
	}


z<-length(control[,1])
y<-NULL
for(i in 1:z){
	w<-mean(as.numeric(target[i,]))/mean(as.numeric(control[i,]))
	x<-t.test(as.numeric(target[i,]),as.numeric(control[i,]))
	q<-c(w,x$p.value)
	y<-rbind(y,q)
	}
fdr<-fdrtool(y[,2], statistic="pvalue")
graphics.off()
w<-cbind(d[,1],e,d[,1],y,fdr$lfdr)
colnames(w)<-c('GeneID',b[1,][2:length(b[1,])],'ID','LogFC','P.Value','adj.P.Val')
write.table(w,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION QPCR
#######################################################################
qPCR=function(file,aggregation,missing_value){
a=read.table(file,sep='\t',header=T)
b=colnames(a)
n=0
m=0
p=0
for (i in 1:length(b)){
	if (b[i]=='Sample'){
		n=i
	}
}
for (i in 1:length(b)){
	if (b[i]=='Detector'){
		m=i
	}
}
for (i in 1:length(b)){
	if (b[i]=='Ct'){
		p=i
	}
}

samples=unique(a[,n])
detector=as.vector(a[,m])
Cts=as.vector(a[,p])
detector_sample=paste(a[,m],a[,n],sep='_')
unique_detector=unique(detector_sample)
calculated_Ct=NULL
unique_detectors=unique(detector)
for(i in 1:length(unique_detector)){
	value=NULL
	for (j in 1:length(detector_sample)){
		if (unique_detector[i]==detector_sample[j]){
			value=append(value,as.double(Cts[j]))
				}
			}
	if (aggregation=='median'){
		calculated_Ct=rbind(calculated_Ct,median(value))
		}
	if (aggregation=='mean'){
		calculated_Ct=rbind(calculated_Ct,mean(value))
		}
	}
result=matrix(calculated_Ct,nrow=length(unique_detectors), ncol=length(samples))
result=cbind(unique_detectors,result)
result[is.na(result)==TRUE]=missing_value
colnames(result)=c('gene_id',as.vector(samples))
write.table(result,'C:\\PG\\processes\\data\\result.txt',sep='\t',col.names=TRUE,row.names=FALSE)
write.table(result,'C:\\PG\\processes\\data\\raw_data.txt',sep='\t',col.names=TRUE,row.names=FALSE)
write.table(result,'C:\\PG\\processes\\data\\custom_result.txt',sep='\t',col.names=TRUE,row.names=FALSE)
}


#######################################################################
###FUNCION PIVOT
#######################################################################
pivot=function(file,columns,rows,values,aggregation){
a=read.table(file,sep='\t',header=T)
b=colnames(a)
n=0
m=0
p=0
for (i in 1:length(b)){
	if (b[i]==columns){
		n=i
	}
}
for (i in 1:length(b)){
	if (b[i]==rows){
		m=i
	}
}
for (i in 1:length(b)){
	if (b[i]==values){
		p=i
	}
}
samples=unique(a[,n])
detector=as.vector(a[,m])
Cts=as.vector(a[,p])
detector_sample=paste(a[,m],a[,n],sep='_')
unique_detector=unique(detector_sample)
calculated_Ct=NULL
unique_detectors=unique(detector)
for(i in 1:length(unique_detector)){
	value=NULL
	for (j in 1:length(detector_sample)){
		if (unique_detector[i]==detector_sample[j]){
			value=append(value,as.double(Cts[j]))
				}
			}
	if (aggregation=='median'){
		calculated_Ct=rbind(calculated_Ct,median(value))
		}
	if (aggregation=='mean'){
		calculated_Ct=rbind(calculated_Ct,mean(value))
		}
	}
result=matrix(calculated_Ct,nrow=length(unique_detectors), ncol=length(samples))
result=cbind(unique_detectors,result)
result[is.na(result)==TRUE]=0
colnames(result)=c('gene_id',as.vector(samples))
write.table(result,'C:\\PG\\processes\\data\\result.txt',sep='\t',col.names=TRUE,row.names=FALSE)
}


#######################################################################
###FUNCION CITOMETRIA
#######################################################################
flowcitometry<-function(file){

if (!("flowCore" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("flowCore")
    }
if (!("flowQ" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("flowQ")
    }
if (!("flowViz" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("flowViz")
    }
library('flowCore')
library('flowQ')
library('flowViz')
flowData <- read.FCS(file)
matrixData<- exprs(flowData)
write.table(matrixData,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(matrixData,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
write.table(matrixData,'c:\\PG\\processes\\data\\raw_data.txt',sep='\t',row.names=FALSE)
}


#######################################################################
###FUNCION CHROMINFO
#######################################################################
chrominfo<-function(file){
X=read.table(file,sep='\t',header=T)
mychrpo <- as.vector(X$GeneID)
chr <- vector()
po <- vector()
num <- 0
for (k in 1:length(mychrpo)) {
	if ( (length(unlist(strsplit(mychrpo[k],":"))) == 2) & (unlist(regexpr('_random',mychrpo[k]))[1] == -1) ) {
		chr_ <- unlist(strsplit(unlist(strsplit(mychrpo[k],":"))[1],"chr"))[2]
		if ( chr_ == "X") {
			chr[k] <- 23
		} else if (chr_ == "Y") {
			chr[k] <- 24
		} else {
			chr[k] <- chr_
		}			 
		po[k] <- unlist(strsplit(unlist(strsplit(mychrpo[k],":"))[2],"-"))[2]
		num <- num + 1
	} else {
		chr[k] <- NA
		po[k] <- NA			
	}
}
X$Chr <- as.numeric(chr)
X$Position <- as.numeric(po)/1000000
write.table(X,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION SEGMENTATION
#######################################################################
segmentation<-function(file,file2,designs,method){
if (!("snapCGH" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("snapCGH")
	}
if (!("strucchange" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("strucchange")
	}
library('snapCGH')
library('strucchange')
X=read.table(file,sep='\t',header=T)
r<-X
Y=read.table(file2,sep='\t',header=T)
Z=read.table(designs,sep='\t',header=T)
b<-grep('normalized',colnames(X))
q<-NULL
q<-cbind(X[1])
q<-cbind(q,X[b])
q<-cbind(q,X[length(X)-1])
q<-cbind(q,X[length(X)])
X<-q
MA<-NULL
MA<-new("MAList")
MA$genes<-as.data.frame(X)
z<-length(X)-2
MA$M<-as.vector(X[,2:z])
MA$A<-as.vector(Y)
MA$design<-as.vector(Z)
if (method=='HomHMM'){
	MA2<-runHomHMM(MA)
	}
if (method=='DNAcopy'){
	MA2<-runDNAcopy(MA)
	}
if (method=='GLAD'){
	MA2<-runGLAD(MA)
	}
a<-sub("normalized", "predicted", colnames(MA2$M.predicted), ignore.case=T) 
f<-sub("normalized", "state", colnames(MA2$state),ignore.case=T)
colnames(MA2$M.predicted)<-a
colnames(MA2$state)<-f
w<-cbind(r,MA2$M.predicted,MA2$state)
write.table(w,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(w,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###GRAFICO DE CORRELACION
#######################################################################
circlecorr <- function(file, col=c("black","white"), bg = "white", 
	cex = 1, order = FALSE, title = "", ...){

	a<-read.table(file,sep='\t',header=T)
	corr<-cor(a[2:length(a)],use="pairwise.complete.obs")
    
	if (is.null(corr)) 
        return(invisible())
    if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 
        6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1)) 
        stop("Need a correlation matrix!")
    n <- nrow(corr)
    m <- ncol(corr)
    
    ## reorder the variables using principal component analysis
    if (order) {
    	if(!n==m){
    		stop("The matrix must be squre if order is TRUE!")
    	}
      x.eigen <- eigen(corr)$vectors[, 1:2]
      e1 <- x.eigen[, 1]
      e2 <- x.eigen[, 2]
      alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
      corr <- corr[order(alpha), order(alpha)]
    }
    
    ## set up variable names
    rname <- rownames(corr)
    cname <- colnames(corr)
    if (is.null(rname)) 
        rname <- 1:n
    if (is.null(cname)) 
        cname <- 1:m
    rname <- as.character(rname)
    cname <- as.character(cname)

    ## calculate label-text width approximately
    par(mar = c(0, 0, 2, 0), bg = "white")
    plot.new()
    plot.window(c(0, m), c(0, n), asp = 1)
    xlabwidth <- max(strwidth(rname, cex = cex))
    ylabwidth <- max(strwidth(cname, cex = cex))

    ## set up an empty plot with the appropriate dimensions
    plot.window(c(-xlabwidth + 0.5, m + 0.5), c(0, n + 1 + ylabwidth),
                asp = 1, xlab="", ylab="")
    rect(0.5, 0.5, m + 0.5, n + 0.5, col = bg)	##background color

    ## add variable names and title
    text(rep(-xlabwidth/2, n), n:1, rname, col = "red", cex = cex)
    text(1:m, rep(n + 1 + ylabwidth/2, m), cname, srt = 90, col = "red", 
        cex = cex)
    title(title)

    ## add grid
    segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1), 
        0.5 + 0:n, col = "gray")
    segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5, 
        m), col = "gray")

    ## assign circles' fill color
    nc <- length(col)
    if(nc==1)
        	bg <- rep(col, n*m)
    else{
        ff <- seq(-1,1, length=nc+1)
        bg2 = rep(0, n * m)
        for (i in 1:(n * m)){
            bg2[i] <- rank(c(ff[2:nc], as.vector(corr)[i]), 
                            ties.method = "random")[nc]
        }
        bg <- (col[nc:1])[bg2]
    }

    ## plot n*m circles using vector language, suggested by Yihui Xie
    ## the area of circles denotes the absolute value of coefficient
    symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F, 
        circles = as.vector(sqrt(abs(corr))/2), bg = as.vector(bg))
}



#######################################################################
###GRAFICO BOXPLOT
#######################################################################
boxview<-function(file,boxcolor,font_size){
	a<-read.table(file,sep='\t',header=T)
	boxplot(a[2:length(a)],col=boxcolor,notch=TRUE,cex.axis=(font_size/18),las=2)
}


######################################################################
###GRAFICO PERFILES
######################################################################
profile<- function(file,font_size){
datos<-read.table(file,sep='\t',header=T)
maximo<-max(datos[2:length(datos)])
samples<-colnames(datos)
if (samples[1]=='probeset_id'){
	x<-c(1:length(datos$probeset_id))
	}
if (samples[1]=='hyb_id'){
	x<-c(1:length(datos$hyb_id))
	}
maximo<-ceiling(maximo)
colors = rainbow(ncol(datos)-1)
plot(x,datos[,samples[2]],xlab=' ',ylab='Expression levels',ylim=c(0,maximo),type='l',col='blue',xaxt='n',lwd=3)
if (samples[1]=='hyb_id'){
	axis(1,at=1:length(datos[,1]),labels=c(as.character(datos[,1])),las=2,cex.axis=(font_size/18))
	}
if (samples[1]=='probeset_id'){
	axis(1,at=1:length(datos[,1]),labels=c(as.character(datos[,1])),las=2,cex.axis=(font_size/18))
	}
for (i in 3:length(samples)){
	lines(datos[,samples[i]],col=colors[i%%length(colors)],lwd=3)
	}
}

#######################################################################
###GENE EXPRESSION OMNIBUS
#######################################################################

GEO<-function(geodata){
if (!("GEOquery" %in% installed.packages())) {
  	  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("GEOquery")
	}
library('GEOquery')
datos<-getGEO(geodata)
datos<-datos[[1]]
data<-exprs(datos)
write.table(data,'c:\\PG\\processes\\data\\result.txt',sep='\t',row.names=FALSE)
write.table(data,'c:\\PG\\processes\\data\\custom_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###COX MODEL
#######################################################################
coxsurvive<-function(file,design,event,time){
if (!("survival" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("survival")
    }
	library('survival')
	data<-read.table(file,sep='\t',header=TRUE,stringsAsFactor=TRUE)
	datos<-data[2:length(data)]
	info<-read.table(design,sep='\t',header=TRUE,stringsAsFactor=TRUE)
	lista<-list(time=info$time,evento=info$event)
	resultado<-matrix(NA,nrow=nrow(datos),ncol=5)
	for (i in 1:nrow(datos)){
		x<-datos[i,]
		cox<-summary(coxph( Surv(time, evento) ~ as.numeric(x), lista))
		resultado[i,]<-cox$coef
		}
	resultados<-cbind(data[1],resultado[,1:2],resultado[,4:5])
	colnames(resultados)<-c('Gene ID', 'coef', 'exp(coef)', 'z', 'P.Value')
	write.table(resultados,'c:\\PG\\processes\\data\\survival_result.txt',sep='\t',row.names=FALSE)
}



#######################################################################
###FUNCION DISTRIBUCION INTENSIDADES
#######################################################################
histogram<-function(file,column,color){
a<-read.table(file,sep='\t',header=T)
f<-nchar(column,type='char')
for (i in 1:f){
	b<-sub("-", ".", column, ignore.case=T) 
	column<-b
	b<-sub(" ", ".", column, ignore.case=T) 
	column<-b
	}
q<-0
for(i in 1:length(names(a))){
	if (names(a[i])==column){
		q<-i
		}
	}
hist(a[q][,1],xlab='Distribution',breaks=100,col=color,main = paste('Histogram of ' , names(a[q]),sep=''))
}



#######################################################################
###GRAFICO TIDY_UP
#######################################################################
tidyview<-function(file){
a<-read.table(file,sep='\t',header=T)
b<-grep('Chr',colnames(a))
y<-cbind(a[b],a[b+1])
q<-max(y[1],na.rm=TRUE)
p<-q-2
plot(y,xaxt='n',col='#c9c9f1',xlab='Chromosome',main='Clone position along chromosomes')
axis(1,at=1:p,cex.axis=0.5)
axis(1,at=as.double(p+1):q,labels=c('X','Y'),cex.axis=0.5)
}

#######################################################################
###GRAFICO SEGMENTATION
#######################################################################
predictedview<-function(file,column){
a<-read.table(file,sep='\t',header=T)
b<-names(a)
d<-a$Position
for(i in 1:length(b)){
	if (column==b[i]){
		d<-cbind(d,a[i])
		}
	}
plot(d,type='h',xlab='Position',ylab='M Predicted',main=paste('Gain and Loss positions of ',column),ylim=c(-2,2),xlim=c(0,max(a$Position,na.rm=TRUE)),col='#c9c9f1')
e<-NULL
f<-NULL
for (i in 1:length(d[,1])){
	if (is.na(d[i,2])==FALSE){
		if (d[i,2]>=1){
			e<-rbind(e,d[i,])
			}
		if (d[i,2]<=-1){
			f<-rbind(f,d[i,])
			}
		}
	}
par(new=T)
plot(e,type='h',xlab='Position',ylab='M Predicted',ylim=c(-2,2),xlim=c(0,max(a$Position,na.rm=TRUE)),col='#ff0000')
par(new=T)
plot(f,type='h',xlab='Position',ylab='M Predicted',ylim=c(-2,2),xlim=c(0,max(a$Position,na.rm=TRUE)),col='#00ff00')
}


#######################################################################
###VENN DIAGRAM
####################################################################### 
myvenn <- function (file1, file2, file3){
if (!("limma" %in% installed.packages())) {
		  source("http://www.bioconductor.org/getBioC.R")
		  biocLite("limma")
	}
	library("limma")
	a<- read.table(file1)
	set1 <- t(a)
	b<- read.table(file2)
	set2 <- t(b)
	# Venn diagram for 2 sets
	if (file3=='NULL'){
		extra <- c("x", "y")
		universe <- sort( unique( c(set1, set2, extra) ) )
		Counts <- matrix(0, nrow=length(universe), ncol=2)
		colnames(Counts) <- c("gene list 1", "gene list 2")
		for (i in 1:length(universe))
		{
			Counts[i,1] <- universe[i] %in% set1
			Counts[i,2] <- universe[i] %in% set2			
		}
	}
	if (file3!='NULL'){
		d<- read.table(file3)
		set3 <- t(d)
		extra <- c("x", "y", "z")
		universe <- sort( unique( c(set1, set2, set3, extra) ) )
		Counts <- matrix(0, nrow=length(universe), ncol=3)
		colnames(Counts) <- c("gene list 1", "gene list 2", "gene list 3")
		for (i in 1:length(universe))
		{
			Counts[i,1] <- universe[i] %in% set1
			Counts[i,2] <- universe[i] %in% set2
			Counts[i,3] <- universe[i] %in% set3
		}
	}
	if (file3=='NULL'){
	vennDiagram(Counts,counts.col=c('red','blue'))
	}
	if (file3!='NULL'){
	vennDiagram(Counts,counts.col=c('red','blue','green'))
	}
}



