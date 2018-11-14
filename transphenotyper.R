#no explicit evolutionary model
#no model of wothin host or macro evolution purely a distance based metric- incoprorates natural selection?
library(data.table)
library(reshape2)

readMat<-function(mat)
{
	mat<-as.data.frame(fread(mat, header=T))
	return(mat)
}

getDistance<-function(mat)
{
	dists<-dist(t(mat[, -c(1,2,3,4)]), method='manhattan')
	df<-melt(as.matrix(dists), varnames=c('row','col'))
	return(dist)
}

readDistance<-function(dist)
{
	dist<-read.table(dist, header=T)
	return(dist)
}

#function to return cluster table of cluster size and character vector of sample names from distance table
getClusters<-function(distance_table, distance_cutoff)
{
	pairs<-distance_table[which(distance_table[,3]<=distance_cutoff),]
	remove_row<-c()
	count<-0
	for(row in 1:nrow(pairs))
	{
		if(as.character(pairs[row, 1]==as.character(pairs[row, 2])))
		{
			count<-count+1
			remove_row[count]<-row
		}
	}
	pairs<-pairs[-remove_row,]
	unique<-levels(as.factor(c(as.character(pairs[,1]), as.character(pairs[,2]))))
	cluster_size<-c()
	cluster_col<-c()
	for(u in 1:length(unique))
	{
		primary=unique[u]
		related_1<-pairs[which(pairs[,1]==primary),2] 
		related_2<-pairs[which(pairs[,2]==primary),1]
		related<-c(as.character(related_1), as.character(related_2))
		related<-levels(as.factor(related))
		cluster<-levels(as.factor(c(as.character(primary), as.character(related))))
		cluster_col[u]<-as.character(paste(as.character(cluster),collapse=','))
		cluster_size[u]<-length(cluster)
	} 
	df<-as.data.frame(cbind(cluster_col, cluster_size))
	df[!duplicated(df),]

}

#function to return phenotype data of all sample names and binary or numeric cluster phenotype
getPhenotype<-function(mat, cluster_table, cluster_size_cutoff=2, null_phenotype=0)
{
	clusters<-cluster_table[which(as.numeric(cluster_table$cluster_size)>=cluster_size_cutoff),]
	df<-data.frame()
	for(i in 1:nrow(clusters))
	{
		samples<-unlist(strsplit(as.character(clusters$cluster_col[i]), split=','))
		newrows<-cbind(Sample=c(unlist(samples)), Phenotype=rep(as.numeric(clusters$cluster_size[i]), length=length(samples)))
		df<-rbind(df, newrows)
	}

	othersamples<-colnames(mat)[is.na(match(colnames(mat), df$Sample))]
	newrows<-cbind(Sample=othersamples, Phenotype=rep(null_phenotype,length(othersamples)))
	alldf<-rbind(df, newrows)
	return(alldf)
}


#################################################################################################################################
####Distance cutoff against numebr of clusters##
##could make a video visual using Jody's visualiser##
plotClusterSize<-function(distance_table,  max_distance, bin_size, out)
{
	distances<-seq(from=1, to=max_distance, by=bin_size)
	df<-matrix(data=NA, nrow=length(distances), ncol=2)
	count=1
	for(i in 1:length(distances))
	{
		a<-getClusters(distance_table, distances[i])
		df[count,1]<-distances[i]
		df[count,2]<-dim(a)[1]
		count=count+1
	}
	jpeg(out)
	plot(df[,2], df[,1], xlab='Distance Cutoff', ylab='Number of Clusters')
	dev.off()
}
#################################################################################################################################
##example##
matbin<-'TW20_minDP10.mat.bin'
distance_cutoff<-10
cluster_size_cutoff<-2
out<-'TW20_DC10_CC2.txt'

#################################################################################################################################

mat<-readMat(matbin)
distance_table<-getDistance(mat)
cluster_table<-getClusters(distance_table, distance_cutoff)
phenotype_table<-getPhenotype(mat, cluster_table, cluster_size_cutoff)
write.table(phenotype_table, out, col.name=T, row.name=F, quote=F)

#################################################################################################################################
