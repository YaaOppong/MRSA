#################################################################################################
library(data.table)
library(seqinr)
library(phylotools)

#preprocessing 
getmat2fastainput<-function(nucmersnps, out_prefix)
{
	snps<-read.delim(nucmersnps, header=F)
	snps2<-as.data.frame(cbind(refpos=as.character(snps[,1]), ref=as.character(snps[,2]), alt=as.character(snps[,3])))
	#deal with duplicate positions?
	write.table(snps2, paste(out_prefix, '_refpos_ref_alt.txt', sep=''), quote=F, row.name=F, col.name=T)
}

#functiom to read IUPAC matrix and store as data frame###########################################
readMat<-function(mat)
{
	mat<-as.data.frame(fread(mat, header=T))
	return(mat)
}

#function to read file of additional genome- refpos, ref, alt####################################
readsnps<-function(snps)
{
	m<-read.table(snps, header=T)
	m<-cbind(refpos=m$refpos,alt=as.character(m$alt))
	return(m)
}
#function to return union of postions from IUPAC matrix and additional genomes###################
unionpos<-function(mat, nucmersnps)
{
	matpos<-mat$POS
	nucmersnpspos=nucmersnps[,1]
	union=sort(unique(c(matpos,nucmersnpspos)))
	return(union)
}

#iterate unionpos over multiple additional genomes (snps only)###################################
iterateunionpos<-function(mat, additional_snps)
{
	union=c()
	for(additional_snp in 1:length(additional_snps))
	{
		additional_snp<-readsnps(additional_snps[additional_snp])
		union<-c(union,unionpos(mat, additional_snp))
	}
	union<-sort(unique(union))
	return(union)
}

#function to generate new IUPAC row for monomorphic snps#########################################
lookupnewrowIUPAC<-function(mat, reference, pos)
{
	ref<-reference[pos]
	alt<-'NA'
	newrow<-as.vector(unlist(c(pos, ref, alt, 'SNP', rep(ref, (length(colnames(mat)) - 4)))))
	return(newrow)
}

#function to generate new binary row for monomorphic snps#########################################
lookupnewrowBIN<-function(mat, reference, pos)
{
	ref<-reference[pos]
	alt<-'NA'
	newrow<-as.vector(unlist(c(pos, ref, alt, 'SNP', rep(0, (length(colnames(mat)) - 4)))))
	return(newrow)
}

#add refernce isolate to IUPAC matrix using REF column
addReferenceIUPAC<-function(mat, reference)
{
	out<-paste(mat, reference, sep='.')
	#mat<-readMat(mat)
	name<-reference
	mat<-cbind(mat, name=mat$REF)
	colnames(mat)[ncol(mat)]<-name
	return(mat)
	#write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)

}

#add refernce isolate to bin matrix using REF column
addReferenceBIN<-function(mat, reference)
{
	out<-paste(mat, reference, sep='.')
	#mat<-readMat(mat)
	name<-reference
	mat<-cbind(mat, name=rep(0, nrow(mat)))
	colnames(mat)[ncol(mat)]<-name
	return(mat)
	#write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)
}

#merge matrices with same columns but differing positions together
mergeMats<-function(mat1, mat2, out)
{
	m1<-readMat(mat1)
	m2<-readMat(mat2)
	m<-rbind(m1, m2)
	new<-m[order(m$POS),]
	dim(new)
	new[1:20,1:20]
	write.table(new, out, quote=F, row.name=F, col.name=T)
}

##########################################################################################################
#function to return IUPAC matrix adding new positions#####################################################
unionposmatIUPAC<-function(matfile, reference, additionals_snps)
{
	mat<-readMat(matfile)
	reference<-read.fasta(reference, seqonly=TRUE)
	reference<-strsplit(unlist(reference), '')[[1]]
	unionp<-iterateunionpos(mat, additionals_snps)
	mpos<-mat$POS
	newpos<-unionp[is.na(match(unionp, mpos))]
	newpos<-unique(newpos)
	matrix<-matrix(, nrow=length(newpos), ncol= ncol(mat))
	df<-data.frame(matrix)
	colnames(df)<-colnames(mat)
	for(i in 1:length(newpos))
	{
		newrow<-lookupnewrowIUPAC(mat, reference, as.numeric(newpos[i]))
		df[i,]<-newrow
	}
	newmat<-rbind(mat, df)
	newmat<-newmat[order(as.numeric(newmat$POS)),]
	write.table(newmat, paste(matfile, 'extended.temp', sep='_'), col.name=TRUE, row.name=FALSE, quote=FALSE)
	return(newmat)
}

#function to return binary matrix adding new positions#####################################################
unionposmatBIN<-function(matfile, reference, additionals_snps)
{
	mat<-readMat(matfile)
	reference<-read.fasta(reference, seqonly=TRUE)
	reference<-strsplit(unlist(reference), '')[[1]]
	unionp<-iterateunionpos(mat, additionals_snps)
	mpos<-mat$POS
	newpos<-unionp[is.na(match(unionp, mpos))]
	newpos<-unique(newpos)
	headers<-colnames(mat)
	matrix<-matrix(, nrow=length(newpos), ncol= ncol(mat))
	df<-data.frame(matrix)
	colnames(df)<-colnames(mat)
	for(i in 1:length(newpos))
	{
		newrow<-lookupnewrowBIN(mat, reference, as.numeric(newpos[i]))
		df[i,]<-newrow
	}
	mat<-rbind(mat, df)
	mat<-mat[order(df$POS),]
	colnames(mat)<-headers
	write.table(mat, paste(matfile, 'extended.temp', sep='_'), col.name=TRUE, row.name=FALSE, quote=FALSE)
	return(mat)
}

#add isolates to IUPAC matrix containing all reference positions plus more
addIsolatesIUPAC<-function(extended_mat, reference_file, additionals_snps, add_reference=TRUE)
{
	out<-paste(extended_mat, '.additionals', sep='')
	mat<-readMat(extended_mat)
	reference<-read.fasta(reference_file, seqonly=TRUE)
	reference<-strsplit(unlist(reference), '')[[1]]
	for(i in 1:length(additionals_snps))
	{
		name<-additionals_snps[i]
		newisolate<-as.data.frame(readsnps(additionals_snps[i]))
		#next step removes duplicated positions (presumably insertion positions)
		newisolate<-newisolate[-which(duplicated(newisolate$refpos)==TRUE),]
		positions<-mat$POS
		missing_positions<-positions[(is.na(match(positions, newisolate$refpos)))]
		newpositions<-c()
		IUPACcol<-c()
		for(j in 1:length(missing_positions))
		{
			newpositions[j]<-missing_positions[j]
			IUPACcol[j]<-reference[missing_positions[j]]
		}
		newposdf<-as.data.frame(cbind(refpos=as.vector(newpositions), alt=as.character(IUPACcol)))
		origposdf<-as.data.frame(cbind(refpos=as.vector(newisolate$refpos), alt=as.character(newisolate$alt)))
		newdf<-rbind(newposdf, origposdf)
		newdf<-newdf[order(as.numeric(as.character(newdf$refpos))),]
		mat<-cbind(mat, name=newdf$alt)
		colnames(mat)[ncol(mat)]<-name
	}
	if(add_reference==TRUE)
	{
		mat<-addReferenceIUPAC(mat, reference_file)
	}
	write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)
}

#add isolates to bin matrix containing all reference positions plus more
addIsolatesBIN<-function(extended_mat,reference_file, additionals_snps, add_reference=TRUE)
{
	out<-paste(extended_mat, '.additionals', sep='')
	mat<-readMat(extended_mat)
	reference<-read.fasta(reference_file, seqonly=TRUE)
	reference<-strsplit(unlist(reference), '')[[1]]
	for(i in 1:length(additionals_snps))
	{
		name<-additionals_snps[i]
		newisolate<-as.data.frame(readsnps(additionals_snps[i]))
		positions<-mat$POS
		missing_positions<-positions[is.na(match(positions, newisolate$refpos))]
		newpositions<-c()
		BINcol<-c()
		for(i in 1:length(missing_positions))
		{
			newpositions[i]<-missing_positions[i]
			BINcol[i]<-0
		}
		newposdf<-as.data.frame(cbind(refpos=as.vector(newpositions), alt=as.character(BINcol)))
		origposdf<-as.data.frame(cbind(refpos=as.vector(newisolate$refpos), alt=as.character(newisolate$alt)))
		newdf<-rbind(newposdf, origposdf)
		newdf<-newdf[order(as.numeric(as.character(newdf$refpos))),]
		mat<-cbind(mat, name=newdf$bin)
		colnames(mat)[ncol(mat)]<-name
	}
	if(add_reference==TRUE)
	{
		mat<-addReferenceBIN(mat, reference_file)
	}
	write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)
}


IUPAC2fasta<-function(additionals_mat)
{
	mat<-readMat(additionals_mat)
	mat_snp<-mat[which(mat$TYPE=='SNP'),]
	out<-paste(additionals_mat, '.fasta', sep='')
	for(i in 5:ncol(mat_snp))
	{
		seq<-mat_snp[,i]
		name<-colnames(mat_snp)[i]
		write.fasta(seq,name, out, open='a')
	}
}

FASTA2phylip<-function(fasta_file)
{
	dat<-read.fasta(fasta_file)
	dat2phylip(dat, paste(fasta_file, '.phy',sep=''))
}
	


########################################################################################################################
#### DECIDE AT WHICH POINT TO FILTER OUT INDELS FROM MAT- MAYBE NOT JUST AT LST FASTA GENERATION STEP?
#### FIGURE OUT WHY THERE IS AN X IN THE FASTA FILES
# write shell script to run nucmer to include in repository
