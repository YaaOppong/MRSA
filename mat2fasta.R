##################################################################################
library(data.table)
library(seqinr)


#functiom to read IUPAC matrix and store as data frame############################
readMat<-function(mat)
{
	mat<-as.data.frame(fread(mat, header=T))
	return(mat)
}

#function to read file of additional genome- refpos, ref, alt#####################

readsnps<-function(snps)
{
	m<-read.table(snps, header=T)
	m<-cbind(refpos=m$refpos,alt=as.character(m$alt))
	return(m)
}
#function to return union of postions from IUPAC matrix and additional genomes####

unionpos<-function(mat, nucmersnps)
{
	matpos<-mat$POS
	nucmersnpspos=nucmersnps[,1]
	union=sort(unique(c(matpos,nucmersnpspos)))
	return(union)
}
#iterate unionpos over multiple additional genomes (snps only)####################

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
#function to generate new IUPAC row for monomorphic snps##########################

lookupnewrowIUPAC<-function(mat, reference, pos)
{
	ref<-reference[pos]
	alt<-'NA'
	newrow<-as.vector(unlist(c(pos, ref, alt, 'SNP', rep(ref, (length(colnames(mat)) - 4)))))
	return(newrow)
}

#function to generate new binary row for monomorphic snps#########################
lookupnewrowBIN<-function(mat, reference, pos)
{
	ref<-reference[pos]
	alt<-'NA'
	newrow<-as.vector(unlist(c(pos, ref, alt, 'SNP', rep(0, (length(colnames(mat)) - 4)))))
	return(newrow)
}

#function to return IUPAC matrix adding new positions#############################
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
		#write.table(newrow, temp, append=TRUE, col.name=FALSE, row.name=FALSE, quote=FALSE)
	}
	mat<-df[order(df$POS),]
	write.table(mat, paste(matfile, 'extended.temp', sep='_'), col.name=TRUE, row.name=FALSE, quote=FALSE)
	return(mat)
}

#function to return binary matrix adding new positions############################
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
	mat<-df[order(df$POS),]
	colnames(mat)<-headers
	write.table(mat, paste(matfile, 'extended.temp', sep='_'), col.name=TRUE, row.name=FALSE, quote=FALSE)
	return(mat)
}

###following requires debugging
addIsolatesIUPAC<-function(extended_mat, reference, additionals_snps)
{
	out<-paste(extended_mat, '.additionals', sep='')
	mat<-readMat(extended_mat)
	positions<-mat$POS
	reference<-read.fasta(reference, seqonly=TRUE)
	reference<-strsplit(unlist(reference), '')[[1]]
	for(i in 1:length(additionals_snps))
	{
		name<-additionals_snps[i]
		newisolate<-readsnps(additionals_snps[i])
		missing_postions<-positions[is.na(match(positions, newisolate$refpos))]
		newpositions<-c()
		for(i in 1:length(missing_positions))
		{
			newpositions[i]<-missing_positions[i]
			IUPACcol[i]<-reference[missing_positions[i]]
		}
		newposdf<-cbind(refpos=newpositions, alt=IUPACcol)
		origposdf<-cbind(newisolate$refpos, newisolate$alt)
		newdf<-rbind(newposdf, origposdf)
		newdf<-newdf[order(newdf$refpos),]
		mat<-cbind(mat, name=newdf$alt)
		colname(mat$name)<-name
		write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)
	}
}

addIsolatesBIN<-function(extended_mat, additionals_snps)
{
	out<-paste(extended_mat, '.additionals', sep='')
	mat<-readMat(extended_mat)
	positions<-mat$POS
	reference<-read.fasta(reference, seqonly=TRUE)
	reference<-strsplit(unlist(reference), '')[[1]]
	for(i in 1:length(additionals_snps))
	{
		name<-additionals_snps[i]
		newisolate<-readsnps(additionals_snps[i])
		missing_postions<-positions[is.na(match(positions, newisolate$refpos))]
		newpositions<-c()
		for(i in 1:length(missing_positions))
		{
			newpositions[i]<-missing_positions[i]
			BINcol[i]<-0
		}
		newposdf<-cbind(refpos=newpositions, bin=BINcol)
		origposdf<-cbind(newisolate$refpos, bin=rep(1,nrow(newisolate)))
		newdf<-rbind(newposdf, origposdf)
		newdf<-newdf[order(newdf$refpos),]
		mat<-cbind(mat, name=newdf$bin)
		colname(mat$name)<-name
		write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)
	}
}

addReferenceIUPAC<-function(mat, reference)
{
	out<-paste(mat, reference, sep='.')
	mat<-readMat(extended_mat)
	name<-reference
	mat<-cbind(mat, name=mat$REF)

}

addReferenceBIN<-function(mat, reference)
{
	out<-paste(mat, reference, sep='.')
	mat<-readMat(extended_mat)
	name<-reference
	mat<-cbind(mat, name=rep(0, nrow(mat)))
	colname(mat$name)<-name
	write.table(mat, out, col.name=TRUE, row.name=FALSE, quote=FALSE)
}

#####BELOW REQUIRES DEBUGGING
#functions to add snps to iupac and add new IUPAC to matrix as new column##########
#def additionalsIUPAC(mat, additional_snps):
#	for additional_snp in additional_snps:
#		name=additional_snp
#		additional_snp=readsnps(additional_snp)
#		mpos=mat['POS']
#		npos=additional_snp.iloc[:,0]
#		newpos=set(mpos).difference(set(npos))
#		mydict={}
#		for i in  range(len(mat['POS'])):
#			mydict[mat.loc[i, 'POS']]= mat.loc[i, 'REF']
#		newdata=[]
#		for pos in newpos:
#			ref=mydict[pos]
#			newrow=[pos,ref]
#			newdata.append(newrow)
#		newdata=pd.DataFrame(newdata)
#		newdata.columns=['refpos','alt']
#		additional_snp=pd.concat([additional_snp, newdata], axis=0)
#		additional_snp=additional_snp.sort_values(by=['refpos'])
#		newmatcol=pd.DataFrame(additional_snp.iloc[:,1])
#		newmatcol.columns=[name]
#		#check mat['POS']==additional_SNPS['refpos']### STILL TO DO###
#		mat=pd.concat([mat, newmatcol], axis=1)
#	return(mat)


#functions to add additional references to BINARY matrix############################
#def additionalsBIN(mat, additional_snps):
#	for additional_snp in additional_snps:
#		name=additional_snp
#		additional_snp=readsnps(additional_snp)
#		mpos=mat['POS']
#		npos=additional_snp.iloc[:,0]
#		newpos=set(mpos).difference(set(npos))
#		mydict={}
#		for i in  range(len(mat['POS'])):
#			mydict[mat.loc[i, 'POS']]= mat.loc[i, 'REF']
#		newcol=[]
#		for pos in additional_snp.iloc[:,0]:
#			newcol.append(1)
#		additional_snp=pd.concat([additional_snp, pd.DataFrame(newcol)], axis=1)
#		additional_snp.columns=['refpos','alt', 'bin' ]
#		newdata=[]
#		for pos in newpos:
#			ref=mydict[pos]
#			newrow=[pos,ref,0]
#			newdata.append(newrow)
#		newdata=pd.DataFrame(newdata)
#		newdata.columns=['refpos','alt', 'bin' ]
#		additional_snp=pd.concat([additional_snp, newdata], axis=0)
#		additional_snp=additional_snp.sort_values(by=['refpos'])
#		newmatcol=pd.DataFrame(additional_snp.iloc[:,2])
#		newmatcol.columns=[name]
#		#check mat['POS']==additional_SNPS['refpos']### STILL TO DO####
#		mat=pd.concat([mat, newmatcol], ignore_index=True, axis=1)
#	return(mat)



#def additionalsIUPAC(mat, additional_snps):
#	for additional_snp in additional_snps:
#		name=additional_snp
#		additional_snp=readsnps(additional_snp)
#		for pos in mat.loc['POS']:
#			newcol.iloc[pos]=additional_snp.iloc[np.where(additional_snp.loc['REFPOS'] == pos, 'ALT')]
#			except
#				np.where(additional_snp.loc['REFPOS'] == pos)= NA:
#					newcol.iloc[pos]=mat.loc[pos, 'REF']

#	return(mat)
#generate multifasta from IUPAC matrix##############################################
#def matfasta(mat, out):
#	outfasta=open(out, 'w')
#	for col in mat.columns()[3:]:
#		a=mat[col]
#		Bio.SeqIO.write(fasta, outfasta, "fasta")
#	outfasta.close()


###################################################################################
###arg input##
#parser=argparse.ArgumentParser(description='')
#parser.add_argument('--out_type',help='output type required; f (fasta), b \
#	(matbin), i (matiupac)')
#parser.add_argument('--mat', help="Variant matrix in IUPAC format (with 4 d\
#	escriptor columns;'POS' 'REF' 'AL' 'TYPE')")
#parser.add_argument('--matbin', help="Binary variant matrix (with 4 descriptor \
#	columns;'POS' 'REF' 'AL' 'TYPE')")
#parser.add_argument('--out_prefix', help='Prefix for output files')
#parser.add_argument('--additional_snps', help='List of SNP only nucmer files of \
#	additional genomes generated using show-snps', nargs='+')
#parser.add_argument('--reference', help='Reference fasta file')
#args=parser.parse_args()

##pipeline execution###############################################################
#if args.out_type=='i':
#	outmatIUPAC=unionposmatIUPAC(args.mat, args.reference, args.additional_snps)
#	outmat=additionalsIUPAC(outmatIUPAC, additional_snps)
#	a.write(outmat)
#	a.close()m


#if args.out_type=='b':
#	outmatBIN=unionposmatBIN(args.matbin, args.reference, args.additional_snps)
#	outmat=additionalsBIN(outmatBIN, additional_snps)
#	b=open(args.out_prefix + 'mat.bin', 'w')
#	b.write(outmat)
#	b.close()


#if args.out_type=='f':
#	outmat=unionposmatIUPAC(args.mat, args.reference, args.additional_snps)
#	a.write(outmat)
#	a.close()
#	matfasta(args.outmat, args.out_prefix + '.fasta')

##################################################################################
##############	outmatIUPAC=unionposmatIUPAC(args.mat, args.reference, args.additional_snps)
#outmatBIN=unionposmatBIN(matbin, reference, additional_snps)
#outmat=additionalsIUPAC(outmatIUPAC, additional_snps)

#####################################################################
