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
	alt<-'Na'
	newrow<-as.vector(unlist(c(pos, ref, alt, 'SNP', rep(ref, (length(colnames(mat)) - 4)))))
	return(newrow)
}

#function to generate new binary row for monomorphic snps#########################
lookupnewrowBIN<-function(mat, reference, pos)
{
	ref<-reference[pos]
	alt<-'Na'
	newrow<-as.vector(unlist(c(pos, ref, alt, 'SNP', rep(0, (length(colnames(mat)) - 4)))))
	return(newrow)
}

#function to return IUPAC matrix adding new positions#############################
unionposmatIUPAC<-function(matfile, reference, additionals_snps)
{
	mat<-readMat(matfile)
	reference<-read.fasta(reference)
	reference<-strsplit(unlist(reference), '')[[1]]
	unionp<-iterateunionpos(mat, additionals_snps)
	mpos<-mat$POS
	newpos<-unionp[is.na(match(mpos, unionp))]
	newpos<-unique(newpos)
	temp<-'tempBIN.txt'
	write.table(mat, temp, quote=FALSE, row.name=FALSE, col.name=TRUE)
	for(pos in 1:length(newpos))
	{
		newrow<-lookupnewrowBIN(mat, reference, newpos[pos])
		write.table(newrow, temp, append=TRUE, col.name=FALSE, row.name=FALSE, quote=FALSE)
	}
	mat<-read.delim('tempBIN.txt', header=TRUE)
	mat<-mat[sort(mat$POS),]
	write.table(mat, paste(matfile, 'extended.temp', sep='_'), col.name=TRUE, row.name=FALSE, quote=FALSE)
	return(mat)
}

#function to return binary matrix adding new positions############################
unionposmatBIN<-function(matfile, reference, additionals_snps)
{
	mat<-readMat(matfile)
	reference<-read.fasta(reference)
	reference<-strsplit(unlist(reference), '')[[1]]
	unionp<-iterateunionpos(mat, additionals_snps)
	mpos<-mat$POS
	newpos<-unionp[is.na(match(mpos, unionp))]
	newpos<-unique(newpos)
	temp<-'tempIUPAC.txt'
	write.table(mat, temp, quote=FALSE, row.name=FALSE, col.name=TRUE)
	for(pos in 1:length(newpos))
	{
		newrow<-lookupnewrowIUPAC(mat, reference, newpos[pos])
		write.table(newrow, temp, append=TRUE, col.name=FALSE, row.name=FALSE, quote=FALSE)
	}
	mat<-read.delim('tempIUPAC.txt', header=TRUE)
	mat<-mat[sort(mat$POS),]
	write.table(mat, paste(matfile, 'extended.temp', sep='_'), col.name=TRUE, row.name=FALSE, quote=FALSE)
	return(mat)
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
outmatBIN=unionposmatBIN(matbin, reference, additional_snps)
outmat=additionalsIUPAC(outmatIUPAC, additional_snps)

#####################################################################
matbin='TW20_minDP10.mat.bin'
reference='TW20.fasta'
additional_snps=c('references/BX_refpos_ref_alt.txt','references/ST22_refpos_ref_alt.txt','references/SE_refpos_ref_alt.txt')