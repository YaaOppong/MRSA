##################################################################################
import pandas as pd
import Bio.SeqIO
import argparse

#functiom to read IUPAC matrix and store as data frame############################
def readmat(mat):
	m=pd.read_csv(mat, delimiter=' ')
	return(m)

#function to read nucmer file of additional genome and store snps and their pos###
def readnucmer(nucmersnps):
	m=pd.read_csv(nucmersnps, delimiter='\t')
	m=m.iloc[:,[0,2]]
	return(m)

#function to return union of postions from IUPAC matrix and additional genomes####
def unionpos(mat, nucmersnps):
	matpos=mat.ix[:,'POS'].tolist()
	nucmersnpspos=nucmersnps.iloc[:,0].tolist()
	union=matpos+nucmersnpspos
	union=sorted(set(union))
	return(union)

#iterate unionpos over multiple additional genomes (snps only)####################
def iterateunionpos(mat, additional_snps):
	union=[]
	for additional_snp in additional_snps:
		additional_snp=readnucmer(additional_snp)
		union=union+unionpos(mat, additional_snp)
	union=sorted(set(union))
	return(union)

#function to generate new IUPAC row for monomorphic snps##########################
def lookupnewrowIUPAC(mat, reference, pos):
	ref=reference.seq[pos]
	alt='NaN'
	newrow=[pos, ref, alt, 'SNP'] + [ref] * (len(mat.columns) - 4)
	newrow=pd.DataFrame(newrow)
	newrow=newrow.transpose()
	return(newrow)

#function to generate new binary row for monomorphic snps#########################
def lookupnewrowBIN(mat, reference, pos):
	ref=reference.seq[pos]
	alt='NaN'
	newrow=[pos, ref, alt, 'SNP'] + [0] * (len(mat.columns) - 4)
	newrow=pd.DataFrame(newrow)
	newrow=newrow.transpose()
	return(newrow)

#function to return IUPAC matrix adding new positions#############################
def unionposmatIUPAC(mat, reference, additionals_snps):
	mat=readmat(mat)
	reference=Bio.SeqIO.read(reference, 'fasta')
	unionp=iterateunionpos(mat, additionals_snps)
	mpos=mat['POS']
	newpos=set(unionp).difference(set(mpos))
	for pos in newpos:
		newrow=lookupnewrowIUPAC(mat, reference, pos)
		mat=pd.concat([mat,newrow])
	mat=mat.sort(mat['POS'])
	return(mat)

#function to return binary matrix adding new positions############################
def unionposmatBIN(mat, reference, additionals_snps):
	mat=readmat(mat)
	reference=Bio.SeqIO.read(reference, 'fasta')
	unionp=iterateunionpos(mat, additionals_snps)
	mpos=mat['POS']
	newpos=set(unionp).difference(set(mpos))
	temp=open('temp.txt', 'w')
	mat.to_csv(temp, mode='w', index=False)
	for pos in newpos:
		newrow=lookupnewrowBIN(mat, reference, pos)
		newrow.to_csv(temp, mode='a', header=False, index=False)
	temp.close()
	mat=pd.read_csv('temp.txt', low_memory=False)
	mat=mat.sort_values(by=['POS'])
	return(mat)

#functions to add snps to iupac and add new IUPAC to matrix as new column##########
def additionalsIUPAC(mat, additional_snps):
	for additional_snp in additional_snps:
		additional_snp=readnucmer(additional_snp)
		unionp=unionpos(mat, additional_snp)
		npos=additional_snp.iloc[:,0]
		newpos=unionp.difference(npos)
		for pos in newpos:
			ref=mat[pos, 'REF']
			newrow=[pos, ref ]
			newdata=newdata.append(newrow)
		additional_snp=additional_snp.concat(pd.DataFrame(newdata))
		additional_snp=additional_snp.sort(additional_snp.iloc[:,0])
		mat=mat.concat(pd.DataFrame(additional_snp.iloc[:,1]), axis=1)
	return(mat)


#functions to add snps to new iupac and add new IUPAC to matrix as new column######
def additionalsBIN(mat, additional_snps):
	for additional_snp in additional_snps:
		additional_snp=readnucmer(additional_snp)
		mpos=mat['POS']
		npos=additional_snp.iloc[:,0]
		newpos=set(mpos).difference(set(npos))
		newcol=[]
		matdict=mat.loc[:,['POS','REF']].to_dict()
		for pos in additional_snp:
			posindex=which
			ref=mat.loc[pos, 'REF']
			if additional_snp[pos,1]==ref:
				newcol=newcol.append(0)
			else:
				newcol=newcol.append(1)
		additional_snp=additional_snp.concat(pd.DataFrame(newcol), axis=1)
		for pos in newpos:
			ref=mat[pos, 'REF']
			newrow=[pos,ref,0]
			newdata=newdata.append(newrow)
		additional_snp=additional_snp.concat(pd.DataFrame(newdata))
		additional_snp=additional_snp.sort(additional_snp.iloc[:,0])
		mat=mat.concat(pd.DataFrame(additional_snp.iloc[:,2]), axis=1)
	return(mat)

#functions to add additional references to BINARY matrix############################

#generate multifasta from IUPAC matrix##############################################
def matfasta(mat, out):
	outfasta=open(out, 'w')
	for col in mat.columns()[3:]:
		a=mat[col]
		Bio.SeqIO.write(fasta, outfasta, "fasta")
	outfasta.close()


###################################################################################
###arg input##
parser=argparse.ArgumentParser(description='')
parser.add_argument('--out_type',help='output type required; f (fasta), b \
	(matbin), i (matiupac)')
parser.add_argument('--mat', help="Variant matrix in IUPAC format (with 4 d\
	escriptor columns;'POS' 'REF' 'AL' 'TYPE')")
parser.add_argument('--matbin', help="Binary variant matrix (with 4 descriptor \
	columns;'POS' 'REF' 'AL' 'TYPE')")
parser.add_argument('--out_prefix', help='Prefix for output files')
parser.add_argument('--additional_snps', help='List of SNP only nucmer files of \
	additional genomes generated using show-snps', nargs='+')
parser.add_argument('--reference', help='Reference fasta file')
args=parser.parse_args()

##pipeline execution###############################################################
if args.out_type=='i':
	outmat=unionposmatIUPAC(args.mat, args.reference, args.additional_snps)
	a.write(outmat)
	a.close()


if args.out_type=='b':
	outmatbin=unionposmatBIN(args.matbin, args.reference, args.additional_snps)
	b=open(args.out_prefix + 'mat.bin', 'w')
	b.write(outmatbin)
	b.close()


if args.out_type=='f':
	outmat=unionposmatIUPAC(args.mat, args.reference, args.additional_snps)
	a.write(outmat)
	a.close()
	matfasta(args.outmat, args.out_prefix + '.fasta')

##################################################################################
##################################################################################
