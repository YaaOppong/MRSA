################################################################################
import pandas as pd
import Bio.SeqIO
import argparse

#functiom to read IUPAC matrix and store as data frame##########################
def readmat(mat):
	m=pd.read_csv(mat, delimiter=' ')
	return(m)

#function to read nucmer file of additional genome and store snps and their pos#
def readnucmer(nucmersnps):
	m=pd.read_csv(nucmersnps, delimiter='\t')
	m=m.iloc[:,[0,2]]
	return(m)

#function to return union of postions from IUPAC matrix and additional genomes##
def unionpos(mat, nucmersnps):
	matpos=mat.ix[:,'POS'].to_dict()
	nucmersnpspos=nucmersnps.iloc[:,0].to_dict()
	union=set().union(matpos,nucmersnpspos)
	union=sorted(union)
	return(union)

##iterate unionpos over multiple additional genomes (snps only)#################
#def iterateunionpos(mat, additional_snps):
#	for additional_snp in additional_snps:
#		additional_snp=readnucmer(additional_snp)
#		union=unionpos(mat, additional_snp)
#	return(union)

#function to generate new IUPAC row for monomorphic snps########################
def lookupnewrowIUPAC(mat, reference, pos):
	#reference=SeqIO.read(reference, 'fasta')
	ref=reference.seq[pos]
	alt='NaN'
	newrow=[pos, ref, alt, 'SNP'] + [ref] * (len(mat.columns) - 4)
	newrow=pd.DataFrame(newrow)
	newrow=newrow.transpose()
	return(newrow)

#function to generate new binary row for monomorphic snps#######################
def lookupnewrowBIN(mat, reference, pos):
	reference=SeqIO.read(reference, 'fasta')
	ref=reference.seq[pos]
	alt='NaN'
	newrow=[pos, ref, alt, 'SNP'] + [0] * (len(mat.columns) - 4)
	newrow=pd.DataFrame(newrow)
	newrow=newrow.transpose()
	return(newrow)

#function to return IUPAC matrix adding new positions###########################
def unionposmatIUPAC(mat, reference, additional_snps):
	mat=readmat(mat)
	reference=SeqIO.read(reference, 'fasta')
	#union=iterateunionpos(mat, additionals_snps)
	additionals_snps=readnucmer(additionals_snps)
	unionp=unionpos(mat, additionals_snps)
	mpos=mat['POS']
	newpos=set(unionp).difference(set(mpos))
	for pos in newpos:
		newrow=lookupnewrowIUPAC(mat, reference, pos)
		mat=mat.append(newrow)
	mat=mat.sort(mat['POS'])
	return(mat)

#function to return binary matrix adding new positions##########################
def unionposmatBIN(mat, reference, additionals_snps):
	mat=readmat(mat)
	reference=SeqIO.read(reference, 'fasta')
	#union=iterateunionpos(mat, additionals_snps)
	additionals_snps=readnucmer(additionals_snps)
	unionp=unionpos(mat, additionals_snps)
	mpos=mat['POS']
	newpos=set(unionp).difference(set(mpos))
	for pos in newpos:
		newrow=lookupnewrowBIN(mat, reference, pos)
		mat=mat.append(newrow)
	mat=mat.sort(mat['POS'])
	return(mat)

#functions to add additional references to IUPAC matrix#########################
def addadditionalsIUPAC(mat, additional_snps):
	for additional_snp in additional_snps:
		additional_snp=readnucmer(additional_snp)
		union=unionpos(mat, additional_snp)
		npos=additional_snp[0,]
		newpos=union.difference(npos)
		for pos in newpos:
			ref=mat[pos, 'REF']
			newrow=pd.DataFrame([pos, ref ])	
			additional_snp=additional_snp.append(newrow)
		additional_snp=additional_snp.sort(additional_snp[,0])
		mat=mat.appendcol(additional_snp[,1])
	return(mat)

#functions to add additional references to BINARY matrix########################

#generate multifasta from IUPAC matrix##########################################
def matfasta(mat, out):
	outfasta=open(out, 'w')
	for col in cols(mat)[3:]:
		a=as.fasta(mat[col])
		SeqIO.write(fasta, outfasta, "fasta")
	outfasta.close()

################################################################################
###arg input##
parser=argparse.ArgumentParser(description='')
parser.add_argument('out_type',help='output type required; -f (fasta) -b \
	(matbin) -i (matiupac)')
parser.add_argument('mat', help="Variant matrix in IUPAC format (with 4 d\
	escriptor columns;'POS' 'REF' 'AL' 'TYPE')")
parser.add_argument('matbin', help="Binary variant matrix (with 4 descriptor \
	columns;'POS' 'REF' 'AL' 'TYPE')")
parser.add_argument('additional_snps', help='List of SNP only nucmer files of \
	additional genomes generated using show-snps')
parser.add_argument('out_prefix', help='Prefix for output files')
args=parser.parse_args()

## pipeline execution###########################################################
outmat=unionposmat(args.mat, args.additional_snps)
outmatbin=unionposmatbin(args.matbin, args.additional_snps)
matfasta(args.outmat, outfasta)

a=open(args.out_prefix + 'mat', 'w')
b=open(args.out_prefix + 'mat.bin', 'w')

a.write(outmat)
b.write(outmatbin)

a.close()
b.close()
###############################################################################
###############################################################################
