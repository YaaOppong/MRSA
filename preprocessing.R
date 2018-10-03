getmat2fastainput<-function(nucmersnps, out_prefix)
{
	snps<-read.delim(nucmersnps, header=F)
	snps2<-as.data.frame(cbind(refpos=as.character(snps[,1]), ref=as.character(snps[,2]), alt=as.character(snps[,3])))
	#deal with duplicate positions?
	write.table(snps2, paste(out_prefix, '_refpos_ref_alt.txt', sep=''), quote=F, row.name=F, col.name=T)
}