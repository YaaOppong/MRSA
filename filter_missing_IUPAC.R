samples_to_remove<-read.table('removed_samples_0.1.txt')
positions_to_remove<-read.table('removed_samples_0.1.txt')


#filter by missingness on sample
columns<-match(samples_to_remove, colnames(mat))
columns<-na.omit(columns)
mat3<-mat2[,-columns]

#filter by missingness on variant
rows<-match(mat$POS, positions_to_remove)
rows<-na.omit(rows)
mat2<-mat[-rows,]


write.table(mat3, '', col.names=T, row.names=F, quote=F)