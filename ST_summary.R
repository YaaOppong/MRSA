a<-read.delim('merged_mlsts.txt', header=T)
jpeg('ST_plot.jpeg', width = 2000, height = 480)
plot(table(a$ST), xlab='ST Type', ylab='Frequency', xaxt='n')
par(las=2)
axis(1,at=1:length(levels(a$ST)), labels=levels(a$ST), cex.axis=0.9)
dev.off()
