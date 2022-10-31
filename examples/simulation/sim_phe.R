# sim_phe.R

args <- commandArgs(trailingOnly=TRUE)
ii = args[1]
hsq = as.numeric(args[2])
print(paste("h2=",hsq,sep=""))

gv = read.table(paste(ii,".gv.csv",sep=""),sep=",")
n = nrow(gv)
vg = var(gv[,2])
ve = vg*(1-hsq)/hsq
print(paste("VarG=",vg,sep=""))

err = rnorm(n) * sqrt(ve)
phe = matrix(nrow=n, ncol=3)
phe[,1] = 0
phe[,2] = gv[,1]
phe[,3] = gv[,2] + err
print(paste("VarP=",var(phe[,3]),sep=""))

colnames(phe) = c("FID","IID","QT")
write.table(phe,paste(ii,".bolt.txt",sep=""),row.names = FALSE, quote=FALSE)
phe = phe[,2:ncol(phe)]
write.csv(phe,paste(ii,".slemm.csv",sep=""),row.names = FALSE, quote=FALSE)
