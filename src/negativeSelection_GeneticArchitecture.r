setwd("prslev/result/")

source('d://tools/qqunif.r')

################# liability threshold

a = read.table('probit_a/h2p_0.1_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)
b = read.table('probit_b/h2p_0.1_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)
c = read.table('probit_c/h2p_0.1_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)


qqunif(a$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), cex.axis=2.0, cex.lab=2.0, cex=1.8)
par(new=T)
qqunif(b$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), col="red", cex.axis=2.0, cex.lab=2.0, cex=1.8)
par(new=T)
qqunif(c$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), col="blue", cex.axis=2.0, cex.lab=2.0, cex=1.8)


legend('topleft', c("Negative selection","LDAK","Polygenic"),
             col=c('red','blue','black'),lty=c(1,1,1), cex=1.8)



################# logit

a = read.table('probit_a/h2p_0.2_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)
b = read.table('probit_b/h2p_0.2_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)
c = read.table('probit_c/h2p_0.2_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)


qqunif(a$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), cex.axis=2.0, cex.lab=2.0)
par(new=T)
qqunif(b$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), col="red", cex.axis=2.0, cex.lab=2.0)
par(new=T)
qqunif(c$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), col="blue", cex.axis=2.0, cex.lab=2.0)


legend('topleft', c("Negative selection","LDAK","Polygenic"),
             col=c('red','blue','black'),lty=c(1,1,1), cex=1.8)


################# logit vs probit, negative selection

m = read.table('probit_b/h2p_0.1_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)
n = read.table('logit_b/h2p_0.1_nsnp_all_snps_vr_0.01_freqr_0.001_k_0.01_Nsample_10000_pi0_0',header=T)

qqunif(m$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), col="red", cex.axis=2.0, cex.lab=2.0)
par(new=T)
qqunif(n$case_inverse_wilcox_p, BH=F, CI=F, BF=F, xlim=c(0,3), ylim=c(0,7), col="blue", cex.axis=2.0, cex.lab=2.0)


legend('topleft', c("Liability threshold","Logit"),
             col=c('red','blue'),lty=c(1,1,1), cex=1.8)


