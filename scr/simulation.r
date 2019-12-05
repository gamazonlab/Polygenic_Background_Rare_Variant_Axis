args=as.numeric(commandArgs(TRUE))

#setting
h2_p=c(0.1,0.2,0.3,0.4,0.5)[3]
N_causal_SNPs=c(100,200,300)[2]  #N of causal SNP
vr=c(0.001,0.002,0.003,0.005,0.01)[5]
freq_r=c(0.001,0.005,0.01)[2]
k=c(0.02,0.05,0.10)[2]  #prevalence
N_samples=c(1000,2000,3000,4000,5000)[3]

run_i=1
run_array=list()
for (h2_p in c(0.1,0.2,0.3,0.4,0.5)){
  for (N_causal_SNPs in c(100)){
    for (vr in c(0.01,0.02,0.03,0.05,0.10)){
      for (freq_r in c(0.005,0.01,0.02)){
        for (k in c(0.02,0.05,0.10)){
          for (N_samples in c(1000,2000,3000,5000,10000)){
            run_array[[run_i]]=c(h2_p,N_causal_SNPs,vr,freq_r,k,N_samples)
            run_i=run_i+1
          }
        }
      }
    }
  }
}

run_array<-run_array[[args]]
h2_p=run_array[1]
N_causal_SNPs=run_array[2]
vr=run_array[3]
freq_r=run_array[4]
k=run_array[5]
N_samples=run_array[6]
#e_h2_r=0.2

#mkdir
cmd<-paste0('mkdir /data/coxvgi/zhoud2/projects/prslev/pheno/y_hat/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples)
system(cmd,wait = T)


#extract common and rare variant
#recode genotype to dosage file

#load biovu snp list (ld-pruned)
biovu<-read.table('/data/coxvgi/zhoud2/projects/prslev/geno/test/biovu.bim',stringsAsFactors = F)


#snp
snp<-biovu[c(grep("rs",biovu[,2])),]

#rare variant
rare_df<-read.table(paste0('/data/coxvgi/zhoud2/data/biovu/geno/chr/rare/maf_',freq_r,'.txt'),header = T,stringsAsFactors = F)
rare=biovu[which(biovu$V2 %in% rare_df$SNP),]
freq_snp=sapply(rare$V2,function(x) rare_df[which(rare_df$SNP==x),5])
rare$effect=(vr/2/freq_snp/(1-freq_snp))^0.5
set.seed(2019)
rare<-rare[sample(nrow(rare),100,replace = F),]

#output
output<-data.frame(h2_p=h2_p,N_causal_SNPs=N_causal_SNPs,vr=vr,freq_r=freq_r,k=k,N_samples=N_samples,seed=seq(1,100),median_mut=NA,median_ref=NA,p_wilcox=NA,h2_rare=NA,penetrance=NA,h2_poly_est=NA,var_y=NA,mean_prs=NA,mean_prs_cases=NA)


#sampling snps
i=1
for (i in 1:100){
  #sample causal snps
  set.seed(i)
  df<-snp[sample(seq(1,nrow(snp)),N_causal_SNPs),]
  
  #generate PG beta
  #Beta~(0,(h2-vr)/N)
  set.seed(i)
  df$effect=rnorm(N_causal_SNPs,0,(h2_p/N_causal_SNPs)^0.5)
  #hist(beta_p)
  
  #merge rare and common variants 
  df<-rbind(rare[i,],df)
  df<-df[,c(2,7)]
  write.table(df,paste0('/data/coxvgi/zhoud2/projects/prslev/snplist/tmp/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_','seed_',i,'.txt'),quote = F,row.names = F,col.names = F,sep = '\t')
  
  #dosage file
  cmd=paste0('plink --bfile /data/coxvgi/zhoud2/projects/prslev/geno/test/biovu --extract /data/coxvgi/zhoud2/projects/prslev/snplist/tmp/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_','seed_',i,'.txt --recodeA --out /data/coxvgi/zhoud2/projects/prslev/geno/dosage/tmp/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_','seed_',i)
  system(cmd,wait = T)
  
  #load dosage
  dosage<-read.table(paste0('/data/coxvgi/zhoud2/projects/prslev/geno/dosage/tmp/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_','seed_',i,'.raw'),header = T,stringsAsFactors = F)
  #rm
  cmd=paste0('rm /data/coxvgi/zhoud2/projects/prslev/geno/dosage/tmp/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_','seed_',i,'*; rm /data/coxvgi/zhoud2/projects/prslev/snplist/tmp/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_','seed_',i,'*');system(cmd)
  
  set.seed(i)
  dosage<-dosage[sample(nrow(dosage),nrow(dosage),replace = F),]
  dosage<-dosage[1:N_samples,-c(2,3,4,5,6)]
  
  #
  df[1,1]<-paste0('X',strsplit(df[1,1],"[:]")[[1]][1],'.',strsplit(df[1,1],"[:]")[[1]][2])
  colnames(dosage)[-1]<-sapply(colnames(dosage)[-1], function(x) strsplit(x,"[_]")[[1]][1])
  lev_pos<-which(colnames(dosage)==df[1,1])
  case_pos<-which(dosage[,lev_pos]!=0)
  
  #scale dosage
  dosage[,-c(1,lev_pos)]<-sapply(dosage[,-c(1,lev_pos)], function(x) scale(x))
  
  for (j in 1:nrow(df)){
    dosage[,j+1]=df[which(df[,1]==colnames(dosage)[j+1]),2]*dosage[,j+1]
  }
  
  set.seed(i)
  dosage$y=rowSums(dosage[,-1])+rnorm(N_samples,0,(1-h2_p-vr)^0.5)
  dosage$y=as.numeric(scale(dosage$y))  #scale
  dosage$prs=rowSums(dosage[,-c(1,which(colnames(dosage)==df[1,1]),ncol(dosage))])
  
  dosage<-dosage[,c(1,which(colnames(dosage)==df[1,1]),ncol(dosage)-1,ncol(dosage))]
  dosage$y_binary<-ifelse(dosage$y>quantile(dosage$y,1-k)[[1]],1,0)
  
  cases<-dosage[dosage$y_binary==1,]
  cases[,2]<-ifelse(cases[,2]!=0,1,0)
  
  mut<-dosage[case_pos,]
  
  output[i,8]<-median(cases[which(cases[,2]==1),4])
  output[i,9]<-median(cases[which(cases[,2]==0),4])
  if (length(which(cases[,2]==0))>0 & length(which(cases[,2]==1))>0){
    fit<-wilcox.test(cases[which(cases[,2]==1),4],cases[which(cases[,2]==0),4])
    output[i,10]<-fit$p.value
  }else if(length(which(cases[,2]==1))==0){
    output[i,10]<-1
  }
  
  fit<-summary(lm(dosage$y~dosage$prs))
  r2_poly<-fit$r.squared
  beta_poly<-fit$coefficients[2,1]
  h2_poly_est<-r2_poly/beta_poly^2
  
  output[i,11]<-vr/var(dosage$y)   #h2_rare
  output[i,12]<-sum(mut$y_binary)/nrow(mut)  #penetrance
  output[i,13]<-h2_poly_est
  output[i,14]<-var(dosage$y)
  output[i,15]<-mean(dosage$prs)
  output[i,16]<-mean(cases$prs)
}


write.table(output,paste0('/data/coxvgi/zhoud2/projects/prslev/result/simulation/h2p_',h2_p,'_nsnp_',N_causal_SNPs,'_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples),quote = F,row.names = F,sep='\t')












