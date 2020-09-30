
#comparison of cSEV-OR between LEV carriers and non-carriers 

#disease risk model setting
#LTM (liability threshold model)
#logit model

#genetric architecture setting
# a=Infinitesimal architecture
# b=with negative selection
# c=LD-adjusted kinship

args=commandArgs(TRUE)  #subjob id

subjob_id=as.numeric(args[1]) #subjob id 

N_simulation=500  #simulation times

#h2_p=0.1;N_causal_SNPs=95593;vr=0.01;freq_r=0.001;k=0.01;genetic_model=1;pi0=0


#subjobs 
run_i=1
run_array=list()
for (h2_p in c(0.1,0.2,0.3,0.4,0.5)){ #h2 poly
  for (N_causal_SNPs in c(95593)){ #N of causal SNPs
    #    for (vr in c(0.01,0.02,0.03,0.05,0.10)){ #h2 rare (var rare)
    for (vr in c(0.01,0.02,0.03,0.05,0.10)){ #h2 rare (var rare)
      for (freq_r in c(0.001,0.005,0.01)){ #freq rare
        for (k in c(0.001,0.005,0.01,0.02,0.05)){ #prevalence
          for (genetic_model in c(1,2,3)){ #genetic architecture
            for (pi0 in c(0,0.1,0.3,0.5,0.7,0.9)[c(1,2,3,4,5,6)]){  #pi0, proportion of causal variants pi0=0: utility; pi0!=0: power
              run_array[[run_i]]=c(h2_p,N_causal_SNPs,vr,freq_r,k,genetic_model,pi0)
              run_i=run_i+1
            }
          }
        }
      }
    }
  }
}

run_array<-run_array[[subjob_id]]
h2_p=run_array[1]
#N_causal_SNPs=run_array[2]
vr=run_array[3]
freq_r=run_array[4]
k=run_array[5]
N_samples=run_array[6]
genetic_model=switch(run_array[6],'a','b','c')
# a=Infinitesimal architecture
# b=with negative selection
# c=LD-adjusted kinship
pi0=run_array[7]

#load LDsc
ldsc<-read.table('/data/c***/z***/projects/prslev/geno/ldsc/ldsc.simulation_all_snps',header = T,stringsAsFactors = F)

#load common snp table
common_all<-read.table('/data/c***/z***/projects/prslev/common/common_all.bim',header = F,stringsAsFactors = F)
common_all<-data.frame(rsid=common_all$V2,counted_allele=common_all$V5,effect=NA,stringsAsFactors = F)

#load MAF for all the common variants (used for maf-dependent genetic architectures)
maf_all<-read.table('/data/c***/z***/projects/prslev/geno/geno_copies/common_all.frq',header = T,stringsAsFactors = F)

#n of causal common variants
N_causal_SNPs=nrow(common_all)

#load rare snps
rare_dosage<-readRDS('/data/c***/z***/data/biovu/geno/23k/prslev_simulation_rare.rds')
colnames(rare_dosage)[-1]<-sapply(colnames(rare_dosage)[-1], function(x) strsplit(x,"[_]")[[1]][1])
rare_all<-data.frame(rsid=colnames(rare_dosage)[-1],effect=NA,stringsAsFactors = F)

#load rare variant info (freq)
rare_info<-read.table(paste0('/data/c***/z***/data/biovu/geno/chr/rare/rare.list'),header = T,stringsAsFactors = F)
rare_info$SNP=paste0('X',sub(':','.',rare_info$SNP)) #match variants name with dosage data
rare_info<-rare_info[rare_info$maf==freq_r,]

rare=rare_all[which(rare_all$rsid %in% rare_info$SNP),]

#allele frequency for rare variants (the actual number of freq, close to the setting)
freq_snp=sapply(rare$rsid,function(x) rare_info[which(rare_info$SNP==x),'maf'])
#calculate the effect size for rare vatiants
rare$effect=(vr/2/freq_snp/(1-freq_snp))^0.5
#sample N rare variants for simulation (N= number of simulation times)
set.seed(2019)
rare<-rare[sample(nrow(rare),N_simulation,replace = T),]


#simulation
N_samples=10000;i=1 #for testing

for (N_samples in c(10000)){  
  print(N_samples)
  
  #output dataframe (liability-threshold model (LTM), probit)
  output_p<-data.frame(h2_p=h2_p,N_causal_SNPs=N_causal_SNPs,vr=vr,freq_r=freq_r,N_samples=N_samples,k=k,pi0=pi0,seed=seq(1,N_simulation),model=genetic_model)
  
  #output dataframe (logit)
  output_l<-data.frame(h2_p=h2_p,N_causal_SNPs=N_causal_SNPs,vr=vr,freq_r=freq_r,N_samples=N_samples,k=k,pi0=pi0,seed=seq(1,N_simulation),model=genetic_model)
  
  for (i in 1:N_simulation){
    print(i)
    
    #---generate effect sizes for common variants according to different genetic architectures---
    df<-common_all
    # a Infinitesimal architecture
    #   beta_poly~(0,h2_p/N_caucal_SNPs)
    if(genetic_model=='a'){ 
      set.seed(i)
      df$effect=(1-pi0)*rnorm(N_causal_SNPs,0,(h2_p/N_causal_SNPs)^0.5)
      
      # b Negative selection
      #   beta_poly~N(O,k_constant([f(1-f)]^(1+alpha)))
    }else if(genetic_model=='b'){ 
      set.seed(i)
      df$effect=df[,'effect']<-sapply(maf_all$MAF,function(f) rnorm(1,mean=0,sd=((f*(1-f))^0.63)^0.5))
      #find k constant
      k_constant=(h2_p/N_causal_SNPs/var(df[,'effect']))
      #multiply by k constant
      df[,'effect']<-df[,'effect']*k_constant^0.5*(1-pi0)
      
      # c LD-adjusted kinship
      #   beta_poly~N(O,k_constant([f(1-f)]^(1+alpha)*(1/(1+ldsc))))
    }else if(genetic_model=='c'){ 
      set.seed(i)
      df$effect=df[,'effect']<-mapply(function(f,ld) rnorm(1,mean=0,sd=(((f*(1-f))^0.75)*1/(1+ld))^0.5), maf_all$MAF, ldsc$ldscore)
      #find k constant
      k_constant=(h2_p/N_causal_SNPs/var(df[,'effect']))
      #multiply by k constant
      df[,'effect']<-df[,'effect']*k_constant^0.5*(1-pi0)
    }
    
    #write weight file
    write.table(df,paste0('/data/g***/z***/prslev/tmp/',subjob_id,'.weight'),row.names = F,col.names = F,sep='\t',quote = F)
    
    #use plink to calculate the PB (PRS of common variants)
    #if(i %in% c(1,26,51,76,101,201,301,401)){
    cmd=paste0('/home/zhoud2/tools/plink2/plink2 --bfile /data/c***/z***/projects/prslev/common/common_all --score /data/g***/z***/prslev/tmp/',subjob_id,'.weight variance-standardize  --out /data/g***/z***/prslev/tmp/',subjob_id,'.score')
    system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
    #}
    
    #load PB results (A)
    PB=read.table(paste0('/data/g***/z***/prslev/tmp/',subjob_id,'.score.sscore'),stringsAsFactors = F)
    PB$IID=PB$V1;PB$A=PB$V5*2*nrow(common_all)
    PB<-PB[,c('IID','A')]
    
    #LEV (R)
    rare_id=rare[i,'rsid']
    rare_effect=rare[i,'effect']
    R=as.numeric(rare_dosage[,rare_id]*rare_effect)
    LEV<-data.frame(IID=rare_dosage$FID,R=R)
    
    #merge A and R
    combined<-merge(PB,LEV,by='IID')
    
    #---sampling subjects---
    set.seed(i)
    combined<-combined[sample(seq(1,nrow(combined)),N_samples,replace = F),]
    
    #x=A+R
    combined$x=combined$A+combined$R
    
    #h2_e
    h2_e=1-h2_p-vr
    
    #probability of disease (same with LTM)
    set.seed(i)
    combined$L=combined$x+rnorm(nrow(combined),0,h2_e^0.5)
    t=quantile(combined$L,1-k)[[1]]
    combined$prob=pnorm(-(t-combined$x)/h2_e^0.5)
    
    #------logit model---start----
    set.seed(i)
    combined$y_binary<-sapply(combined$prob,function(x) rbinom(1,1,x))
    combined$R<-ifelse(combined$R!=0,1,0)
    
    #OR in carrier and non-carrier
    carrier=combined[combined$R==1,]
    noncarrier=combined[combined$R==0,]
    #scale A
    carrier$A=scale(carrier$A)
    noncarrier$A=scale(noncarrier$A)
    
    if(length(which(carrier$y_binary==1))>=3 & length(which(noncarrier$y_binary==1))>=3 & length(which(carrier$y_binary==0))>=3 & length(which(noncarrier$y_binary==0))>=3){
      ans_carr<-summary(glm(carrier$y_binary~carrier$A,family = 'binomial'))
      ans_noncarr<-summary(glm(noncarrier$y_binary~noncarrier$A,family = 'binomial'))
      
      output_l[i,'carrier_PRS_OR']<-2.718^ans_carr$coefficients['carrier$A','Estimate']
      #output_l[i,'carrier_PRS_se']<-ans_carr$coefficients['carrier$A','Std. Error']
      output_l[i,'noncarrier_PRS_OR']<-2.718^ans_noncarr$coefficients['noncarrier$A','Estimate']
      #output_l[i,'noncarrier_PRS_se']<-ans_noncarr$coefficients['noncarrier$A','Std. Error']
      
      #A binary
      #carrier$A_binary=ifelse(carrier$A>median(carrier$A,na.rm = T),1,0)
      #noncarrier$A_binary=ifelse(noncarrier$A>median(noncarrier$A,na.rm = T),1,0)
      #ans_carr<-summary(glm(carrier$y_binary~carrier$A_binary,family = 'binomial'))
      #ans_noncarr<-summary(glm(noncarrier$y_binary~noncarrier$A_binary,family = 'binomial'))
      
      #output_l[i,'carrier_bi_PRS_OR']<-2.718^ans_carr$coefficients['carrier$A_binary','Estimate']
      #output_l[i,'noncarrier_bi_PRS_OR']<-2.718^ans_noncarr$coefficients['noncarrier$A_binary','Estimate']
      
    }
    #------logit model---end----
    
    
    #------liability-threshold model (probit) ---start----
    combined$y_binary<-ifelse(combined$L>t,1,0)
    table_check<-table(combined$y_binary,combined$R)
    
    #OR in carrier and non-carrier
    carrier=combined[combined$R==1,]
    noncarrier=combined[combined$R==0,]
    #scale A
    carrier$A=scale(carrier$A)
    noncarrier$A=scale(noncarrier$A)
    
    
    if(length(which(carrier$y_binary==1))>=3 & length(which(noncarrier$y_binary==1))>=3 & length(which(carrier$y_binary==0))>=3 & length(which(noncarrier$y_binary==0))>=3){
      ans_carr<-summary(glm(carrier$y_binary~carrier$A,family = 'binomial'))
      ans_noncarr<-summary(glm(noncarrier$y_binary~noncarrier$A,family = 'binomial'))
      
      output_p[i,'carrier_PRS_OR']<-2.718^ans_carr$coefficients['carrier$A','Estimate']
      #output_p[i,'carrier_PRS_se']<-ans_carr$coefficients['carrier$A','Std. Error']
      output_p[i,'noncarrier_PRS_OR']<-2.718^ans_noncarr$coefficients['noncarrier$A','Estimate']
      #output_p[i,'noncarrier_PRS_se']<-ans_noncarr$coefficients['noncarrier$A','Std. Error']
    }
    
    #------liability-threshold model (probit) ---end----
    
    
  }
  
  #write result liability-threshold (probit)
  write.table(output_p,paste0('/data/g***/z***/prslev/result/carrier_OR/probit_',genetic_model,'/h2p_',h2_p,'_nsnp_','all_snps','_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_pi0_',pi0),quote = F,row.names = F,sep='\t')
  
  #write result logit
  write.table(output_l,paste0('/data/g***/z***/prslev/result/carrier_OR/logit_',genetic_model,'/h2p_',h2_p,'_nsnp_','all_snps','_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_pi0_',pi0),quote = F,row.names = F,sep='\t')
  
}














