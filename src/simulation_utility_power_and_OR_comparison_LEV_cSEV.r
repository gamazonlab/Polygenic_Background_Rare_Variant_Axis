#The simulation script includes three aims.
#when pi=0, to estimate the utility of detecting invers-axis in cases only.
#when pi!=0, to estimate the power of detecting invers-axis in cases only.
#simulations to compare the OR of disease risk between LEV and cSEV-PB

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

#subjobs 
run_i=1
run_array=list()
for (h2_p in c(0.1,0.2,0.3,0.4,0.5)){ #h2 poly
  for (N_causal_SNPs in c(95593)){ #N of causal SNPs
    for (vr in c(0.01,0.02,0.03,0.05,0.10)){ #h2 rare (var rare)
      for (freq_r in c(1e-4,0.001,0.005,0.01,0.05)){ #freq rare
        for (k in c(0.001,0.005,0.01,0.02,0.05)){ #prevalence
          for (genetic_model in c(1,2,3)){ #genetic architecture
            for (pi0 in c(0,0.1,0.3,0.5,0.7,0.9)[1]){  #pi0, proportion of causal variants pi0=0: utility; pi0!=0: power
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

for (N_samples in c(500,1000,2000,3000,5000,10000)){  
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
    #if(i %in% c(1,101,201,301,401)){
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
    
    #--------logit model---start-------
    set.seed(i)
    combined$y_binary<-sapply(combined$prob,function(x) rbinom(1,1,x))
    combined$R<-ifelse(combined$R!=0,1,0)
    
    #test inverse axis in cases only
    if(length(which(combined$y_binary==1 & combined$R==1))>0 & length(which(combined$y_binary==1 & combined$R==0))>0){
      cases<-combined[combined$y_binary==1,]
      #test inverse axis
      ans_inverse<-wilcox.test(cases[cases$R==1,'A'],cases[cases$R==0,'A'],alternative = c("less"))
      output_l[i,'case_inverse_wilcox_p']<-ans_inverse$p.value
    }else{
      output_l[i,'case_inverse_wilcox_p']<-1
    }
    
    #test inverse axis in controls only
    if(length(which(combined$y_binary==0 & combined$R==1))>0 & length(which(combined$y_binary==0 & combined$R==0))>0){
      controls<-combined[combined$y_binary==0,]
      #test inverse axis
      ans_inverse<-wilcox.test(controls[controls$R==1,'A'],controls[controls$R==0,'A'],alternative = c("less"))
      output_l[i,'control_inverse_wilcox_p']<-ans_inverse$p.value
    }else{
      output_l[i,'control_inverse_wilcox_p']<-1
    }
    
    
    #--OR comparision between LEV and PB--
    for (top_p in c(0.01,0.05,0.1)){
      cut_off=quantile(combined$A,1-top_p)[[1]]
      combined$A_binary<-ifelse(combined$A>cut_off,1,0)
      
      A_check=table(combined$y_binary,combined$A_binary)
      R_check=table(combined$y_binary,combined$R)
      
      if(length(A_check)==4 & min(A_check)>0 & length(R_check)==4 & min(R_check)>0){
        ans_A<-summary(glm(combined$y_binary~combined$A_binary,family = 'binomial'))
        ans_R<-summary(glm(combined$y_binary~combined$R,family = 'binomial'))
        
        output_l[i,paste0('A_top_',top_p,'_beta')]<-ans_A$coefficients['combined$A_binary','Estimate']
        output_l[i,paste0('A_top_',top_p,'_se')]<-ans_A$coefficients['combined$A_binary','Std. Error']
        output_l[i,paste0('R_beta')]<-ans_R$coefficients['combined$R','Estimate']
        output_l[i,paste0('R_se')]<-ans_R$coefficients['combined$R','Std. Error']
      }
    }
    #--------logit model---end-------
    
    #------liability-threshold model (probit) ---start-----
    combined$y_binary<-ifelse(combined$L>t,1,0)
    table_check<-table(combined$y_binary,combined$R)
    
    #test inverse axis in cases only
    if(length(which(combined$y_binary==1 & combined$R==1))>0 & length(which(combined$y_binary==1 & combined$R==0))>0){
      cases<-combined[combined$y_binary==1,]
      #test inverse axis
      ans_inverse<-wilcox.test(cases[cases$R==1,'A'],cases[cases$R==0,'A'],alternative = c("less"))
      output_p[i,'case_inverse_wilcox_p']<-ans_inverse$p.value
    }else{
      output_p[i,'case_inverse_wilcox_p']<-1
    }
    
    #test inverse axis in controls only
    if(length(which(combined$y_binary==0 & combined$R==1))>0 & length(which(combined$y_binary==0 & combined$R==0))>0){
      controls<-combined[combined$y_binary==0,]
      #test inverse axis
      ans_inverse<-wilcox.test(controls[controls$R==1,'A'],controls[controls$R==0,'A'],alternative = c("less"))
      output_p[i,'control_inverse_wilcox_p']<-ans_inverse$p.value
    }else{
      output_p[i,'control_inverse_wilcox_p']<-1
    }
    
    
    #--OR comparision between LEV and PB--
    for (top_p in c(0.01,0.05,0.1)){
      cut_off=quantile(combined$A,1-top_p)[[1]]
      combined$A_binary<-ifelse(combined$A>cut_off,1,0)
      
      A_check=table(combined$y_binary,combined$A_binary)
      R_check=table(combined$y_binary,combined$R)
      
      if(length(A_check)==4 & min(A_check)>0 & length(R_check)==4 & min(R_check)>0){
        ans_A<-summary(glm(combined$y_binary~combined$A_binary,family = 'binomial'))
        ans_R<-summary(glm(combined$y_binary~combined$R,family = 'binomial'))
        
        output_p[i,paste0('A_top_',top_p,'_beta')]<-ans_A$coefficients['combined$A_binary','Estimate']
        output_p[i,paste0('A_top_',top_p,'_se')]<-ans_A$coefficients['combined$A_binary','Std. Error']
        output_p[i,paste0('R_beta')]<-ans_R$coefficients['combined$R','Estimate']
        output_p[i,paste0('R_se')]<-ans_R$coefficients['combined$R','Std. Error']
      }
    }
    
    #------liability-threshold model (probit) ---end-----
    
  }
  
  #write result liability-threshold (probit)
  write.table(output_p,paste0('/data/g***/z***/prslev/result/probit_',genetic_model,'/h2p_',h2_p,'_nsnp_','all_snps','_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_pi0_',pi0),quote = F,row.names = F,sep='\t')
  
  #write result logit model
  write.table(output_l,paste0('/data/g***/z***/prslev/result/logit_',genetic_model,'/h2p_',h2_p,'_nsnp_','all_snps','_vr_',vr,'_freqr_',freq_r,'_k_',k,'_Nsample_',N_samples,'_pi0_',pi0),quote = F,row.names = F,sep='\t')
  
}













