
cat('\n---Polygenic Burden (PB) - Large Effect Variant (LEV) SCAN (PB-LEV-SCAN)---\n')
cat('---version 1.0---\n')

library('ggplot2')
library("optparse")

options(scipen = 10) 

option_list = list(
  
  #path
  make_option("--plink_path", action="store", default="plink", type='character',
              help="Path to plink1.9 [%default]"),
  make_option("--plink2_path", action="store", default="plink2", type='character',
              help="Path to plink2 [%default]"),
  make_option("--maf_distribution_path", action="store", default=NA, type='character',
              help="Path to the maf distribution file (estimated from empirical data) [required]"),
  make_option("--ldsc_path", action="store", default=NA, type='character',
              help="Path to the ldsc file (estimated from empirical data) [required]"),
  make_option("--tmp_folder", action="store", default=NA, type='character',
              help="tmp folder for intermediate files [required]. Please vary the tmp folder names (e.g., tmp_1; tmp_2; tmp_3; etc.) when you have multiple jobs to run in parallel. Otherwise it will overwrite each other"),
  make_option("--out_folder", action="store", default=NA, type='character',
              help="path for output folder [required]"),
  make_option("--out_prefix", action="store", default=NA, type='character',
              help="prefix for output [required]"),
  
  #model
  make_option("--genetic_architecture", action="store", default='polygenic', type='character',
              help="polygenic=polygenic architecture \nNegativeSelection=with negative selection \nLDAK=LD-adjusted kinship"),
  make_option("--disease_model", action="store", default='LTM', type='character',
              help="LTM=liability threshold model\nlogit=logit model"),
  
  #parameter
  make_option("--n_simu", action="store", default=500, type='integer',
              help="simulation times"),
  make_option("--sample_size", action="store", default=10000, type='integer',
              help="number of sample size"),
  make_option("--h2_PB", action="store", default=0.3, type='double',
              help="heritability due to common variant based polygenicity"),
  make_option("--h2_LEV", action="store", default=NA, type='double',
              help="heritability due to large effect variant(s). Use h2_LEV or OR_LEV. If both h2_LEV and OR_LEV are provided, OR_LEV will be ignored"),
  make_option("--OR_LEV", action="store", default=NA, type='double',
              help="OR of large effect variant. Use h2_LEV or OR_LEV. If both h2_LEV and OR_LEV are provided, OR_LEV will be ignored"),
  make_option("--freq_LEV", action="store", default=0.01, type='double',
              help="allele frequency of large effect variant(s)"),
  make_option("--prevalence", action="store", default=0.01, type='double',
              help="disease prevalence in population"),
  make_option("--pi0", action="store", default=0, type='double',
              help="proportion of non-causal common, small effect variants"),
  make_option("--seed", action="store", default=1, type='integer',
              help="seed for sampling"),
  
  #other options
  make_option("--generate_a_figure", action="store_true", default=TRUE,
              help="Generate a figure comparsion the number of LEV-carriers among patients with different polygenic burdens [default: %default]"),
  make_option("--clean_tmp", action="store_true", default=TRUE,
              help="Delete the tmp folder and all temporary files [default: %default]")
  
)

opt = parse_args(OptionParser(option_list=option_list))

#input
plink_path = opt$plink_path
plink2_path = opt$plink2_path
maf_distribution_path = opt$maf_distribution_path
ldsc_path = opt$ldsc_path
tmp_folder = opt$tmp_folder
out_folder = opt$out_folder
out_prefix = opt$out_prefix

genetic_architecture = opt$genetic_architecture
disease_model = opt$disease_model

N_simulation = opt$n_simu
N_samples = opt$sample_size
h2_PB = opt$h2_PB
freq_LEV = opt$freq_LEV
k = opt$prevalence
seed = opt$seed
pi0 = opt$pi0
generate_a_figure = opt$generate_a_figure
clean_tmp = opt$clean_tmp

if(is.na(opt$h2_LEV) & is.na(opt$OR_LEV)){
  h2_LEV = 0.01
  beta_LEV = sqrt(h2_LEV/2/freq_LEV/(1-freq_LEV))
}else if(!is.na(opt$h2_LEV) & is.na(opt$OR_LEV)){
  h2_LEV = opt$h2_LEV
  beta_LEV = sqrt(h2_LEV/2/freq_LEV/(1-freq_LEV))
}else if(is.na(opt$h2_LEV) & !is.na(opt$OR_LEV)){
  #r2=var(ln(OR)*G)/(var(ln(OR)*G)+3.29)  ref Hong S. Lee, Gen Epi
  h2_LEV = (log(opt$OR_LEV,2.718)^2 * 2 * freq_LEV * (1-freq_LEV)) / (log(opt$OR_LEV,2.718)^2 * 2 * freq_LEV * (1-freq_LEV) + 3.29)
  beta_LEV = sqrt(h2_LEV/2/freq_LEV/(1-freq_LEV))
}else if(!is.na(opt$h2_LEV) & !is.na(opt$OR_LEV)){
  h2_LEV = opt$h2_LEV
  beta_LEV = sqrt(h2_LEV/2/freq_LEV/(1-freq_LEV))
  warning('Please only provide h2_LEV or OR_LEV, or OR_LEV will be ignored')
}


#load empirical ldsc
cat('INFO loading ldsc from empirical dataset\n')
ldsc_raw<-read.table(ldsc_path,header = T,stringsAsFactors = F)

#make tmp dir
if(!dir.exists(tmp_folder)){dir.create(tmp_folder)}
#make output dir
if(!dir.exists(out_folder)){dir.create(out_folder)}

#-----------------------------------------------

#df for simulation results
simu<-data.frame(n_simu=seq(1,N_simulation),wilcox_p=NA,p1=NA,p2=NA,p3=NA,p4=NA,p5=NA,p6=NA,p7=NA,p8=NA,p9=NA,p10=NA)

cat('INFO start simulation\n')

for (i in 1:N_simulation){
  #if(i %in% c(1,seq(1,1000)*10+1)){
    cat(paste0('INFO       ',i,' / ',N_simulation,' \n'))
  #}
  
  ###------simulate genotype for common causal variants------
  
  #use plink to simulate genotype data. Note: N_sample/2 doesn't means half cases and half control. See plink website for details.
  cmd=paste0('plink --simulate ',maf_distribution_path,' --simulate-ncases ',N_samples/2,' --simulate-ncontrols ',N_samples/2,' --seed ',seed,' --silent --make-bed --out ',tmp_folder,'/tmp_geno; plink --bfile ',tmp_folder,'/tmp_geno --silent --freq --out ',tmp_folder,'/tmp_freq')
  system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
  
  #10k samples ~30 secs
  #100k samples ~3 mins
  
  #df of effect size for small effect common variants
  df<-read.table(paste0(tmp_folder,'/tmp_freq.frq'),stringsAsFactors = F,header = T)
  N_causal_SNPs=nrow(df)
  
  #ldsc, matching the correlation between ldsc and MAF from empirical data
  df$id=seq(1,nrow(df))
  df<-df[order(df$MAF),]
  ldsc<-ldsc_raw[sample(nrow(ldsc_raw),nrow(df),replace = T),]
  ldsc<-ldsc[order(ldsc$MAF),]
  df$ldsc=ldsc$ldscore
  
  # a Infinitesimal architecture
  #   beta_poly~(0,h2_PB/N_caucal_SNPs)
  if(genetic_architecture=='polygenic'){
    set.seed(i)
    df$effect=(1-pi0)*rnorm(N_causal_SNPs,0,(h2_PB/N_causal_SNPs)^0.5)
    
    # b Negative selection
    #   beta_poly~N(O,k_constant([f(1-f)]^(1+alpha)))
  }else if(genetic_architecture=='NegativeSelection'){
    set.seed(i)
    df$effect=df[,'effect']<-sapply(df$MAF,function(f) rnorm(1,mean=0,sd=((f*(1-f))^0.63)^0.5))
    #find k constant
    k_constant=(h2_PB/N_causal_SNPs/var(df[,'effect']))
    #multiply by k constant
    df[,'effect']<-df[,'effect']*k_constant^0.5*(1-pi0)
    
    # c LD-adjusted kinship
    #   beta_poly~N(O,k_constant([f(1-f)]^(1+alpha)*(1/(1+ldsc))))
  }else if(genetic_architecture=='LDAK'){
    set.seed(i)
    df$effect=df[,'effect']<-mapply(function(f,ld) rnorm(1,mean=0,sd=(((f*(1-f))^0.75)*1/(1+ld))^0.5), df$MAF, df$ldsc)
    #find k constant
    k_constant=(h2_PB/N_causal_SNPs/var(df[,'effect']))
    #multiply by k constant
    df[,'effect']<-df[,'effect']*k_constant^0.5*(1-pi0)
  }
  
  
  #write weight file
  write.table(df[,c('SNP','A1','effect')],paste0(tmp_folder,'/tmp.weight'),row.names = F,col.names = F,sep='\t',quote = F)
  
  #use plink to calculate the PB (polygenic burden, i.e.,PRS of common variants)
  cmd=paste0('plink2 --bfile ',tmp_folder,'/tmp_geno --score ',tmp_folder,'/tmp.weight variance-standardize  --out ',tmp_folder,'/tmp.score')
  system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
  
  #load PB results (A)
  PB=read.table(paste0(tmp_folder,'/tmp.score.sscore'),stringsAsFactors = F)
  PB$IID=PB$V1
  PB$A=PB$V6*2*N_causal_SNPs
  PB<-PB[,c('IID','A')]
  
  
  ###------simulate genotype for the large effect variants (LEV)------
  set.seed(seed+10)
  R_raw_genotype = rbinom(n = N_samples, size = 2, prob = freq_LEV)
  #genetic risk due to LEV (R)
  R<- R_raw_genotype * beta_LEV
  
  #merge A and R
  LEV<-data.frame(IID=PB$IID,R=R,R_raw_genotype=R_raw_genotype)
  combined<-merge(PB,LEV,by='IID')
  
  #x=A+R  total genetic risk
  combined$x=combined$A+combined$R
  
  #h2_e  error/environment
  h2_e=max(1-h2_PB-var(combined$R),0)
  
  # #
  # print(h2_PB)
  # print(h2_e)
  # print(sum(combined$x))
  # print(nrow(combined))
  
  #disease liability (L)  L=A+R+e
  set.seed(i+100)
  combined$L=combined$x+rnorm(nrow(combined),0,h2_e^0.5)
  
  #threshold for liability-threshold model (LTM)
  t=quantile(combined$L,1-k)[[1]]
  
  #disease probability
  combined$prob=pnorm(-(t-combined$x)/h2_e^0.5)
  
  #generate binary trait
  if(disease_model=='LTM'){
    combined$y_binary<-ifelse(combined$L>t,1,0)
  }else if (disease_model=='logit'){
    set.seed(i*10)
    combined$y_binary<-sapply(combined$prob,function(x) rbinom(1,1,x))
  }else{
    stop('unknown disease model, please provide "LTM" or "logit"')
  }
  
  #LEV carrier
  combined$carrier<-ifelse(combined$R_raw_genotype!=0,1,0)
  
  #check the 2 by 2 table
  table_check<-table(combined$y_binary,combined$carrier)
  
  #simulated cases (patients)
  cases<-combined[combined$y_binary==1,]
  
  #test PB-LEV in cases only
  if(length(which(combined$y_binary==1 & combined$carrier==1))>0 & length(which(combined$y_binary==1 & combined$carrier==0))>0){
    #test PB-LEV correlation, one side
    ans_inverse<-wilcox.test(cases[cases$carrier==1,'A'],cases[cases$carrier==0,'A'],alternative = c("less"))
    simu[i,'wilcox_p']<-ans_inverse$p.value
  }else{
    #not enough samples in the 2 by 2 table
    simu[i,'wilcox_p']<-1
  }
  
  #group cases in to 10 equally sized bins based on their PB risk
  cases$group <- as.numeric(cut(cases$A, 10))
  
  #pool cases from each simulation
  if(i==1){
    cases_pool=cases
  }else{
    cases_pool=rbind(cases_pool,cases)
  }
  
  #number of LEV-carriers for each bin
  tmp_c<-data.frame(Group.1=seq(1,10))
  tmp_c<-merge(tmp_c,as.data.frame(aggregate(cases$carrier,list(cases$group),sum)),all=T)
  #number of samples in each bin
  tmp_n<-as.data.frame(table(cases$group))
  #merge
  tmp_c<-merge(tmp_c,tmp_n,by=1,all=T)
  tmp_c[is.na(tmp_c)]<-0
  colnames(tmp_c)<-c('group','x','n')  #'group'; # of carriers; # of cases
  #number of LEV-carriers per 1000 cases
  simu[i,3:12]<-round(tmp_c$x/tmp_c$n*1000,2)
  
}

cat(paste0('INFO simulation finished\n'))

cat(paste0('INFO saving results\n'))

#save cases results
write.table(cases_pool,paste0(out_folder,'/',out_prefix,'_cases.txt'),quote = F,sep='\t',row.names = F)

#save simu results
write.table(simu,paste0(out_folder,'/',out_prefix,'_simu.txt'),quote = F,sep='\t',row.names = F)

#result df
output<-data.frame(genetic_architecture=genetic_architecture,
                   disease_model=disease_model,
                   h2_PB=h2_PB,
                   OR_LEV=2.718^beta_LEV,
                   freq_LEV=freq_LEV,
                   prevalence=k,
                   pi0=pi0,
                   N_samples=N_samples,
                   N_simulation=N_simulation,
                   seed=seed,
                   utility=length(which(simu[,'wilcox_p']<0.05))/nrow(simu)
)

#mean and sd for the number of LEV-carriers per 1000 cases
for(j in 1:10){
  output[1,paste0('mean_p',j)]<-mean(simu[,paste0('p',j)],na.rm = T)
  output[1,paste0('sd_p',j)]<-sd(simu[,paste0('p',j)],na.rm = T)
}

#save result
write.table(output,paste0(out_folder,'/',out_prefix,'.txt'),quote = F,sep='\t',row.names = F)



#---figure generator (LEV_carriers_per_1000_cases_comparison)---


if(generate_a_figure){
  
  cat(paste0('INFO figure generating\n'))
  
  #get the x-coordinates for each cSEV-PB bins
  p_tmp<-ggplot(data=cases_pool,aes(x=A))+  #use the A from last simulation
    geom_histogram(bins = 10)
  
  #df for # carriers in each bin
  carrier_info<-data.frame(pos=ggplot_build(p_tmp)$data[[1]]$x,
                           mean=as.numeric(output[1,paste0('mean_p',seq(1,10))]),
                           sd=as.numeric(output[1,paste0('sd_p',seq(1,10))]))
  carrier_info$lower=carrier_info$mean-carrier_info$sd
  carrier_info$upper=carrier_info$mean+carrier_info$sd
  carrier_info$lower=ifelse(carrier_info$lower<0,0,carrier_info$lower)
  carrier_info$upper=ifelse(carrier_info$upper>1000,1000,carrier_info$upper)
  
  #find ylim (the highest density)
  y_density_max<-max(ggplot_build(p_tmp)$data[[1]]$density)*1.2
  
  #coeff for adjusting dual y-axis
  coeff=max(carrier_info$upper,na.rm = T)/y_density_max
  
  #plot
  p<-ggplot(data=carrier_info)+
    
    #cSEV-PB distribution
    geom_histogram(data=cases_pool,aes(x=A,y=..density..), colour="black", fill="white",bins=10)+
    geom_density(data=cases_pool,aes(x=A),alpha=.2, fill="#FF6666",bw=0.2)+
    
    #mean of 'number of LEV-carriers per 1000 cases' 
    geom_point(aes(x=carrier_info$pos,y=carrier_info$mean/coeff),color='dodgerblue3')+
    
    #error bar for 'number of LEV-carriers per 1000 cases' (+- 1sd)
    geom_pointrange(aes(x=carrier_info$pos,y=carrier_info$mean/coeff,ymin=carrier_info$lower/coeff, ymax=carrier_info$upper/coeff),color='dodgerblue3')+
    
    scale_y_continuous(
      # Features of the first axis
      name = "density",
      #limits for density
      limits = c(0,y_density_max),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="number of LEV-carriers per 1000 cases")
    )+
    
    xlab('polygenic burden (PB)')+
    
    theme(legend.title = element_blank(),
          panel.grid =element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = 'none',
          axis.title.y.right = element_text(color = 'dodgerblue3'),
          axis.text.y.right = element_text(color = 'dodgerblue3'))
  
  #figure
  pdf(paste0(out_folder,'/',out_prefix,'_LEV_carriers_per_1000_cases_comparison_figure.pdf'),height = 5,width = 5)
  suppressWarnings(print(p))
  dev.off()
  
}

if(clean_tmp){
  cat(paste0('INFO cleaning up tmp folder\n'))
  cmd=paste0('rm -r ',tmp_folder)
  system(cmd,wait = F)
}

cat(paste0('INFO finished.\n'))













