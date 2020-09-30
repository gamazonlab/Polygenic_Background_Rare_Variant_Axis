#subtype analysis

args=1

subjob_id=as.numeric(args[1]) #subjob id 

N_simulation=500  #simulation times

#subjobs 
run_i=1
run_array=list()
for (h2_p in c(0.1,0.2,0.3,0.4,0.5)[3]){ #h2 poly
  for (N_causal_SNPs in c(95593)){ #N of causal SNPs
    for (vr in c(0.01,0.02,0.03,0.05,0.10)[1]){ #h2 rare (var rare)
      for (freq_r in c(1e-4,0.001,0.005,0.01,0.05)[4]){ #freq rare
        for (k in c(0.001,0.005,0.01,0.02,0.05)[3]){ #prevalence
          for (genetic_model in c(1,2,3)[1]){ #genetic architecture
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

print(N_samples)

#output dataframe (liability-threshold model (LTM), probit)
output_p<-data.frame()

#output dataframe (logit)
output_l<-data.frame()

i=1

#for (i in 1:N_simulation){
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
write.table(df,paste0('/data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.weight'),row.names = F,col.names = F,sep='\t',quote = F)

#use plink to calculate the PB (PRS of common variants)
#if(i %in% c(1,101,201,301,401)){
cmd=paste0('/home/zhoud2/tools/plink2/plink2 --bfile /data/c***/z***/projects/prslev/common/common_all --score /data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.weight variance-standardize  --out /data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.score')
system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
#}

#load PB results (A)
PB=read.table(paste0('/data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.score.sscore'),stringsAsFactors = F)
PB$IID=PB$V1;PB$A=PB$V5*2*nrow(common_all)
PB<-PB[,c('IID','A')]

#LEV (R)
rare_id=rare[i,'rsid']
rare_effect=rare[i,'effect']
R=as.numeric(rare_dosage[,rare_id]*rare_effect)
LEV<-data.frame(IID=rare_dosage$FID,R=R)

#combine 10 R

R_sum=R
for(r in 1:9){
  set.seed(r)
  R_sum=R_sum+R[sample(length(R),length(R),replace = F)]
}

LEV$R=R_sum


#merge A and R
combined<-merge(PB,LEV,by='IID')

#---sampling subjects---
set.seed(i)
combined<-combined[sample(seq(1,nrow(combined)),N_samples,replace = F),]


prop_subtype=0.2
effect_subtype=0.5

library(ggplot2)

i_plot=1
list_plot=list()

for(prop_subtype in c(0,0.1,0.2,0.3,0.5)){
  for(effect_subtype in c(0,0.5,1)){
    
    #subtyped trait
    set.seed(i+100+i_plot)
    combined$sub=rbinom(n = nrow(combined),size = 1,prob = prop_subtype)
    
    combined$x=ifelse(combined$sub==1,(1+effect_subtype)*combined$A+combined$R,combined$A+combined$R)
    
    #x=A+R
    #combined$x=combined$A+combined$R
    
    #h2_e
    h2_e=1-h2_p-vr
    
    #probability of disease (same with LTM)
    set.seed(i+i_plot)
    combined$L=combined$x+rnorm(nrow(combined),0,h2_e^0.5)
    t=quantile(combined$L,1-k)[[1]]
    combined$prob=pnorm(-(t-combined$x)/h2_e^0.5)
    
    #------logit model-----
    set.seed(i+i_plot)
    #combined$y_binary<-sapply(combined$prob,function(x) rbinom(1,1,x))
    combined$y_binary<-ifelse(combined$L>t,1,0)
    #combined$R<-ifelse(combined$R!=0,1,0)
    
    case<-combined[combined$y_binary==1,]
    
    list_plot[[letters[i_plot]]]<-ggplot(case,aes(x=R,y=A))+
      geom_point(size=0.5)+
      geom_smooth(method = "lm", se = T,data =case, color='red', fill='red')+
      xlab('LEV burden')+
      ylab('cSEV-Polygentic burden')+
      xlim(0,3)+
      ylim(-0.5,2)+
      ggtitle(paste0('proportion of subtype = ',prop_subtype,'\nadditional effect = ',effect_subtype))+
      theme(legend.position = 'none',
            panel.grid =element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            title = element_text(size=8))+
      labs(tag=letters[i_plot])
    
    i_plot=i_plot+1
    
  }
}




library(grid)
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 



png(paste0('/data/c***/z***/projects/prslev/plot/subtype/subtype_10_LEV_probit.png'),width = 2400,height = 1300,res = 150)
lay_out(list(list_plot$a,1,1),
        list(list_plot$b,2,1),
        list(list_plot$c,3,1),
        list(list_plot$d,1,2),
        list(list_plot$e,2,2),
        list(list_plot$f,3,2),
        list(list_plot$g,1,3),
        list(list_plot$h,2,3),
        list(list_plot$i,3,3),
        list(list_plot$j,1,4),
        list(list_plot$k,2,4),
        list(list_plot$l,3,4),
        list(list_plot$m,1,5),
        list(list_plot$n,2,5),
        list(list_plot$o,3,5))
dev.off()























































args=1

subjob_id=as.numeric(args[1]) #subjob id 

N_simulation=500  #simulation times

#subjobs 
run_i=1
run_array=list()
for (h2_p in c(0.1,0.2,0.3,0.4,0.5)[3]){ #h2 poly
  for (N_causal_SNPs in c(95593)){ #N of causal SNPs
    for (vr in c(0.01,0.02,0.03,0.05,0.10)[4]){ #h2 rare (var rare)
      for (freq_r in c(1e-4,0.001,0.005,0.01,0.05)[4]){ #freq rare
        for (k in c(0.001,0.005,0.01,0.02,0.05)[3]){ #prevalence
          for (genetic_model in c(1,2,3)[1]){ #genetic architecture
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

print(N_samples)

#output dataframe (liability-threshold model (LTM), probit)
output_p<-data.frame()

#output dataframe (logit)
output_l<-data.frame()

i=1

#for (i in 1:N_simulation){
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
write.table(df,paste0('/data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.weight'),row.names = F,col.names = F,sep='\t',quote = F)

#use plink to calculate the PB (PRS of common variants)
#if(i %in% c(1,101,201,301,401)){
cmd=paste0('/home/zhoud2/tools/plink2/plink2 --bfile /data/c***/z***/projects/prslev/common/common_all --score /data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.weight variance-standardize  --out /data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.score')
system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
#}

#load PB results (A)
PB=read.table(paste0('/data/g_gamazon_lab/zhoud2/prslev/tmp/',subjob_id,'.score.sscore'),stringsAsFactors = F)
PB$IID=PB$V1;PB$A=PB$V5*2*nrow(common_all)
PB<-PB[,c('IID','A')]

#LEV (R)
rare_id=rare[i,'rsid']
rare_effect=rare[i,'effect']
R=as.numeric(rare_dosage[,rare_id]*rare_effect)
LEV<-data.frame(IID=rare_dosage$FID,R=R)

#combine 10 R

R_sum=R
# for(r in 1:9){
#   set.seed(r)
#   R_sum=R_sum+R[sample(length(R),length(R),replace = F)]
# }

LEV$R=R_sum


#merge A and R
combined<-merge(PB,LEV,by='IID')

#---sampling subjects---
set.seed(i)
combined<-combined[sample(seq(1,nrow(combined)),N_samples,replace = F),]


prop_subtype=0.2
effect_subtype=0.5

library(ggplot2)

i_plot=1
list_plot=list()

for(prop_subtype in c(0,0.1,0.2,0.3,0.5)){
  for(effect_subtype in c(0,0.5,1)){
    
    #subtyped trait
    set.seed(i+100+i_plot)
    combined$sub=rbinom(n = nrow(combined),size = 1,prob = prop_subtype)
    
    combined$x=ifelse(combined$sub==1,(1+effect_subtype)*combined$A+combined$R,combined$A+combined$R)
    
    #x=A+R
    #combined$x=combined$A+combined$R
    
    #h2_e
    h2_e=1-h2_p-vr
    
    #probability of disease (same with LTM)
    set.seed(i+i_plot)
    combined$L=combined$x+rnorm(nrow(combined),0,h2_e^0.5)
    t=quantile(combined$L,1-k)[[1]]
    combined$prob=pnorm(-(t-combined$x)/h2_e^0.5)
    
    #------logit model-----
    set.seed(i+i_plot)
    #combined$y_binary<-sapply(combined$prob,function(x) rbinom(1,1,x))
    combined$y_binary<-ifelse(combined$L>t,1,0)
    #combined$R<-ifelse(combined$R!=0,1,0)
    
    case<-combined[combined$y_binary==1,]
    
    list_plot[[letters[i_plot]]]<-ggplot(case,aes(x=R,y=A))+
      geom_point(size=0.5)+
      geom_smooth(method = "lm", se = T,data =case, color='red', fill='red')+
      xlab('LEV burden')+
      ylab('cSEV-Polygentic burden')+
      xlim(0,3)+
      ylim(-0.5,2)+
      ggtitle(paste0('proportion of subtype = ',prop_subtype,'\nadditional effect = ',effect_subtype))+
      theme(legend.position = 'none',
            panel.grid =element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            title = element_text(size=8))+
      labs(tag=letters[i_plot])
    
    i_plot=i_plot+1
    
  }
}




library(grid)
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 


png(paste0('/data/c***/z***/projects/prslev/plot/subtype/subtype_1_LEV_probit.png'),width = 2400,height = 1300,res = 150)
lay_out(list(list_plot$a,1,1),
        list(list_plot$b,2,1),
        list(list_plot$c,3,1),
        list(list_plot$d,1,2),
        list(list_plot$e,2,2),
        list(list_plot$f,3,2),
        list(list_plot$g,1,3),
        list(list_plot$h,2,3),
        list(list_plot$i,3,3),
        list(list_plot$j,1,4),
        list(list_plot$k,2,4),
        list(list_plot$l,3,4),
        list(list_plot$m,1,5),
        list(list_plot$n,2,5),
        list(list_plot$o,3,5))
dev.off()







