#find, merge rare variants to common variant genotype files----

#calculate freq
for (i in 1:22){
  cmd=paste0('plink --bfile /data/coxvgi/zhoud2/data/biovu/geno/chr/chr',i,' --freq --keep /data/coxvgi/zhoud2/data/biovu/geno/rsid_v37/biovu_rsid_v37.fam --out /data/coxvgi/zhoud2/data/biovu/geno/chr/rare/chr',i,'_freq')
  system(cmd,wait = T)
}

#load freq
freq<-read.table('/data/coxvgi/zhoud2/data/biovu/geno/chr/rare/chr1_freq.frq',header = T,stringsAsFactors = F)
for (i in 2:22){
  print(i)
  freq_tmp<-read.table(paste0('/data/coxvgi/zhoud2/data/biovu/geno/chr/rare/chr',i,'_freq.frq'),header = T,stringsAsFactors = F)
  freq<-rbind(freq,freq_tmp)
}

rare_all<-freq[1,]
for (maf in c(0.001,0.002,0.003,0.005,0.01,0.02,0.03,0.04,0.05)){
  rare_df<-freq[abs(freq$MAF-maf)<maf*0.1,]
  
  if (nrow(rare_df)>50000){
    rare_df<-rare_df[sample(seq(1,nrow(rare_df)),50000,replace = F),]
  }
  
  rare_all<-rbind(rare_all[-1,],rare_df)
  write.table(rare_df,paste0('/data/coxvgi/zhoud2/data/biovu/geno/chr/rare/maf_',maf,'.txt'),quote = F,row.names = F,sep='\t')
}



#extract
write.table(rare_all,'/data/coxvgi/zhoud2/data/biovu/geno/chr/rare/rare_list.txt',quote = F,sep='\t',row.names = F)

for (i in 1:22){
  cmd=paste0('plink --bfile /data/coxvgi/zhoud2/data/biovu/geno/chr/chr',i,' --extract /data/coxvgi/zhoud2/data/biovu/geno/chr/rare/rare_list.txt --keep /data/coxvgi/zhoud2/data/biovu/geno/rsid_v37/biovu_rsid_v37.fam --make-bed --out /data/coxvgi/zhoud2/data/biovu/geno/chr/extract/chr',i)
  system(cmd,wait = T)
}

#get 3 allele snp list
cmd=paste0('plink --bfile /data/coxvgi/zhoud2/data/biovu/geno/rsid_v37/prune/ld0.1 --merge-list /data/coxvgi/zhoud2/data/biovu/geno/chr/extract/chr_list.txt  --make-bed --out /data/coxvgi/zhoud2/projects/prslev/geno/test/biovu')
system(cmd,wait = T)

#rm 3 allele snps
for (i in 1:22){
  cmd=paste0('plink --bfile /data/coxvgi/zhoud2/data/biovu/geno/chr/extract/chr',i,' --exclude /data/coxvgi/zhoud2/projects/prslev/geno/test/biovu-merge.missnp --make-bed --out /data/coxvgi/zhoud2/data/biovu/geno/chr/extract/ex_chr',i)
  system(cmd,wait = T)
}

#merge
cmd=paste0('plink --bfile /data/coxvgi/zhoud2/data/biovu/geno/rsid_v37/prune/ld0.1 --merge-list /data/coxvgi/zhoud2/data/biovu/geno/chr/extract/chr.list  --make-bed --out /data/coxvgi/zhoud2/projects/prslev/geno/test/biovu')
system(cmd,wait = T)
























