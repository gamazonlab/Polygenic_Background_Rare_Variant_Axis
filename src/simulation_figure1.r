
library(ggplot2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL,by_row=T) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = by_row)
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


h2_poly=0.5
h2_levs=0.1
n_levs=5  #number of LEVs
freq_snp=0.05 #freq
h2_levs_each=h2_levs/n_levs #h2 lev for each 
k=0.01 #prevalence

#generate the LEV burden for each LEV
R_single_func<-function(){
  dosage=rbinom(n = 1, size = 1, prob = freq_snp) + rbinom(n = 1, size = 1, prob = freq_snp)
  effect=(h2_levs_each/2/freq_snp/(1-freq_snp))^0.5
  #effect=rnorm(1,(h2_levs_each/2/freq_snp/(1-freq_snp))^0.5,sd=0.005) #make a bit noise to sightly shuffle the dots up and down
  return(effect*dosage)
}

#sum up the LEV burden to get the R
R_sum_func<-function(){
  sum(replicate(n_levs,R_single_func()))
}

R<-replicate(10000,R_sum_func()) #R for 10000 subjects
A<-rnorm(10000, mean = 0, sd=h2_poly^0.5) #A for 10000 subjects

df<-data.frame(A=A,R=R)
df$AR=df$A+df$R

h2_e=1-h2_poly-h2_levs 

df$L=df$AR+rnorm(nrow(df),0,h2_e^0.5)

threshold=as.numeric(quantile(df$L,(1-k)))

p_hist<-ggplot(df, aes(x=AR)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666",bw=0.3) +
  geom_vline(aes(xintercept=3),color="red", linetype="dashed", size=1)+
  xlab('genetic burden (cSEV-PB + LEVs)')+
  ylab('Density')+
  labs(tag='a')



df$color=ifelse((df$A+df$R)>threshold,'case','control')
df$weight=ifelse((df$A+df$R)>threshold,1,0)

case<-df[df$color=='case',]
control<-df[df$color=='control',]

p_scatter<-ggplot(df, aes(x=A, y=R)) + 
  geom_point(data = df, aes(x=A, y=R,color=color),size=0.5)+
  scale_color_manual(values=c('case'='red','control'='black'))+
  geom_point(data = case, aes(x=A, y=R),color='red',size=0.5) + 
  geom_smooth(method = "lm", se = T,data = case,color='red')+
  #geom_smooth(method = "lm", se = T,data = control,color='black')+
  xlab('common, small-effect variant \npolygenic burden (cSEV-PB) ')+
  ylab('rare, pathogenic \nlarge-effect variant burden (LEVs)')+
  theme(legend.position="top",legend.title=element_blank())+
  labs(tag='b')


write.table(df,paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/fig1_new_maf_',freq_snp,'_',n_levs,'_LEV.txt'),quote = F,row.names = F,sep='\t')


#---plot---
png(paste0('~/../Dropbox/POLYGENIC_LARGE_EFFECT/simulation/plots/fig1_new_maf_',freq_snp,'_',n_levs,'_LEV.png'),width = 2000,height = 1200,res=300)
multiplot(p_hist,p_scatter,cols=2,by_row = T)
dev.off()



