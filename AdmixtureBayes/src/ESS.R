args <- commandArgs(TRUE)
filename=args[1]
proportion=as.numeric(args[2])
outname=args[3]
summaries=args[4:length(args)]
print('summaries')
print(summaries)


library('coda')

df=read.csv(filename, header=T)
dfa=subset(df, layer==0)
dfa=as.data.frame(dfa[,summaries])
colnames(dfa) <- summaries
df2=apply(dfa,c(1,2),as.numeric)

all_nums=function(df){
  mcmcobj=mcmc(df[(floor(proportion*nrow(df))):nrow(df)])
  res_df=data.frame(summary=colnames(df))
  res_df$ESS=effectiveSize(mcmcobj)
  return(res_df)
}

esss=all_nums(df2)

write.table(esss, file= outname)
