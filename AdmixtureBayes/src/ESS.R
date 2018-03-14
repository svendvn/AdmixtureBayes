args <- commandArgs(TRUE)
filename=args[1]
proportion=as.numeric(args[2])
outname=args[3]
summaries=args[4:length(args)]
print('summaries')
print(summaries)


library('coda')
library('rwty')


summaries_without_trees=c()
tree_summaries=c()

for(summary in summaries){
	if(grepl('Ntree', summary)){
		tree_summaries=c(tree_summaries, summary)
	}
	else{
		summaries_without_trees=c(summaries_without_trees, summary)
	}
}




df=read.csv(filename, header=T)
dfa=subset(df, layer==0)
dfa=dfa[(floor(proportion*nrow(dfa))):nrow(dfa),]
print(dfa)
dfb=as.data.frame(dfa[,summaries_without_trees])
colnames(dfb) <- summaries_without_trees
df2=apply(dfb,c(1,2),as.numeric)

all_nums=function(df){
  mcmcobj=mcmc(df)
  res_df=data.frame(summary=colnames(df))
  return(effectiveSize(mcmcobj))
}

tree_nums=function(df){
	ids=floor(seq(1,nrow(df), length.out = min(nrow(df),500)))
	df=df[ids,]
	res=c()
	for(tree_summary in tree_summaries){
		write(paste0('\t','tree gen.', ids,' = [&U] ',as.character(df[,tree_summary]),';'),'trees_tmp.txt')
		chain=load.trees('trees_tmp.txt', type='newick')
		 invisible(capture.output(a <- topological.pseudo.ess(chain,n=1)))
		print(a)
		res=c(res, a[1,1])
	}
	return(res)
}

esss=all_nums(df2)

print(colnames(dfa))
print(tree_summaries)
dfc=as.data.frame(dfa[,tree_summaries])
colnames(dfc) <- tree_summaries
tree_esss=tree_nums(dfc)

print(summaries_without_trees)
print(tree_summaries)
print(esss)
print(tree_esss)

to_print=cbind(c(summaries_without_trees,tree_summaries),c(esss,tree_esss))

write.table(to_print, file= outname, quote=FALSE)
