
myArgs <- commandArgs()

draw_to_file=myArgs[1]
leaves.f=myArgs[2]
inner_nodes.f=myArgs[3]
edges.f=myArgs[4]
admixture_events.f=myArgs[5]

leaves.t=read.csv(leaves.f, header=F)
leaves=as.character(leaves.t$V1)

inner_nodes.t=read.csv(inner_nodes.f, header=F)
inner_nodes=c(as.character(inner_nodes.t$V1),'r')

edges.t=apply(read.csv(edges.f,header=F),c(1,2),as.character)
vector_of_edges=c()
for(i in 1:nrow(edges.t)){
  vector_of_edges=c(vector_of_edges, edge(edges.t[i,1],edges.t[i,2]))
}
edges=parent_edges(vector_of_edges)

admixture_events.t=apply(read.csv(admixture_events.f,header=F),c(1,2),as.character)
vector_of_adm=c()
for(i in 1:nrow(admixture_events.t)){
  vector_of_adm=c(vector_of_adm, admix_props(admixture_events.t[i,1],admixture_events.t[i,2],
                                      admixture_events.t[i,3],admixture_events.t[i,4]))
}
admixtures=admixture_proportions(vector_of_adm)

bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
png(draw_to_file)
plot(bears_graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)
dev.off()

