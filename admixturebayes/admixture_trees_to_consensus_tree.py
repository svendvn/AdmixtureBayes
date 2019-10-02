from Rtree_to_covariance_matrix import get_populations
from Treemix_to_AdmixtureBayes import Node
from construct_nodes_choices import read_one_line
from collections import Counter
import pandas as pd
from argparse import ArgumentParser, SUPPRESS
from tree_statistics import identifier_to_tree_clean,generate_predefined_list_string,topological_identifier_to_tree_clean, identifier_to_tree
from copy import deepcopy
from generate_sadmix_trees import effective_number_of_admixes
from Rtree_operations import get_number_of_admixes, node_is_admixture, rename_key, get_admixture_proportion_from_key, get_all_admixture_origins
from tree_plotting import plot_node_structure_as_directed_graph, plot_as_directed_graph #NOTICE THAT THIS IS CALLED ELSEWHERE!! IN THE SCRIPT
import sys
from posterior_quantiles import branch_and_proportion_quantiles

def main(args):
    parser = ArgumentParser(usage='pipeline for plotting posterior distribution summaries.', version='0.3')

    parser.add_argument('--posterior_distribution_file', required=True, type=str, help='The file containing posterior distributions from the "AdmixtureBayes posterior" command. It needs the two columns "pops" and topology.')
    parser.add_argument('--plot', choices=['consensus_trees', 'top_node_trees', 'top_trees','estimates'], required=True,
                        help='The type of plot to make. Choose between: 1) consensus_trees. '
                             'It plots an admixture graph based on all nodes that have a higher (marginal) posterior probability of X. '
                             'Different X\'s can be supplied with the command --consensus_threshold \n'
                             '2) top_node_trees. It plots the X highest posterior combinations of node types '
                             'and creates the corresponding minimal topologies.  X can be supplied through the command --top_node_trees_to_plot'
                             '3) top_trees. It plots the X highest posterior topologies. X can be supplied by the command --top_trees_to_plot.'
                             '4) estimates. It creates a table  with continuous parameters estimated from the posterior sample'
                             'It also plots the concerned topologies with labels. It does this for either the X highest posterior topologies '
                             'or the topologies specified by --estimate_topologies.'
                             'If no --estimate_topologies is set, X can be set by top_trees_to_estimate (default=3). ')
    parser.add_argument('--outgroup', default='outgroup', help='name of the outgroup to plot')
    parser.add_argument('--prefix', default='', type=str, help='string to prepend before each file created by this routine. '
                                                               'That means that any rankings written to a file by setting --write_ranking_to_file or --write_estimates_to_file will have this prefix prepended.')
    parser.add_argument('--consensus_threshold', default=[0.25, 0.5, 0.75, 0.9, 0.95, 0.99], type=float, nargs='+',
                        help='The posterior thresholds for which to draw different consensus trees.')
    parser.add_argument('--top_node_trees_to_plot', type=int, default=3,
                        help='The number of node trees (or minimal topologies) to plot')
    parser.add_argument('--top_trees_to_plot', type=int, default=3,
                        help='The number of trees (or topologies) to plot ')
    parser.add_argument('--top_trees_to_estimate', type=int, default=3,
                        help='The number of trees (or topologies) to plot ')
    parser.add_argument('--estimate_topologies', default=[], type=str, nargs='+',
                        help='The topologies whose conitnuous parameters should be estimated with "--plot estimates"')
    parser.add_argument('--write_estimates_to_file', default=[], type=str, nargs='+',
                        help='The file in which to put the tables when plotting estimates. ')
    parser.add_argument('--write_ranking_to_file', type=str, default='', help='if a file is supplied here, the natural rankings for each of the plots is written here.')
    parser.add_argument('--rankings_to_write_to_file', type=int, default=1000,
                        help='the number of rankings(nodes, min topology or topology depending on --plot) to write to the ranking file.')
    parser.add_argument('--dont_annotate_node_posterior', default=False, action='store_true',
                        help='This will not color the nodes according to their posterior probability.')
    parser.add_argument('--nodes', default='', type=str, help='file where the first line is the leaf nodes')
    parser.add_argument('--suppress_plot', default=False, action='store_true')
    parser.add_argument('--popup', default=False, action='store_true')
    parser.add_argument('--no_sort', default=False, action='store_true', help='often the tree is sorted according to the leaf names. no_sort willl assumed that they are not sorted according to this but sorted according to ')
    parser.add_argument('--sep', default=',', type=str, help='the separator used in the input file')

    #parser.add_argument('--no_header', default=False, action='store_true',help='will assume that there is no header in the file')
    #parser.add_argument('--burn_in_rows', default=0, type=int, help='the number of rows that will be skipped in the input file as burn-in period')
    #parser.add_argument('--burn_in_fraction', default=0.0, type=float, help='the proportion of the rows that are discarded as burn in period')
    #parser.add_argument('--tree_column_name', default='tree', type=str, help='the name in the header of the column with all the trees.')
    parser.add_argument('--consensus_method', choices=['descendant_frequencies'], default='descendant_frequencies', help=SUPPRESS)#'Which method should be used to calculate the consensus tree?')
    #parser.add_argument('--min_w', default=0.0, type=float, help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')









    #parser.add_argument('--plot_tops_file', action='store_true', default=False, help='this will assume that the file is a tops file from downstream_analysis_parser and plot each line numbered.')

    #parser.add_argument('--get_effective_number_of_admixtures', action='store_true', default=False, help='this will cancel all the other analysis and only print the effective number of admixes(tadmixes/sadmixes or admixes) to a a file.')
    #parser.add_argument('--effective_number_of_admixtures_file', type=str, default='no_tadmixes.txt', help='this is the file in which to write the effective number of admixes in the file')
    #parser.add_argument('--type_of_effective_admixtures', type=str, choices=['sadmix','tadmix','admix'], help='this is the type of admixes to write to the file.')

    #parser.add_argument('--node_count_file', default='', type=str, help='if plot_tops option is supplied')
    #parser.add_argument('--node_count_probs', default='', type=str, help='if supplied this will make a new ')
    #parser.add_argument('--test_run', default=False, action='store_true',
    #                    help='will overwrite everything and run a test function')

    options= parser.parse_args(args)



    def combine_nodes(node_structure, new_node, seen_sets):
        candidate=new_node.name
        seen=[]
        for lists_of_fixed_size in seen_sets[::-1]:
            for attached_branch in lists_of_fixed_size:
                if( attached_branch.issubset(candidate) and
                   ((not attached_branch.issubset(seen)) or (not node_structure[attached_branch].has_parent()))):
                    seen.extend(list(attached_branch))
                    new_node.add_child(node_structure[attached_branch])
                    node_structure[attached_branch].add_parent(new_node)
        return node_structure


    def get_number_of_tadmixtures(node_structure):
        total=0
        for key in node_structure:
            total+=max(0,node_structure[key].get_number_of_parents()-1)
        return total

    def node_combinations_to_node_structure(node_combinations):
        length_sorted={}
        for node_combination in node_combinations:
            leaves=frozenset(node_combination.split('.'))
            k=len(leaves)
            if k in length_sorted:
                length_sorted[k].append(leaves)
            else:
                length_sorted[k]=[leaves]
        length_sorted_list=[length_sorted.get(k,[]) for k in range(1,max(length_sorted.keys())+1)]
        #length_sorted_list is of the form [[[A],[B],[C]],[[A,B],[B,C]],...,[[A,B,C]]]
        node_structure={}
        for leaf_node in length_sorted_list[0]:
            node_structure[leaf_node]=Node(leaf_node)
        added_sets=[length_sorted_list[0]]
        for lists_of_fixed_size in length_sorted_list[1:]:
            for branch_set in lists_of_fixed_size:
                new_node=Node(branch_set)
                combine_nodes(node_structure, new_node, added_sets)
                node_structure[branch_set]=new_node
            added_sets.append(lists_of_fixed_size)
        return node_structure

    # if options.node_count_file:
    #     with open(options.node_count_file, 'r') as f:
    #         node_count_dic={}
    #         for lin in f.readlines():
    #             key,freq=lin.rstrip().split()
    #             node_count_dic[frozenset(key.split('.'))]=float(freq)
    # else:
    #     node_count_dic=None

    if options.plot=='consensus_trees' or options.plot=='top_node_trees':
        df = pd.read_csv(options.posterior_distribution_file, sep=options.sep, usecols=['pops'])
        nodes_list = df['pops'].tolist()
        #print(nodes_list)
        seen_combinations = {}
        for nodes in nodes_list:
            #print(nodes)
            for node in nodes.split('-'):
                #print(node)
                seen_combinations[node] = seen_combinations.get(node, 0) + 1
        N = len(nodes_list)
        #print(seen_combinations)
        if options.plot=='consensus_trees':
            node_combinations = []
            for threshold in options.consensus_threshold:
                total_threshold = int(N * threshold)
                final_node_combinations = [k for k, v in seen_combinations.items() if v > total_threshold]
                node_combinations.append(final_node_combinations)
            if not options.dont_annotate_node_posterior:
                node_count_dic={frozenset(k.split('.')):float(v)/N for k,v in seen_combinations.items()}
            else:
                node_count_dic=None
            for i, final_node_combinations in enumerate(node_combinations):
                #print(final_node_combinations)
                final_node_structure = node_combinations_to_node_structure(final_node_combinations)
                if not options.suppress_plot:
                    from tree_plotting import plot_node_structure_as_directed_graph
                    plot_node_structure_as_directed_graph(final_node_structure, drawing_name=options.prefix+'consensus_'+str(int(100*options.consensus_threshold[i]))+'.png', node_dic=node_count_dic,  popup=options.popup)
            if options.write_ranking_to_file:
                with open(options.prefix+options.write_ranking_to_file, 'w') as f:
                    c = Counter(seen_combinations)
                    to_write = c.most_common(options.rankings_to_write_to_file)
                    for node, frequency in to_write:
                        f.write(node+','+str(float(frequency)/N)+'\n')
        elif options.plot=='top_node_trees':
            c=Counter(nodes_list)
            to_plots=c.most_common(options.top_node_trees_to_plot)
            if options.write_ranking_to_file:
                with open(options.prefix+options.write_ranking_to_file, 'w') as f:
                    for tree, frequency in c.most_common(options.rankings_to_write_to_file):
                        f.write(tree + ',' + str(float(frequency) / N) + '\n')
            if not options.dont_annotate_node_posterior:
                c=Counter(seen_combinations)
                node_count_dic={frozenset(key.split('.')):float(count)/N for key,count in c.most_common(1000)}
            else:
                node_count_dic=None
            if not options.suppress_plot:
                from tree_plotting import plot_node_structure_as_directed_graph
                for i, (to_plot,count) in enumerate(to_plots):
                    node_structure = node_combinations_to_node_structure(to_plot.split('-'))
                    plot_node_structure_as_directed_graph(node_structure, drawing_name=options.prefix+'minimal_topology_' +str(i+1)+'.png',
                                                          node_dic=node_count_dic,  popup=options.popup)
    elif options.plot=='top_trees':
        df = pd.read_csv(options.posterior_distribution_file, sep=options.sep, usecols=['pops','topology'])
        trees_list = df['topology'].tolist()
        no_leaves=len(trees_list[0].split('-')[0].split('.'))
        N=len(trees_list)
        c = Counter(trees_list)
        to_plots = c.most_common(options.top_trees_to_plot)

        #obtaining nodes:
        if not options.nodes:
            nodes=df['pops'].tolist()[0].split('-')
            leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
            if len(leaves)==no_leaves:
                pass #everything is good
            elif len(leaves)==no_leaves-1:
                #adding outgroup
                leaves.append(options.outgroup)
            else:
                assert False, 'The number of leaves could not be obtained'
            assert not options.no_sort, 'When nodes are not specified, they will always be sorted'
            leaves=sorted(leaves)
        else:
            leaves=read_one_line(options.nodes)
            if not options.no_sort:
                leaves=sorted(leaves)

        if options.write_ranking_to_file:
            with open(options.prefix+options.write_ranking_to_file, 'w') as f:
                for tree, frequency in c.most_common(options.rankings_to_write_to_file):
                    f.write(tree + ',' + str(float(frequency) / N) + '\n')

        if not options.suppress_plot:
            from tree_plotting import plot_as_directed_graph
            for i, (to_plot, count) in enumerate(to_plots):
                tree=topological_identifier_to_tree_clean(to_plot, leaves=generate_predefined_list_string(deepcopy(leaves)))
                plot_as_directed_graph(tree,drawing_name=options.prefix+'topology_' + str(i + 1) + '.png', popup=options.popup)
    elif options.plot=='estimates':
        try:
            df = pd.read_csv(options.posterior_distribution_file, sep=options.sep, usecols=['string_tree', 'topology', 'pops'])
        except ValueError as e:
            raise Exception('Unexpected columns in the posterior_distribution file. Did you turn on the --faster flag in AdmixtureBayes posterior?')

        topologies_list = df['topology'].tolist()
        string_tree_list=df['string_tree'].tolist()
        if options.estimate_topologies:
            cleaned_topology_list=[s.split('=')[-1].split(';')[0] for s in options.estimate_topologies]
        else:
            c = Counter(topologies_list)
            to_plots = c.most_common(options.top_trees_to_estimate)
            cleaned_topology_list=[d[0] for d in to_plots]
        no_leaves = len(topologies_list[0].split('-')[0].split('.'))

        # obtaining nodes:
        if not options.nodes:
            nodes=df['pops'].tolist()[0].split('-')
            leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
            if len(leaves)==no_leaves:
                pass #everything is good
            elif len(leaves)==no_leaves-1:
                #adding outgroup
                leaves.append(options.outgroup)
            else:
                assert False, 'The number of leaves could not be obtained'
            assert not options.no_sort, 'When nodes are not specified, they will always be sorted'
            leaves=sorted(leaves)
        else:
            leaves=read_one_line(options.nodes)
            if not options.no_sort:
                leaves=sorted(leaves)


        if not options.suppress_plot:
            from tree_plotting import plot_as_directed_graph
        for i, to_plot in enumerate(cleaned_topology_list):
            relevant_string_trees=[]
            for string_tree, topology in zip(string_tree_list, topologies_list):
                if topology==to_plot:
                    relevant_string_trees.append(string_tree)
            branches_intervals, admixture_proportion_intervals=branch_and_proportion_quantiles(relevant_string_trees)
            branch_names=[branches_interval[0] for branches_interval in branches_intervals]

            admixture_names=[ad[0] for ad in admixture_proportion_intervals]
            tree = identifier_to_tree(to_plot,
                                      leaves=generate_predefined_list_string(deepcopy(leaves)),
                                      branch_lengths=generate_predefined_list_string(deepcopy(branch_names)),
                                      admixture_proportions=generate_predefined_list_string(deepcopy(admixture_names)))

            org_keys = tree.keys()
            for key in org_keys:
                node = tree[key]
                if node_is_admixture(node):
                    new_name = get_admixture_proportion_from_key(tree, key)
                    tree = rename_key(tree, key, new_name)
            adms=get_all_admixture_origins(tree)
            adm_interpretation={}
            for key, (branch_name, node_destination) in adms.items():
                adm_interpretation[key]='For the lineages that passes through {}, this is the proportion that follows branch {} to node {}'.format(key, branch_name,node_destination)
            if not options.suppress_plot:
                plot_as_directed_graph(tree, drawing_name=options.prefix+'topology_labels_' + str(i + 1) + '.png', plot_edge_lengths=True,  popup=options.popup)
            if options.write_estimates_to_file:
                branch_file=options.write_estimates_to_file[i*2+0]
                admixtures_file=options.write_estimates_to_file[i*2+1]
            else:
                branch_file=options.prefix+'topology_estimates_branches_'+str(i+1)+'.txt'
                admixtures_file=options.prefix+'topology_estimates_admixtures_'+str(i+1)+'.txt'
            with open(branch_file, 'w') as f:
                f.write(','.join(['branch label','lower 95%','mean','upper 95%'])+'\n')
                for v in branches_intervals:
                    f.write(','.join(map(str,v))+'\n')
            with open(admixtures_file, 'w') as f:
                f.write(','.join(['branch label','lower 95%', 'mean', 'upper 95%','interpretation'])+'\n')
                for v in admixture_proportion_intervals:

                    f.write(','.join(map(str,list(v)+[adm_interpretation[v[0]]]))+'\n')



    sys.exit()






    if options.plot_tops_file:
        with open(options.input_file, 'r') as f:
            for n,lin in enumerate(f.readlines()):
                rank, probability, combination=lin.rstrip().split(',')
                all_nodes=[c.split('.') for c in combination.split('_')]
                flattened=[item for sublist in all_nodes for item in sublist]
                a=list(set(flattened))
                code=rank+'_'+str(int(100*round(float(probability),2)))+'_'+'_'.join(a)
                print 'code',code
                node_structure=node_combinations_to_node_structure(combination.split('_'))

                print node_structure
                plot_node_structure_as_directed_graph(node_structure, drawing_name=code+'.png', node_dic=node_count_dic)
        sys.exit()


    if options.test_run:
        from generate_prior_trees import generate_phylogeny
        from tree_statistics import unique_identifier_and_branch_lengths
        from tree_plotting import plot_node_structure_as_directed_graph, plot_as_directed_graph
        N=5
        tree1=generate_phylogeny(N,1)
        plot_as_directed_graph(tree1, drawing_name='tree1.png')
        tree2=generate_phylogeny(N,1)
        plot_as_directed_graph(tree2, drawing_name='tree2.png')
        stree1=unique_identifier_and_branch_lengths(tree1)
        stree2=unique_identifier_and_branch_lengths(tree2)
        with open('tmp_tree.txt','w') as f:
            f.write(' '.join(['s'+str(i) for i in range(1,N+1)])+'\n')
            f.write(stree1)
        with open('trees.txt','w') as f:
            f.write(stree1+'\n'+stree2+'\n'+stree1)

        options.input_file='trees.txt'
        options.nodes='tmp_tree.txt'
        options.no_header=True
        options.posterior_threshold=[0.25,0.5,0.9]

    if options.input_file==options.node_count_file:
        node_combinations=[]
        print 'using population sets from ', options.node_count_file
        for threshold in options.posterior_threshold:
            final_node_combinations=['.'.join(sorted(list(k))) for k,v in node_count_dic.items() if v > threshold]
            node_combinations.append(final_node_combinations)
    else:
        print 'Reading file...'
        #loading trees
        if options.no_header:
            strees=[]
            with open(options.input_file, 'r') as f:
                for lin in f.readlines():
                    strees.append(lin.rstrip())
        else:
            df=pd.read_csv(options.input_file, sep=options.sep, usecols=[options.tree_column_name])
            strees=df[options.tree_column_name].tolist()
        n=len(strees)
        print 'trees read: ',n

        #thinning tree list

        rows_to_remove_from_fraction=int(options.burn_in_fraction*n)
        rows_to_remove=max(rows_to_remove_from_fraction, options.burn_in_rows)
        strees=strees[rows_to_remove:]

        print 'removed burn-in:', rows_to_remove
        print 'In list are now', len(strees),'trees'

        #thinning

        distance_between=max(1,len(strees)//options.max_number_of_trees)
        nstrees=[]
        for a,stree in enumerate(strees):
            if a%distance_between==0 and len(nstrees)<options.max_number_of_trees:
                nstrees.append(stree)
        print 'thinned'
        print 'In list are now', len(nstrees),'trees'

        N=len(nstrees)

        seen_node_combinations={}

        nodes=read_one_line(options.nodes)
        if not options.no_sort:
            nodes=sorted(nodes)

        tenth=len(nstrees)//10
        trees=[]
        for i,stree in enumerate(nstrees):
            if tenth>0 and i%tenth==0:
                print i//tenth*10, '%'
            if ';' in stree:
                tree=identifier_to_tree_clean(stree, leaves=generate_predefined_list_string(deepcopy(nodes)))
            else:
                tree=topological_identifier_to_tree_clean(stree, leaves=generate_predefined_list_string(deepcopy(nodes)))
            trees.append(tree)
            ad=get_populations(tree, min_w=options.min_w)
            for a in ad:
                seen_node_combinations[a]=seen_node_combinations.get(a,0)+1
        node_combinations=[]
        for threshold in options.posterior_threshold:
            total_threshold=int(N*threshold)
            final_node_combinations=[k for k,v in seen_node_combinations.items() if v > total_threshold]
            node_combinations.append(final_node_combinations)

    for i,final_node_combinations in enumerate(node_combinations):
        print 'final_node_combinations', final_node_combinations
        final_node_structure=node_combinations_to_node_structure(final_node_combinations)
        if options.get_effective_number_of_admixtures:
            with open(options.effective_number_of_admixtures_file, 'w') as f:
                if options.type_of_effective_admixtures=='tadmix':
                    effictive_admixtures=get_number_of_tadmixtures(final_node_structure)
                    f.write(str(effictive_admixtures))
                elif options.type_of_effective_admixtures=='sadmix':
                    val=0
                    count=0
                    for tree in trees:
                        val+=effective_number_of_admixes(tree)
                        count+=1
                    if count==1:
                        f.write(str(int(val)))
                    else:
                        f.write(str(float(val)/count))
                elif options.type_of_effective_admixtures=='admix':
                    val=0
                    count=0
                    for tree in trees:
                        val+=get_number_of_admixes(tree)
                        count+=1
                    if count==1:
                        f.write(str(int(val)))
                    else:
                        f.write(str(float(val)/count))
        if not options.suppress_plot:
            from tree_plotting import plot_node_structure_as_directed_graph, plot_as_directed_graph
            plot_node_structure_as_directed_graph(final_node_structure, drawing_name='tmp'+str(i+1)+'.png', node_dic=node_count_dic)




if __name__=='__main__':
    main(sys.argv[1:])






