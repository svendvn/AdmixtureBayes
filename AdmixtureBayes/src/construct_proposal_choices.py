from meta_proposal import simple_adaptive_proposal

def get_proposals(options):
    all_proposals=['deladmix', 'addadmix', 'rescale', 
                   'regraft', 'rescale_add', 'rescale_admixtures',
                   'rescale_admix_correction', 
                   'rescale_constrained', 'rescale_marginally', 
                   'sliding_regraft', 'sliding_rescale']
    all_proportions=[options.deladmix, options.addadmix, options.rescale, 
                     options.regraft, options.rescale_add, options.rescale_admix, 
                     options.rescale_admix_correction,
                     options.rescale_constrained, options.rescale_marginally, 
                     options.sliding_regraft, options.sliding_rescale]
    
    thinned_proportions=[]
    thinned_proposals=[]
    
    for proposal, proportion in zip(all_proposals, all_proportions):
        if proportion> 1e-8:
            thinned_proportions.append(proportion)
            thinned_proposals.append(proposal)
    return thinned_proportions, thinned_proposals

def make_proposal(options):
    if options.cancel_preserve_root_distance:
        extras={'deladmix':{'preserve_root_distance':False}, 'addadmix':{'preserve_root_distance':False}}
    else:
        extras={}
        
    if options.no_add:
        extras['rescale_constrained']={'update_add':False}
    proportions, proposals = get_proposals(options)
    mp= [simple_adaptive_proposal(proposals, proportions, extras) for _ in xrange(options.MCMC_chains)]
    return mp