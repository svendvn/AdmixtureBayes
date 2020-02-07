from meta_proposal import simple_adaptive_proposal
import warnings

def get_proposals(all_proposals, all_proportions):
#    all_proposals=['deladmix', 'addadmix', 'rescale',
#                   'regraft', 'rescale_add', 'rescale_admixtures',
#                   'rescale_admix_correction',
#                   'rescale_constrained', 'rescale_marginally',
#                   'sliding_regraft', 'sliding_rescale']
#    all_proportions=[options.deladmix, options.addadmix, options.rescale,
#                     options.regraft, options.rescale_add, options.rescale_admix,
#                     options.rescale_admix_correction,
#                     options.rescale_constrained, options.rescale_marginally,
#                     options.sliding_regraft, options.sliding_rescale]
    
    thinned_proportions=[]
    thinned_proposals=[]
    
    for proposal, proportion in zip(all_proposals, all_proportions):
        if proportion> 1e-8:
            thinned_proportions.append(proportion)
            thinned_proposals.append(proposal)
    return thinned_proportions, thinned_proposals

def make_proposal(deladmix,
                  addadmix,
                  rescale,
                  regraft,
                  rescale_add,
                  rescale_admix,
                  rescale_admix_correction,
                  rescale_constrained,
                  rescale_marginally,
                  sliding_regraft,
                  sliding_rescale,
                  MCMC_chains,
                  cancel_preserve_root_distance=False,
                  branch_rate_dispersion=None,
                  short_extension_mean=0,
                  short_extension_proportion=0,
                  no_add=False,
                  ):
    extras={}
    if cancel_preserve_root_distance:
        extras.update({'deladmix':{'preserve_root_distance':False}, 'addadmix':{'preserve_root_distance':False}})
    if branch_rate_dispersion is not None:
        if 'deladmix' in extras:
            extras['deladmix'].update({'gamma_branch_rate':branch_rate_dispersion})
        else:
            extras['deladmix']={'gamma_branch_rate':branch_rate_dispersion}
        if 'addadmix' in extras:
            extras['addadmix'].update({'gamma_branch_rate': branch_rate_dispersion})
        else:
            extras['addadmix']={'gamma_branch_rate':branch_rate_dispersion}
    if short_extension_mean>0 and short_extension_proportion>0:
        if 'deladmix' in extras:
            extras['deladmix'].update({'short_extension_proportion':short_extension_proportion,
                                       'short_extension_mean':short_extension_mean})
        else:
            extras['deladmix']={'short_extension_proportion':short_extension_proportion,
                                       'short_extension_mean':short_extension_mean}
        if 'addadmix' in extras:
            extras['addadmix'].update({'short_extension_proportion':short_extension_proportion,
                                       'short_extension_mean':short_extension_mean})
        else:
            extras['addadmix']={'short_extension_proportion':short_extension_proportion,
                                       'short_extension_mean':short_extension_mean}
    if no_add:
        extras['rescale_constrained']={'update_add':False}
        if rescale_add>0:
            rescale_add=0
            warnings.warn('removing the proposal rescale_add because no_add was used')

    all_proposals = ['deladmix', 'addadmix', 'rescale',
                     'regraft', 'rescale_add', 'rescale_admixtures',
                     'rescale_admix_correction',
                     'rescale_constrained', 'rescale_marginally',
                     'sliding_regraft', 'sliding_rescale']
    all_proportions = [deladmix, addadmix, rescale,
                       regraft, rescale_add, rescale_admix,
                       rescale_admix_correction,
                       rescale_constrained, rescale_marginally,
                       sliding_regraft, sliding_rescale]

    proportions, proposals = get_proposals(all_proposals, all_proportions)

    mp= [simple_adaptive_proposal(proposals, proportions, extras) for _ in xrange(MCMC_chains)]
    return mp
