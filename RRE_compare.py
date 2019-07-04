
def parse_results(infile):
    d = {}
    with open(infile) as handle:
        readfirst = False
        for line in handle:
            if not readfirst:
                readfirst = True
                continue
            tabs = line.strip().split('\t')
            d[tabs[1]] = tabs[2:]
    return(d)
    
    

    
def parse_hmm_domtbl(p):
    # Parse per domain found
    # Sort out overlap later
    with open(p) as f:
        outd = {}
        for _ in range(3):
            l = f.readline()
        for l in f:
            if l.startswith('#'):
                continue
            tabs = l.strip().split(' ')
            tabs = [tab for tab in tabs if tab != '']
            protein_name = tabs[3]
            domain_found = tabs[1]
#                print protein_name,domain_found
            if '.' in domain_found:
                domain_found = domain_found.rpartition('.')[0]
            domain_evalue = tabs[12]
            if 'e' in domain_evalue:
                parts = domain_evalue.split('e')
#                    print domain_evalue,parts
                domain_evalue = float(parts[0])*10**int(parts[1])
            else:
                domain_evalue = float(domain_evalue)
#                print(tabs[17],tabs[18])
            seq_start = int(tabs[17])
            seq_end = int(tabs[18])
            if seq_start > seq_end:
#                print 'Seq start after end: %s' %protein_name
                pass
            if protein_name not in outd:
                outd[protein_name] = []
            outd[protein_name].append([domain_found,seq_start,seq_end,domain_evalue])
    return outd
    
def filter_hmm(d,ev=0.001):
    d_out = {}
    for key,value in d.items():
        doms_remain = []
        for dom_found in value:
            if dom_found[-1] <= ev:
                doms_remain.append(dom_found)
        if doms_remain != []:
            d_out[key] = doms_remain
    return(d_out)

def compare_pairs(dicts):
    shared = set(dicts[0].keys())
    for d in dicts[1:]:
        shared = shared & set(d.keys())
    return shared

def get_all_subgroups(groups):
    if len(groups) == 2:
        yield
    else:
        for group in groups:
            new_groups = groups[:]
            new_groups.remove(group)
            yield new_groups
            for gr in get_all_subgroups(new_groups):
                if gr:
                    yield gr
    
if __name__ == '__main__':
    all_res = {}
    # all_res['repr_rre'] = parse_results('output/MIBIG_repr_only_rre_resubmit_addss3/MIBIG_repr_only_rre_resubmit_addss3_RRE_resubmit_results.txt')
 #   all_res['uniclust'] = parse_results('output/MIBIG_uniclust_resubmit_addss3/MIBIG_uniclust_resubmit_addss3_RRE_resubmit_results.txt')

#    all_res['rre'] = parse_results('output/MIBIG_repr_resubmit_addss3/MIBIG_repr_resubmit_addss3_RRE_resubmit_results_p40.txt')    
    # all_res['hmm_uni'] = parse_results('/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/victorc/van_wezel_group/alex/rrefinder/output/MIBIG_hmm_resubmit_uni_addss3_test8/MIBIG_hmm_resubmit_uni_addss3_test8_RRE_resubmit_results.txt')
    # all_res['hmm_rre'] = parse_results('/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/victorc/van_wezel_group/alex/rrefinder/output/MIBIG_hmm_resubmit_custom_db_addss3_test7/MIBIG_hmm_resubmit_custom_db_addss3_test7_RRE_resubmit_results.txt')

#    all_res['hmm_v3_uni'] = parse_results('output/MIBiG_vs_phmm_v3_iter3_resubmit_uni/MIBiG_vs_phmm_v3_iter3_resubmit_uni_RRE_resubmit_results.txt')
#   all_res['hmm_v3_uni_ss'] = parse_results('output/MiBIG_vs_phmm_v3_resubmit_uni_addss3/MiBIG_vs_phmm_v3_resubmit_uni_addss3_RRE_resubmit_results.txt')
#    all_res['hmm_v3_custom_ss'] = parse_results('output/MIBIG_vs_phmm_v3_iter3_resubmit_custom_addss3/MIBIG_vs_phmm_v3_iter3_resubmit_custom_addss3_RRE_resubmit_results.txt')

#    all_res['ripp_hmm_custom'] = parse_results('output/RiPPTIDEtest_vs_phmm_v3_resubmit_custom_addss3/RiPPTIDEtest_vs_phmm_v3_resubmit_custom_addss3_RRE_resubmit_results.txt')

#    all_res['mibig_all_hmm_custom'] = parse_results('output/MIBiG_all_phmm_v3_resubmit_custom_addss3/MIBiG_all_phmm_v3_resubmit_custom_addss3_RRE_resubmit_results.txt')

#    all_res['hmm_v2_i1'] = parse_results('output/MIBiG_vs_phmm_v2_iter2_timetest/MIBiG_vs_phmm_v2_iter2_timetest_RRE_resubmit_results.txt')
#    all_res['hmm2'] = parse_results('output/MIBiG_vs_phmm_v2_iter2/MIBiG_vs_phmm_v2_iter2_results.txt')
#    all_res['hmm_v2_i1'] = parse_results('output/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit_RRE_resubmit_results.txt')
#    all_res['hmm2'] = parse_results('output/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit_results.txt')

#    all_res['hmm2'] = parse_results('output/MIBiG_phmm_iter3_minlen55_ev_min3_resubmit/MIBiG_phmm_iter3_minlen55_ev_min3_resubmit_RRE_resubmit_results.txt')
#    all_res['hmm2'] = parse_results('output/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit_custom_iter1/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit_custom_iter1_RRE_resubmit_results.txt')
#    all_res['hmm2'] = parse_results('output/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit_custom_iter3/MIBiG_phmm_iter3_minlen50_ev_min3_resubmit_custom_iter3_RRE_resubmit_results.txt')

#    all_res['hmmv2_custom'] = parse_results('output/RiPPTIDEtest_vs_phmm_v2_iter3_resubmit_custom/RiPPTIDEtest_vs_phmm_v2_iter3_resubmit_custom_RRE_resubmit_results.txt')
#    all_res['hmmv2_uni'] = parse_results('output/RiPPTIDEtest_vs_phmm_v2_iter3_resubmit_uni/RiPPTIDEtest_vs_phmm_v2_iter3_resubmit_uni_RRE_resubmit_results.txt')
#    all_res['hmmv3_custom'] = parse_results('output/RiPPTIDEtest_vs_phmm_v3_iter3_resubmit_custom/RiPPTIDEtest_vs_phmm_v3_iter3_resubmit_custom_RRE_resubmit_results.txt')
#    all_res['hmmv3_uni'] = parse_results('output/RiPPTIDEtest_vs_phmm_v3_iter3_resubmit_uni/RiPPTIDEtest_vs_phmm_v3_iter3_resubmit_uni_RRE_resubmit_results.txt')

    all_res['custom'] = parse_results('output/RiPPTIDEtest_vs_phmm_v3_resubmit_custom_addss3/RiPPTIDEtest_vs_phmm_v3_resubmit_custom_addss3_RRE_resubmit_results.txt')
    all_res['uni'] = parse_results('output/RiPPTIDEtest_vs_phmm_v3_resubmit_uni_addss3/RiPPTIDEtest_vs_phmm_v3_resubmit_uni_addss3_RRE_resubmit_results.txt')

    for res in all_res:
        print('Group: %s. Number of results: %s' %(res,len(all_res[res])))
    all_shared = compare_pairs(all_res.values())
    print('All groups: %i' %len(all_shared))
    
    for group in get_all_subgroups(all_res.keys()):
        ds = [all_res[i] for i in group]
        shared = compare_pairs(ds)
        print('Comparing group with %s. Shared: %i' %(', '.join(group),len(shared)))
    
    
    


