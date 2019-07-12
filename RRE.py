import os
import sys
import argparse
import time

from Bio import SeqIO


PYTHON_VERSION = sys.version_info[0]
if PYTHON_VERSION == 2:
    import ConfigParser as configparser
elif PYTHON_VERSION == 3:
    import configparser

from subprocess import call

from multiprocessing import Process, Queue


def parse_genbank(path,out={}):
    if type(path) == list:
        for item in path:
            out = parse_genbank(item,out)
    else:
        all_seqs = open_genbank(path)
        seq_dict = gbk_to_dict(all_seqs)
        out.update(seq_dict)
    return out

def open_genbank(in_path):
    all_seqs = []
    for seq_record in SeqIO.parse(in_path, "genbank"):
        all_seqs.append(seq_record)
    return(all_seqs)

def gbk_to_dict(all_seqs):
    seq_dict = {}
    for seq in all_seqs:
        contig_id = seq.id
        if seq.id == 'unknown':
            contig_id = seq.name
        for feature in seq.features:
            if feature.type == 'CDS':
                name = contig_id
                name += '_%s-%s' %(feature.location.start,feature.location.end)
                for item in ['gene','locus_tag','product_id']:
                    if item in feature.qualifiers:
                        name += '_%s' %(feature.qualifiers[item][0])
                if 'translation' in feature.qualifiers:
                    seq = feature.qualifiers['translation'][0]
                else:
                    seq = feature.extract(seq).translate()
                seq_dict[name] = seq
    return seq_dict

def parse_fasta(path,out=None):
    if out == None:
        out = {}
    if type(path) == list:
        print('Parsing multiple')
        for item in path:
            out = parse_fasta(item,out)
    else:
        infile = open(path)
        name = False
        for line in infile:
            line = line.strip()
            if ">" in line:
                if name and name not in out:
                    out[name] = seq
                name_start = line.find(">")
                name = line[name_start+1:]
                seq = ""
            else:
                seq += line
        if name not in out:
            out[name] = seq
    return(out)
    
def make_diamond_database(protein_file,dbfile=False,threads=False,quiet=False):
    if not dbfile:
        name_clean = protein_file[:-6]
        db_file = '%s.dmnd' %name_clean
    commands = ['diamond','makedb','--in',protein_file,'-d',dbfile]
    if threads:
        commands += ['-p',str(threads)]
    if quiet:
        commands += ['--quiet']
    print(' '.join(commands))
    call(commands)
    return(dbfile)

def run_diamond(protein_file,database,outname,replacetabs=False,tmpdir='/dev/shm/',maxhits=100,sens=False,moresens=False,threads=False,evalue=False,quiet=False):
    commands = ['diamond','blastp','--query',protein_file,'--db',database,\
    '--max-target-seqs',str(maxhits),'--out',outname,'--tmpdir',tmpdir,'--outfmt',\
    '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qlen', 'slen', 'evalue', 'bitscore']
    if threads:
        commands += ['-p',str(threads)]
    if evalue:
        commands += ['--evalue',str(evalue)]
    if sens:
        commands.append('--sensitive')
    elif moresens:
        commands.append('--more-sensitive')
    if quiet:
        commands += ['--quiet']

    print(' '.join(commands))
    call(commands)
    
def run_mcl(inf,outf,threads,I=1.5,quiet=False):
    print("MCL START")
    I = str(I)
    cmd = ["mcl",inf,"-o",outf,"--abc","-I", I, '-te', str(threads)]
    if quiet:
        cmd += ['-show-log','n']
    print(' '.join(cmd))
    call(cmd)
    
def read_mcl(mcl_file,min_group_size=0):
    all_groups = []
    with open(mcl_file,'r') as f:
        for l in f:
            group = l.strip().split('\t')
            if len(group) >= min_group_size:
                all_groups.append(group)
    return(all_groups)

def parse_allvall(allvall_file,pident_cutoff):
    score_dict = {}
    with open(allvall_file) as handle:
        for line in handle:
            tabs = line.strip().split('\t')
            qseqid,sseqid,pident = tabs[0:3]
            pident = float(pident)
            if pident < pident_cutoff or qseqid == sseqid:
                continue
            if qseqid not in score_dict:
                score_dict[qseqid] = {}
            score_dict[qseqid][sseqid] = pident
    return score_dict

def average_blast_scores(allvall):
    # Average the pident from hits that gave hits as queries and as subjects
    # Only keep one entry per pair
    out = {}
    for qseqid in allvall:
        hits = allvall[qseqid]
        for hit in hits:
            if qseqid < hit:
                pident = allvall[qseqid][hit]
                try:
                    pident_b = allvall[hit][qseqid]
                    pident = float(pident + pident_b) / 2
                except KeyError:
                    pass
                if qseqid not in out:
                    out[qseqid] = {}
                out[qseqid][hit] = pident
    return out

def write_ssn(score_dict,outfile):
    with open(outfile,'w') as handle:
        for qseqid in score_dict:
            hits = score_dict[qseqid]
            for hit in hits:
                pident = score_dict[qseqid][hit]
                handle.write('%s\t%s\t%.2f\n' %(qseqid,hit,pident))
                
def write_fasta(all_groups):
    for group in all_groups:
        with open(group.fasta_file,'w') as handle:
            handle.write(group.fasta)
            
def dict_to_fasta(d,f=False,mode='w'):
    s = ''
    for key in d:
        seq = d[key]
        s += '>%s\n%s\n' %(key,seq)
    if f:
        with open(f,mode) as handle:
            handle.write(s)
    return s
            
def expand_alignment(all_groups,settings,resubmit=False):
    print('Expanding alignments')
    db_path = settings.expand_database_path
    if resubmit:
        nr_iter = settings.resubmit_initial_hhblits_iter
    else:
        nr_iter = 3
    for group in all_groups:
        if group.group:
            infile = group.a3m_file
        else:
            infile = group.fasta_file
        if not os.path.isfile(group.exp_alignment_file) or settings.overwrite_hhblits:
            a3m_hhblits(infile,group.exp_alignment_file,db_path,settings.cores,n=nr_iter)

def a3m_hhblits(inf,outf,db,threads=1,n=3):
    clean = inf.rpartition('.')[0]
    dumpfile = clean + '_expalign.hhr'
    commands = ['hhblits','-cpu',str(threads),'-d',db,'-i',inf,'-oa3m',outf,'-o',dumpfile,'-v','0','-n',str(n)]
    print(' '.join(commands))
    call(commands)
    
def all_muscle(all_groups):
    print('Aligning groups')
    for group in all_groups:
        if group.group:
            if not os.path.isfile(group.muscle_file):
                muscle(group.fasta_file,group.muscle_file)
            if not os.path.isfile(group.a3m_file):
                reformat(group.muscle_file,group.a3m_file)

def reformat(inf,outf):
    commands = ['reformat.pl','fas','a3m',inf,outf]
    print(' '.join(commands))
    call(commands)

def muscle(inf,outf,clw=False,quiet=True):
    commands = ['muscle','-in',inf,'-out',outf]
    if clw:
        commands.append('-clw')
    if quiet:
        commands.append('-quiet')
    print(' '.join(commands))
    call(commands)
    return outf

def add_ss(group,resubmit=False,remove=False):
    # Only if the alignment was also expanded
    if resubmit:
        infile = group.RRE_expalign_file
    else:
        infile = group.exp_alignment_file
    # Count the number of lines of entries
    # Only add ss if the a3m doesn't consist of a single sequence
    with open(infile) as handle:
        text = handle.read()
    if text.count('>') == 1:
        print('Only one entry found in infile %s. Not adding secondary structure' %(infile))
    else:
        newfile = infile[:-4] + '_ss.a3m'
        if not os.path.isfile(newfile):
            cmds = ['addss.pl',infile,newfile,'-a3m']
            print(' '.join(cmds))
            call(cmds)
        if resubmit:
            group.RRE_expalign_file = newfile
            if remove:
                os.remove(infile)
        else:
            group.exp_alignment_file = newfile
            if remove:
                os.remove(infile)

def add_all_ss(groups,resubmit=False):
    print('Adding secondary structure')
    for group in groups:
        add_ss(group,resubmit=resubmit)
        

def hhsearch_all(all_groups,settings):
    print('Running hhsearch')
    db_path = settings.rre_database_path
    for group in all_groups:
        if settings.rrefinder_primary_mode == 'hhpred':
            infile = group.exp_alignment_file
        elif group.group:
            infile = group.a3m_file
        else:
            infile = group.fasta_file
        outfile = group.results_file
        if not os.path.isfile(outfile) or settings.overwrite_hhblits:
            hhsearch(infile,outfile,db_path,settings.cores)

def hhsearch(inf,outf,db,threads):
    commands = ['hhsearch','-cpu',str(threads),'-d',db,'-i',inf,'-o',outf, '-v','0']
    call(commands)

def find_RRE_hits(results,targets,min_prob = 50.0,min_len = 0):
    sign_hits = {}
    for res in results:
        for target in targets:
            if res.startswith(target):
                prob = results[res][0]
                start,end = (int(i) for i in results[res][6].split('-'))
                if prob >= min_prob and (end-start) > min_len:
                    sign_hits[res] = results[res]
                    break
    return sign_hits

def read_hhr(f):
    results = {}
    with open(f) as inf:
        # skip the header
        for _ in range(9):
            l = inf.readline()
        # Now start reading
        for l in inf:
            if l == '\n':
                return results
            name = l[4:35].strip(' ')
            rest = l[35:]
            tabs = rest.strip('\n').split(' ')
            tabs = [t for t in tabs if t != '']
            prob,ev,pv,score,ss = [float(i) for i in tabs[0:5]]
            colums = int(tabs[5])
            query_hmm = tabs[6]
            template_hmm = tabs[7]
            if name not in results:
                results[name] = [prob,ev,pv,score,ss,colums,query_hmm,template_hmm]  
    return results

def make_gene_objects(seq_dict,fasta_folder,results_folder,settings):
    all_genes = []
    for gene,seq in seq_dict.items():
        gene = gene.replace(' ','_')
        for char in '();:<>':
            gene = gene.replace(char,'')
        fasta_file = os.path.join(fasta_folder,gene + '.fasta')
        results_file = os.path.join(results_folder, '%s.hhr' %(gene))
        fasta = '>%s\n%s\n' %(gene,seq)
        gene_obj = Container()
        gene_obj.setattrs(fasta_file=fasta_file,results_file=results_file,fasta=fasta,name=gene,group=False)
        if settings.rrefinder_primary_mode == 'hhpred':
            exp_alignment_file = os.path.join(fasta_folder,'%s_expalign.a3m' %gene)
            gene_obj.exp_alignment_file = exp_alignment_file
        all_genes.append(gene_obj)
        
        if settings.resubmit:
            # Set extra paths for the files
            RRE_fasta_file = os.path.join(settings.fasta_folder,gene+'_RRE.fasta')
            RRE_results_file = os.path.join(settings.results_folder,gene+'_RRE.hhr')
            RRE_expalign_file = os.path.join(settings.fasta_folder,gene+'_RRE_expalign.a3m')
            group.setattrs(RRE_fasta_file=RRE_fasta_file,RRE_results_file=RRE_results_file,RRE_expalign_file=RRE_expalign_file)
        
    return all_genes

def make_group_objects(groups,seq_dict,fasta_folder,results_folder,settings):
    all_groups = []
    group_nr = 1
    genes_seen = set()
    for group in groups:
        group_name = '%s_group_%i' %(settings.project_name, group_nr)
        seq_dict_group = dict([(gene,seq_dict[gene]) for gene in group])
        fasta_group = dict_to_fasta(seq_dict_group)
        fasta_file = os.path.join(fasta_folder,group_name + '.fasta')
        text = '_expanded' if settings.rrefinder_primary_mode == 'hhpred' else ''
        results_file = os.path.join(results_folder, '%s%s.hhr' %(group_name,text))
        muscle_file = os.path.join(fasta_folder,group_name + '_aligned.fasta')
        a3m_file = os.path.join(fasta_folder,group_name + '.a3m')
        group_obj = Container()
        group_obj.setattrs(name=group_name,genes=group,fasta_file=fasta_file,results_file=results_file,\
                              muscle_file=muscle_file,a3m_file=a3m_file,fasta=fasta_group,group=True)
        if settings.rrefinder_primary_mode == 'hhpred':
            exp_alignment_file = os.path.join(fasta_folder,'%s_expalign.a3m' %group_name)
            setattr(group_obj,'exp_alignment_file',exp_alignment_file)
        all_groups.append(group_obj)
        group_nr += 1
        
        if settings.resubmit:
            # Set extra paths for the files
            RRE_fasta_file = os.path.join(settings.fasta_folder,group_name+'_RRE.fasta')
            RRE_results_file = os.path.join(settings.results_folder,group_name+'_RRE.hhr')
            RRE_expalign_file = os.path.join(settings.fasta_folder,group_name+'_RRE_expalign.a3m')
            group.setattrs(RRE_fasta_file=RRE_fasta_file,RRE_results_file=RRE_results_file,RRE_expalign_file=RRE_expalign_file)
            RRE_alignment_file = os.path.join(fasta_folder,group_name+'_RRE.a3m')
            group.setattrs(RRE_alignment_file=RRE_alignment_file)
        
        for gene in group:
            genes_seen.add(gene)
    return all_groups,genes_seen
        
def determine_RRE_locations(groups,settings,resubmit=False):
    for group in groups:
        RRE_locations = {}
        if resubmit:
            hittype,RRE_data = group.RRE_resubmit_data
        else:
            hittype,RRE_data = group.RRE_data
        if hittype == 'hmm':
            RRE_locations[RRE_data[0]] = RRE_data[1],RRE_data[2]
        elif hittype == 'hhpred':
            for hit in RRE_data:
                prob,ev,pv,score,ss,colums,query_hmm,template_hmm = RRE_data[hit]
                if prob < settings.resubmit_initial_prob:
                    continue
                start,end = query_hmm.split('-')
                RRE_locations[RRE_hit] = start,end
        group.RRE_locations = RRE_locations
        
        
def extract_RRE(group,settings):
    # print(group.name)
    min_start = 10e6
    max_end = 0
    hittype,RRE_data = group.RRE_data
    if hittype == 'hmm':
        min_start = RRE_data[1]
        max_end = RRE_data[2]
    elif hittype == 'hhpred':
        for hit in RRE_data:
            prob,ev,pv,score,ss,colums,query_hmm,template_hmm = RRE_data[hit]
            if prob < settings.resubmit_initial_prob:
                continue
            start,end = query_hmm.split('-')
            start = int(start)
            end = int(end)
            if start < min_start:
                min_start = start
            if end > max_end:
                max_end = end
                
    fasta_out = {}

    if group.group:
        fasta = parse_fasta(group.muscle_file)
    else:
        parts = group.fasta.split('\n')
        fasta = {parts[0]:parts[1]}

    for name,seq in fasta.items():
        # If multiple genes are part of this group, extract each RRE and make a new alignment from them
        newname = '%s_RRE' %(group.name)
        RRE_part = seq[max(0,min_start-settings.extra_left):max_end+settings.extra_right]
        fasta_out[newname] = RRE_part
    dict_to_fasta(fasta_out,group.RRE_fasta_file)
    
def resubmit_group(group,RRE_targets,settings,cores):
    # First extract the RRE
    extract_RRE(group,settings)
    # Convert the alignment to .a3m if necessary and select the relevant infile for MSA expansion
    if group.group:
        if not os.path.isfile(group.RRE_alignment_file):
            reformat(group.RRE_fasta_file,group.RRE_alignment_file)
        infile = group.RRE_alignment_file
    else:
        infile = group.RRE_fasta_file
    # Expand the alignment
    db_path = settings.resubmit_database
    if not os.path.isfile(group.RRE_expalign_file):
        a3m_hhblits(infile,group.RRE_expalign_file,db_path,cores,n=3)
    if settings.addss == 2 or settings.addss == 3:
        add_ss(group,resubmit=True)
    # Run HHsearch
    db_path = settings.rre_database_path
    if not os.path.isfile(RRE_results_file):
        hhsearch(group.RRE_expalign_file,group.RRE_results_file,db_path,settings.cores)
    # Parse the results
    parse_hhpred_res(group,RRE_targets,settings,resubmit=True)
    
def parse_hhpred_res(group,RRE_targets,settings,resubmit=False):
    if resubmit:
        results = read_hhr(group.RRE_results_file)
        RRE_hits = find_RRE_hits(results,RRE_targets,min_prob=settings.min_prob,min_len = settings.min_len_alignment)
        if len(RRE_hits) > 0:
            group.RRE_resubmit_hit = True
            group.RRE_resubmit_data = ['hhpred',RRE_hits]
        else:
            group.RRE_resubmit_hit = False
    else:
        if settings.resubmit:
            prob = settings.resubmit_initial_prob
        else:
            prob = settings.min_prob
        results = read_hhr(group.results_file)
        RRE_hits = find_RRE_hits(results,RRE_targets,min_prob=prob,min_len = settings.min_len_alignment)
        if len(RRE_hits) > 0:
            group.RRE_data = ['hhpred',RRE_hits]
            group.RRE_hit = True
        else:
            group.RRE_hit = False
    
def parse_all_RREs(groups,RRE_targets,settings,resubmit=False):
    print('Parsing results')
    for group in groups:
        parse_hhpred_res(group,RRE_targets,settings,resubmit)

def resubmit_all(groups,RRE_targets,settings):
    print('Resubmitting %i found RREs' %(len([i for i in groups if i.RRE_hit])))
    for group in groups:
        if group.RRE_hit:
            resubmit_group(group,RRE_targets,settings,settings.cores)

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
            domain_found = tabs[0]
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
    
def parse_hmm_domtbl_hmmsearch(p):
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
            protein_name = tabs[0]
            domain_found = tabs[4]
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

def find_RRE_hits_hmm(groups,results,min_len=0):
    for group in groups:
        if group.name in results:
            domains_found = results[group.name]
            # Get the domain with the lowest e-value
            lowest_ev = 100
            best_domain = None
            for domain in domains_found:
                length = domain[2] - domain[1]
                if length < min_len:
                    continue
                ev = domain[-1]                
                if ev < lowest_ev:
                    lowest_ev = ev
                    best_domain = domain
            if best_domain:
                group.RRE_hit = True
                group.RRE_data = ['hmm',best_domain]
            else:
                group.RRE_hit = False
        else:
            group.RRE_hit = False
    
            
def run_hmm(groups,settings):
    tbl_out = os.path.join(settings.results_folder,'RREfinder_hmm_results.tbl')
    hmm_out = os.path.join(settings.results_folder,'RREfinder_hmm_results.txt')
    if not hasattr(settings,'fasta_file_all'):
        settings.fasta_file_all = fasta_file_all = os.path.join(settings.fasta_folder,'fasta_all.fasta')
        print('Rewriting fasta')
        # Write all fasta files
        with open(fasta_file_all,'w') as handle:
            for group in groups:
                handle.write(group.fasta)

    else:
        # Just use the input fasta if it was a singular sequence
        fasta_file_all = settings.fasta_file_all
    
    # Now run hmmer
    if not os.path.isfile(tbl_out): #Reuse old results
        hmmsearch(settings.fasta_file_all,settings.hmm_db,hmm_out,tbl_out,settings.hmm_evalue,settings.cores)
    
    # Parse the results
    results = parse_hmm_domtbl_hmmsearch(tbl_out)
    # Interpret the results and assign RRE hits
    find_RRE_hits_hmm(groups,results,settings.hmm_minlen)
    
    
def hmmsearch(infile,database,outfile,domtblout,evalue=False,cut=False,cores=1):
    commands = ['hmmsearch','--cpu',str(cores),'-o',outfile,'--domtblout',domtblout]
    if evalue:
        commands.append('-E',str(evalue))
    elif cut:
        commands.append('--%s' %cut)
    commands.extend([database,infile])
    database,infile]
    print(' '.join(commands))
    call(commands)
    
def write_results_summary(all_groups,outfile,resubmit=False,hmm=False,regulators=False):
    
    def write_group_res(group,data_attr,handle,gene=None):
        hittype,data = getattr(group,data_attr)
        text = ''
        max_reg_found = 0
        if hittype == 'hhpred':
            for hit in data:
                res = data[hit]
                if gene:
                    out = [group.name,gene,hit] + [str(i) for i in res[0:9]]
                else:
                    # Group is actually a gene
                    out = ['N\\A',group.name,hit] + [str(i) for i in res[0:9]]
                if regulators and len(res) > 9:
                    # Regulator overlap found
                    nr_reg_found = len(res[9:])
                    if nr_reg_found > max_reg_found:
                        max_reg_found = nr_reg_found
                    for reg in res[9:]:
                        out.extend(reg)
                text += '\t'.join(out) + '\n'             
        elif hittype == 'hmm':
            if gene:
                out = [group.name,gene,data[0],data[3],data[1],data[2]]
            else:
                out = ['N\\A',group.name,str(data[0]),str(data[3]),str(data[1]),str(data[2])]
            if regulators and len(data) > 4:
                max_reg_found = data[4:]
                for reg in data[4:]:
                    out.extend(reg)
            text += '\t'.join(out) + '\n'
        return text,max_reg_found

    with open(outfile,'w') as handle:
        if not hmm:
            header = ['Group name','Gene name','Hit name','Probability','E-value','P-value','Score','SS','Columns','Query HMM','Template HMM']
        else:
            header = ['Group name','Gene name','Domain name','E-value','Start','End']
        if regulators:
            header.extend(['Regulators overlapping']
            
        table_text = ''
        max_reg_found = 0
        if resubmit:
            hit_attr = 'RRE_resubmit_hit'
            data_attr = 'RRE_resubmit_data'
        else:
            hit_attr = 'RRE_hit'
            data_attr = 'RRE_data'
        for group in all_groups:
            if hasattr(group,hit_attr) and getattr(group,hit_attr):
                if group.group:
                    for gene in group.genes:
                        text,nr_reg_found = write_group_res(group,data_attr,handle,gene)
                else:
                    text,nr_reg_found = write_group_res(group,data_attr,handle)
            table_text += text
            if nr_reg_found > max_reg_found:
                max_reg_found = nr_reg_found
            
        for i in max_reg_found:
            header.extend(['Regulator name','Regulator start','Regulator end','Regulator e-value'])   
        handle.write('\t'.join(header) + '\n')
        handle.write(table_text)
        

def pipeline_operator(all_groups,settings,worker_function):
    print('Splitting work over %i processes' %settings.cores)
    work_queue = Queue()
    results_queue = Queue()
    
    def put_jobs(groups,queue,nr):
        for group in groups[0:nr]:
            queue.put(group)
        return groups[nr:]
    
    all_groups = put_jobs(all_groups,work_queue,5*settings.cores)
    print('Creating workers')
    workers = []
    data_worker = [settings,work_queue,results_queue]
    for i in range(settings.cores):
        worker = Process(target=worker_function,args=data_worker)
        workers.append(worker)
    results = []
    for worker in workers:
        worker.start()
    
    while any([w.is_alive() for w in workers]):
        print('%i workers still alive' %(len([w.is_alive() for w in workers])))
        print('Work queue: %i; all_groups: %i' %(work_queue.qsize(),len(all_groups)))
        while not results_queue.empty():
            print('Found %i results in queue' %results_queue.qsize())
            res = results_queue.get()
            results.append(res)
        if work_queue.qsize() < settings.cores:
            if len(all_groups) > 0:
                all_groups = put_jobs(all_groups,work_queue,5*settings.cores)
            else:
                for _ in range(settings.cores):
                    work_queue.put(False)
        time.sleep(1)
    while not results_queue.empty():
        res = results_queue.get()
        results.append(res)
    print('Joining workers')
    for worker in workers:
        worker.join()
    return results
    
def pipeline_worker(settings,queue,done_queue):
    print('Starting process id: %s'  %os.getpid())
    RRE_targets = parse_fasta(settings.rre_fasta_path).keys()
    expand_db_path = settings.expand_database_path
    db_path = settings.rre_database_path
    while True:
        try:
            group = queue.get()
        except:
            time.sleep(10)
            continue
        if group == False:
            break
        print('Process id %s starting work on group %s' %(os.getpid(),group.name))
        # Write out fasta files for each gene
        if not os.path.isfile(group.fasta_file):
            with open(group.fasta_file,'w') as handle:
                handle.write(group.fasta)
        # For gene groups, align and write a3m files
        if group.group:
            if not os.path.isfile(group.muscle_file):
                muscle(group.fasta_file,group.muscle_file)
            if not os.path.isfile(group.a3m_file):
                reformat(group.muscle_file,group.a3m_file)
        # If the alignment needs to be expanded, do so here
        if settings.rrefinder_primary_mode == 'hhpred':
            print('Process id %s expanding alignment' %(os.getpid()))
            if group.group:
                infile = group.a3m_file
            else:
                infile = group.fasta_file
            if settings.resubmit: 
                nr_iter = settings.resubmit_initial_hhblits_iter
            else:
                nr_iter = 3
            if not os.path.isfile(group.exp_alignment_file) or settings.overwrite_hhblits:
                a3m_hhblits(infile,group.exp_alignment_file,expand_db_path,1,nr_iter)
            if settings.addss == 1 or settings.addss == 3:
                add_ss(group)
        # Run hhblits
        if settings.rrefinder_primary_mode == 'hhpred':
            infile = group.exp_alignment_file
        elif group.group:
            infile = group.a3m_file
        else:
            infile = group.fasta_file
        outfile = group.results_file
        if not os.path.isfile(outfile) or settings.overwrite_hhblits:
            print('Process id %s running hhsearch' %(os.getpid()))
            hhsearch(infile,outfile,db_path,1)
        # Parse the results
        parse_hhpred_res(group,RRE_targets,settings,resubmit=False)
        # print('%s %s' %(group.name,group.RRE_hit))
        if group.RRE_hit and settings.resubmit:
            resubmit_group(group,RRE_targets,settings,1)
        
        print('Process id %s depositing results' %os.getpid())
        done_queue.put(group)
    print('Process %s exiting' %os.getpid())

def pipeline_resubmit_worker(settings,queue,done_queue):
    print('Starting process id: %s'  %os.getpid())
    RRE_targets = parse_fasta(settings.rre_fasta_path).keys()
    db_path = settings.rre_database_path
    while True:
        try:
            group = queue.get()
        except:
            time.sleep(10)
            continue
        if group == False:
            break
        print('Process id %s starting work on group %s' %(os.getpid(),group.name))
        resubmit_group(group,RRE_targets,settings,1)
        print('Process id %s finished resubmitting' %os.getpid())
        done_queue.put(group)
    print('Process %s exiting' %os.getpid())

def count_alignment(group,resubmit=False):
    if resubmit:
        infile = group.RRE_expalign_file
    else:
        infile = group.exp_alignment_file
    # Two lines per sequence
    c = 0
    with open(infile) as handle:
        for line in handle:
            c += 1
    return(int(c/2))
        

def rrefinder_main(settings,all_groups):
    # RREfinder main function
    
    # For hhpred, it is more efficient to split up all the sequences individually, which is done here. Each thread works on a single sequence 
    # at a time and finishes it completely
    # For hmm, it is less efficient to split up all the jobs individually, since hmmsearch simply takes a fasta file containing all sequences
    if settings.rrefinder_primary_mode == 'hmm':
        run_hmm(all_groups,settings)
        if settings.resubmit:
            if int(settings.cores) < len(all_groups) and int(settings.cores) > 1:
                # Resubmit with pipeline operator if multiple cores are used and the amount of cores is smaller than the amount of jobs
                pos_groups = pipeline_operator([i for i in all_groups if i.RRE_hit],settings,pipeline_resubmit_worker)
                all_groups = [i for i in all_groups if not i.RRE_hit] + pos_groups
            else:
                resubmit_all(all_groups,RRE_targets,settings)
            
    elif settings.rrefinder_primary_mode == 'hhpred':
        if int(settings.cores) < len(all_groups) and int(settings.cores) > 1:
            # Resubmit with pipeline operator if multiple cores are used and the amount of cores is smaller than the amount of jobs
            all_groups = pipeline_operator(all_groups,settings,pipeline_worker)
        else:
            # Run the jobs one step at a time
            # Write out fasta files for each gene
            write_fasta(all_groups)
            # For gene groups, align and write a3m files (OBSOLETE)
            all_muscle(all_groups)
            # If the alignment needs to be expanded, do so here
            if settings.rrefinder_primary_mode == 'hhpred':
                expand_alignment(all_groups,settings)
                # Add secondary structure
                if settings.addss == 1 or settings.addss == 3:
                    add_all_ss(all_groups)
            # Run hhsearch to find RRE hits
            hhsearch_all(all_groups,settings)
            # Parse the results
            parse_all_RREs(all_groups,RRE_targets,settings)
            # Resubmit if the option is given
            if settings.resubmit:
                resubmit_all(all_groups,RRE_targets,settings)
                    
    return all_groups
    
def rrefam_main(settings,all_groups):
    # Write a fasta file and run hmm
    # store the results slightly differently
    
def scan_regulators(settings,all_groups):
    # Get the fasta file of all relevant hits
    fasta_file_hits = os.path.join(settings.fasta_folder,'fasta_RRE_hits.fasta')
    hmm_file_out = os.path.join(settings.results_folder,'RRE_hits_hmmsearch_regulator.txt')
    hmm_file_tbl = os.path.join(settings.results_folder,'RRE_hits_hmmsearch_regulator.tbl')
    hits = []
    for group in all_groups:
        if settings.resubmit:
            if group.RRE_resubmit_hit:
                hits.append(group)
        else:
            if group.RRE_hit:
                hits.append(group)
    with open(fasta_file_hits,'w') as handle:
        for group in hits:
            handle.write(group.fasta)
    
    # Now run hmmer
    if not os.path.isfile(tbl_out): #Reuse old results
        hmmsearch(fasta_file_hits,settings.regulator_database,hmm_file_out,hmm_file_tbl,cores=settings.cores,cut='cut_tc')
    
    # Parse the results
    results = parse_hmm_domtbl_hmmsearch(tbl_out)
    # The RRE locations are necessary to determine the overlap
    determine_RRE_locations(all_groups,settings,resubmit=settings.resubmit):
    # Mark hits that are overlapping
    determine_regulator_overlap(results,all_groups)


def determine_regulator_overlap(hmm_results,all_groups):
    for group in groups:
        if group.name in results:
            domains_found = results[group.name]
            # Determine the overlap with the RRE for each domain
            # Append to hit data
            for hit in group.RRE_locations:
                start_hit,end_hit = group.RRE_locations[hit]
                length_hit = end_hit - start_hit
                
                hittype,RRE_data = group.RRE_data
                if hittype == 'hmm':
                    RRE_data_hit = RRE_data
                else:
                    RRE_data_hit = RRE_data[hit]    
                  
                overlaps_to_append = []
                for domain in domains_found:
                    start_hmm,end_hmm = domain[1],domain[2]
                    if start_hmm <= start_hit and end_hmm >= end_hit:
                        fraction_overlap = 1.0
                    else:
                        if start_hmm <= start_hit and end_hmm < end_hit:
                            overlap = end_hmm - start_hit 
                        elif start_hmm > start_hit and end_hmm >= end_hit:
                            overlap = end_hit - start_hmm
                        fraction_overlap = float(overlap) / length_hit
                    if fraction_overlap >= settings.min_reg_overlap:
                        # Overlap
                        overlap_data = domain + [fraction_overlap]
                        overlaps_to_append.append(overlap_data)
                # Only add up to three domains
                # Sort by the highest fraction of overlap
                overlaps_to_append.sort(key=lambda x: (x[-1],float(x[-2])))
                RRE_data_hit.append(overlaps_to_append[0:3])
                        
                        
def main(settings):
    # Prepwork
    # Make some folders
    if not hasattr(settings,'project_name'):
        settings.project_name = os.path.basename(settings.infile).rpartition('.')[0]
    if not os.path.isdir(settings.outputfolder):
        os.mkdir(settings.outputfolder)
    data_folder = os.path.join(settings.outputfolder,settings.project_name)
    if os.path.isdir(data_folder):
        log('Warning! Output folder with name %s already found - results may be overwritten' %(settings.project_name))
    fasta_folder = os.path.join(data_folder,'fastas')
    results_folder = os.path.join(data_folder,'results')
    settings.setattrs(fasta_folder=fasta_folder,results_folder=results_folder)
        
    for folder in data_folder,fasta_folder,results_folder:
        if not os.path.isdir(folder):
            os.mkdir(folder)

    # Get the names of the targets that are considered RRE hits
    RRE_targets = parse_fasta(settings.rre_fasta_path).keys()


    # Now parse the files
    if os.path.isfile(settings.infile):
        # Singular file
        infile = settings.infile
        log('Reading in file %s' %infile)
    elif os.path.isdir(settings.infile):
        # Get all the relevant files
        files = os.listdir(settings.infile)
        if settings.intype == 'fasta':
            exts = ['.fas','.fasta','.faa']
        elif settings.intype == 'genbank':
            exts = ['.gbk','.gbff']
        else:
            log('Non-legal file type given. Please choose from genbank or fasta')
            exit()
        infile = [os.path.join(settings.infile,f) for f in files if any([f.endswith(ext) for ext in exts]) and os.path.isfile(os.path.join(settings.infile,f))]
        if len(infile) == 0:
            log('No %s files found in folder %s' %(settings.intype,settings.infile))
            exit()
        else:
            log('%i %s files found in folder %s' %(len(infile),settings.intype,settings.infile))
    
    if settings.intype == 'fasta':
        seq_dict = parse_fasta(infile)
    elif settings.intype == 'genbank':
        seq_dict = parse_genbank(infile)
        
    # Make the "groups"; each group concerns a single sequence, and holds some information about the relevant files
    # In earlier versions, genes would be grouped together as a single query, since HHpred can take an alignment as 
    # input as well as a single fasta file. This was done to save time, but is incompatible with HMMER
            
    if settings.group_genes: 
        if settings.rrefinder_primary_mode == 'hmm':
            print('hmm setting as primary mode is incompatible with --group_genes')
            exit()
        blast_folder = os.path.join(data_folder,'blast')
        if not os.path.isdir(blast_folder):
            os.mkdir(blast_folder)
        
        if os.path.isfile(settings.infile) == str and settings.intype == 'fasta':
            # No need to make a new fasta file
            fasta_all = settings.infile
            pass
        else:
            fasta_all = os.path.join(blast_folder,'%s_sequences.fasta' %(settings.project_name))
            _ = dict_to_fasta(seq_dict,fasta_all)
        diamond_database = os.path.join(blast_folder,'%s_sequences.dmnd' %settings.project_name)
        allvall_path = os.path.join(blast_folder,'%s_allvall.txt' %(settings.project_name))
        ssn_file = os.path.join(blast_folder,'%s_pairs.ssn' %(settings.project_name))
        mcl_file = os.path.join(blast_folder,'%s_groups.mcl' %(settings.project_name))
        
        if not os.path.isfile(diamond_database):
            make_diamond_database(fasta_all,diamond_database,threads=settings.cores,quiet=True)
        if not os.path.isfile(allvall_path):
            run_diamond(fasta_all,diamond_database,allvall_path,threads=settings.cores,quiet=True)
        if not os.path.isfile(ssn_file):
            allvall_dict = parse_allvall(allvall_path,settings.gene_pid_cutoff)
            allvall_filtered = average_blast_scores(allvall_dict)
            write_ssn(allvall_filtered,ssn_file)
        if not os.path.isfile(mcl_file):
            run_mcl(ssn_file,mcl_file,settings.cores,quiet=True)
        groups = read_mcl(mcl_file)
        all_groups,genes_seen = make_group_objects(groups,seq_dict,fasta_folder,results_folder,settings)
        remaining_seq_dict = dict([(gene,seq_dict[gene]) for gene in seq_dict if gene not in genes_seen])
        remaining_groups = make_gene_objects(remaining_seq_dict,fasta_folder,results_folder,settings)
        all_groups += remaining_groups
    else:
        all_groups = make_gene_objects(seq_dict,fasta_folder,results_folder,settings)
    print('Continuing with %i queries, %i of which are groups of genes' %(len(all_groups), len([g for g in all_groups if g.group])))
    
    if settings.mode == 'rrefinder' or settings.mode == 'both':
        rrefinder_main(settings,all_groups)
        
    if settings.mode == 'rrefam' or settings.mode == 'both':
        rrefam_main(settings,all_groups)
    
    if settings.regulator_filter:
        scan_regulators(settings,all_groups)
    
    # For RREfinder
    outfile = os.path.join(data_folder,'%s_results.txt' %settings.project_name)
    write_results_summary(all_groups,outfile,hmm=(settings.rrefinder_primary_mode=='hmm'))
    if settings.resubmit:
        outfile_resubmit = os.path.join(data_folder,'%s_RRE_resubmit_results.txt' %(settings.project_name))
        write_results_summary(all_groups,outfile_resubmit,resubmit=True)


class Container():
    # Container for provided settings and to store info on genes
    def __init__(self):
        pass
    def __repr__(self):
        if hasattr(self,'name'):
            return('Container object (%s)' %self.name)
        else:
            return('Container object')
    def setattrs(self,**kwargs):
        for key,value in kwargs.items():
            setattr(self,key,value)
    

def parse_arguments():
    config = configparser.ConfigParser()
    config.read('config.ini')
    settings = Container()
    for section in config.sections():
        items = config.items(section)
        for item,value in items:
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    try:
                        value = config.getboolean(section,item)
                    except ValueError:
                        if ',' in value:
                            value = value.split(',')
            setattr(settings,item,value)
    
    parser = argparse.ArgumentParser()

    parser.add_argument('project_name',metavar='PROJECT NAME',help='A name for your project')
    parser.add_argument('-i','--infile',help='File or folder to be analyzed')
    parser.add_argument('-t','--intype',help='Type of input file to be analyzed (fasta or genbank; default genbank)',default='genbank')
    parser.add_argument('-o','--outputfolder',help='Folder where the output will be generated (default: output)',default='output')
    parser.add_argument('-c','--cores',help='Number of cores to use',type=int)
    parser.add_argument('-m','--mode',help='The mode to run (rrefinder,rrefam,both)', choices=['rrefinder','rrefam','both'])
    parser.add_argument('-v','--verbose',help='Verbosity (0-2)',choices=[0,1,2],type=int)
    parser.add_argument('--regulator_filter',help='Filter out found regulatory/HTH pfams',default=False,action='store_true')
    
    rrefinder = parser.add_argument_group('RREfinder settings')
    rrefinder.add_argument('--rrefinder_primary_mode',help='Choose from either hhpred or hmm for the initial scan (default: hmm)',choices=['hmm','hhpred'],default='hmm')
    rrefinder.add_argument('--no_resubmit',help='Do not resubmit found RRE hits with the resubmit database', default=True, action='store_false',dest='resubmit')
    rrefinder.add_argument('--rrefinder_min_prob',help='The minimum HHpred predicted probability for a hit to be considered significant (reads from config file if none is given)')
    rrefinder.add_argument('--addss',help=argparse.SUPPRESS, type=int)
    rrefinder.add_argument('--group_genes',help=argparse.SUPPRESS, default=False,action='store_true')

    rrefam = parser.add_argument_group('RREfam settings')
    rrefam.add_argument('--rrefam_cutoff',help='RREfam cutoff',type=float)
    
    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            setattr(settings,key,value)
            
    return settings

def log(text,f):
    print(text)
    with open(f,'a') as handle:
        handle.write(text+'\n')

if __name__ == '__main__':
    settings = parse_arguments()
    exit()
    if settings.rrefinder_primary_mode == 'hhpred' and not settings.expand_database_path:
        log('Using HHpred as initial mode for RREfinder requires a database. Please set the path in the config file')
        exit()
    t0 = time.time()
    res,data_folder = main(settings)
    t1 = time.time()
    logfile = os.path.join(data_folder,'log.txt')
    if os.path.isfile(logfile):
        os.remove(logfile)
    log('Finished. Total time: %.2f seconds (on %i cores)' %((t1-t0),settings.cores),logfile)
    log('Hits found: %i out of %i' %(len([i for i in res if i.RRE_hit]),len(res)),logfile)
    if settings.resubmit:
        log('Resubmit hits found: %i out of %i' %(len([i for i in res if i.RRE_hit and i.RRE_resubmit_hit]),len(res)),logfile)





