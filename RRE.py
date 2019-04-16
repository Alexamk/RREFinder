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
        if settings.expand_alignment:
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

def find_RRE_hits(results,targets,min_prob = 50.0):
    sign_hits = {}
    for res in results:
        for target in targets:
            if res.startswith(target):
                prob = results[res][0]
                if prob >= min_prob:
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

def make_gene_objects(seq_dict,fasta_folder,results_folder,expanded):
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
        if expanded:
            exp_alignment_file = os.path.join(fasta_folder,'%s_expalign.a3m' %gene)
            gene_obj.exp_alignment_file = exp_alignment_file
        all_genes.append(gene_obj)
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
        text = '_expanded' if settings.expand_alignment else ''
        results_file = os.path.join(results_folder, '%s%s.hhr' %(group_name,text))
        muscle_file = os.path.join(fasta_folder,group_name + '_aligned.fasta')
        a3m_file = os.path.join(fasta_folder,group_name + '.a3m')
        group_obj = Container()
        group_obj.setattrs(name=group_name,genes=group,fasta_file=fasta_file,results_file=results_file,\
                              muscle_file=muscle_file,a3m_file=a3m_file,fasta=fasta_group,group=True)
        if settings.expand_alignment:
            exp_alignment_file = os.path.join(fasta_folder,'%s_expalign.a3m' %group_name)
            setattr(group_obj,'exp_alignment_file',exp_alignment_file)
        all_groups.append(group_obj)
        group_nr += 1
        for gene in group:
            genes_seen.add(gene)
    return all_groups,genes_seen
        
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
            # print(prob,ev,pv,score,ss,colums,query_hmm,template_hmm)
            if prob < settings.resubmit_initial_prob:
                continue
            start,end = query_hmm.split('-')
            # print(start,end)
            start = int(start)
            end = int(end)
            if start < min_start:
                # print('Replacing start')
                min_start = start
            if end > max_end:
                # print('Replacing end')
                max_end = end
    fasta_out = {}
    
    if group.group:
        fasta = parse_fasta(group.muscle_file)
    else:
        parts = group.fasta.split('\n')
        fasta = {parts[0]:parts[1]}
        
    for name,seq in fasta.items():
        # print('Determining region')
        # print(min_start,max_end)
        # If multiple genes are part of this group, extract each RRE and make a new alignment from them
        newname = '%s_RRE' %(group.name)
        RRE_part = seq[max(0,min_start-settings.extra_left):max_end+settings.extra_right]
        # print('Extracting from %i to %i' %(max(0,min_start-settings.extra_left),max_end+settings.extra_right))
        fasta_out[newname] = RRE_part
    dict_to_fasta(fasta_out,group.RRE_fasta_file)
    
def resubmit_group(group,RRE_targets,settings,cores):
    # Set extra paths for the files
    RRE_fasta_file = os.path.join(settings.fasta_folder,group.name+'_RRE.fasta')
    RRE_results_file = os.path.join(settings.results_folder,group.name+'_RRE.hhr')
    RRE_expalign_file = os.path.join(settings.fasta_folder,group.name+'_RRE_expalign.a3m')
    group.setattrs(RRE_fasta_file=RRE_fasta_file,RRE_results_file=RRE_results_file,RRE_expalign_file=RRE_expalign_file)
    if group.group:
        RRE_alignment_file = os.path.join(fasta_folder,group.name+'_RRE.a3m')
        group.setattrs(RRE_alignment_file=RRE_alignment_file)
    # Now extract the RRE
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
    if not os.path.isfile(RRE_expalign_file):
        a3m_hhblits(infile,RRE_expalign_file,db_path,cores,n=3)
    if settings.addss == 2 or settings.addss == 3:
        add_ss(group,resubmit=True)
    # Run HHsearch
    db_path = settings.rre_database_path
    if not os.path.isfile(RRE_results_file):
        hhsearch(group.RRE_expalign_file,RRE_results_file,db_path,settings.cores)
    # Parse the results
    parse_res(group,RRE_targets,settings,resubmit=True)
    
def parse_res(group,RRE_targets,settings,resubmit=False):
    if resubmit:
        results = read_hhr(group.RRE_results_file)
        RRE_hits = find_RRE_hits(results,RRE_targets,min_prob=settings.min_prob)
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
        RRE_hits = find_RRE_hits(results,RRE_targets,min_prob=prob)
        if len(RRE_hits) > 0:
            group.RRE_data = ['hhpred',RRE_hits]
            group.RRE_hit = True
        else:
            group.RRE_hit = False
    
def parse_all_RREs(groups,RRE_targets,settings,resubmit=False):
    print('Parsing results')
    for group in groups:
        parse_res(group,RRE_targets,settings,resubmit)

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

def find_RRE_hits_hmm(groups,results):
    for group in groups:
        if group.name in results:
            domains_found = results[group.name]
            # Get the domain with the lowest e-value
            lowest_ev = 100
            for domain in domains_found:
                ev = domain[-1]
                if ev < lowest_ev:
                    lowest_ev = ev
                    best_domain = domain
            group.RRE_hit = True
            group.RRE_data = ['hmm',best_domain]
        else:
            group.RRE_hit = False
    
            
def run_hmm(groups,settings):
    tbl_out = os.path.join(settings.results_folder,'hmm_results.tbl')
    hmm_out = os.path.join(settings.results_folder,'hmm_results.txt')
    if not hasattr(settings,'fasta_file_all'):
        fasta_file_all = os.path.join(settings.fasta_folder,'fasta_all.fasta')
        # Write all fasta files
        with open(fasta_file_all,'w') as handle:
            for group in groups:
                handle.write(group.fasta)

    else:
        # Just use the input fasta if it was a singular sequence
        fasta_file_all = settings.fasta_file_all
    
    # Now run hmmer
    commands = ['hmmscan','--cpu',str(settings.cores),'-E',str(settings.hmm_evalue),'-o',hmm_out,'--domtblout',tbl_out,\
                settings.hmm_db,settings.fasta_file_all]
    print(' '.join(commands))
    call(commands)
    
    # Parse the results
    results = parse_hmm_domtbl(tbl_out)
    # Interpret the results and assign RRE hits
    find_RRE_hits_hmm(groups,results)
    
    
def write_results_summary(all_groups,outfile,resubmit=False,hmm=False):
    
    def write_group_res(group,data_attr,handle,gene=None):
        hittype,data = getattr(group,data_attr)
        if hittype == 'hhpred':
            for hit in data:
                res = data[hit]
                if gene:
                    out = [group.name,gene,hit] + [str(i) for i in res]
                else:
                    # Group is actually a gene
                    out = ['N\\A',group.name,hit] + [str(i) for i in res]
                handle.write('\t'.join(out) + '\n')
        elif hittype == 'hmm':
            if gene:
                out = [group.name,gene,data[0],data[3],data[1],data[2]]
            else:
                out = ['N\\A',group.name,str(data[0]),str(data[3]),str(data[1]),str(data[2])]
            handle.write('\t'.join(out) + '\n')

    
    with open(outfile,'w') as handle:
        if not hmm:
            header = 'Group name','Gene name','Hit name','Probability','E-value','P-value','Score','SS','Columns','Query HMM','Template HMM'
        else:
            header = 'Group name','Gene name','Domain name','E-value','Start','End'
        handle.write('\t'.join(header) + '\n')
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
                        write_group_res(group,data_attr,handle,gene)
                else:
                    write_group_res(group,data_attr,handle)

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
        if settings.expand_alignment:
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
        if settings.expand_alignment:
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
        parse_res(group,RRE_targets,settings,resubmit=False)
        # print('%s %s' %(group.name,group.RRE_hit))
        if group.RRE_hit and settings.resubmit:
            resubmit_group(group,RRE_targets,settings,1)
        
        print('Process id %s depositing results' %os.getpid())
        done_queue.put(group)
    print('Process %s exiting' %os.getpid())

def pipeline_resubmit_worker(settings,queue,done_queue):
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
        

def main(settings):
    if os.path.isfile(settings.infile):
        # Singular file
        infile = settings.infile
        if settings.hmm:
            settings.fasta_file_all = settings.infile
        print('Reading in file %s' %infile)
    elif os.path.isdir(settings.infile):
        # Get all the relevant files
        files = os.listdir(settings.infile)
        if settings.intype == 'fasta':
            exts = ['.fas','.fasta','.faa']
        elif settings.intype == 'genbank':
            exts = ['.gbk','.gbff']
        else:
            print('Non-legal file type given. Please choose from genbank or fasta')
            exit()
        infile = [os.path.join(settings.infile,f) for f in files if any([f.endswith(ext) for ext in exts]) and os.path.isfile(os.path.join(settings.infile,f))]
        if len(infile) == 0:
            print('No %s files found in folder %s' %(settings.intype,settings.infile))
            exit()
        else:
            print('%i %s files found in folder %s' %(len(infile),settings.intype,settings.infile))
    
    if settings.intype == 'fasta':
        seq_dict = parse_fasta(infile)
    elif settings.intype == 'genbank':
        seq_dict = parse_genbank(infile)
    if not hasattr(settings,'project_name'):
        settings.project_name = os.path.basename(settings.infile).rpartition('.')[0]
    if not os.path.isdir('output'):
        os.mkdir('output')
    data_folder = os.path.join('output',settings.project_name)
    if os.path.isdir(data_folder):
        print('Warning! Output folder with name %s already found - results may be overwritten' %(settings.project_name))
    fasta_folder = os.path.join(data_folder,'fastas')
    results_folder = os.path.join(data_folder,'results')
    settings.setattrs(fasta_folder=fasta_folder,results_folder=results_folder)
    if settings.hmm:
        settings.fasta_file_all = fasta_file_all = os.path.join(fasta_folder,'fasta_all.fasta')
        
    for folder in data_folder,fasta_folder,results_folder:
        if not os.path.isdir(folder):
            os.mkdir(folder)
            
    if settings.group_genes:
        if settings.hmm:
            print('--hmm setting is incompatible with --group_genes')
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
        remaining_groups = make_gene_objects(remaining_seq_dict,fasta_folder,results_folder,settings.expand_alignment)
        all_groups += remaining_groups
    else:
        all_groups = make_gene_objects(seq_dict,fasta_folder,results_folder,settings.expand_alignment)
    print('Continuing with %i queries, %i of which are groups of genes' %(len(all_groups), len([g for g in all_groups if g.group])))
    
    # Get the names of the targets that are considered RRE hits
    RRE_targets = parse_fasta(settings.rre_fasta_path).keys()

    # For hhpred, it is more efficient to split up all the sequences individually, which is done here. Each thread works on a single sequence 
    # at a time and finishes it completely
    # For hmm, it is less efficient to split up all the jobs individually, since hmmscan simply takes a fasta file containing all sequences
    if settings.hmm:
        run_hmm(all_groups,settings)
        if settings.resubmit:
            if int(settings.cores) < len(all_groups) and int(settings.cores) > 1:
                # Resubmit with pipeline operator
                pos_groups = pipeline_operator([i for i in all_groups if i.RRE_hit],settings,pipeline_resubmit_worker)
                all_groups = [i for i in all_groups if not i.RRE_hit] + pos_groups
            else:
                resubmit_all(all_groups,RRE_targets,settings)
            
    else:
        if int(settings.cores) < len(all_groups) and int(settings.cores) > 1:
            all_groups = pipeline_operator(all_groups,settings,pipeline_worker)
        else:
            if not settings.hmm:
                # Write out fasta files for each gene
                write_fasta(all_groups)
                # For gene groups, align and write a3m files
                all_muscle(all_groups)
                # If the alignment needs to be expanded, do so here
                if settings.expand_alignment:
                    expand_alignment(all_groups,settings)
                    # Add secondary structure
                    if settings.addss == 1 or settings.addss == 3:
                        add_all_ss(all_groups)
                # Run hhsearch
                hhsearch_all(all_groups,settings)
                # Parse the results
                parse_all_RREs(all_groups,RRE_targets,settings)
                # Resubmit if the option is given
                if settings.resubmit:
                    resubmit_all(all_groups,RRE_targets,settings)
                    
    outfile = os.path.join(data_folder,'%s_results.txt' %settings.project_name)
    write_results_summary(all_groups,outfile,hmm=settings.hmm)
    if settings.resubmit:
        outfile_resubmit = os.path.join(data_folder,'%s_RRE_resubmit_results.txt' %(settings.project_name))
        write_results_summary(all_groups,outfile_resubmit,resubmit=True)
    return all_groups

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
    parser.add_argument('-m','--min_prob',help='The minimum probability for a hit to be considered significant (reads from config file if none is given)')
    parser.add_argument('--expand_alignment',help='Indicate whether or not the queries should be expanded',default=False,action='store_true')
    parser.add_argument('--group_genes',help='Group found genes first with Diamond/mcl', default=False,action='store_true')
    parser.add_argument('--resubmit',help='Resubmit found RRE hits with the resubmit database', default=False, action='store_true')
    parser.add_argument('--addss',help='Add secondary structure prediction (improves sensitivity)\n0 = Never\n1 = Only during initial run\n2 = Only during resubmitting\n3 = Always', type=int, default=0)
    parser.add_argument('--hmm',help='Run hmmer to prefilter hits instead of hhpred',default=False,action='store_true')
    parser.add_argument('-c','--cores',help='Number of cores to use',type=int)
    
    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            setattr(settings,key,value)
            
    return settings

if __name__ == '__main__':
    settings = parse_arguments()
    if settings.expand_alignment and not settings.expand_database_path:
        print('Expanding the alignment requires a database. Please set the path in the config file')
        exit()
    t0 = time.time()
    res = main(settings)
    t1 = time.time()
    print('Finished. Total time: %.2f seconds' %(t1-t0))
    print('Hits found: %i out of %i' %(len([i for i in res if i.RRE_hit]),len(res)))
    if settings.resubmit:
        print('Resubmit hits found: %i out of %i' %(len([i for i in res if i.RRE_hit and i.RRE_resubmit_hit]),len(res)))




