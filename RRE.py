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
            
def expand_alignment(all_groups,settings):
    print('Expanding alignments')
    db_path = settings.expand_database_path
    for group in all_groups:
        if group.group:
            infile = group.a3m_file
        else:
            infile = group.fasta_file
        if not os.path.isfile(group.exp_alignment_file) or settings.overwrite_hhblits:
            a3m_hhblits(infile,group.exp_alignment_file,db_path,settings.cores)

def a3m_hhblits(inf,outf,db,threads=1):
    clean = inf.rpartition('.')[0]
    dumpfile = clean + '_expalign.hhr'
    commands = ['hhblits','-cpu',str(threads),'-d',db,'-i',inf,'-oa3m',outf,'-o',dumpfile,'-v','0','-n','3']
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
            hhblits(infile,outfile,db_path,settings.cores)

def hhsearch(inf,outf,db,threads):
    commands = ['hhsearch','-cpu',str(threads),'-d',db,'-i',inf,'-o',outf, '-v','0']
    call(commands)

def parse_all_RREs(genes,RRE_targets,settings):
    for gene in genes:
        results = read_hhr(gene.results_file)
        RRE_hits = find_RRE_hits(results,RRE_targets,min_prob=settings.min_prob)
        if len(RRE_hits) > 0:
            gene.RRE_data = RRE_hits
            gene.RRE_hit = True
        else:
            gene.RRE_hit = False

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
        
    
def write_results_summary(all_groups,outfile):
    with open(outfile,'w') as handle:
        header = 'Group name','Gene name','Hit name','Probability','E-value','P-value','Score','SS','Columns','Query HMM','Template HMM'
        handle.write('\t'.join(header) + '\n')
        for group in all_groups:
            if group.RRE_hit:
                if group.group:
                    for gene in group.genes:
                        for hit in group.RRE_data:
                            res = group.RRE_data[hit]
                            out = [group.name,gene,hit] + [str(i) for i in res]
                            handle.write('\t'.join(out) + '\n')
                else:
                    for hit in group.RRE_data:
                        res = group.RRE_data[hit]
                        out = ['n/a',group.name,hit] + [str(i) for i in res]
                        handle.write('\t'.join(out) + '\n')

def pipeline_operator(all_groups,settings):
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
        worker = Process(target=pipeline_worker,args=data_worker)
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
        time.sleep(10)
    while not results_queue.empty():
        res = results_queue.get()
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
            if not os.path.isfile(group.exp_alignment_file) or settings.overwrite_hhblits:
                a3m_hhblits(infile,group.exp_alignment_file,expand_db_path,1)
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
            hhsearch(infile,outfile,db_path,settings.cores)
        # Parse the results
        results = read_hhr(group.results_file)
        RRE_hits = find_RRE_hits(results,RRE_targets,min_prob=settings.min_prob)
        if len(RRE_hits) > 0:
            group.RRE_data = RRE_hits
            group.RRE_hit = True
        else:
            group.RRE_hit = False
        # print('%s %s' %(group.name,group.RRE_hit))
        print('Process id %s depositing results' %os.getpid())
        done_queue.put(group)
    print('Process %s exiting' %os.getpid())


def main(settings):
    if os.path.isfile(settings.infile):
        # Singular file
        infile = settings.infile
        print('Reading in file %s' %infile)
    elif os.path.isdir(settings.infile):
        # Get all the relevant files
        files = os.listdir(settings.infile)
        if settings.intype == 'fasta':
            exts = ['.fas','.fasta','.faa']
        elif settings.intype == 'genbank':
            exts = ['.gbk','.gbff']
        infile = [os.path.join(settings.infile,f) for f in files if any([f.endswith(ext) for ext in exts]) and os.path.isfile(os.path.join(settings.infile,f))]
        if len(infile) == 0:
            print('No %s files found in folder %s' %(settings.intype,settings.infile))
            exit()
        else:
            print('%i %s files found in folder %s' %(len(infile),settings.intype,settings.infile))
    
    # Temporary: small testset
    # infile = infile[0:10]
    
    if settings.intype == 'fasta':
        seq_dict = parse_fasta(infile)
    elif settings.intype == 'genbank':
        seq_dict = parse_genbank(infile)
    if not hasattr(settings,'project_name'):
        settings.project_name = os.path.basename(settings.infile).rpartition('.')[0]
    data_folder = os.path.join('output',settings.project_name)
    if os.path.isdir(data_folder):
        print('Warning! Output folder with name %s already found - results may be overwritten' %(settings.project_name))
    fasta_folder = os.path.join(data_folder,'fastas')
    results_folder = os.path.join(data_folder,'results')
    for folder in data_folder,fasta_folder,results_folder:
        if not os.path.isdir(folder):
            os.mkdir(folder)
    if settings.group_genes:
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
    
    if int(settings.cores) < len(all_groups):
        all_groups = pipeline_operator(all_groups,settings)
    else:
        # Get the names of the targets that are considered RRE hits
        print('Doing this')
        RRE_targets = parse_fasta(settings.rre_fasta_path).keys()
        # Write out fasta files for each gene
        write_fasta(all_groups)
        # For gene groups, align and write a3m files
        # TODO: Multiprocess muscle
        all_muscle(all_groups)
        # If the alignment needs to be expanded, do so here
        if settings.expand_alignment:
            expand_alignment(all_groups,settings)
        # Run hhsearch
        hhblits_all(all_groups,settings)
        # Parse the results
        parse_all_RREs(all_groups,RRE_targets,settings)
    outfile = os.path.join(data_folder,'%s_results.txt' %settings.project_name)
    write_results_summary(all_groups,settings.outfile)
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

    parser.add_argument('project_name',metavar='PROJECT NAME',help='A name for your project (Uses filename if none is provided')
    parser.add_argument('-i','--infile',help='File or folder to be analyzed')
    parser.add_argument('-t','--intype',help='Type of input file to be analyzed (fasta or genbank; default genbank',default='genbank')
    parser.add_argument('-m','--min_prob',help='The minimum probability for a hit to be considered significant (reads from config file if none is given)')
    parser.add_argument('--expand_alignment',help='Indicate whether or not the queries should be expanded',default=False,action='store_true')
    parser.add_argument('--group_genes',help='Group found genes first with Diamond/mcl', default=False,action='store_true')
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





