import os
import sys
import argparse

from Bio import SeqIO


PYTHON_VERSION = sys.version_info[0]
if PYTHON_VERSION == 2:
    import ConfigParser as configparser
elif PYTHON_VERSION == 3:
    import configparser

from subprocess import call

# from lib import muscle, parse_fasta
# from Genes import GeneCollection, CollectionCollection


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
                    seq = feature.qualifiers['translation']
                else:
                    seq = feature.extract(seq).translate()
                seq_dict[name] = seq
    return seq_dict

def parse_fasta(path):
    infile = open(path)
    out = {}
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

def set_gene_objects(seq_dict,fasta_folder,results_folder,sens):
    all_genes = []
    for gene,seq in seq_dict.items():
        fasta_file = os.path.join(fasta_folder,gene + '.fasta')
        results_file = os.path.join(results_folder, '%s_%s.hhr' %(gene,sens))
        gene_obj = Container()
        gene_obj.setattrs(fasta_file=fasta_file,results_file=results_file,seq=seq,name=gene)
        if sens:
            exp_alignment_file = os.path.join(fasta_folder,'%s_expalign.a3m' %gene)
            gene_obj.exp_alignment_file = exp_alignment_file
        all_genes.append(gene_obj)
    return all_genes

def write_fasta(all_genes):
    for gene in all_genes:
        with open(gene.fasta_file,'w') as handle:
            handle.write('>%s\n%s' %(gene.name,gene.seq))
            
def expand_alignment(all_genes,settings):
    print('Expanding alignments')
    db_path = settings['uniclust_database_path']
    for gene in all_genes:
        a3m_hhblits(gene.fasta_file,gene.exp_alignment_file,db_path,settings.cores)

def a3m_hhblits(inf,outf,db,threads=1):
    clean = inf.rpartition('.')[0]
    dumpfile = clean + '_expalign.hhr'
    commands = ['hhblits','-cpu',str(threads),'-d',db,'-i',inf,'-oa3m',outf,'-o',dumpfile]
    call(commands)

def hhblits_all(all_genes,settings):
    print('Running hhblits')
    db_path = settings.rre_database_path
    for gene in all_genes:
        if settings.sensitivity == 'highsens':
            infile = gene.exp_alignment_file
        else:
            infile = gene.fasta_file
        outfile = gene.results_file
        hhblits(infile,outfile,db_path,settings.cores)

def hhblits(inf,outf,db,threads):
    commands = ['hhblits','-cpu',str(threads),'-d',db,'-i',inf,'-o',outf, '-v','0']
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
            results[name] = [prob,ev,pv,score,ss,colums,query_hmm,template_hmm]

def write_results_summary(all_genes,outfile):
    with open(outfile,'w') as handle:
        header = 'Gene name','Hit name','Probability','E-value','P-value','Score','SS','Columns','Query HMM','Template HMM'
        handle.write('\t'.join(header) + '\n')
        for gene in all_genes:
            if gene.RRE_hit:
                for hit in gene.RRE_data:
                    res = gene.RRE_data[hit]
                    out = [gene.name,hit] + res
                    outfile.write('\t'.join(out) + '\n')
                    
    

def regroup_genes(genes,groups,operons,data_path,results_path,RRE_sens):
    print('Regrouping genes with COG groups')
    all_groups = {}
    all_genes_grouped = set()
    group_basename = 'group'
    counter = 1
    for group in groups:
        group_name = group_basename + '_' + str(counter)
        genes_used = {}
        for gene in group:
            gene_obj = operons[gene]
            if gene_obj and gene_obj.is_core:
                genes_used[gene] = gene_obj
        if len(genes_used) < 2:
            # No added value if the gene can not be grouped
            continue
        fasta_file = os.path.join(data_path, group_name + '.fasta')
        muscle_file = os.path.join(data_path, group_name + '_aligned.fas')
        a3m_file = os.path.join(data_path, group_name + '.a3m')
        results_file = os.path.join(results_path, group_name + '.hhr')
        gene_coll = GeneCollection(genes_used,name=group_name,fasta_file=fasta_file,muscle_file=muscle_file,a3m_file=a3m_file,results_file=results_file,grouped=True)
        if RRE_sens == 'high_sens':
            exp_alignment_file = os.path.join(data_path, group_name + '_expalign.a3m')
            gene_coll.exp_alignment_file = exp_alignment_file
        all_groups[group_name] = gene_coll
        counter += 1
        for gene in genes_used:
            all_genes_grouped.add(gene)
    # Now add the remaining genes
    print('Adding singleton genes')
    for gene in genes:
        if gene not in all_genes_grouped:
            group_name = group_basename + '_' + str(counter)
            fasta_file = os.path.join(data_path, group_name + '.fasta')
            results_file = os.path.join(results_path, group_name + '.hhr')
            gene_coll = GeneCollection(dict([(gene,operons[gene])]), fasta_file = fasta_file, results_file = results_file, grouped=False, name = group_name)
            if RRE_sens == 'high_sens':
                exp_alignment_file = os.path.join(data_path, group_name + '_expalign.a3m')
                gene_coll.exp_alignment_file = exp_alignment_file
            all_groups[group_name] = gene_coll
            counter += 1
    group_collection = CollectionCollection(all_groups)
    return group_collection

def read_groups(folder,basename):
    files = [f for f in os.listdir(folder) if f.startswith(basename) and f.endswith('.fasta')]
    groups_named = {}
    for f in files:
        d = parse_fasta(os.path.join(folder,f))
        groupname = f.split('.')[0]
        groups_named[groupname] = sorted(d.keys())
    return groups_named



def all_muscle(groups):
    print('Aligning groups')
    for group in groups:
        if group.grouped:
            muscle(group.fasta_file,group.muscle_file)
            reformat(group.muscle_file,group.a3m_file)

def reformat(inf,outf):
    commands = ['reformat.pl','fas','a3m',inf,outf]
    print(' '.join(commands))
    call(commands)

def main(settings):
    if settings.intype == 'fasta':
        seq_dict = parse_fasta(settings.infile)
    elif settings.intype == 'genbank':
        all_seqs = open_genbank(settings.infile)
        seq_dict = gbk_to_dict(all_seqs)
    if not hasattr(settings,'project_name'):
        settings.project_name = os.path.basename(settings.infile).rpartition('.')[0]
    data_folder = os.path.join('output',settings.project_name)
    fasta_folder = os.path.join(data_folder,'fastas')
    results_folder = os.path.join(data_folder,'results')
    for folder in data_folder,fasta_folder,results_folder:
        if not os.path.isdir(folder):
            os.mkdir(folder)
    all_genes = set_gene_objects(seq_dict,fasta_folder,results_folder,settings.sensitivity)
    # return all_genes
    # Get the names of the targets that are considered RRE hits
    RRE_targets = parse_fasta(settings.rre_fasta_path).keys()
    # Write out fasta files for each gene
    write_fasta(all_genes)
    # If the alignment needs to be expanded, do so here
    if settings.sensitivity == 'highsens':
        expand_alignment(all_genes,settings)
    # Run hhblits
    hhblits_all(all_genes,settings)
    # Parse the results
    parse_all_RREs(all_genes,RRE_targets,settings)
    write_results_summary(all_genes,settings.outfile)
    


class Container():
    # Container for provided settings
    def __init__(self):
        pass
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

    parser.add_argument('-p','--project_name',help='A name for your project (Uses filename if none is provided')
    parser.add_argument('-i','--infile',help='File or folder to be analyzed')
    parser.add_argument('-t','--intype',help='Type of input file to be analyzed (fasta or genbank; default genbank',default='genbank')
    parser.add_argument('-o','--outfile',help='File where the results will be written to')
    parser.add_argument('-m','--min_prob',help='The minimum probability for a hit to be considered significant (reads from config file if none is given)')
    parser.add_argument('-s','--sensitivity',help='The sensitivity (highsens or lowsens); Alignments are first expanded in the highsens mode, which requires the uniclust database')
    
    
    args = parser.parse_args()
    for key,value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            setattr(settings,key,value)
            
    return settings


if __name__ == '__main__':
    settings = parse_arguments()
    res = main(settings)
    





