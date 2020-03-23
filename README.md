# RREFinder
Bioinformatic application for the detection of RREs in protein sequences of interest

# Download guide:
1) Clone the repository

For exploratory mode:

2) Install Git LFS (https://developer.lsst.io/v/DM-7539/tools/git_lfs.html) and setup git LFS (option 1, for read-only, is sufficient)
3) Pull the databases with git lfs pull

# Installation guide:

Install the included conda environment with "conda env create -f RREfinder.yml"
Or make sure the following packages/programs are installed and in your path

- HMMer
- PSIPRED
- python package Biopython
    
At this point you can already run RREFinder in precision mode. Use python RRE.py with -m precision or --mode precision to do so.

For exploratory mode, please follow the instructions below.

Install the HHsuite V3 or above (https://github.com/soedinglab/hh-suite)
    - Make sure the following parameters are set (which is part of their recommended installation)
    
1) $HHLIB is set to the folder containing the HHsuite tool
2) The binary files and the script files are in your path 
(i.e. if you open a terminal and type "hhblits" or "addss.pl", both should be recognized commands)
Easiest way to do this is to modify the .bashrc file, and add a few extra lines:
export HHLIB="/path/to/HHsuite"
export PATH="$HHLIB/bin:$HHLIB/scripts:$PATH"
Then reboot the terminal or rerun the file (source ~/.bashrc)

- Configure the HHsuite paths:
1) Find the file HHPaths.pm in the scripts folder 
2) You should see a section looking like this:

        ##############################################################################################
        #PLEASE COMPLETE THE PATHS ... TO PSIPRED AND OLD-STYLE BLAST (NOT BLAST+) (NEEDED FOR PSIPRED)
        #our $execdir = ".../psipred/bin";         # path to PSIPRED V2 binaries
        #our $datadir = ".../psipred/data";        # path to PSIPRED V2 data files
        #our $ncbidir = ".../blast/bin";           # path to NCBI binaries (for PSIPRED in addss.pl)

3) Complete the paths
 The $execdir and the $ncbidir need to contain the psipred and the blastpgp binaries, respectively.
 The easiest way is to see where your psipred binary is located (which psipred), and point to that folder.
 In the conda package, the ncbidir is usually the same thing, as all the binaries are usually in the same folder.
 e.g. /path/to/conda/envs/RREfinder/bin

 For the $datadir, you need to find where the psipred data is stored in the conda package.
 This depends on your conda setup. Example locations include: 
 
    /path/to/conda/envs/RREfinder/share/psipred_4.01/data
 
    /path/to/conda/pkgs/psipred_4.01/share/psipred_4.01/data
 
 The data folder should contain seven files all named weights.dat or some variation of that name.

Alternatively, if you manually install PSIPRED, let the execdir point to the folder containing the psipred binary,
the datadir to the folder containing the psipred data and the ncbidir to the folder containing the legacy BLAST binary (blastpgp)
         
RREfinder is now ready for use!

# Usage guide

RREFinder accepts both .fasta files and .genbank files. 
Alternatively, you can use --antismash RiPP and let the --infile point to a .final.gbk file from antiSMASH analysis. 
This way all the antiSMASH RiPP gene clusters will be analyzed.

    python RRE.py -i my_infile.gbk project_name
or

    python RRE.py -i my_infile.fasta -t fasta project_name
or 

    python RRE.py -i antismash_output.final.gbk --antismash RiPP project_name


Output can be found under output/project_name
Both modes are used by default. To use only precision or exploratory mode, use the -m tag. 

    python RRE.py -i my_infile.gbk -m precision project_name
or 

    python RRE.py -i my_infile.fasta -t fasta -m exploratory project_name

You can also specify a range of other options, such as number of cores to use, or the bitscore cutoffs.
Use 'python RRE.py -h'  to see a list of options.

# Use RREFinder to detect RREs with the uniclust30 database (Advanced)
You can also use RREFinder as a way to run the HHPred pipeline in a commandline fashion, just to detect RREs. This provides analogous results as when submitting queries online in the HHPred pipeline. In addition, RRE regions are automatically resubmitted in the same pipeline, to verify them. This requires downloading of the uniclust30 database (https://uniclust.mmseqs.com/).

Warning! This process is very time-consuming, and not typically suitable for large-scale analysis.

HHPred is used first to find low-probability hits (<= 40). The region with flanking regions are extracted, and resubmitted  using the same pipeline, using a high-probability cutoff (>= 90).

Setup:
1) Make sure the pipeline is prepared for exploratory mode as described above.
2) Download the uniclust30 database from https://uniclust.mmseqs.com/. It is recommended to store the database on an SSD drive.
3) Point to the database in the config file, with both the expand_database and the resubmit_database variables

        resubmit_database=path/to/uniclust30
        expand_database=path/to/uniclust30
        
3) Run RRE.py with the --rrefinder-primary-mode flag set to hhpred. E.g.

        python RRE.py -i my_infile.gbk -m exploratory --rrefinder-primary-mode hhpred my_project

To run the same process without resubmitting the hits, use the --no-resubmit flag

        python RRE.py -i my_infile.gbk -m exploratory --rrefinder-primary-mode hhpred --no-resubmit my_project
