# RREFinder
Bioinformatic application for the detection of RREs in protein sequences of interest

## Download guide:
1) Clone the repository

For exploratory mode:

Databases need to be downloaded independently from https://zenodo.org/record/3733240#.XoRkIUEzbCI. 
Please download all the files, except the folder with the seed sequences, and place them into the data/database folder within the downloaded repository.
This can also be done by running ''download_RRE_databases.py''

## Installation guide:

Install the included conda environment with "conda env create -f RREfinder.yml"
    
At this point you can already run RREFinder in precision mode. Use python RRE.py with -m precision or --mode precision to do so.

For exploratory mode, four more steps are necessary to use all functionalities of the HHSuite. These are accomplished by running the ''setup_RRE_exploratory.py'' script in the activated environment for conda-style installations. 

In detail:

1) The environment variables HHLIB needs to be set and point towards the folder containing the HHSuite bin and scripts folder 

2) The scripts folder needs to be added to path

3) The HHPaths.pm file in the scripts folder needs to be updated to point to the right binaries

4) The databases need to be downloaded

Exploratory mode relies on the secondary structure predictions of HHSuite. Although the HHSuite is installed with Conda, the secondary structure prediction script is not. Nevertheless, it is downloaded along with the rest of the package, and can be found in the environment scripts folder (Typically found under ''$CONDA_PREFIX/scripts''

Running the script "setup_RRE_exploratory.py" will create shell files in the /etc/conda/activate.d and /etc/conda/deactivate.d folders that set and unset/reset the HHLIB and PATH variables. It will find the HHPaths.pm file and update it. And it will download the databases from the zenodo repository.
Note that the script has not been tested in a wide variety of conda environments. Use at your own risk, and let me know if it can be improved.


If you prefer to set these yourself (or if you have done a manual installation), you can follow the steps below.


## Manual install

1) Make sure the following packages/programs are installed and in your path

- HMMER version >= 3.3
- PSIPRED version >= 4.01
- python package Biopython version >= 1.76

RREFinder has not been tested for previous versions of these packages.

2) Install the HHsuite version >= 3.3  (https://github.com/soedinglab/hh-suite)

3) Make sure the following parameters are set (which is part of their recommended installation)
    
- $HHLIB is set to the folder containing the HHsuite tool
- The binary files and the script files are in your path 
(i.e. if you open a terminal and type "hhblits" or "addss.pl", both should be recognized commands)
Easiest way to do this is to modify the .bashrc file, and add a few extra lines:


        export HHLIB="/path/to/HHsuite/"
        export PATH="$HHLIB/bin:$HHLIB/scripts:$PATH"
        
**If you compiled the code from source, make sure your $HHLIB points towards the newly built HHsuite folder, e.g. in the /build folder.**
If you followed the instructions as on the GitHub of HHsuite, the skeleton used to build it remains, which also contains a script folder. Pointing to this folder will cause errors.

Then reboot the terminal or rerun the file (source ~/.bashrc)

4) Configure the HHsuite paths
Find the file HHPaths.pm in the HHsuite scripts folder 

        ##############################################################################################
        #PLEASE COMPLETE THE PATHS ... TO PSIPRED AND OLD-STYLE BLAST (NOT BLAST+) (NEEDED FOR PSIPRED)
        #our $execdir = ".../psipred/bin";         # path to PSIPRED V2 binaries
        #our $datadir = ".../psipred/data";        # path to PSIPRED V2 data files
        #our $ncbidir = ".../blast/bin";           # path to NCBI binaries (for PSIPRED in addss.pl)

Complete these paths. The $execdir and the $ncbidir need to contain the psipred and the blastpgp binaries, respectively.
The see where a binary is located use 
        
        which psipred
        which blastpgp
        
within the active environment. Use the folder, not the binary itself.
In the conda package, the two folders are likely the same.

For the $datadir, you need to find where the psipred data is stored within the conda package.
This depends on your conda setup. Example locations include: 
 
    /path/to/conda/envs/RREfinder/share/psipred_4.01/data
 
    /path/to/conda/pkgs/psipred_4.01/share/psipred_4.01/data
 
The data folder should contain seven files all named weights.dat or some variation of that name.

If you manually install PSIPRED, let the execdir point to the folder containing the psipred binary,
the datadir to the folder containing the psipred data and the ncbidir to the folder containing the legacy BLAST binary (blastpgp)
         
RREfinder is now ready for use in both modes!

## Usage guide

RREFinder accepts both protein .fasta files and .genbank files. 
Alternatively, you can use '--antismash ripp' and let the --infile point to a .final.gbk file from antiSMASH analysis. 
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

    python RRE.py -i my_infile.gbk -m exploratory project_name

You can also specify a range of other options, such as number of cores to use, or the bitscore cutoffs.
Use 'python RRE.py -h'  to see a list of options. If no options are given, they are read from the config.ini file.

## Use RREFinder to detect RREs with the uniclust30 database (Advanced)
You can also use RREFinder as a way to run the HHPred pipeline in a commandline fashion, just to detect RREs. This provides analogous results as when submitting queries online in the HHPred pipeline. In addition, RRE regions are automatically resubmitted in the same pipeline, to verify them. Running the pipeline in this way requires downloading of the uniclust30 database (https://uniclust.mmseqs.com/).

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
