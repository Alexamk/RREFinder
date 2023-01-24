# This script is meant to ease the setup of exploratory mode of RREFinder
# Unfortunately, the addss.pl script required from the HHSuite is not supported 
# as a standard part of the conda package. It is downloaded along with it though.
# This script will do three things:

# 1) modify the conda (de/)activate scripts, as found in 
# $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh and 
# $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
# To set the HHLIB variable, and add the scripts part of the HHSuite module to the path
# These are unset again when the environment is deactivates

# 2) Modify the HHPaths.pm scripts so that addss.pl can find the required binaries

# 3) Download the databases, if necessary

# NOTE: These actions might interfere with other custom conda operations. Use at your own risk

import os
import download_RRE_databases
import shutil

# Step 1

conda_path = os.environ['CONDA_PREFIX']
activate_folder = os.path.join(conda_path, 'etc', 'conda', 'activate.d')
deactivate_folder = os.path.join(conda_path, 'etc', 'conda', 'deactivate.d')
env_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'envs')

activate_filename = 'RRE_activate.sh'
deactivate_filename = 'RRE_deactivate.sh'

for folder, filename in zip((activate_folder, deactivate_folder), (activate_filename, deactivate_filename)):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    shutil.copy(os.path.join(env_folder, filename), os.path.join(folder, filename))

# Step 2
hhpaths_file = os.path.join(conda_path, 'scripts', 'HHPaths.pm')
if not os.path.isfile(hhpaths_file):
    raise ValueError(f'HHPaths file {hhpaths_file} not found. Exiting')
# Read through the file, replace the necessary three lines

replacement_dict = {'our $execdir': 'our $execdir = $ENV{"HHLIB"}."/bin";                            # path to PSIPRED V2 binaries\n',
                    'our $datadir': 'our $datadir = $ENV{"HHLIB"}."/share/psipred_4.01/data";        # path to PSIPRED V2 data files\n',
                    'our $ncbidir': 'our $ncbidir = $ENV{"HHLIB"}."/bin";                            # path to NCBI binaries (for PSIPRED in addss.pl)\n'}

text_out = ''
with open(hhpaths_file) as handle:
    for line in handle:
        if line[0:12] in replacement_dict:
            text_out += replacement_dict[line[0:12]]
        else:
            text_out += line

# Make a backup of the old paths file just in case, unless it already exists
backup_file = os.path.join(conda_path, 'scripts', 'HHPaths_backup.pm')
if not os.path.isfile(backup_file):
    os.rename(hhpaths_file, backup_file)
with open(hhpaths_file, 'w') as handle:
    handle.write(text_out)

# Step 3
download_RRE_databases.main()

print('Setup done. Please deactivate and reactivate the environment')
