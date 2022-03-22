# This script is meant to ease the setup of exploratory mode of RREFinder
# Unfortunately, the addss.pl script required from the HHSuite is not supported 
# as a standard part of the conda package. It is downloaded along with it though.
# This script will do three things:

# 1) modify the conda activate script, as found in 
# $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh and 
# $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh
# To set the HHLIB variable, and add the scripts part of the HHSuite module to the path
# These are unset again when the environment is deactivates

# 2) Modify the HHPaths.pm scripts so that addss.pl can find the required binaries

# 3) Download the databases, if necessary

# NOTE: These actions might interfere with other custom conda operations. Use at your own risk

import os
import download_databases

# Step 1

activate_text = '''#!/bin/sh

export HHLIB="$CONDA_PREFIX"
export RRE_OLD_PATH="$PATH"
export PATH="$PATH:$HHLIB/scripts"'''

deactivate_text = '''#!/bin/sh

unset HHLIB
export PATH="$RRE_OLD_PATH"
unset OLD_PATH'''

conda_path = os.environ['CONDA_PREFIX']
activate_folder = os.path.join(conda_path, 'etc', 'conda', 'activate.d')
deactivate_folder = os.path.join(conda_path, 'etc', 'conda', 'deactivate.d')
filename = 'RRE_envs.sh'
for folder, text in zip((activate_folder, deactivate_folder), (activate_text, deactivate_text)):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    with open(os.path.join(folder, filename), 'w') as handle:
        handle.write(text)

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
download_databases.main()

print('Setup done. Please deactivate and reactivate the environment')
