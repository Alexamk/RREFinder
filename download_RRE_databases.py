# Standalone script to download the exploratory databases from the Zenodo database
# Called by setup_exploratory.py

import urllib.request
import os
import json
from itertools import product


def main():
    datafolder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    database_folder = os.path.join(datafolder, 'database')

    # First get the zenodo_record if it is not here
    zenodo_record = '3733240#.Yjn1ejwo_JU'
    url_record = f'https://zenodo.org/api/records/{zenodo_record}'
    record_file = os.path.join(datafolder, 'RRE_zenodo_record.json')
    if not os.path.isfile(record_file):
        with urllib.request.urlopen(url_record) as handle:
            text = handle.read().decode('utf-8')
        with open(record_file, 'w') as handle:
            handle.write(text)

    # Now load the json file
    with open(record_file) as handle:
        json_dict = json.load(handle)

    required_files = set(f'RRE_v7_3_{base}.{ext}' for base, ext in product(['a3m', 'cs219', 'hhm'], ['ffindex', 'ffdata']))
    database_files_already_present = os.listdir(database_folder)
    links_to_download = {}
    for available_file in json_dict['files']:
        filename = available_file['key']
        if filename in required_files and filename not in database_files_already_present:
            links_to_download[filename] = available_file['links']['self']

    for filename, link in links_to_download.items():
        print(f'Downloading {filename} from {link}')
        with urllib.request.urlopen(link) as handle:
            text = handle.read()
        with open(os.path.join(database_folder, filename), 'wb') as handle:
            handle.write(text)

if __name__ == '__main__':
    main()
