import random
import os 
import math
import urllib2

def read_csv(infile,type):
    all_data = {}
    ripp_nr = 1
    with open(infile) as handle:
        # Skip header
        handle.readline()
        for line in handle:
            tabs = line.strip().split('\t')
            ripptype = tabs[5]
            if ripptype != type.lower(): # No hybrid clusters
                continue
            cluster_name = '%s_%i' %(ripptype,ripp_nr)
            all_data[cluster_name] = tabs           
            ripp_nr += 1
    return all_data
            

def download(data,folder):
    for name in data:
        print('Downloading %s' %name)
        url = data[name][-1]
        path = os.path.join(folder,name + '.gbk')
        urlobject = urllib2.urlopen(url)
        text = urlobject.read()
        with open(path,'w') as handle:
            handle.write(text)   
    

if __name__ == '__main__':
    all_data = {}
    infolder = '../antismash_db_ripps/'
    download_folder = os.path.join(infolder,'downloads')
    files = [i for i in os.listdir(infolder) if os.path.isfile(os.path.join(infolder,i))]
    for f in files:
        ripptype = f.partition('.')[0]
        file_path = os.path.join(infolder,f)
        data = read_csv(file_path,ripptype)
        all_data[ripptype] = data
        
    # Always download the smaller classes
    to_download = {}
    max_download = download_remaining = 5000
    classes_remaining = all_data.keys()
    while download_remaining > 0 and len(classes_remaining) > 0:
        download_remaining = max_download - len(to_download)
        classes_remaining.sort(key=lambda x: len(all_data[x]))
        cl = classes_remaining[0]
        rem_per_class = int(math.ceil(float(download_remaining) / len(classes_remaining)))
        if len(all_data[cl]) <= rem_per_class:
            print('Downloading all (%i) for class %s' %(len(all_data[cl]),cl))
            to_download.update(all_data[cl])
        else:
            print('Downloading %i for class %s' %(rem_per_class,cl))
            sampled = random.sample(all_data[cl],rem_per_class)
            for s in sampled:
                to_download[s] = all_data[cl][s]
        if classes_remaining == []:   
            break
        classes_remaining = classes_remaining[1:]
    download(to_download,download_folder)
            