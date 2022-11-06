import pandas as pd
import pysam
import numpy as np
import gzip
import Levenshtein
import multiprocessing 
import time
import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--barcodes", help="Barcodes file", required=True)
parser.add_argument("-i", "--inputFile", help="Rawdata to demultiplex", required=True)
parser.add_argument("-o", "--outputPrefix", help="Output directory", required=True)
parser.add_argument("-t", "--threads", help="Number of threads", required=False, default=10)
args = parser.parse_args()

barcodes = pd.read_csv(args.barcodes ,header=None).values.T[0]
filename = args.inputFile
out_prefix = args.outputPrefix
filepaths  = [out_prefix + barcode + ".fastq.gz" for barcode in barcodes]
thread_num = int(args.threads)


def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

def get_cellbarcode(entry,barcodelength = 8):
    ind = entry.name.split("_")[0]#.split(",")
    ind1 = ind[0:barcodelength]
    #ind1 = ind[0][0:barcodelength]
    #ind2 = rev_compl(ind[1][-barcodelength:])
    barcode_cell = ""

    ind1_distance = np.array([Levenshtein.distance(ind1,i) for i in barcodes])
    #ind2_distance = np.array([Levenshtein.distance(ind2,i) for i in barcodes])
    if barcodes[ind1_distance == 0].size>0:
        barcode_cell = ind1
    #elif barcodes[ind2_distance == 0].size>0:
    #    barcode_cell = ind2
    elif barcodes[ind1_distance == 1].size>0:
        barcode_cell = barcodes[ind1_distance == 1][0]
    #elif barcodes[ind2_distance == 1].size>0:
    #    barcode_cell = barcodes[ind2_distance == 1][0]
    return [entry,barcode_cell]

def write_entry(result):
    barcode_cell = result[1]
    entry = result[0]
    if (barcode_cell!=""):
        file_id = np.where(barcodes==barcode_cell)[0][0]
        all_files[file_id].write(str(entry) + '\n')

all_files = []

for path in filepaths:
    file = gzip.open(path, 'wt',compresslevel=1)
    all_files.append(file)

start = time.time()
print("Process with " + str(thread_num) + " processes!")
with pysam.FastxFile(filename) as fin:
    pool = multiprocessing.Pool(thread_num,maxtasksperchild=10000)
    n = 0 
    for entry in tqdm.tqdm(fin,desc="Demultiplexing: ",mininterval=2):
        last = pool.apply_async(get_cellbarcode, args=(entry,), callback=write_entry)
        if len(pool._cache) > 1e4:
            last.wait()

    pool.close()
    pool.join()
    
for file in all_files:
    file.close()