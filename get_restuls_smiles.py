import getopt
import sys
import psycopg2
import pandas as pd
import time
import numpy as np
import os

SE_FOLDER_NAME = 'search'
SE_FILE_NAME = '_output.txt'

CM_FOLDER_NAME = 'compare'
CM_GET_LOC_NAME = '_compare.txt'
FILE_SMILES = '/data/workspace/out_test/out_smiles'


def save_mols(rand, ids, nq, top_k, table_name, dis, nprobe):
    filename = CM_FOLDER_NAME + '/' + table_name + '_' + str(nprobe) + CM_GET_LOC_NAME
    filenames_smiles = os.listdir(FILE_SMILES)
    filenames_smiles.sort()
    count = 0
    with open(filename, 'w') as f:
        for i in ids:
            loca = int(i[1:5])
            offset = int(i[5:11])
            smiles_list = load_smiles(filenames_smiles[loca])
            smiles = smiles_list[offset]
            f.write(str(rand[count]) + ' ' + str(ids[count]) + ' ' + str(smiles) + ' ' + str(dis[count]) + '\n')
            count += 1
            if count%top_k==0:
                f.write('\n')


def load_smiles(file):
    file = FILE_SMILES + '/' + file
    # print("smiles_file:",file)
    smiles = []
    for line in open(file, 'r'):
        data = line.strip('\n')
        # smiles.append(data.encode())
        smiles.append(data)
    return smiles


def load_search_out(table_name, nprobe,ids=[], rand=[], distance=[]):
    file_name = SE_FOLDER_NAME + '/' + table_name + '_' + nprobe + SE_FILE_NAME
    print(file_name)
    nq = 0
    with open(file_name, 'r') as f:
        for line in f.readlines():
            data = line.split()
            if data:
                rand.append(data[0])
                ids.append(data[1])
                distance.append(data[2])
            else:
                nq += 1
    return rand, ids, distance, nq


def main(argv):
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "ht:n:g",
            ["help", "table=", "nprobe=", "gen"]
        )
    except getopt.GetoptError:
        print("Usage: python get_results_smiles.py --table=<table_name> -n <nprobe> -g")
        sys.exit(2)

    for opt_name, opt_value in opts:
        if opt_name in ("-h", "--help"):
            print("python get_results_smiles.py --table=<table_name> -n <nprobe> -g")
            sys.exit()
        elif opt_name in ("-n", "--nprobe"):
            nprobe = opt_value
        elif opt_name in ("-t", "--table"):
            table_name = opt_value
        elif opt_name in ("-g", "--gen"):
            if not os.path.exists(CM_FOLDER_NAME):
                os.mkdir(CM_FOLDER_NAME)
            rand, ids, dis, nq = load_search_out(table_name, nprobe)
            top_k = int(len(rand) / nq)
            print("nq:", nq, "top_k:", top_k)
            save_mols(rand, ids, nq, top_k, table_name, dis, nprobe)


if __name__ == "__main__":
    main(sys.argv[1:])

