import threading
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import os
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw
import math
import numpy as np
import gc
import linecache


file_length = 1000000
total_length = 60000000
vec_dim = 2048

FILE_NAME = '/mnt/pubchem/pub_1kw.smi'
OUT = '/mnt/pubchem/out'
OUT_SMILES = 'out_smiles'
OUT_IDS = 'out_ids'
OUT_NPY = 'out_npy'


def thread_runner(smiles, ids, filename, cycle):
    #thread_num = math.ceil(len(smiles)/file_length)
    thread_num = math.ceil(total_length/file_length)
    print("thread_num:", thread_num)
    with ProcessPoolExecutor(thread_num) as executor:
        for i in range(thread_num):
            try:
                executor.submit(get_smiles_fp, smiles[i*file_length:i*file_length+file_length], ids[i*file_length:i*file_length+file_length], i+cycle*thread_num)
            except:
                executor.submit(get_smiles_fp, smiles[i*file_length:len(smiles)], ids[i*file_length:len(ids)], i+cycle*thread_num)


def get_smiles_fp(smiles, ids, num):
    hex_fps = []
    new_ids = []
    new_smiles = []
    count = 0
    for smile in smiles:
        m = Chem.MolFromSmiles(smile)
        if m is None:
            count = count + 1
            continue
        # the method for generating fingerprints - morgan_fp or rdk_fp
        # fp2 = AllChem.GetMorganFingerprintAsBitVect(m, 2, vec_dim)
        fp2 = Chem.RDKFingerprint(m, fpSize=vec_dim)
        hex_fp = DataStructs.BitVectToFPSText(fp2)
        hex_fps.append(hex_fp)
        new_ids.append(ids[count])
        new_smiles.append(smile)
        count = count + 1
    # print(len(hex_fps),len(new_ids),len(new_smiles))
    hex_fps = np.array(hex_fps)
    np.save(OUT + '/' + OUT_NPY + '/' + "%05d"%num +'.npy', hex_fps)
    save_file(new_smiles, OUT + '/' + OUT_SMILES + '/' + "%05d"%num +'.smi')
    save_file(new_ids, OUT + '/' + OUT_IDS + '/' + "%05d"%num +'.txt')
    del hex_fps
    del new_smiles
    del new_ids
    gc.collect()


def save_file(news, fname):
    print(len(news),fname)
    with open(fname,'a') as f:
        for new in news:
            f.write(new + '\n')


def get_files_fp(file_name):
    # cycle = 0
    cycle = 11
    num = 0
    smiles=[]
    ids=[]
    with open(file_name, "r") as infile:
        for line in infile:
            num += 1
            if num < 66000000:
                continue
            else:
                num = 1
            parts = line.split()
            try:
                smiles.append(parts[0])
                ids.append(parts[1])
            except :
                print(parts)
            if num >= total_length:
                print("cycle:",cycle, smiles[0], ids[0])
                thread_runner(smiles, ids, file_name, cycle)
                num = 0
                cycle += 1
                smiles=[]
                ids=[]
    if num != 0:
        print("cycle:",cycle, smiles[0], ids[0])
        thread_runner(smiles, ids, file_name, cycle)


def main():
    if not os.path.exists(OUT):
        os.mkdir(OUT)
        os.mkdir(OUT + '/' + OUT_NPY)
        os.mkdir(OUT + '/' + OUT_SMILES)
        os.mkdir(OUT + '/' + OUT_IDS)
        
    time1 = time.time()
    get_files_fp(FILE_NAME)
    time2 = time.time()
    print("total time: ", time2 - time1)


if __name__ == "__main__":
    main()
