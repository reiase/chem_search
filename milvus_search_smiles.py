import sys, getopt
import numpy  as np
import time
import random
import os
from milvus import Milvus, Prepare, IndexType, Status
from rdkit import Chem
from rdkit import DataStructs

MILVUS = Milvus()
SERVER_ADDR = "192.168.1.58"
SERVER_PORT = 19522

PG_HOST = "192.168.1.85"
PG_PORT = 5420
PG_USER = "postgres"
PG_PASSWORD = "postgres"
PG_DATABASE = "postgres"


NQ_FOLDER_NAME = 'nq_npy_2048'
SE_FOLDER_NAME = 'search'
SE_FILE_NAME = '_output.txt'
BASE_FOLDER_NAME = 'test_3m'

VECTOR_DIMENSION = 2048
GT_NQ = 100

NPROBE = 64

def connect_server():
    print("connect to milvus.")
    status = MILVUS.connect(host=SERVER_ADDR, port=SERVER_PORT)
    handle_status(status=status)
    return status


# the status of milvus
def handle_status(status):
    if status.code != Status.SUCCESS:
        print(status)
        sys.exit(2)


def load_hex_vec():
    filenames = os.listdir(NQ_FOLDER_NAME)
    filenames.sort()
    data = []
    for filename in filenames:
        filename = NQ_FOLDER_NAME + '/' + filename
        data = np.load(filename)
        data = data.tolist()
        vectors = []
        for d in data:
            vectors.append(bytes.fromhex(d))
    return vectors


def connect_postgres_server():
    try:
        conn = psycopg2.connect(host=PG_HOST, port=PG_PORT, user=PG_USER, password=PG_PASSWORD, database=PG_DATABASE)
        print("connect the database!")
        return conn
    except:
        print("unable to connect to the database")


def search_smi_in_pg(cur, table_name, ids):
    try:
        sql = "select smiles from " + table_name+ " where ids = '" + str(ids) + "';"
        cur.execute(sql)
        rows = cur.fetchall()
        return str(rows[0][0])
    except:
        print("search faild!")


def search_milsmi_in_pg(cur, table_name, ids):
    try:
        sql = "select smiles from " + table_name+ " where milvus_ids = '" + str(ids) + "';"
        cur.execute(sql)
        rows = cur.fetchall()
        return str(rows[0][0])
    except:
        print("search faild!")


def search_milids_in_pg(cur, table_name, ids):
    try:
        sql = "select ids from " + table_name+ " where milvus_ids = '" + str(ids) + "';"
        cur.execute(sql)
        rows = cur.fetchall()
        return str(rows[0][0])
    except:
        print("search faild!")


def save_re_to_file(table_name, results):
    if not os.path.exists(SE_FOLDER_NAME):
        os.mkdir(SE_FOLDER_NAME)
    file_name = SE_FOLDER_NAME + '/' + table_name + SE_FILE_NAME
    conn = connect_postgres_server()
    cur = conn.cursor()
    # print(file_name)
    with open(file_name, 'w') as f:
        for i in range(len(results)):
            for j in range(len(results[i])):
                ids = search_milids_in_pg(cur, table_naem, results[i][j].id)
                smiles = search_milsmi_in_pg(cur, table_naem, results[i][j].id)
                line = ids + ',' + smiles + ',' + str(results[i][j].id) + ',' + str(results[i][j].distance)
                f.write(line + '\n')
    cur.close()
    conn.close()


def get_smi_in_pg(table_naem, ids):
    conn = connect_postgres_server()
    cur = conn.cursor()
    index = search_smi_in_pg(cur, table_naem, ids)
    cur.close()
    conn.close()
    return index


def search_ids_smi_list(table_name, topk, ids, smiles):
    rand = None
    query_list = []

    if ids:
        smiles = get_smi_in_pg(table_naem, ids)
    mols = Chem.MolFromSmiles(smiles)
    fp = Chem.RDKFingerprint(mols, fpSize=VECTOR_DIMENSION)
    hex_fp = DataStructs.BitVectToFPSText(fp)
    # print(hex_fp)
    vec = bytes.fromhex(hex_fp)
    query_list.append(vec)

    print("table name:", table_name, "query list:", len(query_list), "topk:", topk)
    time_start = time.time()
    status, results = MILVUS.search(collection_name=table_name, query_records=query_list, top_k=topk, params={})
    time_end = time.time()
    time_cost = time_end - time_start
    print("time_search = ", time_cost)
    print(status,results)

    time_start = time.time()
    save_re_to_file(table_name, results)
    time_end = time.time()
    time_cost = time_end - time_start
    print("time_save = ", time_cost)


def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "hst:q:k:n:i:",
            ["help", "search", "table=", "topk=", "nprobe=", "smi=", "id="],
        )
    except getopt.GetoptError:
        print("Usage: test.py --table <tablename> -k <topk> --smi <smiles>")
        sys.exit(2)
    nq = 0
    nprobe = NPROBE
    for opt_name, opt_value in opts:
        if opt_name in ("-h", "--help"):
            print("test.py --table <tablename> -k <topk> --smi <smiles> or test.py --table <tablename> -k <topk> --ids <ids_txt>")
            sys.exit()
        elif opt_name in ("-t", "--table"):
            table_name = opt_value
        elif opt_name in ("-q", "--nq"):
            nq = int(opt_value)
        elif opt_name in ("-k", "--topk"):
            topk = int(opt_value)
        elif opt_name in ("-n", "--nprobe"):
            nprobe = int(opt_value)
        elif opt_name in ("--smi"):
            ids = ''
            smi = opt_value
            connect_server()
            search_ids_smi_list(table_name, topk, ids, smi)  # test.py --table <tablename> -k <topk> --smi <smiles>
        elif opt_name in ("-i", "--id"):
            smi = ''
            ids = opt_value
            connect_server()
            search_ids_smi_list(table_name, topk, ids, smi)  # test.py --table <tablename> -k <topk> --ids <ids_txt>


if __name__ == '__main__':
    main()
