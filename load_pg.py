import getopt
import os
import sys
import psycopg2
import time
from functools import reduce
import numpy as np


FILE_IDS_SMILES = '/mnt/test/zink_data/ids_smiles.csv'

# FILE_SMILES = '/mnt/workspace/data/pub_0_99/can_smiles_smile'
# FILE_SMILES = '/data/workspace/apptec/1B_data/out_test/out_smiles'
FILE_SMILES = '/mnt/out/out_smiles'
FILE_IDS = '/mnt/out/out_ids'

TO_PG = True
PG_HOST = "localhost"
PG_PORT = 5432
PG_USER = "zilliz"
PG_PASSWORD = "zilliz"
PG_DATABASE = "milvus"


def hex_to_pg(PG_TABLE):
    filenames_smiles = os.listdir(FILE_SMILES)
    filenames_smiles.sort()

    filenames_ids = os.listdir(FILE_IDS)
    filenames_ids.sort()

    count = 0
    if TO_PG:
        conn = connect_postgres_server()
        cur = conn.cursor()
        create_pg_table(conn, cur, PG_TABLE)

    for filename in filenames_ids:
        names = load_smiles(filenames_smiles[count])
        ids_vec = load_ids(filenames_ids[count])

        with open(FILE_IDS_SMILES, 'w') as f:
            for i in range(len(names)):
                f.write(str(ids_vec[i]) + ',' + names[i] + '\n')
        if TO_PG:
            copy_data_to_pg(conn, cur, PG_TABLE)
        count += 1
    cur.close()
    conn.close()


def connect_postgres_server():
    try:
        conn = psycopg2.connect(host=PG_HOST, port=PG_PORT, user=PG_USER, password=PG_PASSWORD,
                                database=PG_DATABASE)
        print("connect the database!")
        return conn
    except:
        print("unable to connect to the database")


def create_pg_table(conn, cur, PG_TABLE_NAME):
    sql = "CREATE TABLE " + PG_TABLE_NAME + " (ids bigint, smiles text);"
    print(sql)
    try:
        cur.execute(sql)
        conn.commit()
        print("create postgres table!")
    except:
        print("can't create postgres table")
        # sys.exit()


def copy_data_to_pg(conn, cur, PG_TABLE_NAME):
    sql = "copy " + PG_TABLE_NAME + " from '" + FILE_IDS_SMILES + "' with CSV delimiter ',';"
    print(sql)
    try:
        cur.execute(sql)
        conn.commit()
        print("copy data to pg sucessful!")
    except:
        print("faild  copy!")
        # sys.exit()


def load_smiles(file):
    file = FILE_SMILES + '/' + file
    print("smiles_file:",file)
    smiles = []
    for line in open(file, 'r'):
        data = line.strip('\n')
        # smiles.append(data.encode())
        smiles.append(data)
    return smiles


def load_ids(file):
    file = FILE_IDS + '/' + file
    print("ids_file:",file)
    ids = []
    for line in open(file, 'r'):
        data = line.strip('\n')
        ids.append(int(data[4:]))
    return ids


def main(argv):
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "t:l",
            ["load","table="],
        )
        # print(opts)
    except getopt.GetoptError:
        print("Usage: load_pg.py --table <table_name> -l")
        sys.exit(2)

    for opt_name, opt_value in opts:
        if opt_name in ("-t", "--table"):
            PG_TABLE_NAME = opt_value
        elif opt_name in ("-l", "--load"):
            hex_to_pg(PG_TABLE_NAME)
        else:
            print("wrong parameter")
            sys.exit(2)

if __name__ == "__main__":
    main(sys.argv[1:])
