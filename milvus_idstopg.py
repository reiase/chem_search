import getopt
import os
import sys
import psycopg2
import time
from functools import reduce
import numpy as np


PG_HOST = "192.168.1.85"
PG_PORT = 5420
PG_USER = "postgres"
PG_PASSWORD = "postgres"
PG_DATABASE = "postgres"


FILE_IDS_SMILES = '/data/pubchem/out_data/test_ids_smiles.csv'
FILE_SMILES = '/data/pubchem/out_data/out_smiles'
FILE_IDS = '/data/pubchem/out_data/out_ids'


def connect_postgres_server():
    try:
        conn = psycopg2.connect(host=PG_HOST, port=PG_PORT, user=PG_USER, password=PG_PASSWORD,
                                database=PG_DATABASE)
        print("connect the database!")
        return conn
    except:
        print("unable to connect to the database")


def create_pg_table(conn, cur, PG_TABLE_NAME):
    sql = "CREATE TABLE " + PG_TABLE_NAME + " (milvus_ids bigint, ids text, smiles text);"
    # print(sql)
    try:
        cur.execute(sql)
        conn.commit()
        print("create postgres table!")
    except:
        print("can't create postgres table")
        # sys.exit()


def copy_data_to_pg(conn, cur, PG_TABLE_NAME):
    sql = "copy " + PG_TABLE_NAME + " from '" + FILE_IDS_SMILES + "' with CSV delimiter ',';"
    # print(sql)
    try:
        cur.execute(sql)
        conn.commit()
        print("copy data to pg successfully!")
    except:
        print("faild  copy!")


def create_index_pg(conn, cur, PG_TABLE_NAME):
    sql = "CREATE INDEX milids_idx ON " + PG_TABLE_NAME +"(milvus_ids); CREATE INDEX ids_idx ON " + PG_TABLE_NAME +"(ids);"
    # print(sql)
    try:
        cur.execute(sql)
        conn.commit()
        print("create index with pg successfully!")
    except:
        print("faild create index!")


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
        ids.append(data)
    return ids


def ids_to_pg(conn, cur, table_name):
    filenames_ids = os.listdir(FILE_IDS)
    filenames_ids.sort()

    filenames_smiles = os.listdir(FILE_SMILES)
    filenames_smiles.sort()

    count = 0
    for filename in filenames_ids:
        if count>10:
            break
        ids = load_ids(filename)
        smiles = load_smiles(filenames_smiles[count])
        count += 1

        milvus_ids = []
        for i in range(len(ids)):
            location = '8' + '%05d'%count  + '%07d'%i
            milvus_ids.append(location)

        with open(FILE_IDS_SMILES, 'w') as f:
            for i in range(len(ids)):
                f.write(milvus_ids[i] + ',' + ids[i] + ',' + smiles[i] + '\n')
        copy_data_to_pg(conn, cur, table_name)



def main(argv):
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "t:h",
            ["hex","table="],
        )
        # print(opts)
    except getopt.GetoptError:
        print("Usage: load_vec_to_milvus.py -n <npy>  -c <csv> -f <fvecs> -b <bvecs>")
        sys.exit(2)

    for opt_name, opt_value in opts:
        if opt_name in ("-t", "--table"):
            PG_TABLE_NAME = opt_value

            conn = connect_postgres_server()
            cur = conn.cursor()

            create_pg_table(conn, cur, PG_TABLE_NAME)
            ids_to_pg(conn, cur, PG_TABLE_NAME)
            create_index_pg(conn, cur, PG_TABLE_NAME)

            cur.close()
            conn.close()

        else:
            print("wrong parameter")
            sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])