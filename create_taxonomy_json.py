import json
import os
import requests
from collections import defaultdict
import numpy as np
import pandas as pd
import pickle
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score


# Read local_db_path from environment variable first, then from config.json
local_db_path = os.environ.get("LOCAL_DB_PATH")

if not local_db_path:
    CONFIG = {}
    config_file = os.path.join(os.path.dirname(__file__), "config.json")
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            CONFIG = json.load(f)
        local_db_path = CONFIG.get("local_db_path", None)

if not local_db_path:
    print("[ERROR] 'local_db_path' not set in LOCAL_DB_PATH environment variable or config.json")
    exit()
    
nodes_path = os.path.join(local_db_path, "taxonomy", "nodes.dmp")
if not os.path.isfile(nodes_path):
    print(f"[ERROR] nodes.dmp not found at {nodes_path}. Local database is not correctly set up. Please retry to execute create_local_db.sh")
    exit()
    
names_path = os.path.join(local_db_path, "taxonomy", "names.dmp")
if not os.path.isfile(names_path):
    print(f"[ERROR] names.dmp not found at {nodes_path}. Local database is not correctly set up. Please retry to execute create_local_db.sh")
    exit()    

parent = {}
rank = {}
children = defaultdict(set)
with open(nodes_path, "r") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]      
        if len(parts) < 3:
            continue
        taxid = int(parts[0])
        par = int(parts[1])
        r = parts[2]
        parent[taxid] = par
        rank[taxid] = r
        children[par].add(taxid)

children_serializable = {k: list(v) for k, v in children.items()}

taxid_to_name = {}
with open(names_path, "r") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 4:
            continue
        taxid = int(parts[0])
        name_txt = parts[1]
        name_class = parts[3]
        if name_class == "scientific name":
            taxid_to_name[taxid] = name_txt

parent_path = os.path.join(local_db_path, "taxonomy", "parent.json")
with open(parent_path, 'w') as f:
    json.dump(parent, f)

rank_path = os.path.join(local_db_path, "taxonomy", "rank.json")
with open(rank_path, 'w') as f:
    json.dump(rank, f)

children_path = os.path.join(local_db_path, "taxonomy", "children.json")
with open(children_path, 'w') as f:
    json.dump(children_serializable, f)

taxid_to_name_path = os.path.join(local_db_path, "taxonomy", "taxid2scientific_name.json")
with open(taxid_to_name_path, 'w') as f:
    json.dump(taxid_to_name, f)
