import os
import sys
import pandas as pd
import glob
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def read_log(file_path):
    log_data = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        step = None
        for line in lines:
            line = line.strip()
            if line.startswith("[INFO] Step") and "Searching among" in line:\n                step = line.split(":")[0].replace("[INFO] ", "")
                log_data[step] = {}
                
                if "pending sequences" in line:
                    parts = line.split("with")
                    dbsize = parts[0].split("among")[1].strip().split(" ")[0]
                    nb_query = parts[1].split("pending sequences")[0].strip()
                    log_data[step]['nb_query'] = int(nb_query)
                    log_data[step]['dbsize'] = int(dbsize)
                
            elif step and line.startswith("[INFO] Elapsed time:"):
                elapsed_time = line.split(":")[1].strip().split(" ")[0]
                log_data[step]['elapsed_time'] = float(elapsed_time)
    return log_data

# Find all .log files in ../runs/*/
log_files = []
runs_dir = '../runs'
if os.path.exists(runs_dir):
    for log_path in glob.glob(os.path.join(runs_dir, '*', '*.log')):
        log_files.append(log_path)
        print(f"Found: {log_path}")
else:
    print(f"Warning: {runs_dir} directory not found")

if not log_files:
    print("No log files found. Run some Brownaming analyses first!")
    sys.exit(1)

data_file_data = {
    'nb_query': [],
    'dbsize': [],
    'time': []
}

print(f"\\nProcessing {len(log_files)} log files...")
for log_file in log_files:
    log_data = read_log(log_file)
    last_elapsed_time = 0.0
    for step in log_data:
        if step.startswith("Step"):
            nb_query = log_data[step]['nb_query'] if 'nb_query' in log_data[step] else None
            dbsize = log_data[step]['dbsize'] if 'dbsize' in log_data[step] else None
            elapsed_time = log_data[step]['elapsed_time'] if 'elapsed_time' in log_data[step] else None
            if nb_query and dbsize and elapsed_time:
                runtime = round(float(elapsed_time) - float(last_elapsed_time), 2)
                data_file_data['nb_query'].append(nb_query)
                data_file_data['dbsize'].append(dbsize)
                data_file_data['time'].append(runtime)
                last_elapsed_time = elapsed_time

    
df = pd.DataFrame(data_file_data)
df.to_csv('data_file.tsv', sep='\t', index=False)
print(f"\\nDataset created with {len(df)} entries")
print(f"Saved to: data_file.tsv")
print("\\nYou can now train the model with: python train_model.py")