import canalizing_function_toolbox_v1_9 as can
import json
import numpy as np
import scipy.stats as st

def load_data():
    f = open("./results/SLURM_ID1.txt")
    dataStr = f.read().split("\n")
    dataStr = [line.split("\t") for line in dataStr]
    return json.loads(dataStr[2][1]),json.loads(dataStr[3][1])

def count_total(loops):
    totalCounts = [0 for i in range(len(loops[0]))]
    for loop_length in FBLs:
        for net_num,loop_count in enumerate(loop_length):
            totalCounts[net_num] = totalCounts[net_num] + loop_count
    return totalCounts

FBLs,rndFBLs = load_data()
nsims = len(rndFBLs)

totalCounts = count_total(FBLs)

rndCounts = []
for sim in rndFBLs:
    simCounts = count_total(sim)
    for count in simCounts:
        rndCounts.append(count)

print("Natural: ", st.norm.interval(alpha=0.99, loc=np.mean(totalCounts), scale=st.sem(totalCounts)))
print("Generated: ", st.norm.interval(alpha=0.99, loc=np.mean(rndCounts), scale=st.sem(rndCounts)))
