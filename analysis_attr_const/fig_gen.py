import matplotlib.pyplot as plt
import numpy as np
import load_database11 as db

rewire_string = open("analysis_attr_const/results/results_operation-activation_prop_nsim100_SLURM_ID1.txt").read().split("\n")
activation_string = open("analysis_attr_const/results/results_operation-activation_prop_nsim100_SLURM_ID1.txt").read().split("\n")
canalization_string = open("analysis_attr_const/results/results_operation-canalization_nsim100_SLURM_ID1.txt").read().split("\n")
real_string = open("analysis_attr_const/results/results_operation-real_nsim100_SLURM_ID1.txt").read().split("\n")

def read_data(data_string, nsims):
    data_arr = []
    title_arr = []
    for i,line in enumerate(data_string):
        if i < 4 or i > 7:
            continue
        lineSplit = line.split("\t")
        title_arr.append(lineSplit[0])
        nums = list(map(float,lineSplit[1:]))
        data_arr.append(nums) 
    return data_arr,title_arr
    
        
def plot_data(data_matrix, title, real_data):
    f,ax = plt.subplots(figsize=(12,4))
    for i in range(data_matrix.shape[0]):
        ax.violinplot(data_matrix[i], positions=[i], showextrema=False)
    ax.set_ylabel(title)
    ax.set_title("Total Attractors - Preserved Canalization vs Real Networks")
    ax.plot(real_data, "o", color='k')
    
def sort_p_shuffle(null_data, real_data):
    null_matrix = np.reshape(null_data, (36,100))
    counts = []
    for i,real in enumerate(real_data):
        counts.append(np.sum(null_matrix[i] > real))
    indices = sorted(range(36),key=lambda x: counts[x])
    out_matrix = np.array([null_matrix[indices[i]] for i in range(36)])
    out_data = np.array([real_data[indices[i]] for i in range(len(real_data))])
    return out_matrix,out_data

            

        
real_data = read_data(real_string, 1)[0]
activation_data,activation_titles = read_data(canalization_string, 100)
idx = 1
activation_matrix,real_data_sub = sort_p_shuffle(activation_data[idx], real_data[idx])

plot_data(activation_matrix, activation_titles[idx], real_data_sub)

