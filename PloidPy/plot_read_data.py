import numpy as np
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt


def plot_joint_dist(a, savefile, max_cutoff = 0.95):
    uniq = np.unique(a[:,1], return_counts = True)
    sns.set_context("talk")
    cumdens = np.cumsum(uniq[1] / np.sum(uniq[1]))
    maxv = uniq[0][np.where(cumdens >= max_cutoff)[0][0]]
    a0 = a[a[:,1] < maxv]
    (sns.jointplot(a0[:,0], a0[:,1], kind = "hex", color="blue")
     .set_axis_labels("Minor Allele Coverage", "Total Read Coverage"))
    plt.tight_layout()
    plt.savefig(savefile)
