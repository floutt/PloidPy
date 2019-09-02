import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_read_hexbin(arr, savefile):
    plt.hexbin(arr[:,0], arr[:, 1])
    plt.colorbar()
    plt.xlabel("Minor Allele Coverage")
    plt.ylabel("Total Coverage")
    plt.savefig(savefile)


def plot_read_histogram(arr, savefile, bins = 30):
    weights = np.ones_like(arr[:,0])/ float(len(arr[:,0]))
    plt.subplot(211)
    plt.hist(arr[:,0], bins = bins, weights = weights)
    plt.ylabel("Frequency")
    plt.xlabel("Coverage")
    plt.title("Minor Allele Coverage")
    plt.subplot(212)
    plt.hist(arr[:,1], bins = bins, weights = weights)
    plt.ylabel("Frequency")
    plt.xlabel("Coverage")
    plt.title("Total Coverage")
    plt.savefig(savefile)
