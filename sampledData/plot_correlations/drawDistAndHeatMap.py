import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import numpy as np
import os

def draw_histogram(data, data2, bin):
    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 12}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

    # Create a figure with two subplots (1 row, 2 columns)
    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots

    # Set common properties for both plots
    plotLabels = {'title': '', 'xlabel': 'SFE $(mJ/m^2)$', 'ylabel': 'Probability Density'}
    #colors = sns.color_palette("husl")

    # Plotting on the first subplot
    axs[0].grid(True)
    axs[0].set_title(f"{plotLabels['title']}")
    axs[0].set_xlabel(f"{plotLabels['xlabel']}")
    axs[0].set_ylabel(f"{plotLabels['ylabel']}")
    #axs[0].hist(data, facecolor=colors[0], edgecolor='k', density=True, bins=bin)
    axs[0].hist(data, density=True, bins=bin)
    axs[0].set_title('Original')

    # Plotting on the second subplot
    axs[1].grid(True)
    axs[1].set_title(f"{plotLabels['title']}")
    axs[1].set_xlabel(f"{plotLabels['xlabel']}")
    axs[1].set_ylabel(f"{plotLabels['ylabel']}")
    #axs[1].hist(data2, facecolor=colors[1], edgecolor='k', density=True, bins=bin)  # Changed color for distinction
    axs[1].hist(data2, density=True, bins=bin)  # Changed color for distinction
    axs[1].set_title('Sampled (n=1000)')

    # Save figure 
    figName = "histogram.png"
    figFolder = './figures'
    if not os.path.exists(figFolder):
        os.makedirs(figFolder)  # Changed to os.makedirs for better practice
    plt.savefig(f'{figFolder}/{figName}')
    plt.close()

def drawHeatMap(data, data2):
    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 12}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

    # Create a figure with two subplots (1 row, 2 columns)
    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots
    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs[0].set_xticks(xTickPos)
    axs[0].set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs[0].set_yticks(yTickPos)
    axs[0].set_yticklabels(labels=yTicks)
    axs[0].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[0].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[0].set_title('Original')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = data
    c = axs[0].pcolormesh(X, Y, Z, cmap=cm.viridis)
    fig.colorbar(c, ax=axs[0])

    x = np.arange(data2.shape[1])
    y = np.arange(data2.shape[0])
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs[1].set_xticks(xTickPos)
    axs[1].set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs[1].set_yticks(yTickPos)
    axs[1].set_yticklabels(labels=yTicks)
    axs[1].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs[1].set_title('Sampled (n=1000)')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = data2
    c = axs[1].pcolormesh(X, Y, Z, cmap=cm.viridis)
    fig.colorbar(c, ax=axs[1])

    # Save figure 
    figName = "heatmap.png"
    figFolder = './figures'
    if not os.path.exists(figFolder):
        os.makedirs(figFolder)  # Changed to os.makedirs for better practice
    plt.savefig(f'{figFolder}/{figName}')
    plt.close()

def main():
    # Load data
    data = np.loadtxt('../original_Correlation.txt')
    data2 = np.loadtxt('../sampled_correlation.txt')

    # draw histogram
    bin = 'auto'
    draw_histogram(data, data2, bin)

    # draw heatmap
    data = np.reshape(data, (27, 30))
    data2 = np.reshape(data2, (27, 30))
    drawHeatMap(data, data2)

if __name__ == "__main__":
    main()
