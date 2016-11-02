import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde

def myboxplot (ax,data,colors) :
    bp = ax.boxplot(data, 1, '')
    plt.setp (bp['boxes'], color='black')
    plt.setp (bp['whiskers'], color='black', linewidth=2)
    plt.setp (bp['caps'], color='black', linewidth=2)
    plt.setp (bp['medians'], color='black')
    for i,box in enumerate(bp['boxes']) :
        boxCoords = zip(box.get_xdata(),box.get_ydata())
        boxPolygon = Polygon(boxCoords, facecolor=colors[i])
        ax.add_patch(boxPolygon)

def plot_hic_matrix (ax,H,start,end,resolution,fmt = "%.1f",scale=1000000) :
    # plot the matrix on the given ax
    ax.matshow (1-np.log2(H),cmap=plt.cm.gray,origin='lower',interpolation='none')
    # fix the ticks
    ticklabels = []
    for tick in ax.get_xticks().tolist() :
        newtick = (start + tick*resolution)/scale
        ticklabels.append(fmt%newtick)
    ax.set_xticklabels(ticklabels)
    ax.set_yticklabels(ticklabels)
    ax.text (1.0,1.01,'Mb',transform=ax.transAxes)

def ax_only_y (ax,show_xaxis=False) :
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().tick_left()
    if not show_xaxis :
        ax.get_xaxis().set_visible(False)
        ax.spines['bottom'].set_visible(False)
    else :
        ax.get_xaxis().tick_bottom()

def line_plot (ax,xvals,yvals,N1=None,N2=None,show_xaxis=False) :
    if N1 is not None and N2 is not None :
        ax.set_xlim (N1,N2)
        mask = np.logical_and(xvals>N1,xvals<N2)
    else :
        ax.set_xlim (min(xvals),max(xvals))
        mask = np.ones_like(xvals,dtype=bool)
    # plot values
    for x,y in zip(xvals[mask],yvals[mask]) :
        ax.add_artist(Line2D((x,x),(0,y),color='k',linewidth=1))
    # plot borders
    ymin = min(yvals)
    ymax = max(yvals)
    delta = ymax-ymin
    ax.set_ylim (ymin-0.01*delta,ymax+0.01*delta)
    # plot style
    ax_only_y (ax,show_xaxis)

def color_density_scatter (ax,x,y) :
    """
    Calculate the gaussian kernel density of the points that we want to look at:
    this way we will be able to color-code the points on the plot by the
    density of the neighbouring points. Taken from
    http://stackoverflow.com/a/20107592/2312821
    Then return the scatter plot color-coded correspondingly.
    """
    # first calculate the Gaussian kernel density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    # we need then to sort the output, so that the points with highest density
    # will be plotted last
    idx = z.argsort()
    x,y,z = x[idx],y[idx],z[idx]
    ax.scatter(x,y,c=np.log(z),s=10,edgecolor='')
