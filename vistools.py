import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

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

def plot_hic_matrix (H,ax,N1,N2,hic_res=2000,offset=0) :
    n1 = N1/hic_res
    n2 = N2/hic_res
    ax.matshow (1-np.log2(H[n1:n2,n1:n2]),
             cmap=plt.cm.gray,
             origin='lower',
             extent=[N1+offset,N2+offset,N1+offset,N2+offset],
             interpolation='none')

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
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().tick_left()
    if not show_xaxis :
        ax.get_xaxis().set_visible(False)
        ax.spines['bottom'].set_visible(False)
    else :
        ax.get_xaxis().tick_bottom()
