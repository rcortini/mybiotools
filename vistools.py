import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
import itertools
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties

def myboxplot (ax,data,colors) :
    bp = ax.boxplot(data, 1, '')
    plt.setp (bp['boxes'], color='black')
    plt.setp (bp['whiskers'], color='black', linewidth=1)
    plt.setp (bp['caps'], color='black', linewidth=1)
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

def line_plot (ax,xvals,yvals,N1=None,N2=None,show_xaxis=False,
               color='k',linewidth=1) :
    if N1 is not None and N2 is not None :
        ax.set_xlim (N1,N2)
        mask = np.logical_and(xvals>N1,xvals<N2)
    else :
        ax.set_xlim (min(xvals),max(xvals))
        mask = np.ones_like(xvals,dtype=bool)
    # plot values
    xvals = np.array(xvals)
    yvals = np.array(yvals)
    for x,y in zip(xvals[mask],yvals[mask]) :
        ax.add_artist(Line2D((x,x),(0,y),color=color,linewidth=linewidth))
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

def plot_triangular_matrix (ax,C,cmap=plt.cm.gray):
    n = C.shape[0]
    t = np.array([[1,0.5],[-1,0.5]])
    A = np.dot(np.array([(i[1],i[0]) 
               for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    ax.pcolormesh(A[:,1].reshape(n+1,n+1),A[:,0].reshape(n+1,n+1),np.flipud(C),
                 cmap=cmap)
    ax.set_ylim(0,n)
    ax.set_frame_on(False)
    plt.tick_params(which='both',length=0)
    ax.add_artist(Line2D((0,n/2),(0,n),color='k',linewidth=1))
    ax.add_artist(Line2D((n,n/2),(0,n),color='k',linewidth=1))
    ax.set_yticklabels([])

def letterAt(letter, x, y, yscale=1, ax=None):
    """
    This function is used in the sequence_logo function, and is not exported in
    the module. From https://stackoverflow.com/a/42631740/2312821
    """
    fp = FontProperties(family="DejaVu Sans", weight="bold") 
    globscale = 1.35
    LETTERS = { "T" : TextPath((-0.305, 0), "T", size=1, prop=fp),
                "G" : TextPath((-0.384, 0), "G", size=1, prop=fp),
                "A" : TextPath((-0.35, 0), "A", size=1, prop=fp),
                "C" : TextPath((-0.366, 0), "C", size=1, prop=fp) }
    COLOR_SCHEME = {'G': 'orange', 
                    'A': 'red', 
                    'C': 'blue', 
                    'T': 'darkgreen'}
    text = LETTERS[letter]
    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + \
        mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=COLOR_SCHEME[letter],  transform=t)
    if ax != None:
        ax.add_artist(p)
    return p

def sequence_logo(pwm) :
    """
    Draws a sequence logo starting from a position weight matrix. The position
    weight matrix should be an instance of the
    Bio.motifs.matrix.PositionWeightMatrix class.
    """
    N = pwm.length
    # first, calculate the information content in each position
    ic = np.zeros(N)
    for i in xrange(N) :
        s = 2.0
        for letter in pwm.iterkeys() :
            x = pwm[letter][i]
            if x>0 : s += x*np.log2(x)
        ic[i] = s
    # then initiate the figure
    fig,ax = plt.subplots(figsize=(N,3))
    ax_only_y(ax,show_xaxis=True)
    # fill the figure with the letters
    for x in xrange(pwm.length) :
        y = 0
        I = ic[x]
        for letter in pwm.iterkeys() :
            score = pwm[letter][x]*I
            letterAt(letter,x+1,y,score,ax)
            y += score
    # final adjustments
    plt.xticks(range(1,x+2))
    plt.xlim((0,x+2))
    plt.ylim(0,2)
    # and return
    return fig
