import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

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
