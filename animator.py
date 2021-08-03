import sys
import time
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from graphics import *
from PIL import Image as NewImage
from tqdm import trange

SCALE = 1
WIDTH = 1921
HEIGHT = 1081


def clear(win):
    for item in win.items[:]:
        item.undraw()


def save_frame(win, i):
    # saves the current TKinter object in postscript format
    win.postscript(file="products/frames/.image.eps", colormode='color')
    # Convert from eps format to gif format using PIL
    img = NewImage.open("products/frames/.image.eps")
    img.save("products/frames/{}.jpeg".format(i, "jpeg"), quality=100, subsampling=0)

def animate():
    df = pd.read_csv(sys.argv[1], header=None)
    df.columns = 't', 'id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'
    N = max(df.id + 1)

    win = GraphWin('{} Body Simulator'.format(N), WIDTH, HEIGHT, autoflush=False)
    for i in trange(0, df.shape[0], 1):

        if i % N == 0:
            clear(win)

        x = df.loc[i]['x'] * SCALE + WIDTH/2
        y = -df.loc[i]['y'] * SCALE + HEIGHT/2
        mass = df.loc[i]['mass']
        id = df.loc[i]['id']
        location = Point(x, y)
        body = Circle(Point(x, y), mass/2)
        body.setFill('black')
        body.draw(win)

        if id == N-1:
            win.update()
            save_frame(win, int((i+1)/N))

    win.close()


if __name__ == '__main__':
    animate()
