import sys
import pandas as pd
import numpy as np
from tqdm import tqdm

from vpython import *

FRAME_RATE = 30     # frames per second


def animate():
    df = pd.read_csv(sys.argv[1], header=None, low_memory=False)
    df.columns = 't', 'id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'

    bodies = df.to_numpy()
    N = int(max(df.id + 1))

    canvas(title='{} body simulator'.format(N), width=1800, height=800)
    balls = [sphere(color=color.cyan) for i in range(N)]

    for i in tqdm(range(0, df.shape[0], N)):
        rate(FRAME_RATE)
        x = bodies[i:i+N, 2]
        y = bodies[i:i+N, 3]
        z = bodies[i:i+N, 4]
        mass = bodies[i:i+N, -1]
        for j in range(N):
            balls[j].pos = vector(x[j], y[j], z[j])
            balls[j].radius = mass[j]/2 


if __name__ == '__main__':
    animate()
