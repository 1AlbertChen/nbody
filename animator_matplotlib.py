import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
from PIL import Image as NewImage
from tqdm import tqdm

#SIZE = 800

def add_boxes():
    nodes = np.array(((-75.000000, 75.000000, 50.000000),
                (-25.000000, 75.000000, 50.000000),
                (-75.000000, 25.000000, 50.000000),
                (-25.000000, 25.000000, 50.000000),
                (-37.500000, 37.500000, 25.000000),
                (-12.500000, 37.500000, 25.000000),
                (-37.500000, 12.500000, 25.000000),
                (-12.500000, 12.500000, 25.000000),
                (50.000000, 50.000000, 100.000000),
                (25.000000, 75.000000, 50.000000),
                (75.000000, 75.000000, 50.000000),
                (25.000000, 25.000000, 50.000000),
                (12.500000, 37.500000, 25.000000),
                (37.500000, 37.500000, 25.000000),
                (12.500000, 12.500000, 25.000000),
                (37.500000, 12.500000, 25.000000),
                (75.000000, 25.000000, 50.000000),
                (-50.000000, -50.000000, 100.000000),
                (-75.000000, -25.000000, 50.000000),
                (-25.000000, -25.000000, 50.000000),
                (-75.000000, -75.000000, 50.000000),
                (-25.000000, -75.000000, 50.000000),
                (50.000000, -50.000000, 100.000000),
                (25.000000, -25.000000, 50.000000),
                (12.500000, -12.500000, 25.000000),
                (37.500000, -12.500000, 25.000000),
                (12.500000, -37.500000, 25.000000),
                (37.500000, -37.500000, 25.000000),
                (75.000000, -25.000000, 50.000000),
                (25.000000, -75.000000, 50.000000),
                (75.000000, -75.000000, 50.000000)))

    ax = plt.gca()
    for j in range(nodes.shape[0]):
        side = nodes[j, 2]
        x = nodes[j, 0]
        y = nodes[j, 1]
        # Top left
        rect = patches.Rectangle((x-side/2, y-side/2), side, side,
                                    linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(rect)


def animate():
    df = pd.read_csv(sys.argv[1], header=None)
    df.columns = 't', 'id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'

    bodies = df.to_numpy()
    N = max(df.id + 1)

    plt.figure(figsize=(20, 14))
    for i in tqdm(range(0, df.shape[0], N)):

        x = bodies[i:i+N, 2]
        y = bodies[i:i+N, 3]
        mass = bodies[i:i+N, -1]

        plt.scatter(x, y, s=mass*10)
        # plt.xlim([-2, 2])
        # plt.ylim([-2, 2])
        plt.xlim([-SIZE, SIZE])
        plt.ylim([-SIZE, SIZE])



        plt.axis('off')
        plt.title('{} Body Simulator'.format(N))
        plt.savefig('products/frames/{}.png'.format(i//N + 1), bbox_inches='tight')

        plt.cla()




if __name__ == '__main__':
    animate()
