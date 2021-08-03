# Constants
N=1000
INIT_FILE=products_open/init.csv
OUTPUT_FILE=products_open/output.csv
OUTPUT_MOVIE=products_globular/nbody.mp4
TIME_STEP=0.1     # seconds
DURATION=1000.0  # seconds
OUTPUT_EVERY=10     # time steps between writing to the output file
G=1.0               # gravitational constant
EPS=2.0             # Smoothing constant to avoid singularities in the forces
THETA=0.5           # Barnes-Hut constant to decide when to combine far-away bodies.
                    #   Lower is more accurate and slower.
BRUTE_FORCE=1     # Whether to use the O(n^2) brute force algorithm


python3 initial_open.py $N

gcc-7 -Ofast -Wall -std=c99 -fopenmp \
             -D TIME_STEP=$TIME_STEP \
             -D DURATION=$DURATION \
             -D OUTPUT_EVERY=$OUTPUT_EVERY \
             -D G=$G \
             -D EPS=$EPS \
             -D THETA=$THETA \
             -D BRUTE_FORCE=$BRUTE_FORCE \
             nbody.c -o nbody.o -lm

nice ./nbody.o $INIT_FILE > $OUTPUT_FILE

python3 animator_vpython.py $OUTPUT_FILE
open $OUTPUT_MOVIE


