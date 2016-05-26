#!/usr/bin/env python

import sys

keeplocs = set()

with open(sys.argv[1], "r") as f:
    for l in f:
        bits = l.split()
        keeplocs.add((bits[0], bits[1]))

with open(sys.argv[2], "r") as f:
    for l in f:
        bits = l.split()
        if (bits[0], bits[1]) in keeplocs:
            sys.stdout.write(l)

