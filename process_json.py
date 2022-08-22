#!/usr/bin/env python

import json
import sys

f = open(sys.argv[1])

data = json.load(f)

for i in range(0, 5):
    name, inp = data['benchmarks'][i]['name'].split("/", 1)
    time = data['benchmarks'][i]['real_time']
    mbu = data['benchmarks'][i]['max_bytes_used']

    print('{:>25}, {:>10}, {:15.2f}, {:20d}'.format(name, inp, time, mbu))
