#!/usr/bin/env python

import json
import sys

if len(sys.argv) == 3:
    f = open(sys.argv[1])
    memory_usage = sys.argv[2]

    data = json.load(f)

    i = 0
    name, inp = data['benchmarks'][i]['name'].split("/", 1)
    time = data['benchmarks'][i]['real_time']
    time_unit = data['benchmarks'][i]['time_unit']
    mbu = data['benchmarks'][i].get('max_bytes_used', -1)
    print('{:>50}, {:>10}, {:25f}, {:>4}, {:20d}, {:>20}'.format(name, inp, time, time_unit, mbu, memory_usage))
elif len(sys.argv) == 1:
    name = "Name"
    inp = "Data/Iters"
    time = "Time"
    time_unit = "unit"
    mbu = "Memory (gbench)"
    memory_usage = "Memory(time)"
    print('{:>50}, {:>10}, {:>25}, {:>4}, {:>20}, {:>20}'.format(name, inp, time, time_unit, mbu, memory_usage))


