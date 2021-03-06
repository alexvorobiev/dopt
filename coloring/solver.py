#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from subprocess import Popen, PIPE


def solve_it(input_data):
    # Modify this code to run your optimization algorithm

    # Writes the inputData to a temporay file

    tmp_file_name = 'tmp.data'
    tmp_file = open(tmp_file_name, 'w')
    tmp_file.write(input_data)
    tmp_file.close()

    # Runs the command: java Solver -file=tmp.data

    #process = Popen(['C:/Program Files/Wolfram Research/Mathematica/9.0/math', '-script', 'solver.m',  tmp_file_name], stdout=PIPE)
    #process = Popen(['c:/Program Files/R/R-3.0.3/bin/x64/Rscript', 'solver.R', '-file=' + tmp_file_name], stdout=PIPE)


# def solve_it_dsatur_2(input_data):
#     """Solver for another DSATUR implementation
#     """

#     # Writes the inputData to a temporay file

#     #tmp_file_name = 'tmp.data'
#     #tmp_file = open(tmp_file_name, 'w')
#     tmp_file = tempfile.NamedTemporaryFile(dir='./', prefix='input_graph.',
#                                            delete=False)
#     tmp_file.write(input_data)
#     tmp_file.close()

#     # Run DSATUR

    process = Popen(['./run_coloring.sh', tmp_file.name], stdout=PIPE)


    (stdout, stderr) = process.communicate()

    # removes the temporay file
    os.remove(tmp_file_name)

    return stdout.strip()


import sys

if __name__ == '__main__':
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        input_data_file = open(file_location, 'r')
        input_data = ''.join(input_data_file.readlines())
        input_data_file.close()
        print solve_it(input_data)
    else:
        print 'This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/gc_4_1)'

