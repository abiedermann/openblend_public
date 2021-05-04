#!/usr/bin/env python

from experiment_designer import experiment_designer
import sys

output_filename = sys.argv[1][:-5]+"_reprocessed.xlsx"

foo = experiment_designer(output_filename=output_filename,
        input_filename=sys.argv[1])
print(foo.volumes.sum(axis=0))
print(foo.stocks)
foo.print_design()


