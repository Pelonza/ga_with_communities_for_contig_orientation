"""
By Karl Schmitt, 4/19/17

This script performs basic visualizations of the logbooks from using GA's to 
solve the contig orientation problem.

Dependencies:
    Bokeh
    DEAP
    json
    
"""

import bokeh as bk

from deap import tools
import json
import os


def get_parser():
    """Get parser object for script xy.py."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--ifile",
                        dest="ifilename",
                        type=lambda x: is_valid_file(parser, x),
                        help="Read tallies from FILE",
                        metavar="FILE")
#    parser.add_argument("-o", "--ofile",
#                        dest="oorient",
#                        help="Write best orientation to FILE",
#                        metavar="FILE")
#    parser.add_argument("-s", "--sfile",
#                        dest="output_stats",
#                        help="Write statistics to FILE",
#                        metavar="FILE")
#    parser.add_argument("-n",
#                        dest="ntrials",
#                        default=50,
#                        type=int,
#                        help="how many trials to do per parameter set")
#    parser.add_argument("-q", "--quiet",
#                        action="store_false",
#                        dest="verbose",
#                        default=True,
#                        help="don't print status messages to stdout")
    return parser

def is_valid_file(parser, arg):
    """
    Check if arg is a valid file that already exists on the file system.

    Parameters
    ----------
    parser : argparse object
    arg : str

    Returns
    -------
    arg
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg
# %%
if __name__ == "__main__":
    
    
    args=get_parser().parse_args()
    
    with open(args.ifilename, 'r') as f:
        #logbook=json.load(f)
        
    
    print(logbook)