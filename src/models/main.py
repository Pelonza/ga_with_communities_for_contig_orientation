# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 08:50:14 2017

@author:Karl R. B. Schmitt


"""


import ga_functions as ga_func #Contains self-defined functions for GA
import load_data as ld #Contains the load_data function.
import os
import json
import multiprocessing

def trial(args):
    
    G, Init_mps = ld.load_data(args.ifilename)
    
    mypop, mystats, myhof, mylog = ga_func.gaopt_Uni(G, Init_mps)
    
    with open(args.oorient, 'w+') as f:
        json.dump(myhof[0], f)
    with open(args.output_stats, 'w+') as f:
        json.dump(mylog,f)
        
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
    parser.add_argument("-o", "--ofile",
                        dest="oorient",
                        help="Write best orientation to FILE",
                        metavar="FILE")
    parser.add_argument("-s", "--sfile",
                        dest="output_stats",
                        help="Write statistics to FILE",
                        metavar="FILE")    
    parser.add_argument("-n",
                        dest="n",
                        default=10,
                        type=int,
                        help="how many lines get printed")
    parser.add_argument("-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        default=True,
                        help="don't print status messages to stdout")
    return parser


if __name__ == "__main__":
    """
    Main that calls a single trial with currently default parameters.
    """
    
    args = get_parser().parse_args()
    
    #pool = multiprocessing.Pool()
    trial(args)


    
        
    
        
        
        
        
        









