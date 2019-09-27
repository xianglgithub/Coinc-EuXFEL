#!/usr/bin/env python


import argparse
import os
import importlib as il
from proc.proc import Coinc


import configparser



parser = argparse.ArgumentParser(description='Coinc')
parser.add_argument('-i', '--ini', type=str)
args = parser.parse_args()


config = configparser.ConfigParser()
config.read(args.ini)


proc = Coinc(config)
proc.start(verbose=False)
    
    
    
