#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script calls seg_study as a function to collect efficiency/deadtime statistics for a set of veto segments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import commands, sys, time, os
from optparse import OptionParser
#import h5py
#from MVIC_tools import *
#from glue import segmentsUtils
#from glue.segments import *
from seg_study import *


web_directory = '/home/dhoak/public_html/'

veto_file_path = '/home/dhoak/detchar/S6_segments/'

ifo = 'H1'

#epochs = ['S6A','S6B','S6C','S6D']
epochs = ['S6C']

for epoch in epochs:

    if epoch=='S6A':

        start_time = '930960015'   # Tue Jul 07 00:00:00 GMT 2009
        end_time   = '935798415'   # Tue Sep 01 00:00:00 GMT 2009
        #end_time   = '937526415'   # Mon Sep 21 00:00:00 GMT 2009

    elif epoch=='S6B':

        start_time = '937785615'   # Thu Sep 24 00:00:00 GMT 2009
        end_time   = '947203215'   # Mon Jan 11 00:00:00 GMT 2010

    elif epoch=='S6C': 

        start_time = '947635215'   # Sat Jan 16 00:00:00 GMT 2010
        #end_time   = '961545615'   # Sat Jun 26 00:00:00 GMT 2010
        end_time   = '961535808'   # Fri Jun 25 21:16:33 GMT 2010

    elif epoch=='S6D': 

        start_time = '961545715'   # Sat Jun 26 00:00:00 GMT 2010
        end_time   = '971654415'   # Thu Oct 21 00:00:00 GMT 2010

    #veto_segfile = veto_file_path + ifo + '_' + epoch + '_sciencecat1_cat3veto.txt'
    veto_segfile = veto_file_path + ifo + '_' + epoch + '_sciencecat1_cat2cat3hveto.txt'

    data_path = 'detchar/S6_segments/' + ifo + '/' + epoch + '/cat2cat3hveto/'

    if not os.path.exists(web_directory + data_path):
        os.makedirs(web_directory + data_path)

    efficiency_per_SNR, deadtime, num_inj, num_vetoed_inj, safety_significance = seg_study(ifo,data_path,start_time,end_time,veto_segfile)

