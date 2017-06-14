#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script generates segment files for the S6 epochs from the category definer files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
#import matplotlib
#matplotlib.use("Agg")
#import pylab
import commands, sys, time, os
#import h5py
#from scipy import signal
#from MCIV_tools import *
from glue import segmentsUtils
from glue.segments import *


def getSegmentOverlap(epoch,segmentList):

    segments_in_epoch = []
    segStart = epoch[0][0]
    segStop = epoch[0][1]

    for seg in segmentList:
        start = seg[0]
        stop = seg[1]

        if start > segStop:
            break
        elif stop < segStart:
            continue
        else:
            segment = [max(segStart,start), min(segStop,stop)]
            segments_in_epoch.append(segment)

    return segments_in_epoch


ifo = 'H1'
epoch = 'S6D'

veto_file_path = '/home/dhoak/detchar/S6_segments/'
veto_segfile = veto_file_path + ifo + '_' + epoch + '_sciencecat1_cat2cat3hveto.txt'
vetoes_in_epoch = segmentsUtils.fromsegwizard(open(veto_segfile), coltype = float, strict=False).coalesce()

start = 961720002
stop = 961724624
science_seg = segmentlist(segmentsUtils.segmentlist_range(start,stop,stop-start))

science_file = veto_file_path + ifo + '_' + epoch + '_sciencecat1.txt'

science_segments = segmentsUtils.fromsegwizard(open(science_file), coltype = float, strict=False).coalesce()

#clock3 = time.time()
#vetoes_in_segment1 = science_seg & vetoes_in_epoch
#clock4 = time.time()
#elapsed = clock4-clock3

#print
#print 'Number of veto segments in epoch is ', len(vetoes_in_segment1)
#print 'Time to perform glue.segments overlap is ', elapsed

clock3 = time.time()
vetoes_in_segment2 = getSegmentOverlap(science_seg,vetoes_in_epoch)
clock4 = time.time()
elapsed = clock4-clock3

print
print 'Number of veto segments in epoch is ', len(vetoes_in_segment2)
print 'Time to perform manual overlap is ', elapsed

print shape(vetoes_in_segment2)

"""
i=0
for science_segment in science_segments:
    i+=1
    start = science_segment[0]
    stop = science_segment[1]
    science_seg = segmentlist(segmentsUtils.segmentlist_range(start,stop,stop-start))


    clock3 = time.time()
    vetoes_in_segment2 = getSegmentOverlap(science_seg,vetoes_in_epoch)
    clock4 = time.time()
    elapsed = clock4-clock3

    print
    print 'Segment ', i
    print 'Number of veto segments in epoch is ', len(vetoes_in_segment2)
    print 'Time to perform manual overlap is ', elapsed


for i in range(len(vetoes_in_segment1)):

    v1_start = vetoes_in_segment1[i][0]
    v1_stop = vetoes_in_segment1[i][1]

    v2_start = vetoes_in_segment2[i][0]
    v2_stop = vetoes_in_segment2[i][1]

    if v1_start != v2_start:
        print 'Start', v1_start, v2_start
    if v1_stop != v2_stop:
        print 'Stop', v1_stop, v2_stop
"""
