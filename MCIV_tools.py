#! /usr/bin/env python

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is a collection of subroutines to assist with detchar analyses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from numpy import *
import commands, time, sys
from pylal import Fr
from scipy import signal
from glue import lal
from pylal import frutils
from scipy.stats import poisson
from glue import segmentsUtils
from glue.segments import *


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# channelRange -- Function to return max and min values of a channel over a defined segment list
#
# 'segments' is a space-delimited list of science segments
# 'ifo' is a string (H1, L1, etc)
# 'channel_name' is the channel, i.e. 'ASC-QPDX_P'
# 'end_buffer' is the number of seconds before the end of a segment to ignore
# 'segment_min_length' is the minimum segment length to consider
#
# Note that end_buffer must be less than segment_min_length!
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#def channelRange(ifo,segments,channel_name,end_buffer):
def channelRange(ifo,segments,channel_name,end_buffer,segment_min_length):

    channel_min = ifo + ':' + channel_name + '.min'
    channel_max = ifo + ':' + channel_name + '.max'

    channel_list = [channel_min, channel_max]

    range_min = False
    range_max = False

    seg_min = []
    seg_max = []

    # Loop over segments
    i = 0
    for segment in segments:
        i = i+1
        seg_start = str(int(segment[0]))
        seg_stop = str(int(segment[1]))
        #seg_start,seg_stop = segment.split()
        start = int(seg_start)
        stop = int(seg_stop)
        seg_length = stop-start

        if seg_length < segment_min_length:
            continue

        if seg_length < 1000:

            seg_start = str(start)
            seg_stop = str(stop-end_buffer)

            # frame type T = second trends
            data = getChannelVector(ifo,seg_start,seg_stop,'T',channel_list,rate=1)

            seg_min = append(seg_min,min(data[:,0]))
            seg_max = append(seg_max,max(data[:,1]))


        else:

            # Use the minute trends for the middle of long segments; it saves time
            seg1 = str(start)
            seg2 = str(start+300)
            seg3 = str(stop-300)
            seg4 = str(stop-end_buffer)

            data1 = getChannelVector(ifo,seg1,seg2,'T',channel_list,rate=1)
            data2 = getChannelVector(ifo,seg2,seg3,'M',channel_list,rate=0.0166667)
            data3 = getChannelVector(ifo,seg3,seg4,'T',channel_list,rate=1)

            # Minute trends are a pain...sample rate less than 1 leads to uncertainty in number of datapoints
            data_min = concatenate((data1[:,0],trim_zeros(data2[:,0]),data3[:,0]))
            data_max = concatenate((data1[:,1],trim_zeros(data2[:,1]),data3[:,1]))

            seg_min = append(seg_min,data_min.min())
            seg_max = append(seg_max,data_max.max())

            #seg_min = append(seg_min,median(data_min))
            #seg_max = append(seg_max,median(data_max))

            #print start, data_min.min(), data_max.max()

    if not(any(seg_min)) or not(any(seg_max)):
        print 'No channel range data!'
        sys.exit()

    range_min = min(seg_min)
    range_max = max(seg_max)

    return range_min, range_max


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# getChannelVector - wrapper that returns multi-channel data from frame files with downsampling
#
# All fields required.
#
# ifo -- H1, L1, etc.
# seg_start -- start time in GPS seconds of data segment to return
# seg_stop  -- stop time in GPS seconds
#
# frame_type -- exact frame type for ligo_data_find (H1_RDS_R_L1, T, M, etc)
#
# channel_list -- list of channel names to return, i.e. ['H1:ASC-QPDX_P.min','H1:ASC-QPDX_P.max']
#
# rate -- desired sample rate of returned data.  Must be <= sample rate of channel_list on frame_type, so do your homework!
#         (noninteger sample rates are probably a bad idea) 
# Returns:
#
# data -- NxM array of M channels for N = (seg_stop - seg_start)*rate samples
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getChannelVector(ifo,seg_start,seg_stop,frame_type,channel_list,rate,filter_rate=None):

    start = int(seg_start)
    stop = int(seg_stop)

    seglength = stop-start

    cmd = 'gw_data_find -o ' + ifo[0] + ' -s ' + seg_start + ' -e ' + seg_stop + ' -t ' + frame_type + ' -u file --lal-cache'
    cache = commands.getoutput(cmd)
    lines = cache.rsplit('\n')
    #print cmd

    ifo, frame_type, s, d, path = lines[0].split()
    duration = int(d)
    frame_start = int(s)
    frame = path[16:len(path)]
    
    frame_stop = frame_start + duration
    
    data_start = max(start,frame_start)
    data_stop  = min(stop,frame_stop)
    
    span = data_stop - data_start
    
    channel_rate = zeros(len(channel_list))
    for i in range(len(channel_list)):
        frame_data = Fr.frgetvect(frame,channel_list[i],data_start,span)
        channel_rate[i] = 1./frame_data[3][0]

    if not(all(channel_rate==channel_rate[0])):
        print 'Not all channel rates are the same!'
        sys.exit()

    # Initialize data array
    # use a fixed length array of zeros - zeros necessary for minute trends, fixed length necessary for pre-allocation
    # As for now we assume all channels are being read at the same rate

    data        = zeros((int(around(seglength*rate)),len(channel_list)))
    data_vector = zeros((int(around(seglength*channel_rate[0])),len(channel_list)))
    data_idx = 0

    for line in lines:
        ifo, frame_type, s, d, path = line.split()
        duration = int(d)
        frame_start = int(s)
        frame = path[16:len(path)]

        frame_stop = frame_start + duration

        data_start = max(start,frame_start)
        data_stop  = min(stop,frame_stop)

        span = data_stop - data_start

        for i in range(len(channel_list)):
            frame_data = Fr.frgetvect(frame,channel_list[i],data_start,span)
            sample_rate = 1./frame_data[3][0]

            if sample_rate != channel_rate[i]:
                print 'Sample rate of channel', channel_list[i], 'has changed!'
                sys.exit()

            x = frame_data[0]

            if frame_type == 'M' or frame_type[-2:] == '_M':
                new_data_idx = len(x) + data_idx

            elif abs(shape(x)[0] - (span*channel_rate[i])) > 0.1:
                print 'Length of returned data vector does not match expectation!'
                print 'getChannelVector', shape(x)[0], data_idx, span, rate
                sys.exit()

            else:
                new_data_idx = len(x) + data_idx

            data_vector[data_idx:new_data_idx,i] = x

        data_idx = new_data_idx


        if abs(sample_rate - rate) < 0.1:
            data = data_vector
        else:
            if sample_rate < rate:
                print 'Sample rate for channel ' + channel_list[i] + ' is smaller than requested.  Check frame type!'                
                sys.exit()

            q = int(sample_rate / rate)
            x = frame_data[0]

            for i in range(len(channel_list)):
                data[:,i] = decimate(data_vector[:,i],q)

    
    return data



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scipy.signal.decimate function
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.decimate.html
#
# Inputs:
#
# x (Required) -- is an array to be anti-aliased and downsampled
# q (Required) -- is the downsampling factor (Ex: 16384Hz -> 256Hz is q=64
#
# n (Optional) -- is the order of the anti-aliasing filter
# ftype (Optional) -- is the filter type ('iir' or 'fir')
# axis (Optional) -- is the axis along which to decimate
#
# Returns:
#
# y -- the decimated data
#
# By default an order 3 Chebyshev type I filter is used. A 30 point FIR filter with hamming window is used if ftype is 'fir'.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def decimate(x, q, n=None, ftype='iir', axis=-1):
    if not isinstance(q, int):
        raise TypeError("q must be an integer")
    if n is None:
        if ftype == 'fir':
            n = 30
        elif q>32:
            n = 2
        else:
            n = 3
    if ftype == 'fir':
        b = signal.firwin(n + 1, 1. / q, window='hamming')
        a = 1.
    else:
        b, a = signal.cheby1(n, 0.05, 0.8 / q)
        #b, a = signal.butter(4, 0.8 / q, 'low', analog=False)

    #print b, a
    y = signal.filtfilt(b, a, x)
    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, q)
    return y[sl]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Returns a list of paths to omega trigger files to be searched when looking for glitches between 'start' and 'stop'
# The files for S6 single-detector omega triggers (D. McLeod) are organized by epoch, then by day.
# Time spans that go across boundaries of S6 epochs are currently not allowed
# 'start' and 'stop' must be integer gps times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getGlitchFiles(ifo,seg_start,seg_stop):

    start = int(seg_start)
    stop = int(seg_stop)

    # S6 epochs

    tStartS6A = 930960015
    tEndS6A = 935798415
    #tEndS6A = 937526415
    
    tStartS6B = 937785615
    tEndS6B = 947203215
    
    tStartS6C = 947635215
    tEndS6C = 961545615
    
    tStartS6D = 961545615
    tEndS6D = 971654415
    
    if start < tStartS6A:
        print "GPS start time is before start of S6."
        sys.exit()
    elif (start >= tStartS6A) and (start <= tEndS6A):
        start_epoch = 's6a'
        epoch_start = tStartS6A
    elif (start >= tStartS6B) and (start <= tEndS6B):
        start_epoch = 's6b'
        epoch_start = tStartS6B
    elif (start >= tStartS6C) and (start <= tEndS6C):
        start_epoch = 's6c'
        epoch_start = tStartS6C
    elif (start >= tStartS6D) and (start <= tEndS6D):
        start_epoch = 's6d'
        epoch_start = tStartS6D
    elif start > tEndS6D:
        print "GPS start time is after S6."
        sys.exit()
    else:
        print "GPS start time is not included in an epoch of S6."
        sys.exit()
        
    if stop < tStartS6A:
        print "GPS end time is before start of S6."
        sys.exit()
    elif (stop >= tStartS6A) and (stop <= tEndS6A):
        stop_epoch = 's6a'
        epoch_stop = tEndS6A
    elif (stop >= tStartS6B) and (stop <= tEndS6B):
        stop_epoch = 's6b'
        epoch_stop = tEndS6B
    elif (stop >= tStartS6C) and (stop <= tEndS6C):
        stop_epoch = 's6c'
        epoch_stop = tEndS6C
    elif (stop >= tStartS6D) and (stop <= tEndS6D):
        stop_epoch = 's6d'
        epoch_stop = tEndS6D
    elif stop > tEndS6D:
        print "GPS end time is after S6."
        sys.exit()
    else:
        print "GPS end time is not included in an epoch of S6."
        sys.exit()

    # The gitch filenames are unforgiving for searches across epoch boundaries.  {sigh}
    if start_epoch != stop_epoch:
        print 'Times span two S6 epochs.  Please restrict your analysis to a single epoch.'
        print 'Start time is in ', start_epoch
        print 'Stop time is in ', stop_epoch
        sys.exit()


    # Find number of days in the span
    cmd = 'lalapps_tconvert -f %j ' + seg_start
    day1 = commands.getoutput(cmd)

    cmd = 'lalapps_tconvert -f %j ' + seg_stop
    day2 = commands.getoutput(cmd)

    num_days = int(day2) - int(day1)

    # Find the day's number in the epoch
    start_day_in_epoch = (start - epoch_start) / 86400

    glitch_files = []

    # For each day, find the glitch file
    for i in range(num_days+1):

        # Convert gpstime to MMM/DD/YY
        date = commands.getoutput('lalapps_tconvert -f %D ' + str(start + i*86400))

        # Find gpstime for the start of the day
        day_start = int(commands.getoutput('lalapps_tconvert ' + str(date)))
    
        # The glitch files start at 4 seconds before midnight
        day_file_start = day_start - 4
    
        # Find day of epoch
        day_of_epoch = start_day_in_epoch + i

        # Construct path to glitch file
    
#        path = '/home/dhoak/detchar/omega/' + ifo.lower() + '1/' + start_epoch + '/' + str(epoch_start) + '-' + str(epoch_stop) + '/'
        path = '/home/detchar/omega/' + start_epoch + '/' + str(epoch_start) + '-' + str(epoch_stop) + '/'

        path = path + ifo.upper() + '1-LDAS_STRAIN/segment_' + str(day_of_epoch) + '-' + str(day_file_start) + '-86408/'

        path = path + ifo.upper() + '1-OMEGA_TRIGGERS_DOWNSELECT_CLUSTERED-' + str(day_start) + '-86400.txt'

        glitch_files.append(path)

    return glitch_files


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Returns the date string, 'DDMMMYYYY', for a given gpstime
# 'gpstime' is a gpstime in string format
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getDateString(gpstime,hms=False):

    cmd = 'lalapps_tconvert -d ' + gpstime
    date_return = commands.getoutput(cmd)
    date_array = date_return.split(' ')
    day = date_array[2]
    month = date_array[1]
    year = date_array[5]
    hms = date_array[3]
    hour, minute, second = hms.split(':')

    """
    cmd = 'lalapps_tconvert -f %d ' + gpstime
    day = commands.getoutput(cmd)
    day_string = ''.join(day.rsplit('\n'))

    cmd = 'lalapps_tconvert -f %b ' + gpstime
    month = commands.getoutput(cmd)
    month_string = ''.join(month.rsplit('\n'))

    cmd = 'lalapps_tconvert -f %m ' + gpstime
    month_number = commands.getoutput(cmd)

    cmd = 'lalapps_tconvert -f %Y ' + gpstime
    year = commands.getoutput(cmd)
    year_string = ''.join(year.rsplit('\n'))

    cmd = 'lalapps_tconvert -f %H ' + gpstime
    hour = commands.getoutput(cmd)

    cmd = 'lalapps_tconvert -f %M ' + gpstime
    minute = commands.getoutput(cmd)

    cmd = 'lalapps_tconvert -f %S ' + gpstime
    second = commands.getoutput(cmd)
    """

    if hms:
        #x = dict((v,k) for k,v in enumerate(calendar.month_abbr))
        #month_number = x[month]
        date_string = day + ' ' + month  + ' ' + year
        #return date_string, month_number, day, year, hour, minute, second
        return date_string, month, day, year, hour, minute, second
    else:
        date_string = day + month + year
        return date_string, month, day, year




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A function to return the veto segments within a science segment.
#
# Inputs:  'epoch', a length=1 glue.segments object
#
#          'segmentList', a length=N glue.segments object
#
#
# Output:  'segments_in_epoch', a list of segments from segmentList that fall between within 'epoch'
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to implement hveto safety test
#
# Inputs:
#
# ifo          - str, H1 or L1
# start_time   - str, start time for study
# end_time     - str, end time for study
# veto_segmets - glue.segments list of veto segments to check for safety
# deadtime     - float, deadtime % in science-cat1 time for the veto segments in the epoch [start_time, end_time]
#
# Output:
#
# num_inj        - number of successful hardware injections at ifo in epoch [start_time, end_time]
# num_vetoed_inj - number of injection times +/- 50msec that overlap with vetos in veto_segments
# significance   - -log10 of probability for vetoing num_vetoed_inj given num_inj and deadtime
#
# Note - the significance calculation here can returns a large positive value if the number of vetoed injections
#        is very improbable from chance Poisson statistics alone.  This is true if the number is much larger than
#        expected, AND if the number is much lower.  If 40 vetoed injections are expected, and ony two are vetoed,
#        the significance is 2.6 (which is typically above the threshold for safety).  Check if  
#        num_vetoed_inj < deadtime * num_inj.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hveto_safety(ifo,veto_segments,start_time,end_time,deadtime):

    if ifo=='H1':
        inj_file = '/home/dhoak/detchar/S6_segments/H1_injection_segments_raw.txt'
    elif ifo=='L1':
        inj_file = '/home/dhoak/detchar/S6_segments/L1_injection_segments_raw.txt'
    else:
        print 'IFO not valid!'
        sys.exit()

    epoch_start = int(start_time)
    epoch_stop = int(end_time)
    epoch_seg = segmentlist(segmentsUtils.segmentlist_range(epoch_start,epoch_stop,epoch_stop-epoch_start))

    # Load the successful hardware injection segments for S6
    inj_seglist = segmentsUtils.fromsegwizard(open(inj_file), coltype = float).coalesce()

    # Get the successful hardware injections in the epoch defined by [start_time, end_time]
    inj_in_epoch = inj_seglist & epoch_seg

    num_inj = len(inj_in_epoch)

    # mu is the expected number of vetoed injections from chance coincidences
    mu = num_inj * deadtime

    vetoed_inj_seg = inj_in_epoch & veto_segments

    #num_vetoed_inj = 0
    #for inj in inj_in_epoch:
    #    num_vetoed_inj += len(getSegmentOverlap(inj,veto_segments))

    num_vetoed_inj = len(vetoed_inj_seg.coalesce())

    # apply Eq. 1-3 from hveto paper to calculate significance
    P_sum = 0
    rv = poisson(mu)
    for jj in range(num_vetoed_inj,(num_vetoed_inj+1)*10):
        P_sum += rv.pmf(jj)

    significance = -log10(P_sum)

    return num_inj, num_vetoed_inj, significance



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Given a start time and an end time, returns an S6 epoch for file handling
# Start time and end time must fall within the same epoch
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getS6epoch(start_time,end_time):

    start = int(start_time)
    stop = int(end_time)

    # S6 epochs

    tStartS6A = 930960015
    tEndS6A = 937526415
    
    tStartS6B = 937785615
    tEndS6B = 947203215
    
    tStartS6C = 947635215
    tEndS6C = 961545615
    
    tStartS6D = 961545615
    tEndS6D = 971654415
    
    if start < tStartS6A:
        print "GPS start time is before start of S6."
        sys.exit()
    elif (start >= tStartS6A) and (start <= tEndS6A):
        start_epoch = 'S6A'
        epoch_start = tStartS6A
    elif (start >= tStartS6B) and (start <= tEndS6B):
        start_epoch = 'S6B'
        epoch_start = tStartS6B
    elif (start >= tStartS6C) and (start <= tEndS6C):
        start_epoch = 'S6C'
        epoch_start = tStartS6C
    elif (start >= tStartS6D) and (start <= tEndS6D):
        start_epoch = 'S6D'
        epoch_start = tStartS6D
    elif start > tEndS6D:
        print "GPS start time is after S6."
        sys.exit()
    else:
        print "GPS start time is not included in an epoch of S6."
        sys.exit()
        
    if stop < tStartS6A:
        print "GPS end time is before start of S6."
        sys.exit()
    elif (stop >= tStartS6A) and (stop <= tEndS6A):
        stop_epoch = 'S6A'
        epoch_stop = tEndS6A
    elif (stop >= tStartS6B) and (stop <= tEndS6B):
        stop_epoch = 'S6B'
        epoch_stop = tEndS6B
    elif (stop >= tStartS6C) and (stop <= tEndS6C):
        stop_epoch = 'S6C'
        epoch_stop = tEndS6C
    elif (stop >= tStartS6D) and (stop <= tEndS6D):
        stop_epoch = 'S6D'
        epoch_stop = tEndS6D
    elif stop > tEndS6D:
        print "GPS end time is after S6."
        sys.exit()
    else:
        print "GPS end time is not included in an epoch of S6."
        sys.exit()

    # The gitch filenames are unforgiving for searches across epoch boundaries.  {sigh}
    if start_epoch != stop_epoch:
        print 'Times span two S6 epochs.  Please restrict your analysis to a single epoch.'
        print 'Start time is in ', start_epoch
        print 'Stop time is in ', stop_epoch
        sys.exit()


    return start_epoch



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A function to return the total time in a list of segments (used for lists returned by getSegmentOverlap)
#
# Inputs:  'segmentList', a list of segments, shape Nx2
#
# Output:  'segment_total', the total time included in the segments
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getSegmentSum(segmentList):

    segment_total = 0

    for seg in segmentList:
        start = seg[0]
        stop = seg[1]

        segment_total += stop-start

    return segment_total
