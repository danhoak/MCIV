#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function implements a X vs Y parameter space glitch study
# Called by MCIV_start with 'deriv' option
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
matplotlib.use("Agg",warn=False)
import pylab
import commands, sys, time, os
import h5py
from scipy import signal
from MCIV_tools import *
from glue import segmentsUtils
from glue.segments import *
import mpl_toolkits.mplot3d.axes3d as p3

# Probably should add some comments...

def glitch_map(ifo,channelX_name,channelY_name,start_time,end_time,sample_rate,frame_choice,snr_low,snr_high,freq_low,freq_high,filter_rate,science_flag,code_version):

    observatory = ifo[0]

    clock_big = time.time()

    command_string = 'python MCIV_start.py -i ' + ifo + ' -t chans -s ' + start_time + ' -e ' + end_time + ' -x ' + channelX_name + \
        ' -y ' + channelY_name + ' -f ' + frame_choice + ' -r ' + sample_rate + \
        ' -l ' + str(snr_low) + ' -u ' + str(snr_high) + \
        ' -g ' + str(freq_low) + ' -j ' + str(freq_high) + \
        ' -a ' + str(filter_rate) + ' -m ' + str(science_flag)

    # End-of-segment buffer to avoid big motions associated with lockloss
    # Define minimum segment length to consider (not including end_buffer)
    end_buffer = 20
    segment_min_length = 256

    # Construct the channel names
    channelX = ifo + ':' + channelX_name
    channelY = ifo + ':' + channelY_name

    # Work out the frame type
    if frame_choice == 'R':
        frame_type = frame_choice
    elif frame_choice == 'RDS':
        frame_type = ifo + '_RDS_R_L1'
    else: frame_type = frame_choice

    # channel1 is vertical; build channel_list in [x,y] format
    channel_list = [channelX, channelY]
    rate = int(sample_rate)

    start_date_string, month, day, year = getDateString(start_time)
    end_date_string, month, day, year = getDateString(end_time)

    # Define where to put all the outputs
    output_path = '/home/dhoak/public_html/detchar/MCIV/' + code_version + '/' + ifo + '/' + start_date_string + '_' + end_date_string + '/' + \
        'chans_' + channelX_name + '_' + channelY_name + '/' 

    # Define a unique filename prefix to avoid confusion with other scripts.
    file_prefix = ifo + '_' + channelX_name + '_' + channelY_name + '_' + start_date_string + '_' + end_date_string + '_'

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    #file_name = 'snr' + str(snr_low)
    file_name = 'veto_inputs'
    logfilename = output_path + 'timemap_logfile.txt'

    num_days = 1 + (int(end_time) - int(start_time))/86400

    # Get the current time, start time & end times of the analysis

    cmd = 'lalapps_tconvert now'
    current_gpstime = commands.getoutput(cmd)
    current_gpstime_string = ''.join(current_gpstime.rsplit('\n'))
    cmd = 'lalapps_tconvert ' + current_gpstime_string
    current_date = commands.getoutput(cmd)
    current_date_string = ''.join(current_date)
    
    cmd = 'lalapps_tconvert ' + start_time
    start_date = commands.getoutput(cmd)
    
    cmd = 'lalapps_tconvert ' + end_time
    end_date = commands.getoutput(cmd)
    
    # Start a log file to store output
    logfile = open(logfilename, 'w',0)
    logfile.write('Multi-channel glitch study: X vs Y\n')
    logfile.write(code_version + '\n\n')
    logfile.write(current_date_string + '\n\n')
    logfile.write(command_string + '\n\n')
    logfile.write('Analyzing glitches for ' + ifo + '.\n')
    logfile.write('Glitch map produced for channels ' + channelX + ' and ' + channelY + '.\n\n')

    logfile.write('~'*100 + '\n\n')
    logfile.write('Parameters for glitch map:\n\n')

    logfile.write('Number of days is ' + str(num_days) + '\n')
    logfile.write('Frame type is ' + frame_type + '\n')
    logfile.write('Sample rate is ' + sample_rate + '\n')
    logfile.write('Low-pass filter rate is ' + str(filter_rate) + '\n')
    logfile.write('End-of-segment buffer is ' + str(end_buffer) + ' seconds.\n')
    logfile.write('Minimum segment length is ' + str(segment_min_length) + ' seconds.\n')
    logfile.write('Range of trigger SNR is from ' + str(snr_low) + ' to ' + str(snr_high) + '\n')
    logfile.write('Range of trigger frequency is from ' + str(freq_low) + ' Hz to ' + str(freq_high) + ' Hz\n\n')
    logfile.write('~'*100 + '\n\n')


    clock1 = time.time()
    if(science_flag):

        S6_epoch = getS6epoch(start_time,end_time)

        logfile.write('Collecting science data from ' + S6_epoch + ': ' + ''.join(start_date) + ' to ' + ''.join(end_date) + '\n\n')

        science_cat1_filename = '/home/dhoak/detchar/S6_segments/' + ifo + '_' + S6_epoch + '_sciencecat1.txt'

        logfile.write('Science-minus-CAT1 segment file is ' + science_cat1_filename + '\n')

        #cmd1 = 'ligolw_segment_query -dqa "' + ifo + ':DMT-SCIENCE" -s ' + start_time + ' -e ' + end_time + ' -o ' + file_prefix + 'segments.xml'
        #cmd2 = 'ligolw_print -t segment -c start_time -c end_time ' + file_prefix + 'segments.xml -d " "'

        #commands.getoutput(cmd1)
        #segment_list = commands.getoutput(cmd2)
        #cmd3 = 'rm '+ file_prefix + 'segments.xml'
        #commands.getoutput(cmd3)
    
        #segments = segment_list.rsplit('\n')

        epoch_start = int(start_time)
        epoch_stop = int(end_time)
        epoch_seg = segmentlist(segmentsUtils.segmentlist_range(epoch_start,epoch_stop,epoch_stop-epoch_start))

        epoch_science = segmentsUtils.fromsegwizard(open(science_cat1_filename), coltype = float, strict=False).coalesce()

        segments = epoch_seg & epoch_science

        if not segments:
            print 'No science data in specified epoch.'
            sys.exit()

        total_seg_time = abs(segments)
        #for segment in segments:
        #    seg_start,seg_stop = segment.split()        
        #    start = int(seg_start)
        #    stop = int(seg_stop)
        #    total_seg_time = total_seg_time + stop - start
            
        #logfile.write('...done.\n\n')

    else:

        logfile.write('Collecting all data from ' + ''.join(start_date) + ' to ' + ''.join(end_date) + '\n\n')

        segments = [start_time + ' ' + end_time]
        total_seg_time = int(end_time) - int(start_time)

    clock2 = time.time()
    elapsed = clock2-clock1
    logfile.write('Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')
    logfile.write('Number of segments is ' + str(len(segments)) + '.  Total time in segments is ' + str(total_seg_time) + ' seconds.\n\n')

    clock1 = time.time()
    logfile.write('Getting channel ranges over time span....\n')

    # Find channel ranges; this determines the endpoints for the plots and histogramming
    # Channel 1 is vertical
    xmin, xmax = channelRange(ifo,segments,channelX_name,end_buffer,segment_min_length)
    ymin, ymax = channelRange(ifo,segments,channelY_name,end_buffer,segment_min_length)

    clock2 = time.time()
    elapsed = clock2-clock1
    logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')

    logfile.write('Range for ' + channelX + ' is from ' + str(xmin) + ' to ' + str(xmax) + '\n')
    logfile.write('Range for ' + channelY + ' is from ' + str(ymin) + ' to ' + str(ymax) + '\n\n')

    xspan = (xmax - xmin)
    yspan = (ymax - ymin)

    xmean = (xmax + xmin)/2.0
    ymean = (ymax + ymin)/2.0

    logfile.write('Middle of ' + channelX + ' range is ' + str(xmean) + '.  Normalization is ' + str(xspan) + '.\n')
    logfile.write('Middle of ' + channelY + ' range is ' + str(ymean) + '.  Normalization is ' + str(yspan) + '.\n\n')

    # Need to think about optimizing this?
    num_bins = 250.0
    bin_spacing = 0.52*2 / num_bins

    bin_edges = arange(-0.52,0.52,bin_spacing)
    
    bins = [bin_edges,bin_edges]
    
    H_all = zeros((len(bin_edges)-1,len(bin_edges)-1))
    
    glitch_count = zeros((len(bin_edges)-1,len(bin_edges)-1))
    glitch_probability = zeros((len(bin_edges)-1,len(bin_edges)-1))

    # Loop through the science data to create the 2-D channel map
    # For each science segment, get the list of frame files that have the data

    logfile.write('~'*100 + '\n\n')
    logfile.write('Looping through science segments...\n\n')

    # Initialize an array to keep track of the median channel value in each segment
    # Make the first entries the channel range means
    segmentMedianValues = array([[ xmean,  ymean]])

    i=0
    total_glitches = 0
    for segment in segments:
        i = i+1
        segment_start = str(int(segment[0]))
        segment_stop = str(int(segment[1]))
        #segment_start,segment_stop = segment.split()
        logfile.write('Segment ' + str(i) + ' ' + segment_start + ' ' + segment_stop + ' ' + str(int(segment_stop)-int(segment_start)) + '\n')

        start = int(segment_start)
        stop = int(segment_stop) - end_buffer
        segment_stop = str(stop)

        #if start == 943172645:
        #    logfile.write('Bad segment.')
        #    continue

        # To handle the glitch files we need to consider segments one day at a time
        # If a segment spans two days we'll chop it up and consider them separately

        # Quick way to tell if segment spans a day
        segmentStartGlitchFile = getGlitchFiles(ifo,segment_start,segment_start)
        segmentStopGlitchFile = getGlitchFiles(ifo,segment_stop,segment_stop)

        segmentSplitFlag = 0
        if segmentStartGlitchFile[0] != segmentStopGlitchFile[0]:
            logfile.write('Segment spans two days.  We will consider each day separately.\n')
            stopGlitchFile = segmentStopGlitchFile[0]
            dayStopTime = stopGlitchFile[-19:-10]  # Note this is only good for S6...
            segmentEndPoints = [segment_start + ' ' + dayStopTime,dayStopTime + ' ' + segment_stop]
            segmentSplitFlag = 1

        else:
            segmentEndPoints = [segment_start + ' ' + segment_stop]

        # Most of the time this loop will be redundant
        j=0
        for seg in segmentEndPoints:
            j = j+1
            seg_start,seg_stop = seg.split()
            start = int(seg_start)
            stop = int(seg_stop)
            
            if stop <= start+segment_min_length:
                logfile.write('Segment is less than ' + str(segment_min_length) + ' seconds.  Skipping segment.\n')
                continue

            logfile.write(str(j) + ' ' + seg_start + ' ' + seg_stop + ' ' + str(stop-start) + '\n')

            glitchFile = getGlitchFiles(ifo[0],seg_start,seg_start)

            if size(glitchFile) > 1:
                logfile.write('\nToo many glitch files!\n')
                print 'Too many glitch files!'
                sys.exit()

            # Open the glitch file and load the events
            f = open(glitchFile[0])
            events = f.readlines()
            f.close()

            clock1 = time.time()
        
            # Get the data for this seg
            # data should have been downsampled to rate
            vector = getChannelVector(ifo,seg_start,seg_stop,frame_type,channel_list,rate)

            xx = ( vector[:,0] - xmean ) / xspan
            yy = ( vector[:,1] - ymean ) / yspan

            # Filter the data

            q = int(rate)/filter_rate
            if q > 1:
                n = 3
                b, a = signal.cheby1(n, 0.05, 0.8 / q)
                x = signal.filtfilt(b, a, xx)
                y = signal.filtfilt(b, a, yy)
            elif q < 1:
                print 'Low pass filter rate is larger than sample rate!'
                sys.exit()
            else:
                x = xx
                y = yy

            m = array([[ median(vector[:,0]),  median(vector[:,1])]])
            segmentMedianValues = vstack((segmentMedianValues,m))

            # data_points should be a Nx2 array with each column corresponding to channel_list
            # have to swap y and x in the order so that the plotting comes out right
            data_points = transpose(vstack((y,x)))

            rows, columns = shape(data_points)

            data_time_vector = arange(start,stop,1.0/rate)

            # Warning if the amount of returned data is not what we expect
            if len(data_time_vector) != len(x):
                logfile.write('\nWrong number of data points\n')
                sys.exit()

            logfile.write('Normalized span of X channel in this segment is ' + str("{0:.3f}".format(min(x))) + ' to ' + str("{0:.3f}".format(max(x))) + '\n')
            logfile.write('Normalized span of Y channel in this segment is ' + str("{0:.3f}".format(min(y))) + ' to ' + str("{0:.3f}".format(max(y))) + '\n')
            logfile.write('Number of data points in this segment is ' + str(rows) + '\n')

            # Bin the data

            Hist, edges = histogramdd(data_points, bins)
            H_all = H_all + Hist

            clock2 = time.time()
            elapsed = clock2-clock1
            logfile.write('...done.  Time to bin data is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n')

            clock1 = time.time()

            num_glitches = 0
            glitch_idx = 0

            # Count the glitches
            for event in events:
                if event[0] == '#':
                    continue

                t, f, d, b, amp, ignore1, ignore2, ignore3 = event.split(' ',7)
                event_time = float(t)
                frequency = float(f)
                amplitude = float(amp)
                snr = sqrt(2*amplitude)

                if (event_time > stop):
                    break

                elif event_time<start or snr<snr_low or snr>snr_high or frequency<freq_low or frequency>freq_high:
                    continue

                else:
                    num_glitches += 1

                    # Find the index of the event in data_points
                    # For simplicity, choose the data point closest in time to the event
                    # This should be ok provided sample rate is large enough

                    # Since glitches are ordered in time, we can use glitch_idx to keep track of the last glitch's position 
                    # in the time vector, and only scan from that point onward?
                    event_idx = (abs(data_time_vector[glitch_idx:]-event_time)).argmin() + glitch_idx
                    glitch_idx = event_idx

                    event_vector = data_points[event_idx,:]

                    # Get the event vector into the right shape to use histogramdd, find the indices of the pixel
                    event_vector.shape = ((2,1))
                    hist_test, edges = histogramdd(event_vector.T, bins)
                    biny, binx = unravel_index(hist_test.argmax(), hist_test.shape)

                    # Add the glitch to the total count - note that x and y are flipped

                    if Hist[biny,binx] == 0:
                        # Sometimes the channel exceeds the limits of the map; in this case the map values at binx, biny will be zero.
                        # This seems to happen very rarely -- for now we just discard these glitches.
                        # Check that the number of glitches within the SNR/frequency band matches the number of glitches in the map!

                        continue

                    else:
                        glitch_count[biny,binx] += 1

            
            # Close the seg loop (not the segment loop)
            logfile.write('Number of glitches in this segment is ' + str(num_glitches) + '\n')

            total_glitches = total_glitches + num_glitches

            clock2 = time.time()
            elapsed = clock2-clock1
            logfile.write('...done.  Time to count glitches is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n')


        logfile.write('\n')


    # Generate histogram of 2-channel parameter space with units of seconds
    time_map = H_all/int(rate)

    # Calculate the glitch rate: number of glitches per bin divided by the time in each bin
    glitch_rate = divide(glitch_count,time_map)

    # Calculate the glitch probability - this will be a more robust measure than glitch rate

    # There are three possible states:
    #
    # 1 - The time map is nonzero, and the glitch count is nonzero -- in this case the rate is well-defined
    #
    # 2 - The time map is nonzero, and the glitch count is zero -- in this case the glitch rate is smaller than we can measure.  These places
    #     should not be vetoed.  We set the glitch count to 0.5 and divide by the time map.  This means that places where the IFO has spent 1 second
    #     but had no glitches will have a glitch rate of 0.5 Hz.  Maybe not optimal, but 1 second is a lot of time when you're sampling at 32Hz or more.
    #
    #     In another approach, we set an artificial glitch count equal to twice the sample frequency (to reflect our lack of knowledge),
    #     and divide by the time map.  (If there is one sample in the bin, the rate will be 2Hz, which will probably be enough for a veto.)
    #
    # 3 - The time map is zero -- in this case we have no knowledge, set to invalid.
    
    # Note that the third case covers places where the IFO hasn't been.  If it goes there in future segments (like when we define vetoes for
    # segments we haven't studied), should we veto that time, or keep it?  Probably better to keep it, but difficult to implement.
    # For now we assume that map is generated using a large, representative sample of IFO time, and previously unseen IFO states are bad.

    for ii in range(len(bin_edges)-1):
        for jj in range(len(bin_edges)-1):

            if glitch_count[ii,jj] > 0 and H_all[ii,jj] > 0:
                glitch_probability[ii,jj] = glitch_count[ii,jj]/H_all[ii,jj]

            elif glitch_count[ii,jj] == 0 and H_all[ii,jj] > 0:
                glitch_probability[ii,jj] = 0.5/H_all[ii,jj]

            else:
                glitch_probability[ii,jj] = NaN  # These invalid numbers will be set for veto in generate_mask


    # Save the data

    datafile = output_path + file_name + '.hdf5'

    logfile.write('Saving data...\n')
    logfile.write('Output file name is ' + datafile + '\n\n')

    channels = [channelX, channelY]
    gpstimes = [int(start_time), int(end_time)]
    snr_limits = [snr_low, snr_high]
    freq_limits = [freq_low, freq_high]

    spans = [xspan, yspan]
    means = [xmean, ymean]

    f = h5py.File(datafile,'w')
    timemap_dset = f.create_dataset('time_map', data=time_map)
    glitchcount_dset = f.create_dataset('glitch_count', data=glitch_count)
    glitchrate_dset = f.create_dataset('glitch_rate', data=glitch_rate)
    glitchprob_dset = f.create_dataset('glitch_probability', data=glitch_probability)
    bins_dset = f.create_dataset('bin_edges', data=bin_edges)
    channels_dset = f.create_dataset('channels', data=channels)
    map_epoch_dset = f.create_dataset('map_epoch', data=gpstimes)
    snr_limits_dset = f.create_dataset('snr_limits', data=snr_limits)
    freq_limits_dset = f.create_dataset('freq_limits', data=freq_limits)
    segMedian_dset = f.create_dataset('segmentMedianValues', data=segmentMedianValues)

    spans_dset = f.create_dataset('spans', data=spans)
    means_dset = f.create_dataset('means', data=means)

    ifo_dset = f.create_dataset('ifo', data=ifo)
    frame_type_dset = f.create_dataset('frame_type', data=frame_type)
    rate_dset = f.create_dataset('rate', data=rate)
    frate_dset = f.create_dataset('filter_rate', data=filter_rate)
    scan_dset = f.create_dataset('scan_type', data='chans')

    # Save the glue.segments object as a numpy array
    segments_dset = f.create_dataset('segments', data=array(segments))

    f.close()

    # generate grid points for surface plots
    X = tile(bin_edges,(num_bins,1))
    Y = tile(bin_edges,(num_bins,1)).T

    logfile.write('~'*100 + '\n\n')
    logfile.write('Making plots...\n')
    
    matplotlib.rcParams.update({'savefig.dpi':250})

    fignum=0
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(log10(time_map),cmap=pylab.cm.jet,origin='lower',interpolation='nearest')
    cbar = pylab.colorbar()
    cbar.set_label('log[seconds]')
    pylab.grid(True)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Two-channel Time Map')
    pylab.savefig(output_path + 'time_map.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(log10(glitch_count),cmap=pylab.cm.jet,origin='lower',interpolation='nearest')
    cbar = pylab.colorbar()
    cbar.set_label('log[#]')
    pylab.grid(True)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Glitch Count')
    pylab.savefig(output_path + 'log_glitch_count.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(log10(glitch_rate),cmap=pylab.cm.jet,origin='lower',interpolation='nearest')
    cbar = pylab.colorbar()
    cbar.set_label('log[Hz]')
    pylab.grid(True)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Log Glitch Rate')
    pylab.savefig(output_path + 'log_glitch_rate.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(log10(glitch_probability),cmap=pylab.cm.jet,origin='lower',interpolation='nearest')
    cbar = pylab.colorbar()
    cbar.set_label('log[#] glitches per bin sample')
    pylab.grid(True)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Two-channel Glitch Probability Map')
    pylab.savefig(output_path + 'glitch_probability_map.png')    
    pylab.close()

    fignum=fignum+1
    fig=pylab.figure(fignum)
    ax = p3.Axes3D(fig)
    surf = ax.plot_surface(X,Y,log10(time_map),cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
    pylab.xlabel(channelX,fontsize=10)
    pylab.ylabel(channelY,fontsize=10)
    ax.set_zlabel('log[seconds]',fontsize=10)
    pylab.title('Time map, ' + channelX + ' vs ' + channelY,fontsize=12)
    pylab.savefig(output_path + 'time_map_surface.png')
    pylab.close()

    fignum=fignum+1
    fig=pylab.figure(fignum)
    ax = p3.Axes3D(fig)
    surf = ax.plot_surface(X,Y,time_map,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
    pylab.xlabel(channelX,fontsize=10)
    pylab.ylabel(channelY,fontsize=10)
    ax.set_zlabel('seconds',fontsize=10)
    pylab.title('Time map, ' + channelX + ' vs ' + channelY,fontsize=12)
    pylab.savefig(output_path + 'time_map_surface_linear.png')
    pylab.close()

    fignum=fignum+1
    fig=pylab.figure(fignum)
    ax = p3.Axes3D(fig)
    surf = ax.plot_surface(X,Y,glitch_count,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
    pylab.xlabel(channelX,fontsize=10)
    pylab.ylabel(channelY,fontsize=10)
    ax.set_zlabel('#',fontsize=10)
    pylab.title('Glitch Count, ' + channelX + ' vs ' + channelY,fontsize=12)
    pylab.savefig(output_path + 'glitch_count_surface_linear.png')
    pylab.close()

    fignum=fignum+1
    fig=pylab.figure(fignum)
    ax = p3.Axes3D(fig)
    surf = ax.plot_surface(X,Y,glitch_rate,cmap=pylab.cm.jet,rstride=8, cstride=8,antialiased=False,linewidth=0.1,vmin=-2.0,vmax=1.0)
    pylab.xlabel(channelX,fontsize=10)
    pylab.ylabel(channelY,fontsize=10)
    ax.set_zlabel('Hz',fontsize=10)
    ax.set_zlim3d(0,40)
    pylab.title('Glitch Rate, ' + channelX + ' vs ' + channelY,fontsize=12)
    pylab.savefig(output_path + 'glitch_rate_surface_linear.png')
    pylab.close()
       
    logfile.write('~'*100 + '\n\n')    

    ttotal = sum(time_map)
    gtotal = sum(glitch_count)
    rtotal = nansum(glitch_rate)
    
    logfile.write('Sum of all time represented in map is ' + str(ttotal) + ' seconds.\n')
    logfile.write('We neglected ' + str(total_seg_time-ttotal) + ' seconds of science data.\n\n')
    
    logfile.write('Number of glitches within SNR threshold is ' + str(total_glitches) +'\n')
    logfile.write('Number of glitch positions included in map is ' + str(gtotal) +'\n\n')
    
    logfile.write('The average rate over this epoch of glitches with ' + str(snr_low) + \
                      ' < SNR < ' + str(snr_high) + ' is ' + str("{0:.2f}".format(gtotal/ttotal)) + '\n\n')
    
    logfile.write('~'*100 + '\n\n')
    
    cmd = 'lalapps_tconvert now'
    current_gpstime = commands.getoutput(cmd)
    current_gpstime_string = ''.join(current_gpstime.rsplit('\n'))
    cmd = 'lalapps_tconvert ' + current_gpstime_string
    current_date = commands.getoutput(cmd)
    current_date_string = ''.join(current_date)
    
    clock3 = time.time()
    elapsed = clock3-clock_big
    logfile.write('Total elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')
    
    logfile.write(current_date_string + '\n')
    logfile.close()

    return output_path
