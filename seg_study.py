#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function reads a segment list and studies the glitches vetoed in science-cat1 time.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
matplotlib.use("Agg", warn=False)
from matplotlib.font_manager import FontProperties
#from matplotlib.dates import DAILY, HOURLY
import pylab
import commands, sys, time, os
import datetime
#from optparse import OptionParser
import h5py
from MCIV_tools import *
from glue import segmentsUtils
from glue.segments import *


def seg_study(ifo,data_path,start_time,end_time,veto_segfile):

    root_directory = '/home/dhoak/public_html/'
    web_directory = root_directory + data_path
    datafile = web_directory + 'veto_outputs.hdf5'

    logfilename = web_directory + 'segstudy_logfile.txt'

    num_days = 1 + (int(end_time) - int(start_time))/86400

    command_string = 'python study_vetoes.py -p ' + data_path + ' -s ' + start_time + ' -e ' + end_time + ' -i ' + ifo + ' -f ' + veto_segfile 

    observatory = ifo[0]

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

    start_date_string, month, day, year = getDateString(start_time)
    end_date_string, month, day, year = getDateString(end_time)


    # Start a log file to store output
    logfile = open(logfilename, 'w',0)
    logfile.write(current_date_string + '\n\n')
    logfile.write(command_string + '\n\n')
    logfile.write('Counting vetoed glitches for ' + ifo + '.\n\n')
    logfile.write('~'*100 + '\n\n')
    logfile.write('Parameters for vetoes and glitches:\n\n')
    logfile.write('Number of days is ' + str(num_days) + '\n')
    logfile.write('Veto segment file is ' + veto_segfile + '\n')
    logfile.write('~'*100 + '\n\n')


    # Get science segments in the specified interval

    S6_epoch = getS6epoch(start_time,end_time)

    logfile.write('Collecting science data from ' + S6_epoch + ': ' + ''.join(start_date) + ' to ' + ''.join(end_date) + '\n\n')

    science_cat1_filename = '/home/dhoak/detchar/S6_segments/' + ifo + '_' + S6_epoch + '_sciencecat1.txt'

    logfile.write('Science-minus-CAT1 segment file is ' + science_cat1_filename + '\n')

    epoch_start = int(start_time)
    epoch_stop = int(end_time)
    epoch_seg = segmentlist(segmentsUtils.segmentlist_range(epoch_start,epoch_stop,epoch_stop-epoch_start))

    epoch_science = segmentsUtils.fromsegwizard(open(science_cat1_filename), coltype = float, strict=False).coalesce()

    segments = epoch_seg & epoch_science

    if not segments:
        print 'No science data in specified epoch.'
        sys.exit()

    total_seg_time = abs(segments)

    logfile.write('...done.\n\n')
    logfile.write('Number of segments is ' + str(len(segments)) + '.  Total time in segments is ' + str(total_seg_time) + ' seconds.\n\n')
    
    logfile.write('Loading veto segments...')

    # Retrieve the segments
    veto_seglist = segmentsUtils.fromsegwizard(open(veto_segfile), coltype = float, strict=False).coalesce()

    # Get durations of veto segments for histogram - note this requires a coalesced veto segment file, not sure how reliable that will be...
    data=genfromtxt(veto_segfile)
    durations=data[:,3]

    # We're only interested in the vetoes that intersect the epoch in question
    vetoes_in_epoch = veto_seglist & epoch_seg

    veto_in_epoch_time = abs(vetoes_in_epoch)

    # Count veto segments for used percentage statistics
    num_vetoes_in_epoch = len(vetoes_in_epoch)
    used_vetoes_in_epoch = 0.0

    logfile.write('done.  Total time in veto segments is ' + str("{0:.2f}".format(abs(veto_seglist))) + ' seconds.\n')
    logfile.write('Time in veto segments in specified epoch is ' + str("{0:.2f}".format(veto_in_epoch_time)) + ' seconds.\n\n')

    logfile.write('~'*100 + '\n\n')
    logfile.write('Looping through science segments...\n\n')


    # Initialize list to hold triggers
    triggers = []

    # initialize lists to hold segment statistics
    segment_efficiency = zeros(len(segments))
    segment_deadtime = zeros(len(segments))
    segment_eff_dead = zeros(len(segments))
    segment_UPV = zeros(len(segments))
    segment_start_gps = [0]*len(segments)   # this is used for plotting


    clock_big = time.time()

    j=0
    for segment in segments:
        j += 1
        seg_start = str(int(segment[0]))
        seg_stop = str(int(segment[1]))
        start = int(seg_start)
        stop = int(seg_stop)

        segment_start_gps[j-1] = start

        # Find the daily glitch files that correspond to this segment
        glitchFiles = getGlitchFiles(ifo[0],seg_start,seg_stop)

        logfile.write('Segment ' + str(j) + ' ' + seg_start + ' ' + seg_stop + ' ' + str(int(seg_stop)-int(seg_start)) + '\n')
        clock1 = time.time()

        # We're only interested in the intersection of the veto segments and the science segment
        science_seg = segmentlist(segmentsUtils.segmentlist_range(start,stop,stop-start))

        #vetoes_in_segment = science_seg & vetoes_in_epoch
        vetoes_in_segment = getSegmentOverlap(science_seg,vetoes_in_epoch)
        if vetoes_in_segment:
            num_veto_segments, columns = shape(vetoes_in_segment)
        else:
            num_veto_segments = 0

        # Count the glitches
        num_glitches = 0
        num_vetoed_glitches = 0

        used_veto_segments = zeros((num_veto_segments))

        for gfile in glitchFiles:

            # Open the glitch file and load the events
            f = open(gfile)
            events = f.readlines()
            f.close()

            # Start a counter for the index of the veto segment; this saves time
            veto_index = 0

            for event in events:
                if event[0] == '#':
                    continue

                t, f, d, b, amp, ignore1, ignore2, ignore3 = event.split(' ',7)
                event_time = float(t)
                amplitude = float(amp)
                snr = sqrt(2*amplitude)
    
                if (event_time > stop):
                    break

                elif (event_time < start) or (snr < 5):
                    continue

                else:
                    num_glitches += 1

                    event_vector = [event_time, float(f), float(d), float(b), float(snr), False]

                    if num_veto_segments:

                        # We shouldn't need to scan through the entire veto list; can assume it's coalesced and chronologically ordered?

                        for vidx in range(num_veto_segments):
                        #for vidx in range(veto_index,num_veto_segments):

                            if event_time >= vetoes_in_segment[vidx][0] and event_time <= vetoes_in_segment[vidx][1]:

                                num_vetoed_glitches = num_vetoed_glitches + 1
                                event_vector[5] = True

                                used_veto_segments[vidx] = 1

                                veto_index = vidx

                    triggers.append(event_vector)


        # Calculate the veto statistics for this science segment
        # index starts at zero, but we incremented the counter above, so subtract 1 from the index
                    
        if num_glitches == 0:
            segment_efficiency[j-1] = 0.0
        else:
            segment_efficiency[j-1] = num_vetoed_glitches/float(num_glitches)

        segment_deadtime[j-1] = getSegmentSum(vetoes_in_segment)/float(stop-start)

        if num_veto_segments == 0:
            segment_UPV[j-1] = 0
            segment_eff_dead[j-1] = 0
        else:
            segment_UPV[j-1] = sum(used_veto_segments)/float(num_veto_segments)
            segment_eff_dead[j-1] = segment_efficiency[j-1]/segment_deadtime[j-1]

        used_vetoes_in_epoch += sum(used_veto_segments)

        # Finished with this science segment
        logfile.write('Number of veto segments in this segment is ' + str(num_veto_segments) + '\n')
        logfile.write('Number of triggers in this segment is ' + str(num_glitches) + '\n')
        logfile.write('Number of vetoed triggers in this segment is ' + str(num_vetoed_glitches) + '\n')
        logfile.write('Number of veto segments that contain a glitch: ' + str(sum(used_veto_segments)) + '\n')
        
        logfile.write('Veto efficiency for this segment is ' + str("{0:.2f}".format(segment_efficiency[j-1])) + '\n')
        logfile.write('Veto deadtime for this segment is ' + str("{0:.2f}".format(segment_deadtime[j-1])) + '\n')
        logfile.write('Efficiency / deadtime is ' + str("{0:.2f}".format(segment_eff_dead[j-1])) + '\n')
        logfile.write('Use percentage for veto segments in this segment is ' + str("{0:.2f}".format(segment_UPV[j-1])) + '\n')
        
        clock2 = time.time()
        elapsed = clock2-clock1
        logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')



    # Finished with all science segments; collect epoch statistics

    logfile.write('~'*100 + '\n\n')

    cmd = 'lalapps_tconvert now'
    current_gpstime = commands.getoutput(cmd)
    current_gpstime_string = ''.join(current_gpstime.rsplit('\n'))
    cmd = 'lalapps_tconvert ' + current_gpstime_string
    current_date = commands.getoutput(cmd)
    current_date_string = ''.join(current_date)
    
    logfile.write(current_date_string + '\n\n')


    glitch_list = vstack((triggers))

    freq = glitch_list[:,1]
    duration = glitch_list[:,2]
    bandwidth = glitch_list[:,3]
    SNR = glitch_list[:,4]
    vetoed = glitch_list[:,5]
    
    vetoed_freq = freq[where(vetoed)]
    vetoed_duration = duration[where(vetoed)]
    vetoed_bandwidth = bandwidth[where(vetoed)]
    vetoed_SNR = SNR[where(vetoed)]
    
    glitch_counter = ones(shape(vetoed))

    total_glitches = len(glitch_list)
    total_vetoed_glitches = len(vetoed_SNR)

    efficiency = total_vetoed_glitches/float(total_glitches)
    deadtime = veto_in_epoch_time/total_seg_time


    logfile.write('Calculating per-SNR efficiency...\n')
    clock1 = time.time()

    # Calculate efficiency for various SNR thresholds

    snr_level = [5, 6, 8, 10, 20, 50, 100]

    glitch_count = zeros(shape(snr_level))
    vetoed_glitch_count = zeros(shape(snr_level))

    for i in range(len(snr_level)):

        glitch_count[i] = sum(glitch_counter[where(SNR>=snr_level[i])])
        vetoed_glitch_count[i] = sum(glitch_counter[where(vetoed_SNR>=snr_level[i])])
    
    efficiency_per_SNR = vetoed_glitch_count / glitch_count

    EFF_DED = efficiency_per_SNR / deadtime

    clock2 = time.time()
    elapsed = clock2-clock1
    logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')


    # Quantify safety

    logfile.write('Calculating safety statistics...\n')
    clock1 = time.time()

    num_inj, num_vetoed_inj, safety_significance = hveto_safety(ifo,vetoes_in_epoch,start_time,end_time,deadtime)

    clock2 = time.time()
    elapsed = clock2-clock1
    logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')

    logfile.write('Collecting segment start times...\n')
    clock1 = time.time()

    study_epoch = [int(start_time), int(end_time)]

    # Convert the segment start times into datetime format so we can plot things
    segment_start_dates = []
    for i in range(len(segment_start_gps)):

        date_string, month, day, year, hour, minute, second = getDateString(str(segment_start_gps[i]),hms=True)

        date_object = datetime.datetime.strptime(date_string + ' ' + hour + ' ' + minute + ' ' + second,"%d %b %Y %H %M %S")
        segment_start_dates.append(matplotlib.dates.date2num(date_object))


    clock2 = time.time()
    elapsed = clock2-clock1
    logfile.write('...done.  Elapsed time is ' + str("{0:.2f}".format(elapsed)) + ' seconds.\n\n')


    # Save the glitch tables and statistics to the data file

    f = h5py.File(datafile,'a')

    #if 'glitch_list' in f:
    #    del f['glitch_list']

    if 'study_epoch' in f:
        del f['study_epoch']

    if 'study_veto_file' in f:
        del f['study_veto_file']

    #if 'veto_in_epoch_time' in f:
    #    del f['veto_in_epoch_time']

    if 'study_deadtime' in f:
        del f['study_deadtime']

    if 'segment_efficiency' in f:
        del f['segment_efficiency']

    if 'segment_deadtime' in f:
        del f['segment_deadtime']

    if 'segment_UPV' in f:
        del f['segment_UPV']

    if 'study_segments' in f:
        del f['study_segments']

    # Note that this is a boolean array of length study_segments, identifying which segments contain a glitch
    if 'used_vetoes_in_epoch' in f:
        del f['used_vetoes_in_epoch']

    if 'num_vetoes_in_epoch' in f:
        del f['num_vetoes_in_epoch']

    if 'num_study_glitches' in f:
        del f['num_study_glitches']

    if 'efficiency_per_SNR' in f:
        del f['efficiency_per_SNR']

    #if 'num_glitches' in f:
    #    del f['num_glitches']

    if 'num_vetoed_glitches' in f:
        del f['num_study_vetoed_glitches']

    if 'num_inj' in f:
        del f['num_inj']

    if 'num_vetoed_inj' in f:
        del f['num_vetoed_inj']

    if 'safety_significance' in f:
        del f['safety_significance']

    #glist_dset = f.create_dataset('glitch_list', data=glitch_list)
    vetofile_dset = f.create_dataset('study_veto_file', data=veto_segfile)
    glitchtime_dset = f.create_dataset('study_epoch', data=study_epoch)
    vetotime_dset = f.create_dataset('study_deadtime', data=veto_in_epoch_time)
    #deadtime_dset = f.create_dataset('study_deadtime', data=deadtime)
    segment_eff_dset = f.create_dataset('segment_efficiency', data=segment_efficiency)
    segment_dtime_dset = f.create_dataset('segment_deadtime', data=segment_deadtime)
    segment_UPV_dset = f.create_dataset('segment_UPV', data=segment_UPV)
    segment_dset = f.create_dataset('study_segments', data=segments)
    
    used_vetoes_dset = f.create_dataset('used_vetoes_in_epoch', data=used_vetoes_in_epoch)
    num_vetoes_dset = f.create_dataset('num_vetoes_in_epoch', data=num_vetoes_in_epoch)
    num_glitches_dset = f.create_dataset('num_study_glitches', data=total_glitches)
    num_vetoed_glitches_dset = f.create_dataset('num_study_vetoed_glitches', data=total_vetoed_glitches)

    SNR_efficiency_dset = f.create_dataset('efficiency_per_SNR', data=efficiency_per_SNR)
    num_inj_dset = f.create_dataset('num_inj', data=num_inj)
    num_vetoed_inj_dset = f.create_dataset('num_vetoed_inj', data=num_vetoed_inj)
    safety_dset = f.create_dataset('safety_significance', data=safety_significance)

    f.close()

    logfile.write('~'*100 + '\n\n')
    logfile.write('Number of veto segments in this epoch is ' + str(num_vetoes_in_epoch) + '.\n\n')
    logfile.write('Veto deadtime in this epoch is ' + str("{0:.2f}".format(deadtime)) + '.\n\n')

    logfile.write('Number of events in this epoch is ' + str(total_glitches) + '.\n')
    logfile.write('Number of vetoed events in this epoch is ' + str(total_vetoed_glitches) + '.\n')
    logfile.write('Number of veto segments with a glitch is ' + str(int(used_vetoes_in_epoch)) + '.\n\n')
    logfile.write('Veto efficiency in this epoch is ' + str("{0:.2f}".format(efficiency)) + '.\n\n')

    logfile.write('Veto [efficiency/deadtime] in this epoch is ' + str("{0:.2f}".format(efficiency/deadtime)) + '.\n\n')
    logfile.write('Used percentage (% of veto segments with a glitch) is ' + str("{0:.2f}".format(used_vetoes_in_epoch/num_vetoes_in_epoch)) + '.\n\n')

    logfile.write('Number of hardware injections in this epoch is ' + str(num_inj) + '.\n')
    logfile.write('Number of vetoed hardware injections is ' + str(num_vetoed_inj) + '.\n')
    logfile.write('Safety significance is ' + str("{0:.2f}".format(safety_significance)) + '.\n')

    logfile.write('~'*100 + '\n\n')
    logfile.write('Generating plots...')

    fontP = FontProperties()
    fontP.set_size('small')

    matplotlib.rcParams.update({'savefig.dpi':250})

    fignum=0
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.loglog(freq,SNR,'b.')
    pylab.loglog(vetoed_freq,vetoed_SNR,'r.')
    pylab.xlabel('Frequency [Hz]',fontsize=14)
    pylab.ylabel('SNR',fontsize=14)
    pylab.grid(True, which='major', linestyle=':')
    pylab.grid(True, which='minor', linestyle=':')
    pylab.axis([40,2050,4.0,max(SNR)*1.1])
    pylab.legend(('All Events','Vetoed Events'),loc=1,prop=fontP,fancybox=True)
    pylab.title('SNR-Frequency Scatterplot')
    pylab.savefig(web_directory + 'strengthFreq.png')
    pylab.close()


    fignum=fignum+1
    pylab.figure(fignum)

    ax1 = pylab.subplot(3,1,1)
    pylab.plot_date(segment_start_dates,segment_efficiency,'r*',markersize=10)
    pylab.plot_date(segment_start_dates,segment_deadtime,'ks')
    # pylab.plot_date(segment_start_dates,segment_UPV,'b^',markersize=6)
    pylab.grid(True, which='major', linestyle=':')
    # pylab.grid(True, which='minor', linestyle=':')
    pylab.legend(('Efficiency','Deadtime','Use Percentage'),loc=1,prop=fontP,bbox_to_anchor=(1.05, 1.2),fancybox=True)
    # pylab.yscale('log')
    pylab.ylim(0.0,1.0)
    pylab.ylabel('Fraction Vetoed',fontsize=12)

    pylab.title('Per-segment veto statistics',fontsize=13)

    # loc1 = ax1.xaxis.get_major_locator()
    # loc1.maxticks[DAILY] = 12
    # loc1.maxticks[HOURLY] = 8
    ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b%d\n%Y'))
    pylab.xticks(fontsize=10)
    pylab.yticks(fontsize=11)
    pylab.tick_params(axis='x',labelbottom='off')

    ax2 = pylab.subplot(3,1,2)
    pylab.plot_date(segment_start_dates,segment_UPV,'b^',markersize=6)
    pylab.grid(True, which='major', linestyle=':')
    # pylab.grid(True, which='minor', linestyle=':')
    # pylab.legend(('Efficiency','Deadtime','Use Percentage'),loc=1,prop=fontP,bbox_to_anchor=(1.05, 1.1),fancybox=True)
    # pylab.yscale('log')
    pylab.ylim(0.0,1.0)
    pylab.ylabel('Use Percentage',fontsize=12)

    # loc2 = ax2.xaxis.get_major_locator()
    # loc2.maxticks[DAILY] = 12
    # loc2.maxticks[HOURLY] = 8
    ax2.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b%d\n%Y'))
    pylab.xticks(fontsize=10)
    pylab.yticks(fontsize=11)
    pylab.tick_params(axis='x',labelbottom='off')
    
    ax3 = pylab.subplot(3,1,3)
    pylab.plot_date(segment_start_dates,segment_efficiency/segment_deadtime,'bs')
    pylab.grid(True, which='major', linestyle=':')
    # pylab.grid(True, which='minor', linestyle=':')
    # pylab.yscale('log')
    pylab.ylabel('Eff/dead',fontsize=12)

    # loc3 = ax3.xaxis.get_major_locator()
    # loc3.maxticks[DAILY] = 12
    # loc3.maxticks[HOURLY] = 8
    ax3.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b%d\n%Y'))
    pylab.xticks(fontsize=10)
    pylab.yticks(fontsize=11)
    pylab.xlabel('Segment Start Time',fontsize=12)
    pylab.savefig(web_directory + 'segment_eff_dead_stats.png')
    pylab.close()


    fignum=fignum+1
    pylab.figure(fignum)

    dmin = 0.01
    dmax = 60.0

    nbins = 70
    logdmin=log10(dmin)
    logdmax=log10(dmax)
    dbinwidth = (logdmax-logdmin)/nbins
    dbins=[0.0]*(nbins+1)
    for i in range(len(dbins)):
        dbins[i] = pow(10,logdmin+i*dbinwidth)

    n, b = histogram(durations, dbins)

    plotbins = zeros((len(b)*2)) + 0.001
    plotcounts = zeros((len(b)*2)) + 0.001

    for i in range(len(b)):
        plotbins[2*i] = b[i] - b[i]*0.0003
        plotbins[2*i+1] = b[i] + b[i]*0.0003

    for i in range(len(n)-1):
        plotcounts[2*i+1] = n[i] + 0.001
        plotcounts[2*i+2] = n[i] + 0.001

    pylab.plot(plotbins,plotcounts,'k-',linewidth=1.0)

    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Duration [sec]',fontsize=14)
    pylab.grid(True, which='major', linestyle=':')
    pylab.grid(True, which='minor', linestyle=':')
    pylab.xlim(dmin,dmax)
    pylab.ylim(0.7,1.2*max(plotcounts))
    pylab.title('Veto Segment Duration')
    pylab.savefig(web_directory + 'studyvetoDurationHist.png')
    pylab.close()



    fignum=fignum+1
    pylab.figure(fignum)

    fmin = 40
    fmax = 2040

    # Construct logspace frequency bins
    nbins = 120
    logfmin=log10(fmin)
    logfmax=log10(fmax)
    fbinwidth = (logfmax-logfmin)/nbins
    fbins=[0.0]*(nbins+1)
    for i in range(len(fbins)):
        fbins[i] = pow(10,logfmin+i*fbinwidth)

    # Histogram total glitches
    n, b = histogram(freq, fbins)

    plotbins = zeros((len(b)*2)) + 0.001
    plotcounts = zeros((len(b)*2)) + 0.001

    for i in range(len(b)):
        plotbins[2*i] = b[i] - b[i]*0.0003
        plotbins[2*i+1] = b[i] + b[i]*0.0003

    for i in range(len(n)-1):
        plotcounts[2*i+1] = n[i] + 0.001
        plotcounts[2*i+2] = n[i] + 0.001

    pylab.plot(plotbins,plotcounts,'k-',linewidth=1.0)

    # Histogram vetoed glitches
    vn, vb = histogram(vetoed_freq, fbins)

    vplotbins = zeros((len(vb)*2)) + 0.001
    vplotcounts = zeros((len(vb)*2)) + 0.001

    for i in range(len(vb)):
        vplotbins[2*i] = vb[i] - vb[i]*0.0003
        vplotbins[2*i+1] = vb[i] + vb[i]*0.0003

    for i in range(len(vn)-1):
        vplotcounts[2*i+1] = vn[i] + 0.001
        vplotcounts[2*i+2] = vn[i] + 0.001

    pylab.plot(vplotbins,vplotcounts,'r-',linewidth=1.0)
        
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Frequency [Hz]',fontsize=14)
    pylab.grid(True, which='major', linestyle=':')
    pylab.grid(True, which='minor', linestyle=':')
    pylab.xlim(fmin,fmax)
    pylab.ylim(0.7,1.2*max(plotcounts))
    pylab.legend(('All Events','Vetoed Events'),loc=1,prop=fontP,fancybox=True)
    pylab.title('Central Frequency Histogram')
    pylab.savefig(web_directory + 'frequencyHist.png')
    pylab.close()


    freq_percent_vetoed = vplotcounts / plotcounts
    for i in range(len(freq_percent_vetoed)):
        if plotcounts[i]==0.001:
            freq_percent_vetoed[i] = 0.0

    fignum=fignum+1
    pylab.figure(fignum)

    pylab.plot(plotbins,freq_percent_vetoed,'k-',linewidth=1.0)

    pylab.xscale('log')
    pylab.xlabel('Frequency [Hz]',fontsize=14)
    pylab.grid(True, which='major', linestyle=':')
    pylab.grid(True, which='minor', linestyle=':')
    pylab.axis([fmin,fmax,0.0,1.1])
    pylab.title('Fraction Events Vetoed per Frequency Bin')
    pylab.savefig(web_directory + 'vetoPercentFreqBin.png')
    pylab.close()


    fignum=fignum+1
    pylab.figure(fignum)

    zmin = 5
    zmax = 5000

    # Construct logspace SNR bins
    nbins = 120
    logzmin=log10(zmin)
    logzmax=log10(zmax)
    zbinwidth = (logzmax-logzmin)/nbins
    zbins=[0.0]*(nbins+1)
    for i in range(len(zbins)):
        zbins[i] = pow(10,logzmin+i*zbinwidth)

    # Histogram total glitches
    n, b = histogram(SNR, zbins)

    plotbins = zeros((len(b)*2)) + 0.001
    plotcounts = zeros((len(b)*2)) + 0.001

    for i in range(len(b)):
        plotbins[2*i] = b[i] - b[i]*0.0003
        plotbins[2*i+1] = b[i] + b[i]*0.0003

    for i in range(len(n)-1):
        plotcounts[2*i+1] = n[i] + 0.001
        plotcounts[2*i+2] = n[i] + 0.001

    pylab.plot(plotbins,plotcounts,'k-',linewidth=1.0)

    # Histogram vetoed glitches
    vn, vb = histogram(vetoed_SNR, zbins)

    vplotbins = zeros((len(vb)*2)) + 0.001
    vplotcounts = zeros((len(vb)*2)) + 0.001

    for i in range(len(vb)):
        vplotbins[2*i] = vb[i] - vb[i]*0.0003
        vplotbins[2*i+1] = vb[i] + vb[i]*0.0003

    for i in range(len(vn)-1):
        vplotcounts[2*i+1] = vn[i] + 0.001
        vplotcounts[2*i+2] = vn[i] + 0.001

    pylab.plot(vplotbins,vplotcounts,'r-',linewidth=1.0)

    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('SNR',fontsize=14)
    pylab.grid(True, which='major', linestyle=':')
    pylab.grid(True, which='minor', linestyle=':')
    pylab.xlim(zmin,zmax*1.2)
    pylab.ylim(0.7,1.2*max(plotcounts))
    pylab.legend(('All Events','Vetoed Events'),loc=1,prop=fontP,fancybox=True)
    pylab.title('SNR Histogram')
    pylab.savefig(web_directory + 'snrHist.png')
    pylab.close()


    logfile.write('done.\n\n')

    logfile.write('~'*100 + '\n\n')

    cmd = 'lalapps_tconvert now'
    current_gpstime = commands.getoutput(cmd)
    current_gpstime_string = ''.join(current_gpstime.rsplit('\n'))
    cmd = 'lalapps_tconvert ' + current_gpstime_string
    current_date = commands.getoutput(cmd)
    current_date_string = ''.join(current_date)
    
    logfile.write(current_date_string + '\n')
    logfile.close()



    return efficiency_per_SNR, deadtime, num_inj, num_vetoed_inj, safety_significance
