#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the gateway wrapper function for the MCIV/dMCIV channel mapping routines.
# Calls deriv_map or glitch_map as functions, with appropriate options.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys, commands
#from optparse import OptionParser

ifo = 'H1'

start_time = '937785625'
end_time = '940032000'

#start_time = '937785625'
#end_time = '937885625'

cfile = open('/home/dhoak/detchar/MCIV/channel_list/H1/H1_S6B_chans_prelim.txt')
channels = cfile.readlines()


subfile = open('MCIV_chans.sub','w')
subfile.write('universe = vanilla\n')
subfile.write('executable = /home/dhoak/detchar/MCIV/v0.1/MCIV_start.py\n\n')

#subfile.write('arguments = " -i H1 -t chans -s 937785625 -e 937885625 -x ASC-WFS2_QP -y OMC-QPD3_P_OUT_DAQ -f RDS -r 128 -l 6 -u 1000 -g 40 -j 2048 -m 1 "\n')
subfile.write('arguments = " $(macroarguments) "\n')
#subfile.write('environment = USER=$ENV(USER);\n')
subfile.write('getenv = true\n')
subfile.write('priority = 0\n')
#subfile.write('Requirements = Memory >= 1000\n')
subfile.write('request_memory = 1000\n')
subfile.write('log = /usr1/dhoak/log/MCIV_test_0123456789\n\n')

subfile.write('error = logs/MCIV-chans-$(cluster)-$(process).err\n')
subfile.write('output = logs/MCIV-chans-$(cluster)-$(process).out\n')
subfile.write('notification = never\n')
subfile.write('queue 1\n')
subfile.close()


dagfile = open('H1_test.dag','w')

i=0
for channel_pair in channels:

    i += 1
    cpair = channel_pair.rsplit('\n')[0]
    [channelX, channelY, rate] = cpair.split(' ')

    #cmd = 'condor_run python MCIV_start.py -i ' + ifo + ' -t chans -s ' + start_time + ' -e ' + end_time + \
    #    ' -x ' + channelX + ' -y ' + channelY + ' -f RDS -r 128 -l 6 -u 1000 -g 40 -j 2048 -m 1'

    inputs = '-i ' + ifo + ' -t chans -s ' + start_time + ' -e ' + end_time + ' -x ' + channelX + ' -y ' + channelY + \
        ' -f RDS -r ' + rate + ' -l 6 -u 1000 -g 40 -j 2048 -m 1'

    #print channelX, channelY
    #print cmd

    dagfile.write('JOB mciv_chans_' + ifo + '_' + channelX + '_' + channelY + ' MCIV_chans.sub\n')
    dagfile.write('RETRY mciv_chans_' + ifo + '_' + channelX + '_' + channelY + ' 1\n')
    dagfile.write('VARS mciv_chans_' + ifo + '_' + channelX + '_' + channelY + ' macroarguments="' + inputs + '"\n')


dagfile.close()
