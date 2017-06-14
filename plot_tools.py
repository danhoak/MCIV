#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script has some plotting routines for MCIV visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
matplotlib.use("Agg")
import pylab
#import commands, sys, time, os
from optparse import OptionParser
import h5py
from MCIV_tools import *
#from glue import segmentsUtils
#from glue.segments import *
import mpl_toolkits.mplot3d.axes3d as p3
#from mpl_toolkits.mplot3d import Axes3D

usage = """usage: %prog [options]

This script will load a two-channel veto mask and generate some plots.

Required options:

--path

Example command line: 

python plot_tools.py -p detchar/dMCIV/beta/H1/01Oct2009_01Dec2009/SUS-ETMX_COIL_UL/

"""

parser = OptionParser(usage=usage)

parser.add_option("-p", "--path", action="store", type="str",\
                      help="path to directory with veto mask file (assumes /home/dhoak/public_html root directory)")


(options, args) = parser.parse_args()

data_path = options.path

root_directory = '/home/dhoak/public_html/'
web_directory = root_directory + data_path
output_path = web_directory
datafile = web_directory + 'snr6_data.hdf5'


# Load the necessary data fields from the data file in /path

g = h5py.File(datafile,'r')

dataset = g['ifo']
ifo = dataset[...]

dataset = g['bin_edges']
bin_edges = dataset[...]

dataset = g['time_map']
time_map = dataset[...]

dataset = g['glitch_count']
glitch_count = dataset[...]

dataset = g['glitch_rate']
glitch_rate = dataset[...]

dataset = g['channels']
channels = dataset[...]
channelX = channels[0]
channelY = channels[1]

dataset = g['spans']
spans = dataset[...]
xspan = spans[0]
yspan = spans[1]

dataset = g['means']
means = dataset[...]
xmean = means[0]
ymean = means[1]

dataset = g['frame_type']
frame_type = dataset[...]

g.close()

Nbins = shape(bin_edges)[0]

log_map = log10(time_map)

X = tile(bin_edges,(Nbins,1))
Y = tile(bin_edges,(Nbins,1)).T

fignum=0
fignum=fignum+1
fig=pylab.figure(fignum)
ax = p3.Axes3D(fig)
surf = ax.plot_surface(X,Y,log10(time_map),cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
pylab.xlabel(channelX,fontsize=10)
pylab.ylabel(channelY,fontsize=10)
ax.set_zlabel('log[seconds]',fontsize=10)
pylab.title('Time map, ' + channelX + ' vs ' + channelY,fontsize=12)
pylab.savefig(output_path + 'time_map_surface.png')

fignum=fignum+1
fig=pylab.figure(fignum)
ax = p3.Axes3D(fig)
surf = ax.plot_surface(X,Y,time_map,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
pylab.xlabel(channelX,fontsize=10)
pylab.ylabel(channelY,fontsize=10)
ax.set_zlabel('log[seconds]',fontsize=10)
pylab.title('Time map, ' + channelX + ' vs ' + channelY,fontsize=12)
pylab.savefig(output_path + 'time_map_surface_linear.png')

fignum=fignum+1
fig=pylab.figure(fignum)
ax = p3.Axes3D(fig)
surf = ax.plot_surface(X,Y,log10(glitch_count),cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
#surf = ax.plot_surface(X,Y,glitch_count,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
pylab.xlabel(channelX,fontsize=10)
pylab.ylabel(channelY,fontsize=10)
ax.set_zlabel('log[#]',fontsize=10)
pylab.title('Glitch Count, ' + channelX + ' vs ' + channelY,fontsize=12)
pylab.savefig(output_path + 'glitch_count_surface.png')

fignum=fignum+1
fig=pylab.figure(fignum)
ax = p3.Axes3D(fig)
surf = ax.plot_surface(X,Y,glitch_count,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
#surf = ax.plot_surface(X,Y,glitch_count,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1,vmin=0.0)
pylab.xlabel(channelX,fontsize=10)
pylab.ylabel(channelY,fontsize=10)
ax.set_zlabel('#',fontsize=10)
pylab.title('Glitch Count, ' + channelX + ' vs ' + channelY,fontsize=12)
pylab.savefig(output_path + 'glitch_count_surface_linear.png')

fignum=fignum+1
fig=pylab.figure(fignum)
ax = p3.Axes3D(fig)
surf = ax.plot_surface(X,Y,log10(glitch_rate),cmap=pylab.cm.jet,rstride=8, cstride=8,antialiased=False,linewidth=0.1,vmin=-2.0,vmax=1.0)
#surf = ax.plot_surface(X,Y,glitch_rate,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1)
pylab.xlabel(channelX,fontsize=10)
pylab.ylabel(channelY,fontsize=10)
ax.set_zlabel('log[Hz]',fontsize=10)
pylab.title('Glitch Rate, ' + channelX + ' vs ' + channelY,fontsize=12)
pylab.savefig(output_path + 'glitch_rate_surface.png')

fignum=fignum+1
fig=pylab.figure(fignum)
ax = p3.Axes3D(fig)
surf = ax.plot_surface(X,Y,glitch_rate,cmap=pylab.cm.jet,rstride=8, cstride=8,antialiased=False,linewidth=0.1,vmin=-2.0,vmax=1.0)
#surf = ax.plot_surface(X,Y,glitch_rate,cmap=pylab.cm.jet,rstride=4, cstride=4,antialiased=False,linewidth=0.1)
pylab.xlabel(channelX,fontsize=10)
pylab.ylabel(channelY,fontsize=10)
ax.set_zlabel('Hz',fontsize=10)
ax.set_zlim3d(0,40)
pylab.title('Glitch Rate, ' + channelX + ' vs ' + channelY,fontsize=12)
pylab.savefig(output_path + 'glitch_rate_surface_linear.png')

fignum=fignum+1
pylab.figure(fignum)
pylab.imshow(time_map,cmap=pylab.cm.jet,origin='lower',interpolation='nearest')
cbar = pylab.colorbar()
cbar.set_label('seconds')
pylab.grid(True)
pylab.xticks(visible=False)
pylab.yticks(visible=False)
pylab.xlabel(channelX,fontsize=14)
pylab.ylabel(channelY,fontsize=14)
pylab.title('Two-channel Time Map')
pylab.savefig(output_path + 'time_map_linear_scale.png')

