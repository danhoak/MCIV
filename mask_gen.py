#! /usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to generate veto masks from glitch map output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg",warn=False)
import pylab
import h5py
import commands, sys, os
#from optparse import OptionParser


def mask_gen(output_path,fixed_threshold):
             
    user_threshold = float(fixed_threshold)


    # a function to change the values of the mask based on neighbors
    # from numpy game of life, http://www.loria.fr/~rougier/teaching/numpy/numpy.html
    def mask_iterate(Z):
        # Count neighbours
        # this includes diagonals; note structure ;)
        N = (Z[0:-2,0:-2] + Z[0:-2,1:-1] + Z[0:-2,2:] +
             Z[1:-1,0:-2]                + Z[1:-1,2:] +
             Z[2:  ,0:-2] + Z[2:  ,1:-1] + Z[2:  ,2:])

        birth = (N>4) & (Z[1:-1,1:-1]==0)
        survive = (N>3) & (Z[1:-1,1:-1]==1)

        # Count neighbours, neglecting diagonals
        #N = (Z[0:-2,1:-1] + Z[1:-1,0:-2] + Z[1:-1,2:] + Z[2:  ,1:-1] )
        
        #birth = (N>2) & (Z[1:-1,1:-1]==0)
        #survive = (N>1) & (Z[1:-1,1:-1]==1)
    
        Z[...] = 0
        Z[1:-1,1:-1][(birth | survive)] = 1
        return Z


    # Build the location of the data file we will use to generate the masks

    web_directory = output_path
    datafile = web_directory + 'veto_inputs.hdf5'

    # Copy the data file for safety...

    cmd = 'cp ' + datafile + ' ' + web_directory + 'veto_inputs_copy.hdf5'
    commands.getoutput(cmd)

    # Read in the necessary data products to generate masks / plots
    
    g = h5py.File(datafile,'r')

    dataset = g['time_map']
    time_map = dataset[...]
    
    dataset = g['glitch_count']
    glitch_count = dataset[...]
    
    dataset = g['glitch_rate']
    glitch_rate = dataset[...]
    
    dataset = g['glitch_probability']
    glitch_probability = dataset[...]
    
    dataset = g['channels']
    channels = dataset[...]
    channelX = channels[0]
    channelY = channels[1]
    
    dataset = g['segmentMedianValues']
    segmentMedianValues = dataset[...]
    
    dataset = g['spans']
    spans = dataset[...]
    xspan = spans[0]
    yspan = spans[1]
    
    dataset = g['means']
    means = dataset[...]
    xmean = means[0]
    ymean = means[1]
    
    dataset = g['rate']
    rate = dataset[...]
    
    g.close()

    segmentMedians = segmentMedianValues[1:]

    xlim = [xmean - xspan/2.0,xmean+xspan/2.0]
    ylim = [ymean - yspan/2.0,ymean+yspan/2.0]

    # Build an array of glitch probabilties to quantify efficiency and deadtime
    # For now the limits are hard-coded, with a floor of 1 in 10,000

    dmin = 1e-4
    dmax = 10.0

    nbins = 200
    logdmin=log10(dmin)
    logdmax=log10(dmax)
    dbinwidth = (logdmax-logdmin)/nbins
    dbins=[0.0]*(nbins+1)
    for i in range(len(dbins)):
        dbins[i] = pow(10,logdmin+i*dbinwidth)

    probability_threshold = dbins

    # Calculate the average glitch rate, and the probability threshold that corresponds to 10x the average glitch rate

    total_time = sum(sum(time_map))
    total_glitches = sum(sum(glitch_count))
    
    mean_glitch_rate = total_glitches/total_time
    
    ten_rate = 10.0*mean_glitch_rate
    ten_probability = ten_rate / rate
    
    # Initialize arrays for veto threshold metrics
    
    deadtime = empty(shape(probability_threshold))
    efficiency = empty(shape(probability_threshold))
    EFF_DED = empty(shape(probability_threshold))

    for idx in range(len(probability_threshold)):

        # First, set any NaN's in the glitch rate matrix to zeros.  These places will be vetoed.
        fixed_glitch_prob = ma.fix_invalid(glitch_probability,fill_value = 0)

        # Now, generate a mask where any probability less than 1e-8 is true.  These places will be vetoed.
        # This step is intended to catch the places that were invalid in the probability map -- these are places in the parameter space where the IFO
        # has not been observed.  We assume that we have observed a representative sample of IFO states, and new states are bad.
        # In the future we may find that 1e-8 is too high for very long, glitch-free datasets.
        glitch_mask_nonzero = ma.masked_less(fixed_glitch_prob.data,1e-8)
        
        # ...and another mask where any probability above a threshold is true.  These places will be vetoed.
        glitch_mask_toohigh = ma.masked_greater(fixed_glitch_prob.data,probability_threshold[idx])
        
        # Now initialize an array to all ones (True) - we start with a mask that doesn't veto anything
        veto_mask = ones(shape(glitch_probability))

        # Find the indices where we want to set the mask to zero
        mdx1 = where(glitch_mask_nonzero.mask)
        mdx2 = where(glitch_mask_toohigh.mask)
        
        # Set the mask to zero at those places
        veto_mask[mdx1] = 0
        veto_mask[mdx2] = 0
        
        # Run the smoothing -- this may not always work...
        step_sum = sum(veto_mask)
        sum_diff = step_sum
        i=0
        while(sum_diff > 2):
            mask_iterate(veto_mask)
            sum_diff = abs(step_sum - sum(veto_mask))
            step_sum = sum(veto_mask)

        # Get the indices where it is zero (should be a union of mdx1 and mdx2)
        tdx = where(veto_mask==0)

        # Figure out how much time / how many glitches are vetoed (yes, this step really works -- may not need the second sum step 
        # on the veto arrays, but it doesn't hurt)
        veto_time = time_map[tdx]
        veto_glitch = glitch_count[tdx]
        
        vetoed_time = sum(sum(veto_time))
        vetoed_glitches = sum(sum(veto_glitch))
        
        deadtime[idx] = vetoed_time / total_time
        efficiency[idx] = vetoed_glitches / total_glitches


    EFF_DED = efficiency / deadtime

    # Find the point where efficiency/deadtime is maximized
    odx = argmax(EFF_DED)

    f = h5py.File(datafile,'a')
    
    if 'peak_veto_mask' in f:
        del f['peak_veto_mask']

    if 'rate_veto_mask' in f:
        del f['rate_veto_mask']

    if 'user_veto_mask' in f:
        del f['user_veto_mask']

    if 'peak_veto_mask_simple' in f:
        del f['peak_veto_mask_simple']

    if 'rate_veto_mask_simple' in f:
        del f['rate_veto_mask_simple']

    if 'user_veto_mask_simple' in f:
        del f['user_veto_mask_simple']

    if 'user_threshold' in f:
        del f['user_threshold']

    if 'ersatz_veto_mask' in f:
        del f['ersatz_veto_mask']

    # Build and save four masks; one at the EFF/DED peak, one at 10x the mean glitch rate, one at the user-specified threshold,
    # and one phony mask for testing

    fixed_glitch_probability = ma.fix_invalid(glitch_probability,fill_value = 0)
    glitch_mask_nonzero = ma.masked_less(fixed_glitch_probability.data,0.000001)
    glitch_mask_toohigh = ma.masked_greater(fixed_glitch_probability.data,probability_threshold[odx])
    peak_veto_mask = ones(shape(glitch_probability))
    mdx1 = where(glitch_mask_nonzero.mask)
    mdx2 = where(glitch_mask_toohigh.mask)
    peak_veto_mask[mdx1] = 0
    peak_veto_mask[mdx2] = 0

    peak_veto_mask_simple = peak_veto_mask.copy()

    step_sum = sum(peak_veto_mask_simple)
    sum_diff = step_sum
    i=0
    while(sum_diff > 2):
        mask_iterate(peak_veto_mask_simple)
        sum_diff = abs(step_sum - sum(peak_veto_mask_simple))
        step_sum = sum(peak_veto_mask_simple)


    peak_mask_dset = f.create_dataset('peak_veto_mask_simple', data=peak_veto_mask_simple)

    fixed_glitch_probability = ma.fix_invalid(glitch_probability,fill_value = 0)
    glitch_mask_nonzero = ma.masked_less(fixed_glitch_probability.data,0.000001)
    glitch_mask_toohigh = ma.masked_greater(fixed_glitch_probability.data,ten_probability)
    rate_veto_mask = ones(shape(glitch_probability))
    mdx1 = where(glitch_mask_nonzero.mask)
    mdx2 = where(glitch_mask_toohigh.mask)
    rate_veto_mask[mdx1] = 0
    rate_veto_mask[mdx2] = 0

    rate_veto_mask_simple = rate_veto_mask.copy()

    step_sum = sum(rate_veto_mask_simple)
    sum_diff = step_sum
    i=0
    while(sum_diff > 2):
        mask_iterate(rate_veto_mask_simple)
        sum_diff = abs(step_sum - sum(rate_veto_mask_simple))
        step_sum = sum(rate_veto_mask_simple)

    rate_mask_dset = f.create_dataset('rate_veto_mask_simple', data=rate_veto_mask_simple)

    fixed_glitch_probability = ma.fix_invalid(glitch_probability,fill_value = 0)
    glitch_mask_nonzero = ma.masked_less(fixed_glitch_probability.data,0.000001)
    glitch_mask_toohigh = ma.masked_greater(fixed_glitch_probability.data,user_threshold)
    user_veto_mask = ones(shape(glitch_probability))
    mdx1 = where(glitch_mask_nonzero.mask)
    mdx2 = where(glitch_mask_toohigh.mask)
    user_veto_mask[mdx1] = 0
    user_veto_mask[mdx2] = 0

    user_veto_mask_simple = user_veto_mask.copy()

    step_sum = sum(user_veto_mask_simple)
    sum_diff = step_sum
    i=0
    while(sum_diff > 2):
        mask_iterate(user_veto_mask_simple)
        sum_diff = abs(step_sum - sum(user_veto_mask_simple))
        step_sum = sum(user_veto_mask_simple)

    user_mask_dset = f.create_dataset('user_veto_mask_simple', data=user_veto_mask_simple)
    user_thresh_dset = f.create_dataset('user_threshold', data=user_threshold)

    ersatz_veto_mask = ones(shape(glitch_probability))
    rows, columns = shape(ersatz_veto_mask)
    for i in range(rows):
        for j in range(columns):
            if i < 50 or i > 200 or j < 50 or j > 200:
                ersatz_veto_mask[i,j] = 0

    ersatz_mask_dset = f.create_dataset('ersatz_veto_mask', data=ersatz_veto_mask)

    f.close()



    # Generate some plots

    fontP = FontProperties()
    fontP.set_size('small')

    matplotlib.rcParams.update({'savefig.dpi':250})

    fignum=0

    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(peak_veto_mask,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Veto Mask -- Glitch Probability > ' + '%.4f' % probability_threshold[odx])
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'peak_veto_mask.png')
    pylab.close()

    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(peak_veto_mask_simple,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Veto Mask -- Glitch Probability > ' + '%.4f' % probability_threshold[odx])
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'peak_veto_mask_simple.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(rate_veto_mask,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Veto Mask -- Glitch Probability > ' + '%.4f' % ten_probability)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'rate_veto_mask.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(rate_veto_mask_simple,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Veto Mask -- Glitch Probability > ' + '%.4f' % ten_probability)
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'rate_veto_mask_simple.png')
    pylab.close()    

    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(user_veto_mask,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Veto Mask -- Glitch Probability > ' + str(user_threshold))
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'user_veto_mask.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(user_veto_mask_simple,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Veto Mask -- Glitch Probability > ' + str(user_threshold))
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'user_veto_mask_simple.png')
    pylab.close()
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.imshow(ersatz_veto_mask,cmap=pylab.cm.gray,origin='lower',interpolation='nearest')
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Ersatz Veto Mask')
    pylab.xticks(visible=False)
    pylab.yticks(visible=False)
    pylab.savefig(web_directory + 'ersatz_veto_mask.png')
    pylab.close()    
    
    fignum=fignum+1
    pylab.figure(fignum)
    pylab.plot((segmentMedians[:,0] - xmean)/xspan,(segmentMedians[:,1] - ymean)/yspan,'r.')
    pylab.axis([-0.52, 0.52, -0.52, 0.52])
    pylab.grid(True)
    pylab.xlabel(channelX,fontsize=14)
    pylab.ylabel(channelY,fontsize=14)
    pylab.title('Median Positions of Science Segments')
    pylab.savefig(web_directory + 'segment_medians.png')
    pylab.close()    
    
    
    fignum=fignum+1
    pylab.figure(fignum)
    
    pylab.subplot(2,1,1)
    
    pylab.loglog(probability_threshold,efficiency,'b-')
    pylab.loglog(probability_threshold,deadtime,'r-')
    pylab.ylabel('Fraction Vetoed',fontsize=12)
    pylab.title('Veto effectiveness')
    pylab.grid(True)
    ylim = pylab.ylim()
    pylab.plot([user_threshold, user_threshold],ylim,'c-.',linewidth=1.1)
    pylab.plot([ten_probability, ten_probability],ylim,'m-.',linewidth=1.1)
    pylab.xticks(visible=False)
    pylab.legend(('Efficiency','Deadtime (preliminary)','User-Defined Threshold','10x Mean Glitch Rate'),loc=1,prop=fontP,bbox_to_anchor=(1.1, 1.2),fancybox=True)
    
    pylab.subplot(2,1,2)
    
    pylab.loglog(probability_threshold,EFF_DED,'k-')
    pylab.ylabel('Eff / dead',fontsize=12)
    pylab.grid(True)
    ylim = pylab.ylim()
    pylab.plot([user_threshold, user_threshold],ylim,'c-.',linewidth=1.1)
    pylab.plot([ten_probability, ten_probability],ylim,'m-.',linewidth=1.1)
    pylab.xlabel('Pixel Glitch Probability Threshold',fontsize=12)
    pylab.savefig(web_directory + 'eff_dead.png')
    pylab.close()
    

    return None
