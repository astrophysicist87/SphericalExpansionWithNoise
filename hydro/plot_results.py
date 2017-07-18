#!/usr/bin/env python

from numpy import *
from numpy.core import numeric as _nx
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.patches import Ellipse
import sys, os

mpl.rcParams['pdf.fonttype'] = 42
labelsize = 20
mpl.rcParams['xtick.labelsize'] = labelsize
mpl.rcParams['ytick.labelsize'] = labelsize

#################################################################
# Parameters to run script with
#################################################################
#file information
#filename = sys.argv[1]

#grid information
nDY = 51
nTrajectories = 3
nTauD = 3

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

panelLabels = ['(a)', '(b)']
panelCounter = 0

chosenParticle = 'pion'
chosenTauDs = ['0.0','0.5','1.0']

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def stack(arrays, axis=0):
	arrays = [np.asanyarray(arr) for arr in arrays]
	if not arrays:
		raise ValueError('need at least one array to stack')
	
	shapes = set(arr.shape for arr in arrays)
	if len(shapes) != 1:
		raise ValueError('all input arrays must have the same shape')
	
	result_ndim = arrays[0].ndim + 1
	if not -result_ndim <= axis < result_ndim:
		msg = 'axis {0} out of bounds [-{1}, {1})'.format(axis, result_ndim)
		raise IndexError(msg)
	if axis < 0:
		axis += result_ndim
	
	sl = (slice(None),) * axis + (_nx.newaxis,)
	expanded_arrays = [arr[sl] for arr in arrays]
	return _nx.concatenate(expanded_arrays, axis=axis)

#################################################################
# Plot spectra
def plotSpectra(pathname, scaleFactor = 1.0):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0

	filename1 = pathname + '/trajectory_1.dat'
	filename2 = pathname + '/trajectory_2.dat'
	filename3 = pathname + '/trajectory_3.dat'
	nCols = 2
	chosenCols = [0, 1]
	dims = [nDY, nCols]

	# read in file
	data1 = loadtxt(filename1, usecols=tuple(chosenCols)).reshape(dims)
	data2 = loadtxt(filename2, usecols=tuple(chosenCols)).reshape(dims)
	data3 = loadtxt(filename3, usecols=tuple(chosenCols)).reshape(dims)
	data = stack([data1, data2, data3])

	#print 'Loading data file...'
	ax.plot(data[0,:,0], scaleFactor * data[0,:,1], linestyle=':', color='red', linewidth=lw, label='Trajectory 1')
	ax.plot(data[1,:,0], scaleFactor * data[1,:,1], linestyle='--', color='green', linewidth=lw, label='Trajectory 2')
	ax.plot(data[2,:,0], scaleFactor * data[2,:,1], linestyle='-', color='blue', linewidth=lw, label='Trajectory 3')
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$C(\Delta y)$', fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	if scaleFactor > 1.1:
		text(0.15, 0.9, r'$\times 10^{-3}$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
		text(0.9, 0.15, r'(a)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	else:
		text(0.9, 0.15, r'(b)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	#plt.title(pathname)
	new_ylimits = asarray(ax.get_ylim())
	ax.set_ylim(bottom=0.99*new_ylimits[0], top=new_ylimits[1])
	
	#plt.show(block=False)
	#plt.show()
	outfilename = pathname + '/' + pathname.replace('/', '_') + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


#################################################################
# Plot HBT correlations
def plotHBT(pathname, scaleFactor = 1.0):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0

	filename1 = pathname + '/trajectory_1.dat'
	filename2 = pathname + '/trajectory_2.dat'
	filename3 = pathname + '/trajectory_3.dat'
	nCols = 2
	chosenCols = [0, 6]
	dims = [nDY, nCols]

	# read in file
	data1 = loadtxt(filename1, usecols=tuple(chosenCols)).reshape(dims)
	data2 = loadtxt(filename2, usecols=tuple(chosenCols)).reshape(dims)
	data3 = loadtxt(filename3, usecols=tuple(chosenCols)).reshape(dims)
	data = stack([data1, data2, data3])

	#print 'Loading data file...'
	ax.plot(data[0,:,0], scaleFactor * data[0,:,1], linestyle=':', color='red', linewidth=lw, label='Trajectory 1')
	ax.plot(data[1,:,0], scaleFactor * data[1,:,1], linestyle='--', color='green', linewidth=lw, label='Trajectory 2')
	ax.plot(data[2,:,0], scaleFactor * data[2,:,1], linestyle='-', color='blue', linewidth=lw, label='Trajectory 3')
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$C_{\mathrm{HBT}}(\Delta y)$', fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	if scaleFactor > 1.1:
		text(0.15, 0.9, r'$\times 10^{-3}$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
		text(0.875, 0.15, r'(a)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	else:
		text(0.875, 0.15, r'(b)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	#plt.title(pathname)
	new_ylimits = asarray(ax.get_ylim())
	ax.set_ylim(bottom=0.99*new_ylimits[0], top=new_ylimits[1])
	
	plt.show(block=False)
	#plt.show()
	outfilename = pathname + '/' + pathname.replace('/', '_') + '.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


def generate_all_plots():
	#plotSpectra('results_1sp/pion/w_delta_muB')
	plotHBT('results_HBT/pion/w_delta_muB')
	#plotSpectra('results_1sp/pion/wo_delta_muB', 1000.0)
	plotHBT('results_HBT/pion/wo_delta_muB', 1000.0)
	#plotSpectra('results_1sp/proton')
	plotHBT('results_HBT/proton')
	pause()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
