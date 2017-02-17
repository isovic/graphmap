#! /usr/bin/python

import os;
import sys;
import math;

import numpy as np;
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from scipy.optimize import curve_fit

USE_MATPLOTLIB = True;
try:
	# import matplotlib;
	# matplotlib.use('Agg')
	import matplotlib.pyplot as plt;
	from matplotlib.font_manager import FontProperties;
	import seaborn as sns;
except Exception, e:
	USE_MATPLOTLIB = False;
	print e;

HIGH_DPI_PLOT = False;
# HIGH_DPI_PLOT = True;



def LineFunction(x, b):
    return (1*x + b);

def PlotMedianLine(ax, min_x, max_x, l_median, color='r'):
	x0 = 0;		y0 = 1*x0 + l_median;
	x1 = max_x;	y1 = 1*x1 + l_median;
	ax.plot([x0, x1], [y0, y1], color, lw=1);

def PlotLines(ax, min_x, max_x, l_median, threshold, color='purple'):
	threshold_l = threshold * 2.0 / (math.sqrt(2.0));
	l_min = l_median - threshold_l;
	l_max = l_median + threshold_l;
	
	x0 = 0;		y0 = 1*x0 + l_min;
	x1 = max_x;	y1 = 1*x1 + l_min;
	ax.plot([x0, x1], [y0, y1], color, lw=1);
	
	x0 = 0;		y0 = 1*x0 + l_max;
	x1 = max_x;	y1 = 1*x1 + l_max;
	ax.plot([x0, x1], [y0, y1], color, lw=1);
	


def load_csv(csv_path):
	fp = open(csv_path, 'r');
	lines = fp.readlines();
	fp.close()

	x = [];
	y = [];
	c = [];
	i = 0;
	for line in lines:
		split_line = line.split('\t');
		if (i == 0):
			if (len(split_line) == 6):
				query_header = split_line[0];
				query_id = int(split_line[1]);
				query_length = int(split_line[2]);
				l1_used = True if (int(split_line[3]) == 1) else 0;
				l_median = float(split_line[4]);
				l_maximum_allowed = float(split_line[5]);
			elif (len(split_line) == 7):
				query_header = split_line[0];
				query_id = int(split_line[1]);
				query_length = int(split_line[2]);
				ref_header = split_line[3];
				# ref_id = int(split_line[4]);
				# ref_length = int(split_line[5]);

				l1_used = True if (int(split_line[4]) == 1) else 0;
				l_median = float(split_line[5]);
				l_maximum_allowed = float(split_line[6]);
			else:
				sys.stderr.write('ERROR: Heading line does not contain a valid number of parameters!\n');
				exit(1);

			i += 1;
			continue;
		x.append(float(split_line[0].strip()));
		y.append(float(split_line[1].strip()));
		if (len(split_line) > 2):
			c.append(int(split_line[2].strip()));
		else:
			c.append(0);
		i += 1;
	return [x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed];

def plot_data(fig, ax, subplot_coords, x, y, c, query_length, ymin, ymax, l_median, threshold_L1_under_max, plot_mode, plot_title, out_png_path=''):
	if USE_MATPLOTLIB == True:
		# plt.figure();
		# plt.clf();
		# plt.subplot(subplot_coords);

		ax.grid();

		
		# if (subplot_coords == 223 or subplot_coords == 224):
		# 	plt.xlabel('Read coordinates');
		# if (subplot_coords == 221 or subplot_coords == 223):
		# 	plt.ylabel('Region coordinates');

		# plt.text(0.5, 1.08, 'Walks from the GraphMap output - %s' % plot_title,
		ax.text(0.5, 1.02, plot_title,
			horizontalalignment='center',
			fontsize=12,
			transform = ax.transAxes)

		# if (subplot_coords == 221 or subplot_coords == 223):
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.xlim(0, query_length);
		plt.ylim(ymin, ymax);
		plt.xticks(np.arange(0, query_length, query_length/5))

		all_colors = 'bgrcmyk';

		colors1 = [all_colors[val%len(all_colors)] for val in c[0::2]];
		colors2 = [all_colors[val%len(all_colors)] for val in c[1::2]];

		i = 0;
		while (i < len(x)):
			ax.plot(x[i:(i+2)], y[i:(i+2)], color=all_colors[(c[i]) % len(all_colors)]);
			i += 2;

		# ax.plot(x, y, 'o');
		# ax.scatter(x[0::2], y[0::2], s=10, facecolor='b', lw = 0.2)
		# ax.scatter(x[1::2], y[1::2], s=10, facecolor='c', lw = 0.2)
		ax.scatter(x[0::2], y[0::2], s=10, edgecolor=colors1, facecolor=colors1, lw = 0.2)
		ax.scatter(x[1::2], y[1::2], s=10, edgecolor=colors2, facecolor=colors2, lw = 0.2)

		if (plot_mode == 1):
			try:
				PlotMedianLine(ax, min(x), max(x), l_median, 'r');
				PlotLines(ax, min(x), max(x), l_median, threshold_L1_under_max, 'purple');
				PlotMedianLine(ax, 0, query_length, l_median, 'r');
				PlotLines(ax, 0, query_length, l_median, threshold_L1_under_max, 'purple');
			except Exception, e:
				sys.stderr.write(str(e) + '\n');

		# plt.xlabel('Read coordinates');
		# plt.ylabel('Region coordinates');

		if (subplot_coords == 223 or subplot_coords == 224):
			ax.set_xlabel('Query coordinates');
		if (subplot_coords == 221 or subplot_coords == 223):
			ax.set_ylabel('Reference coordinates');
		# sns.despine(offset=10, trim=True);
		# plt.setp([a.get_xticklabels() for a in fig.axes[:]], visible=False)
		# if (subplot_coords == 221 or subplot_coords == 222):
		# 	plt.setp(ax.get_xticklabels(), visible=False)

		# if (out_png_path != ''):
		# 	if (HIGH_DPI_PLOT == False):
		# 		plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
		# 	else:
		# 		plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);
		
		# print '';



if __name__ == "__main__":
	if (len(sys.argv) < 3):
		print 'Plots intermediate results from the LCSk-L1 step of the GraphMap algorithm.'
		print '';
		print 'Usage:';
		print '\t%s <path_to_results> local_scores_id_1 [local_scores_id_2 local_scores_id_3 ...]' % sys.argv[0];
		print '';
		exit(1);
	
	results_path = sys.argv[1];

	i = 2;
	while (i < len(sys.argv)):
		local_scores_id = (sys.argv[i]);
		scores_path = '%s/scores-%s' % (results_path, local_scores_id);
		lcs_path = '%s/LCS-%s' % (results_path, local_scores_id);
		lcsl1_path = '%s/LCSL1-%s' % (results_path, local_scores_id);
		l1_path = '%s/double_LCS-%s' % (results_path, local_scores_id);
		# boundedl1_path = '%s/boundedl1-%d' % (results_path, local_scores_id);

		# data_path = lcs_path;
		# print data_path;
		# [x, y] = load_csv(data_path + '.csv');
		# [l_median, threshold_L1_under_max] = FindHoughLine(x, y, error_rate);


		# fig = plt.figure();

		sns.set_style("darkgrid");
		sns.set_style("white")
		# sns.set_style("ticks");
		[fig, ((ax1, ax2), (ax3, ax4))] = plt.subplots(2, 2, sharex=True, sharey=True)

		# plt.clf();

		# fig.subplots_adjust(hspace=-0.5);
		fig.subplots_adjust(wspace=0.1);
		# fig.subplots_adjust(hspace=.5);
		# fig.subplots_adjust(wspace=.5);

		# f, ax = plt.subplots(4, sharex=True, sharey=True)

		# plt.gca().spines['top'].set_visible(False)
		# plt.gca().spines['right'].set_visible(False)
		# plt.gca().get_xaxis().tick_bottom()
		# plt.gca().get_yaxis().tick_left()

		# plt.setp([axarr[0].get_xticklabels(), axarr[1].get_xticklabels(), axarr[1].get_yticklabels(), axarr[3].get_yticklabels()], visible=False)
		# plt.setp([axarr[1].get_yticklabels(), axarr[3].get_yticklabels()], visible=False)

		data_path = scores_path;
		sys.stderr.write('Reading: %s\n' % data_path);
		[x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed] = load_csv(data_path + '.csv');
		ymin = min(y);
		ymax = max(y);
		# FindHoughLine(x, y, error_rate);
		plot_data(fig, ax1, 221, x, y, c, query_length, ymin, ymax, l_median, l_maximum_allowed, 0, 'Anchors', data_path + '.png');

		data_path = lcs_path;
		sys.stderr.write('Reading: %s\n' % data_path);
		[x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed] = load_csv(data_path + '.csv');
		# FindHoughLine(x, y, error_rate);
		plot_data(fig, ax2, 222, x, y, c, query_length, ymin, ymax, l_median, l_maximum_allowed, l1_used, 'LCSk', data_path + '.png');

		data_path = lcsl1_path;
		sys.stderr.write('Reading: %s\n' % data_path);
		[x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed] = load_csv(data_path + '.csv');
		# FindHoughLine(x, y, error_rate);
		plot_data(fig, ax3, 223, x, y, c, query_length, ymin, ymax, l_median, l_maximum_allowed, l1_used, 'LCSk L1 filtered', data_path + '.png');

		data_path = l1_path;
		sys.stderr.write('Reading: %s\n' % data_path);
		[x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed] = load_csv(data_path + '.csv');
		# FindHoughLine(x, y, error_rate);
		plot_data(fig, ax4, 224, x, y, c, query_length, ymin, ymax, l_median, l_maximum_allowed, l1_used, 'Second LCSk after L1', data_path + '.png');

		# data_path = boundedl1_path;
		# print data_path;
		# [x, y] = load_csv(data_path + '.csv');
		# # FindHoughLine(x, y, error_rate);
		# plot_data(224, x, y, l_median, threshold_L1_under_max, 1 if (error_rate > 0.0) else 0, 'boundedl1-%d' % local_scores_id, data_path + '.png');

		out_png_path = '%s/all-%s-qid_%d.png' % (results_path, local_scores_id, query_id);
		if (out_png_path != ''):
			if (HIGH_DPI_PLOT == False):
				sys.stderr.write('Writing image to file: %s\n\n' % out_png_path);
				plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
			else:
				sys.stderr.write('Writing image to file: %s\n\n' % out_png_path);
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);

# scripts/test-scatterplot3.py temp/local_scores/LCS-315.csv temp/local_scores/LCS-314.csv
# scripts/test-scatterplot4.py temp/local_scores/scores-104.csv 0 temp/local_scores/LCS-104.csv 1 temp/local_scores/L1-104.csv 1

			#for par in fitpars:
				#x1 = (1.0 - par) / 1;
				#plt.plot ([0.0, x1], [par, 1.0], 'g');

			#plt.xlim([0.0, 1.0]);
			#plt.ylim([0.0, 1.0]);
			
		i += 1;

	plt.show();
