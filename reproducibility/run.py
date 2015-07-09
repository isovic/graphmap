#! /usr/bin/python

import os;
import sys;
import subprocess;

def execute_command(command):
	sys.stderr.write('Executing command: %s\n' % (command));
	subprocess.call(command, shell=True);

def run_simulations():
	sys.stderr.write('Starting the alignment process on simulated data.\n');
	sys.stderr.write('Note that this might take a very long time.\n');
	sys.stderr.write('E.g. BLAST took 110670 CPU secs in our tests on hg19_chr3 Oxford Nanopore 2d simulated dataset.\n');
	execute_command('aligneval/run-alignment.py');
	sys.stderr.write('\n');

	sys.stderr.write('Alignment script returned.\n');
	sys.stderr.write('\n');

	sys.stderr.write('Running the evaluation script.\n');
	execute_command('aligneval/run-evaluation.py');
	sys.stderr.write('\n');

	sys.stderr.write('Copying the results to reproducibility/results-simulated folder.\n');
	execute_command('cp aligneval/results/*.csv results-simulated');

	sys.stderr.write('Done!\n');
	sys.stderr.write('\n');

def main():
	if (os.path.exists('samscripts') == False or os.path.exists('aligneval') == False):
		sys.stderr.write('Please run setup.py first, to install all dependencies. Exiting.\n');
		exit(1);

	if (len(sys.argv) < 2):
		sys.stderr.write('Run the alignment and evaluation processes from the GraphMap preprint paper.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\tsim - Runs alignment on all simulation datasets. This might take quite a while to execute.\n');
		exit(0);

	if (sys.argv[1] == 'sim'):
		if (len(sys.argv) != 2):
			sys.stderr.write('Runs alignment on all simulation datasets. This might take quite a while to execute.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		run_simulations();
		exit(0);
		
	else:
		sys.stderr.write('ERROR: Unknown subcommand!\n');
		exit(0);

if __name__ == "__main__":
	main();
