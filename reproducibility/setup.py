#! /usr/bin/python

import os;
import sys;
import subprocess;

def execute_command(command):
	sys.stderr.write('Executing command: %s\n' % (command));
	subprocess.call(command, shell=True);

def main():
	if (not os.path.exists('samscripts')):
		execute_command('git clone https://github.com/isovic/samscripts.git');

	if (not os.path.exists('aligneval')):
		execute_command('git clone https://github.com/isovic/aligneval.git');
		execute_command('cd aligneval; ./setup.py all');

	folders_to_generate = ['data/reads', 'data/reference', 'results-simulated', 'results-real'];
	for folder_to_generate in folders_to_generate:
		if (not os.path.exists(folder_to_generate)):
			os.makedirs(folder_to_generate);

if __name__ == "__main__":
	main();
