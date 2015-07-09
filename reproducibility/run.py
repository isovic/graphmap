#! /usr/bin/python

import os;
import sys;
import subprocess;

def execute_command(command):
	sys.stderr.write('Executing command: %s\n' % (command));
	subprocess.call(command, shell=True);

def main():
	if (os.path.exists('samscripts') == False or os.path.exists('aligneval') == False):
		sys.stderr.write('Please run setup.py first, to install all dependencies. Exiting.\n');
		exit(1);

	


if __name__ == "__main__":
	main();
