# -*- coding: utf-8 -*-

import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
import multif
from .. import SU2


def CheckConvergence ( config ) :
	
	# filenames
	plot_format      = config['OUTPUT_FORMAT']
	plot_extension   = SU2.io.get_extension(plot_format)
	history_filename = config['CONV_FILENAME'] + plot_extension
	special_cases    = SU2.io.get_specialCases(config)
	
	history      = SU2.io.read_history( history_filename )
	
	plot = SU2.io.read_plot(history_filename);
	
	RhoRes = history['Res_Flow[0]'];
	NbrIte = len(RhoRes);
	
	IniRes = RhoRes[0];
	FinRes = RhoRes[NbrIte-1];
	
	ResDif = FinRes - IniRes;
		
	print "Initial res = %le, Final res = %lf, DIFF = %lf\n" % (IniRes, FinRes, ResDif);
	
	