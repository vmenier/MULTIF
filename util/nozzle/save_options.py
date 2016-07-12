# -*- coding: utf-8 -*-


import os, time, sys, shutil, copy
from optparse import OptionParser
import textwrap
from .. import SU2
import util



def NozzleSetup( config_name, flevel ):
	
	if not os.path.isfile(config_name) :
		msg = "  ## ERROR : could not find configuration file %s\n\ns" % config_name;
		sys.stdout.write(msg);
		sys.exit(1);
	
	nozzle = util.options.nozzle.Nozzle()
	
	config = SU2.io.Config(config_name)
	
	fidelity_tags = config['FIDELITY_LEVELS_TAGS'].strip('()');
	fidelity_tags = fidelity_tags.split(",");
	
	NbrFidLev = len(fidelity_tags);
	
	if NbrFidLev < 1 :
		print "  ## ERROR : No fidelity level was defined.\n"
		sys.exit(1)
	
	sys.stdout.write('\n%d fidelity level(s) defined. Summary :\n' % (NbrFidLev));
	
	sys.stdout.write('-' * 90);
	sys.stdout.write('\n%s | %s | %s\n' % ("Level #".ljust(10), "Tag".ljust(10),"Description.".ljust(70)));
	sys.stdout.write('-' * 90);
	sys.stdout.write('\n'); 
	
	for i in range(NbrFidLev) :
		tag = fidelity_tags[i];
		kwd = "DEF_%s" % tag;
		
		if kwd not in config :
			print "\n  ## ERROR : The fidelity level tagged %s is not defined.\n" % kwd;
			sys.exit(1)
	
		cfgLvl = config[kwd].strip('()');
		cfgLvl = cfgLvl.split(",");
		
		method = cfgLvl[0];
	
		description = "";
	
		if method == 'NONIDEALNOZZLE' :
			tol = float(cfgLvl[1]);
			if tol < 1e-30 :
				print "  ## ERROR : Wrong tolerance for fidelity level %d (tagged %s)\n" % (i,tag);
			description = "Quasi 1D non ideal nozzle with tolerance set to %lf." % (tol);	
	
		elif method == 'RANS' or method == 'EULER':
			dim = cfgLvl[1];
			if dim != '2D' and dim != '3D' :
				print "  ## ERROR : Wrong dimension for fidelity level %d (tagged %s) : only 2D or 3D simulations" % (i,tag);
			meshsize = cfgLvl[2];	
			if meshsize != 'COARSE' and meshsize != 'MEDIUM' and meshsize != 'FINE' :
				print "  ## ERROR : Wrong mesh level for fidelity level %d (tagged %s) : must be set to either COARSE, MEDIUM or FINE" % (i,tag);
			description = "%s %s CFD simulation using the %s mesh level." % (dim, method, meshsize);
	
		else :
			print "  ## ERROR : Unknown governing method (%s) for fidelity level %s." % (method, tag);
			print "  Note: it must be either NONIDEALNOZZLE, EULER, or RANS\n";
	
	
		sys.stdout.write("   %s | %s | %s \n" % ( ("%d" % i).ljust(7), tag.ljust(10), textwrap.fill(description,70, subsequent_indent="".ljust(26))) );
	
	sys.stdout.write('-' * 90);
	sys.stdout.write('\n\n');
	
	if flevel >= NbrFidLev :
		sys.stderr.write("\n ## ERROR : Level %d not defined !! \n\n" % flevel);
		sys.exit(0);
	
	sys.stdout.write('  -- Info : Fidelity level to be run : %d\n\n' % flevel);
	
	nozzle.level = "low";
	
	return nozzle;
	