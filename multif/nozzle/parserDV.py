import sys


def ParseDesignVariables_Plain (filename):
	
	try:
		fil = open(filename, 'r');
	except:
		sys.stderr.write("  ## ERROR : Could not open %s\n\n" % filename);
		sys.exit(0);
	
	fil = open(filename, 'r');
	
	DV_List = [];
	
	id_line = 0;
	
	for line in fil:
		
		id_line = id_line+1;
		
		hdl = line.split();
				
		if len(hdl) != 1 :
			sys.stderr.write("  ## ERROR : Unexpected input design variables format (line %d of %s).\n\n" % (id_line, filename));
			sys.exit(0);
			
		DV_List.append(float(hdl[0]));
		
	return DV_List;
	


