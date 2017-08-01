import sys, os, copy
import numpy as np
import multiprocessing

import multif


def SampleGetOutput(nozzle):
	Output = [];
	
	prt_item = [];
	prt_comp = [];
	prt_val  = [];
	
	for i in range(0, len(nozzle.outputTags)):
		tag = nozzle.outputTags[i];			
		
		# 6 Get Hessian, gradient, and value Get Hessian and gradient
		# 5 Get Hessian and value
		# 4 Get Hessian
		# 3 Get gradient and value
		# 2 Get gradient
		# 1 Get value
		# 0 No data required, function is inactive
		code = nozzle.outputCode[i];
		
		if isinstance(nozzle.responses[tag],list):
			for i in range(len(nozzle.responses[tag])):
				if isinstance(nozzle.responses[tag][i],list): # i.e. nested list
					for j in range(len(nozzle.responses[tag][i])):
						prt_item.append('%s_%i_%i' % (tag,i,j));
						prt_comp.append(" ");
						prt_val.append(nozzle.responses[tag][i][j]);
				else:
					prt_item.append('%s %i' % (nozzle.responses[tag][i],i));
					prt_comp.append('%s' % nozzle.prefixLabels[i]);
					prt_val.append(nozzle.responses[tag][i]);
		elif isinstance(nozzle.responses[tag],np.ndarray):
			arrayShape = nozzle.responses[tag].shape;
			if len(arrayShape) == 1:
				nr = arrayShape[0];
				for i in range(nr):
					prt_item.append('%s loc %i' % (tag,i));
					prt_comp.append('');
					prt_val.append(nozzle.responses[tag][i]);
			elif len(arrayShape) == 2:
				nr, nc = arrayShape;	  
				for i in range(nr):
					for j in range(nc):
						prt_item.append('%s loc %i %i' % (tag,i,j));
						prt_comp.append('');
						prt_val.append(nozzle.responses[tag][i,j]);								  
		else:
			prt_item.append('%s' % tag);
			prt_comp.append('');
			prt_val.append(nozzle.responses[tag]);				
		
	return prt_val;
	

def RunOneSample(run_id, nozzle, output='verbose'):
	
	dirNam = 'RUN_%d' % run_id;
	
	if not os.path.exists(dirNam):
		os.makedirs(dirNam);
	os.chdir(dirNam);
	
	if output == 'verbose':
		sys.stdout.write('Running sample in %s\n' % dirNam);	
	   
	sav_stdout = sys.stdout;
	
	sys.stdout = open('log_%d'%run_id, 'w');
	
	if nozzle.method == 'NONIDEALNOZZLE' :
		multif.LOWF.Run(nozzle, output);
	elif nozzle.dim == '2D': # nozzle method should be Euler or RANS
		multif.MEDIUMF.Run(nozzle, output);
	elif nozzle.dim == '3D': # nozzle.method should be RANS
		multif.HIGHF.Run(nozzle, output);
	else:
		sys.stderr.write(" ## ERROR runSample: Wrong fidelity level defined.\n");
		sys.exit(1);
		
	sys.stdout = sav_stdout;
	
	# Exit directory
	os.chdir('..');
	
	# --- Outputs
	
	prt_val = SampleGetOutput(nozzle);
	
	return prt_val;
		

def RunSamples(nozzle, samples_filename, beg, end):
		
	samples_hdl = np.loadtxt(samples_filename);
	NbrSam = len(samples_hdl);
	
	# --- No bounds provided? Run all the samples
	
	if not beg:
		beg = 0;
	if not end:
		end = NbrSam-1;
	
	if beg < 0 or beg > NbrSam-1 \
		or end < 0 or end > NbrSam-1\
		or beg >= end :
		
		sys.stderr.write(" ## ERROR RunSamples: Wrong sample id bounds : %d %d\n" % (beg, end));
		sys.exit(1);
	
	nozzleEval = [];
	
	outputs = [];
	
	# --- Start python's multiprocessing pool
	
	if nozzle.partitions >= 1:
		pool = multiprocessing.Pool(processes=nozzle.partitions);
	
	# --- Load data
	
	for i in range(end):
		outputs.append([]);
		nozzleEval.append([]);
	
	for iSam in range(beg,end):
		
		dv_list = [];
		for j in range(len(samples_hdl[iSam])):
			dv_list.append(samples_hdl[iSam,j]);
		
		if len(nozzle.dvList) != len(dv_list):
			sys.stderr.write("  ## ERROR sample %d : Wrong number of DV. SKIP.\n" % (iSam));
			sys.stderr.write("nozzle has %d DV, current sample line has %d DV\n" % (len(nozzle.dvList),len(dv_list)));
			continue;
		
		nozzleEval[iSam] = copy.deepcopy(nozzle);
		
		# --- Update DV using current sample values
		nozzleEval[iSam].dvList = dv_list;		
		nozzleEval[iSam].UpdateDV(output='quiet');
		nozzleEval[iSam].SetupWall(output='quiet');
		

	
	# --- Run the analysis
		
	if nozzle.partitions == 1:
		for iSam in range(beg,end):
			outputs[iSam] = RunOneSample(iSam, nozzleEval[-1]);
	else:
		
		mEval = [];
		for iSam in range(end):
			mEval.append(-1);
				
		for iSam in range(beg,end):
			mEval[iSam] = pool.apply_async(RunOneSample,(iSam,nozzleEval[iSam]))
		pool.close();
		pool.join();
		
		for iSam in range(beg,end):
			outputs[iSam] = mEval[iSam].get();
	        
	# --- Write results in file
	
	resNam = "samples_results.dat";
	
	fil = open(resNam, "w");
	
	if fil:
		sys.stdout.write("--- Writing results file. %s opened.\n" % resNam);
	else:
		sys.stderr.write(" ## ERROR : Can't open result file %s.\n" % resNam);
	
	for i in range(beg,end):
		fil.write("%d " %  i)
		for j in range(len(samples_hdl[i])):
			fil.write("%lf " %  samples_hdl[i][j]);
		for j in range(len(outputs[i])):
			fil.write("%s " %  outputs[i][j]);
		fil.write("\n");
	fil.close();
	
		