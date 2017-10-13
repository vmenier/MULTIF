import time, os, shutil, subprocess, datetime, sys, math
import numpy as np
import multif
from multif import _meshutils_module

try :
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
except:
    sys.stderr.write("## WARNING: Unable to load mpl_toolkits capabilities. Visualization functions won't work properly.\n");    

try:    
    import matplotlib.pyplot as plt
    from   matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.tri as tri
    import matplotlib.cm as cm
    from matplotlib.ticker import NullFormatter  # useful for `logit` scale
    
except:
    sys.stderr.write("## WARNING: Unable to load matplotlib capabilities. Visualization functions won't work properly.\n");
    

def Viz_GetMeshAndSolution(MshNam, SolNam):
    
    sys.stdout.write("Opening %s and %s\n" % (MshNam, SolNam));
        
    pyVer = [];
    pyTri = [];
    pyTet = [];
    pyEdg = [];
    pySol = [];
    
    _meshutils_module.py_ReadMesh (MshNam,SolNam, pyVer, pyTri, pyTet, pyEdg, pySol);
    
    NbrVer = len(pyVer)/3;
    Ver = np.array(pyVer).reshape(NbrVer,3);
    
    NbrTri = len(pyTri)/4;
    Tri = np.array(pyTri).reshape(NbrTri,4);
    
    SolSiz = len(pySol)/NbrVer;
    Sol = np.array(pySol).reshape(NbrVer,SolSiz);
    
    return Ver, Tri, Sol



def Viz_MapColorValue(zval, colormap, vmin, vmax):
    # Map the normalized value zval to a corresponding color in the colormap
    if vmin>vmax:
        raise ValueError('incorrect relation between vmin and vmax')
    t=(zval-vmin)/float((vmax-vmin))
    R, G, B, alpha=colormap(t)
    return R, G, B, alpha
    #return 'rgb('+'{:d}'.format(int(R*255+0.5))+','+'{:d}'.format(int(G*255+0.5))+\
    #       ','+'{:d}'.format(int(B*255+0.5))+')' 


def Viz_PlotMeshAndSolution2D(ax, Ver, Tri, Sol):
    #---
    print "PLOT 2D NO AVAILABLE"

def Viz_PlotMeshAndSolution3D(ax, Ver, Tri, Sol):
    
    NbrTri = len(Tri);
    
    crdMin=[1e6,1e6,1e6];
    crdMax=[-1e6,-1e6,-1e6];    
    
    solLoc = Sol[:,4]/Sol[:,0]; # pressure
    
    solMin = 1e30;
    solMax = -1e30;
    for i in range(NbrTri):
        if int(Tri[i][3]) != 9 and int(Tri[i][3]) != 10:
            continue;
        for j in range(3):
            iVer = int(Tri[i][j])-1
            solMin = min(solMin,solLoc[iVer]);
            solMax = max(solMax,solLoc[iVer]);
    
    
    cpt = 0;
    for i in range(NbrTri):
        
        if int(Tri[i][3]) != 9 and int(Tri[i][3]) != 10:
            continue;
            
        cpt += 1;
        
        ids = [int(Tri[i][0])-1,int(Tri[i][1])-1,int(Tri[i][2])-1]
        
        x = [Ver[ids[0]][0],Ver[ids[1]][0],Ver[ids[2]][0]]
        y = [Ver[ids[0]][1],Ver[ids[1]][1],Ver[ids[2]][1]]
        z = [Ver[ids[0]][2],Ver[ids[1]][2],Ver[ids[2]][2]]
        
        for j in range(3):
            for k in range(3):
                crdMin[k] = min(crdMin[k],Ver[ids[j]][k]);
                crdMax[k] = max(crdMax[k],Ver[ids[j]][k]);
        
        sol = 0.333*(solLoc[ids[0]]+solLoc[ids[1]]+solLoc[ids[2]]);
        
        #colormap=cm.RdBu;
        #colormap=cm.rainbow;
        colormap=cm.Spectral;
        col = Viz_MapColorValue(sol, colormap, solMin, solMax);
        
        verts = [zip(x,y,z)]
        ax.add_collection3d(Poly3DCollection(verts, color=col))
        
        #if cpt > 10 :
        #    break;
    
    
    ax.set_xlim(crdMin[0], crdMax[0])
    ax.set_ylim(crdMin[1], crdMax[1])
    ax.set_zlim(crdMin[2], crdMax[2])
    
    #for ii in xrange(0,360,10):
    #        ax.view_init(elev=5., azim=ii)
    #        plt.savefig("rot%d.png" % ii)
    #
    
    #for ii in xrange(-200,200,20):
    #        ax.view_init(elev=ii, azim=200)
    #        ax.set_aspect('equal')
    #        plt.savefig("rot90_elev%d.png" % ii,bbox_inches='tight')
    
    #elev=5;
    #rot = 90;
    #ax.view_init(elev=elev, azim=rot)
    #ax.set_aspect('equal')
    #plt.savefig("pres_rot%d_elev%d.png" %(rot,elev),bbox_inches='tight')
    #
    #elev=5;
    #rot = 290;
    #ax.view_init(elev=2, azim=rot)
    ax.set_aspect('equal')
    #plt.savefig("pres_rot%d_elev%d.png" %(rot,elev),bbox_inches='tight')
    
    ax.view_init(elev=5., azim=90)

    return ax;





def Viz_NozzleVisu(nozzle, **kwargs):
    
    visu_prefix = "";
    visu_path   = ".";
    
    if 'visu_prefix' in kwargs:
        visu_prefix = kwargs['visu_prefix'];
        
    if 'visu_path' in kwargs:
        visu_path = kwargs['visu_path'];
    
    print "Viz_NozzleVisu"
    print "Dimension : %s" % nozzle.dim;
    print "visu_path %s visu_prefix %s\n" % (visu_path,visu_prefix)
    
    mshNam = "nozzle.su2"
    solNam = "nozzle.dat"
    
    
    #if nozzle.dim == '3D':
    #    Ver, Tri, Sol = Viz_GetMeshAndSolution(mshNam, solNam);
    #    Viz_PlotMeshAndSolution3D(Ver, Tri, Sol);
        
    pdfNam = os.path.join(visu_path,"%ssummary.pdf" % visu_prefix)
    
    with PdfPages(pdfNam) as pdf:
        
        np.random.seed(19680801)
        
        # make up some data in the interval ]0, 1[
        y = np.random.normal(loc=0.5, scale=0.4, size=1000)
        y = y[(y > 0) & (y < 1)]
        y.sort()
        x = np.arange(len(y))
        
        fig = plt.figure(figsize=(8.27,11.69));
        
        #--- Set title
        
        plt.suptitle('%s %s Run Summary' % (nozzle.dim, nozzle.method), fontsize=16)
        
        ax = fig.add_subplot(3, 2, 3)
        
        #--- Plot residual
        
        
        hdl = [];
        if( os.path.isfile('history.csv') ):
            hdl = np.loadtxt('history.csv',skiprows=1,delimiter=',');   
        elif( os.path.isfile('history.dat') ):
            hdl = np.loadtxt('history.dat',skiprows=3,delimiter=',');
        
        ax.set_title('Solver convergence')
        ax.set_xlabel('#ITE')
        ax.set_ylabel('Density residual')
        ax.grid(True)
        
        if len(hdl)>1:
            ax.plot(hdl[:,0], hdl[:,11]) 

        
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig('solver_convergence.pdf', bbox_inches=extent.expanded(1.5, 1.3))
        
        #--- Plot flow solution
        
        if nozzle.dim == '3D':  
            Ver, Tri, Sol = Viz_GetMeshAndSolution(mshNam, solNam);
            ax = fig.add_subplot(3, 1, 1, projection='3d', aspect='equal')
            Viz_PlotMeshAndSolution3D(ax, Ver, Tri, Sol);
            ax.grid(False)
            ax.set_title('Flow solution (Pressure)')
            extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            fig.savefig('flow_solution.pdf', bbox_inches='tight')

        #--- Plot text summary
        
        ax = fig.add_subplot(3, 2, 4)
        ax.set_title('QoIs')
        ax.axis('off')
        
        t="\n";
        qoiNam = "results.out";        
        
        try:
            hdl_qoi=open(qoiNam,"r")
            lines=hdl_qoi.readlines()
            for l in lines:
                print l.split()
                
                tag = l.split()[1];
                qoi = float(l.split()[0]);
                t = "%s%s = %lf\n" % (t,tag, qoi)
            hdl_qoi.close()
        except:
            
            sys.stderr.write(" ## ERROR: Unable to open %s\n" % qoiNam)
            raise;
        
        ax.text(0, 1,t, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, wrap=True)
        
        pdf.savefig()
        plt.close()

        d = pdf.infodict()
        d['Title'] = 'Multipage PDF Example'
        d['Author'] = u'Victorien Menier'
        d['Subject'] = 'MULTI-F Post-processing'
        d['Keywords'] = ''
        d['CreationDate'] = datetime.datetime.today()
        d['ModDate'] = datetime.datetime.today()

    
    
    