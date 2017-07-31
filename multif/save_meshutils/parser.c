#include "meshutils.h"

/*
Victorien Menier Feb 2016
*/

int ParseCommandLine (Options *mshopt,int argc, char *argv[]) 
{
  int  i, j;
  
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
	
			case 'b':
			
			if ( !strcmp(argv[i],"-box") ) {   
				
        for (j=0; j<6; j++) {
          i++;
					
          if ( i < argc ) {
            mshopt->Box[j] = atof(argv[i]);
          }
          else {
            printf("  ## ERROR parser : Wrong box definition : 6 double required.\n");
            printf("                    Ex : -box [xmin,ymin,zmin,xmax,ymax,zmax] \n");
            exit(1);
          }
        }
      }
			
			break;

      case 'h':  
        if ( !strcmp(argv[i],"-header") ) {
          if ( ++i < argc ) 
            strcpy(mshopt->HeaderNam, argv[i]);
	        else
	          i--;
        }  
        break;

	    case 'i':  
	      if ( !strcmp(argv[i],"-in") ) {
	        if ( ++i < argc ) 
	          strcpy(mshopt->InpNam, argv[i]);
		      else
		        i--;
	      }  
	      break;
		
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc ) 
            strcpy(mshopt->OutNam, argv[i]);
	        else
	          i--;
        }
				break;
			
			case 'O':
	      if ( !strcmp(argv[i],"-O") ) {
	        if ( ++i < argc ) 
	          mshopt->Mod = atoi(argv[i]);
		      else
		        i--;
	      }
				break;

			case 's':
	      if ( !strcmp(argv[i],"-sol") ) {
	        if ( ++i < argc ) 
	          strcpy(mshopt->SolNam, argv[i]);
		      else
		        i--;
	      }
				break;
      
      default:
    	  printf("  Unrecognized option %s\n",argv[i]);
      }
    }


    
    else {
      printf("  Argument %s ignored\n",argv[i]);
    }
    i++;
  }    

  return 1;
}


int GetInputFileType (char *FilNam) 
{
  char *ext=NULL;  
	
  if ( !FilNam ) {
    printf("\n ## ERROR GetInputFileType : No file name provided.\n\n");
   	exit(1);
  }
  
  ext = strrchr(FilNam, '.');
  
  if (!ext) {
    return 0;    
  } else {      
      ext = ext+1;
      if ( strcmp(ext,"su2") == 0  ) {
        return FILE_SU2;
      }
      else if ( strcmp(ext,"dat") == 0  ) {
        return FILE_DAT;
      }
      else if ( strcmp(ext,"mesh") == 0 || strcmp(ext,"meshb") == 0 ) {
        return FILE_GMF;
      }
      else if ( strcmp(ext,"sol") == 0 || strcmp(ext,"solb") == 0 ) {
        return FILE_GMFSOL;
      }
      else if ( strcmp(ext,"geo") == 0  ) {
        return FILE_GEO;
      }
      else {
				return 0;
      }
  }
}


void Usage()
{
	
	printf("\n----- meshutils: usage ----- \n\n");
	
	printf("-- Inline options:\n");
	printf("   -O      [n] : Module choice\n");
	printf("   -in     [c] : Input mesh file (*.mesh or *.su2) \n");
	printf("   -out    [c] : Output file name (default: [InpNam].o) \n");
	printf("   -sol    [c] : Input solution file (*.sol or *.dat)\n");
	printf("   -header [c] : Reference *.dat file name. \n");
	
	printf("\n-- Modules:\n");
	printf("   [1] : Output Mach (mach.meshb and mach.solb, GMF format) number from input SU2 files.\n");
	printf("   [2] : Convert SU2 mesh+solution files to the GMF format.\n");
	printf("   [3] : Convert GMF mesh+solution files to the SU2 format. Use -header option to copy header information.\n");
	printf("   [4] : Solution extraction.\n");
	printf("   [5] : CAD operations.\n");
	printf("   [6] : Clean mesh.\n");
	
	printf("\n-- Examples:\n");
	printf("    meshutils -O 1 -in restart_flow.su2 -sol restart_flow.dat \n");
	printf("    meshutils -O 2 -in restart_flow.su2 -sol restart_flow.dat -out restart_flow\n");
	printf("    meshutils -O 3 -in restart_flow.meshb -sol restart_flow.solb -header restart_flow.dat -out restart.new\n");
	
	printf("\n");
}



