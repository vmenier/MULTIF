#include "meshutils.h"

/*
Victorien Menier March 2016
*/

Cad* AllocGeo (int *CadSiz)
{
	int i;
	
	int NbrVer    = CadSiz[CadVertices];
	int NbrVerMax = CadSiz[CadVerMax];
	
	int NbrLin    = CadSiz[CadLines];
	int NbrLinMax = CadSiz[CadLineMax];
	
	int NbrSrf    = CadSiz[CadSurfMax];
	
	int NbrLoo       = CadSiz[CadLoops];
	int NbrLooMax    = CadSiz[CadLoopMax];
	
	int NbrCtrVer = CadSiz[CadCtrVer];
	
	Cad *cad = NULL;
	
	cad = (Cad*) malloc(sizeof(struct S_Cad));
	
	
	cad->NbrVer = 0;
	cad->NbrLin = 0;
	cad->NbrSrf = 0;
	cad->NbrLoo = 0;
	
	cad->VerMax = 0;
	cad->LinMax = 0;
	cad->SrfMax = 0;
	cad->LooMax = 0;
	
	cad->Ver = NULL;
	cad->Lin = NULL;
	cad->Srf = NULL;
	cad->Loo = NULL;
	

	if ( NbrVer > 0 ) {
		cad->Ver = (CadVertex*) malloc(sizeof(CadVertex)*(NbrVer+1));
		cad->VerMax = NbrVerMax;
		
		cad->VidNew = (int*) malloc(sizeof(int)*(NbrVerMax+1));
		memset(cad->VidNew, 0, sizeof(int)*(NbrVerMax+1));
		
		cad->VidOld = (int*) malloc(sizeof(int)*(NbrVer+1));
		memset(cad->VidOld, 0, sizeof(int)*(NbrVer+1));
				
		for (i=0; i<=NbrVer; i++)
			cad->Ver[i].active = 0;
			
		cad->NbrVer = NbrVer;
	}
	
	if ( NbrLin > 0 ) {
		cad->Lin = (CadLine*) malloc(sizeof(CadLine)*(NbrLin+1));
		cad->LinMax = NbrLinMax;
		
		cad->LidNew = (int*) malloc(sizeof(int)*(NbrLinMax+1));
		memset(cad->LidNew, 0, sizeof(int)*(NbrLinMax+1));
		
		cad->LidOld = (int*) malloc(sizeof(int)*(NbrLin+1));
		memset(cad->LidOld, 0, sizeof(int)*(NbrLin+1));
		
		cad->NbrLin = NbrLin;
	}
	
	if ( NbrCtrVer > 0 ) {
		cad->CtrVer = (int*) malloc(sizeof(int)*(NbrCtrVer+1));
		cad->NbrCtrVer = NbrCtrVer;
		memset(cad->CtrVer, 0, sizeof(int)*(NbrCtrVer+1));
	}
	
	if ( NbrSrf > 0 ) {		
		cad->Srf = (CadSurface*) malloc(sizeof(CadSurface)*(NbrSrf+1));
		cad->SrfMax = NbrSrf;
	}
	
	if ( NbrLoo > 0 ) {
		cad->Loo = (CadLoop*) malloc(sizeof(CadLoop)*(NbrLoo+1));
		cad->LooMax = NbrLoo;
		
		cad->NbrLoo = NbrLoo;
	}
	
	return cad; 
}


Cad * FreeGeo (Cad *cad)
{
	if ( !cad )
		return NULL;
	
	
	//--- Vertices
	
	if ( cad->Ver )
		free(cad->Ver);
	
	if ( cad->VidNew )
		free(cad->VidNew);
		
	if ( cad->VidOld )
		free(cad->VidOld);
	
	//--- Lines
	
	if ( cad->Lin ) 	
		free(cad->Lin);
	
	if ( cad->CtrVer )
		free(cad->CtrVer);
	
	if ( cad->LidNew ) 
		free(cad->LidNew);
		
	if ( cad->LidOld ) 
		free(cad->LidOld);
		
	//--- Surf

	if ( cad->Srf ) 	
		free(cad->Srf);
		
	if ( cad->Loo ) 	
		free(cad->Loo);

	free(cad);
	
	return NULL;
}

int GetGeoSize (char *GeoNam, int *CadSiz)
{
  int  i, cpt, typ;
  double res;
  char str[1024], str2[1024], *tok=NULL, *tok1=NULL, buf[1024];

  int cptVer=0, cptBezVer=0, vidMax=0, vidBezMax=0;
	int NbrCtrVer = 0, lidMax = 0, NbrLin = 0, LoopIdMax = 0, NbrLoop=0;

	for (i=0; i<CadKwdSize; i++) 	
		CadSiz[i] = 0;

  FILE *GeoHdl = NULL;

  GeoHdl = fopen(GeoNam, "r");

  if ( !GeoHdl ) {
    printf("  ## ERROR : Couldn't open %s\n", GeoNam);
		return 0;
  }

  printf(" -- Info : %%%% %s opened.\n", GeoNam);

  //-----------------
  //--- Count Points
  //-----------------

  do
	{
		res = fscanf(GeoHdl, "%s", str);
		
    strncpy(str2,str,5);
    if ( strncmp (str2,"Point",5) == 0 ) {

      tok = strtok(str,"()");
      tok = strtok(NULL,"()");
      if ( tok == NULL ) {
        printf("  ## ERROR : Wrong syntax GetGeoSize: %s\n", str);
        exit(1);
      }

      strncpy(str2,tok,1);
      if ( strncmp (str2,"p",1) == 0 ) {
        tok1 = strtok(tok,"p");
        if ( tok1 == NULL ) {
          printf("  ## ERROR : Wrong syntax GetGeoSize : %s\n", str);
          exit(1);
        }
        vidBezMax = max(vidBezMax,atoi(tok1));
        cptBezVer++;
      }
      else {
        vidMax = max(vidMax, atoi(tok));
        cptVer++;
      }   
    }

	}while( (res != EOF) );

	//-----------------
  //--- Count Splines
  //-----------------
  rewind(GeoHdl);
	
	
	
	do
	{
		res = fscanf(GeoHdl, "%s", str);
		
    strncpy(str2,str,13);
		

		typ = 0;
		
		if ( strncmp (str2,"BSpline",6) == 0 )
			typ = CADBSPLINE;
		else if ( strncmp (str2,"Line",4) == 0 )
			typ = CADLINE;
		else if ( strncmp (str2,"Loop",4) == 0 )
			typ = CADLOOP;
		
		if ( typ > 0 ) {
      tok = strtok(str,"()");
      tok = strtok(NULL,"()");
      if ( tok == NULL ) {
        //printf("  ## ERROR GetGeoSize : Wrong syntax here : %s\n", str);
        continue;
        exit(1);
      }
			
			if ( typ == CADBSPLINE || typ == CADLINE ) {
			
				lidMax = max(lidMax, atoi(tok));
				NbrLin++;
				
				printf("lin %d\n", NbrLin);
				
      	strcpy(str2, "");
      	for (i=0; i<1000; i++) {
      	  fread(buf, 1, 1,GeoHdl);
      	  if ( strncmp (buf,";",1) == 0 ) {
      	    break;
      	  }
					
      	  strncat (str2, buf, 1);
printf("str2 = %s\n buf = %s\n", str2, buf);
      	}
      	
      	strcpy(str, "");
      	cpt = 0;
      	tok = strtok(str2,"={};, ");
				
				NbrCtrVer++;
      	
      	tok = strtok(NULL,"={};, ");
      	while ( tok != NULL ) {
      	  strcpy(str, tok);
					NbrCtrVer++;
      	  tok = strtok(NULL,"={};, ");
      	}
			}
			else if ( typ == CADLOOP ) {
				
				LoopIdMax = max(LoopIdMax, atoi(tok));
				NbrLoop++;
				
			}
			
    }
    
    strcpy(str, "");
	}while( (res != EOF) );

	
	CadSiz[CadVertices] = cptVer;
	CadSiz[CadLines]    = NbrLin;
	CadSiz[CadVertices] = cptVer;
	CadSiz[CadCtrVer]   = NbrCtrVer;
	
	CadSiz[CadLoops]   = NbrLoop;
	CadSiz[CadLoopMax] = LoopIdMax;
	
	printf("CadSiz[CadLoops] =%d, CadSiz[CadLoopMax]=%d\n", CadSiz[CadLoops], CadSiz[CadLoopMax]);
	
	CadSiz[CadVerMax]  = vidMax;
	CadSiz[CadLineMax] = lidMax;
	
  fclose(GeoHdl);

	return 1;
}



int ReadGeo (char *GeoNam, Cad *cad)
{
  
  int  i, cpt, typ;
  double res;
  char str[1024], str2[1024], *tok=NULL, buf[1024];
  int NbrVer = 0, NbrCtrVer=0, NbrLin=0, NbrLoo=0;

	int Vid, Lid;
	CadVertex *ver = NULL;
	
	CadLine *lin = NULL;
	CadLoop *loo = NULL;

  FILE *GeoHdl = fopen(GeoNam, "r");
  
  if ( !GeoHdl ) {
    printf("  ## ERROR : Couldn't open %s\n", GeoNam);
    exit(1);
  }
  
  printf(" -- Info : %%%% %s opened.\n", GeoNam);
  
	
  //-------------------
  //--- Read Points
  //-------------------
  
  rewind(GeoHdl);
  
  NbrVer = 0;
  do
	{
		res = fscanf(GeoHdl, "%s", str);
		
    strncpy(str2,str,5);
    if ( strncmp (str2,"Point",5) == 0 ) {
      
      tok = strtok(str,"()");
      tok = strtok(NULL,"()");
      if ( tok == NULL ) {
        printf("  ## ERROR ReadGeo: Wrong syntax : %s\n", str);
        exit(1);
      }
      
      strncpy(str2,tok,1);

			Vid = atoi(tok);
			NbrVer++;
			
			if ( Vid > cad->VerMax || NbrVer > cad->VerMax ) {
        printf("  ## ERROR : inconsistent number of vertices read!\n");
				printf(" Vid %d , NbrVer %d, vermax %d\n", Vid, NbrVer, cad->VerMax);
        exit(1);
      }
			
			ver = &cad->Ver[NbrVer];
			ver->active = 1;
			
			cad->VidNew[Vid]    = NbrVer;
			cad->VidOld[NbrVer] = Vid;
						
      strcpy(str2, "");
      for (i=0; i<1000; i++) {
        fread(buf, 1, 1,GeoHdl);
				
        if ( strncmp (buf,";",1) == 0 ) {
          break;
        }
        
        strncat (str2, buf, 1);
      }
      
      cpt = 0;
      tok = strtok(str2,"{};, ");
      tok = strtok(NULL,"{};, ");
      while ( tok != NULL && cpt < 3 ) {
        ver->Crd[cpt] = atof(tok);
        tok = strtok(NULL,"{};, ");
        cpt++;
      }
  
    }
    		
	}while( (res != EOF) );
  
  //-----------------
  //--- Read Splines
  //-----------------
	
  rewind(GeoHdl);
  
	do
	{
		res = fscanf(GeoHdl, "%s", str);
		
    strncpy(str2,str,6);
    
		typ = 0;
		
		if ( strncmp (str2,"BSpline",6) == 0 )
			typ = CADBSPLINE;
		else if ( strncmp (str2,"Line",4) == 0 )
			typ = CADLINE;
		else if ( strncmp (str2,"Loop",4) == 0 )
			typ = CADLOOP;
		
		if ( typ > 0 ) {
			
      tok = strtok(str,"()");
      tok = strtok(NULL,"()");

      if ( tok == NULL ) {
        printf("  ## ERROR ReadGeo: Wrong syntax here : %s\n", str);
        continue;
        exit(1);
      }
			
			if ( typ == CADBSPLINE || typ == CADLINE ){
				
				Lid = atoi(tok);
				NbrLin++;
				
				if ( NbrLin > cad->LinMax ||  Lid > cad->LinMax ) {
					printf("  ## ERROR : Inconsistent number of lines.\n");
					printf(" NbrLin %d, Lid %d, LinMax = %d\n", NbrLin, Lid, cad->LinMax);
					exit(1);
				}
				
				lin = &cad->Lin[NbrLin];
				lin->Typ = typ;
				
				cad->LidOld[NbrLin] = Lid;
				cad->LidNew[Lid]    = NbrLin;
				
				printf("cad->LidNew[%d]    = %d\n", Lid,	cad->LidNew[Lid] );
				
      	strcpy(str2, "");
      	for (i=0; i<1000; i++) {
      	  fread(buf, 1, 1,GeoHdl);
      	  if ( strncmp (buf,";",1) == 0 ) {
      	    break;
      	  }
      	  strncat (str2, buf, 1);
      	}
      	
      	strcpy(str, "");
      	cpt = 0;
      	tok = strtok(str2,"={};, ");
				
				NbrCtrVer++;
				if ( NbrCtrVer > cad->NbrCtrVer ) {
					printf("  ## ERROR : Inconsistent number of control points.\n");
					printf("NbrCtrVer %d cad->NbrCtrVer %d\n", NbrCtrVer, cad->NbrCtrVer);
					exit(1);
				}
      	
				Vid = atoi(tok);
				cad->CtrVer[NbrCtrVer] = cad->VidNew[Vid];
				
				lin->HeadCtr = NbrCtrVer;
				
				cpt=1;
      	
      	tok = strtok(NULL,"={};, ");
      	
      	while ( tok != NULL ) {
					cpt++;
      	  strcpy(str, tok);
					Vid = atoi(tok);
					NbrCtrVer++;
					cad->CtrVer[NbrCtrVer] = cad->VidNew[Vid];
					
      	  tok = strtok(NULL,"={};, ");
      	}
      	
				lin->NbrCtr = cpt;
			}
      else if ( typ == CADLOOP ) {
				NbrLoo++;
				loo = &cad->Loo[NbrLoo];
				
					strcpy(str2, "");
	      	for (i=0; i<1000; i++) {
	      	  fread(buf, 1, 1,GeoHdl);
	      	  if ( strncmp (buf,";",1) == 0 ) {
	      	    break;
	      	  }
	      	  strncat (str2, buf, 1);
	      	}

	      	strcpy(str, "");
	      	cpt = 0;
	      	tok = strtok(str2,"={};, ");


					Vid = atoi(tok);
					
					loo->Lin[0] = Vid;
					
					cpt=1;

	      	tok = strtok(NULL,"={};, ");

	      	while ( tok != NULL ) {
	      	  strcpy(str, tok);
						Vid = atoi(tok);
						loo->Lin[cpt] = Vid;
						cpt++;
						if ( cpt > 10 ) {
							printf("  ## ERROR ReadGeo: too many lines in loop %d\n", NbrLoo);
							exit(1);
						}
	      	  tok = strtok(NULL,"={};, ");
	      	}
	
					loo->NbrLin = cpt;				
			}
    }
    
    strcpy(str, "");
	}while( (res != EOF) );
  
  
  fclose(GeoHdl);
	

  return 1;
  
}



