#include "meshutils.h"

/*
Victorien Menier March 2016
*/


int SolutionExtraction(Options *mshopt,Mesh *Msh)
{
	int iVer, i, idxVer;
	double x, y;
	
	double eps = 1e-15;
	
	FILE *FilHdl = fopen ("yextract.dat", "wb");
	
	if ( !FilHdl ) {
		printf("  ## ERROR SolutionExtraction: Could not open file.\n");
		return 0;
	}
	
	//--- Write header
	fprintf(FilHdl, "# x y ");
	for (i=0; i<Msh->SolSiz; i++) {
		fprintf(FilHdl, "%s ", Msh->SolTag[i]);
	}
	fprintf(FilHdl, "\n");
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		x = Msh->Ver[iVer][0];
		y = Msh->Ver[iVer][1];
		
		if ( fabs(y) > eps ) continue;
		
		fprintf(FilHdl, "%lf %lf ", x, y);
		
		idxVer = iVer*Msh->SolSiz;
		
		for (i=0; i<Msh->SolSiz; i++) {
			fprintf(FilHdl, "%lf ", Msh->Sol[idxVer+i]);
		}
		
		fprintf(FilHdl, "\n");
		
	}
	
	if ( FilHdl )
		fclose(FilHdl);
	
	return 1;
	
}



//--- Structure for the qsort function
typedef struct cvt_SortInterface
{
	double y;
	double x;
  int    IdxSol;
} SortInterface, *pSortInterface;


int CmpInterface (const void *a, const void *b)
{
  SortInterface *ia = (SortInterface*) a;
  SortInterface *ib = (SortInterface*) b;

	double eps = 1e-20;


  if ( fabs(ia->x-ib->x) < DBL_EPSILON ) {
    if ( fabs(ia->y-ib->y) < DBL_EPSILON ) {
			return 0;
    }
    else if ( ia->y > ib->y ) return 1;
    else return -1;
  }
  else if ( ia->x > ib->x ) return 1;
  else return -1;

	

	//if ( fabs(ia->y-ib->y) < eps )
	//	return 0;
	//else if ( ia->y > ib->y  )
	//	return 1;
	//else return -1;
	
	//if ( fabs(ia->y-ib->y) < eps )
	//	return 0;
	//else if ( ia->y > ib->y  )
	//	return 1;
	//else return -1;

}


double * ExtractAlongLine (Options *mshopt, Mesh *Msh, double *box, int *NbrRes, int *Siz)
{
	
	int iTri, e, i, v, d, idx, iVar;
	int t2e[3][2] = { {1,2}, {2,0}, {0,1} };
	
	int flag=0;
	
	double x, y, alp=0, crd[2], xi, xj;
	
	int    edgVid[2] = {0,0};
	double edgCrd[2][2];
	
	double xExt = box[0];
	double yMin = box[2];
	double yMax = box[3];
	
	
	int NbvMax = 3*Msh->NbrTri;
	double *Sol = NULL;
	
	double *result = NULL;
	
	SortInterface *SrtInt = (SortInterface*) malloc (sizeof(struct cvt_SortInterface)*(NbvMax));	
	Sol = (double*) malloc(sizeof(double)*NbvMax*Msh->SolSiz);
	
	int NbrSort = 0, NbrSortFin=0;
	
	//printf("  -- Info : Box = [%lf, %lf] [%lf, %lf]\n", box[0], box[1], box[2], box[3]);
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		
		for (e=0; e<3; e++) {
			
			for (v=0; v<2; v++){
				edgVid[v] = Msh->Tri[iTri][t2e[e][v]];
				
				for (d=0; d<2; d++) 
					edgCrd[v][d] = Msh->Ver[edgVid[v]][d];
			}
			
			//--- Check if nodes are not both out of bound
			
			flag = 0;
			
			for (v=0; v<2; v++) {
				y = edgCrd[v][1];
				
				if ( y <= yMax && y >= yMin ) {
					flag = 1;
					break;
				}
			}
			
			
			if ( flag == 0 )
				continue;
			
			
			//--- Check if nodes are on the line
			
			flag = 0;
			
			for (v=0; v<2; v++) {
				
				x = edgCrd[v][0];
				y = edgCrd[v][1];
				
				if ( fabs(x-xExt) > 1e-20 )
					continue;
					
				if ( y > yMax || y < yMin )
					continue;
				
				flag = 1;
				
				//--- This point is on the line, add it to the list
				idx = NbrSort * Msh->SolSiz;
				for (iVar=0; iVar<Msh->SolSiz; iVar++) { 
					Sol[idx+iVar] = Msh->Sol[edgVid[v]*Msh->SolSiz+iVar];
					//printf("SOL %d = %lf\n", iVar, Sol[idx+iVar]);
				}
				//exit(1);
				
				SrtInt[NbrSort].x      = x;
				SrtInt[NbrSort].y      = y;
				SrtInt[NbrSort].IdxSol = NbrSort;
				
				NbrSort++;
			}
			
			/* If at least one node is on the line, continue */
			if ( flag == 1 ) 
				continue;
				
			xi = edgCrd[0][0];
			xj = edgCrd[1][0];
			
			flag = (xi<xExt?-1:1)*(xj<xExt?-1:1);
			
			if ( flag == 1 ) continue;
			
			if ( xi < xj )
				alp = (xExt-xi)/(xj-xi);
			else 
				alp = 1.-(xExt-xj)/(xi-xj);
			
			//for (i=0; i<2; i++)
			//	crd[i] = alp*node[jPoint]->GetCoord(i)+(1.-alp)*node[iPoint]->GetCoord(i);
			for (d=0; d<2; d++) 
				crd[d] = alp*edgCrd[1][d]+(1.-alp)*edgCrd[0][d];
			
			SrtInt[NbrSort].x      = crd[0];
			SrtInt[NbrSort].y      = crd[1];
			SrtInt[NbrSort].IdxSol = NbrSort;
			
			
			idx = NbrSort * Msh->SolSiz;
			for (iVar=0; iVar<Msh->SolSiz; iVar++) { 
				Sol[idx+iVar] = 0.;
				Sol[idx+iVar] += alp * Msh->Sol[edgVid[1]*Msh->SolSiz+iVar];
				Sol[idx+iVar] += (1.-alp) * Msh->Sol[edgVid[0]*Msh->SolSiz+iVar];
			}
			
			NbrSort++;
			
		}
		
	}
	
	
	qsort(SrtInt,NbrSort,sizeof(struct cvt_SortInterface),CmpInterface);
	
	//--- Remove duplicated entries
	
	result = (double *) malloc(sizeof(double)*NbrSort*(2+Msh->SolSiz));
	
	NbrSortFin = 0;
	result[0] = SrtInt[0].x;
	result[1] = SrtInt[0].y;
	for (iVar=0; iVar<Msh->SolSiz; iVar++) 
		result[2+iVar] = Sol[SrtInt[0].IdxSol*Msh->SolSiz+iVar];
	NbrSortFin++;
	
	for (i=1; i<NbrSort; i++) {
		
		if ( fabs(SrtInt[i].y-SrtInt[i-1].y) < 1e-20 )
			continue;
		
		idx = NbrSortFin*(Msh->SolSiz+2);
		
		result[idx+0] = SrtInt[i].x;
		result[idx+1] = SrtInt[i].y;
		for (iVar=0; iVar<Msh->SolSiz; iVar++) 
			result[idx+2+iVar] = Sol[SrtInt[i].IdxSol*Msh->SolSiz+iVar];
			
		NbrSortFin++;
	}
	
	
	//for (i=0; i<NbrSortFin; i++) {
	//	idx = i*(Msh->SolSiz+2);
	//	printf("FINAL: (%lf %lf) : ", result[idx+0], result[idx+1]);
	//	
	//	for (iVar=0; iVar<Msh->SolSiz; iVar++) {
	//		printf (" %lf " , result[idx+2+iVar]);
	//	}
	//	printf("\n");
	//}
	
	//if (result)
	//	free(result);
	
	if ( Sol )
		free(Sol);
	
	if ( SrtInt )
		free(SrtInt);
	
	*NbrRes = NbrSortFin;
	*Siz    = Msh->SolSiz+2;
	
	return result;
	
}


double * ExtractSolutionAtRef (Options *mshopt, Mesh *Msh, int *Ref, int NbrRef,  int *NbrRes, int *Siz)
{
	
	int iEfr, iVer, v, cptVer=0, idx=0, iVar, i, flag=0;
	double *result = NULL;
	
	int *Tag = (int*) malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(Msh->NbrVer+1));
	
	SortInterface *SrtInt = (SortInterface*) malloc (sizeof(struct cvt_SortInterface)*(Msh->NbrVer+1));	
	double * Sol = (double*) malloc(sizeof(double)*(Msh->NbrVer+1)*Msh->SolSiz);
	
	//--- Tag vertices that belong to chosen bdr edges
	
	
	for (iEfr=1; iEfr<=Msh->NbrEfr; iEfr++) {
		
		flag = 0;
		for (i=0; i<NbrRef; i++) {
			if ( Msh->Efr[iEfr][2] == Ref[i] ) {
				flag = 1;
				break;
			}
		}
		
		if ( flag == 0 )
			continue;
		
		for (v=0; v<2; v++) {
			iVer = Msh->Efr[iEfr][v];
			Tag[iVer] = 1;
		}	
	}
	
	//--- Count vertices
	
	cptVer = 0;
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		
		if ( Tag[iVer] == 1 ) {
			
			SrtInt[cptVer].x      = Msh->Ver[iVer][0];
			SrtInt[cptVer].y      = Msh->Ver[iVer][1];
			SrtInt[cptVer].IdxSol = cptVer;
			
			idx = cptVer * Msh->SolSiz;
			for (iVar=0; iVar<Msh->SolSiz; iVar++) { 
				Sol[idx+iVar] = Msh->Sol[iVer*Msh->SolSiz+iVar];
			}
			
			//printf("Ver %d, dens = %lf\n", iVer, Sol[idx]);
			
			cptVer++;
		}
	
	}
	
	//--- Sort
	
	qsort(SrtInt,cptVer,sizeof(struct cvt_SortInterface),CmpInterface);
	
	//--- Alloc and output result
	
	result = (double *) malloc(sizeof(double)*cptVer*(2+Msh->SolSiz));
	
	for (i=0; i<cptVer; i++) {
		
		idx = i*(Msh->SolSiz+2);
		
		result[idx+0] = SrtInt[i].x;
		result[idx+1] = SrtInt[i].y;
		
		for (iVar=0; iVar<Msh->SolSiz; iVar++) 
			result[idx+2+iVar] = Sol[SrtInt[i].IdxSol*Msh->SolSiz+iVar];
		
		//if ( i < 10 )
		//	printf(" %d : %lf\n", i, result[idx+2]);
		
	}
	
	//--- Free memory
	
	if ( SrtInt )
		free(SrtInt);
	
	if (Sol)
		free(Sol);
	
	if (Tag)
		free(Tag);
		
	*NbrRes = cptVer;
	*Siz    = Msh->SolSiz+2;
		
	return result;
}



