#include "meshutils.h"

/*
Victorien Menier May 2017
*/

double GetThetaProj(double y, double z);
int GetEllipseNormal (double r1, double r2, double theta, double *nor);
int GetVecCur (double *crd0, double *crd1, double *vec);
int GetVecCur_theta (double *crd0, double r1, double r2, double theta, double *vec);
double NormalEllipseProjection (double r1, double r2, double *Crd, double t0);

int ProjectNozzleWall_Down_DV (double *CrdOld, double *CrdNew, CadNozzle * NozBas)
{
	
	double x=0.0, r1=0.0, r2=0.0, zcenter=0.0, alp=0.0;
	double r1_bas=0.0, r2_bas=0.0, zbas=0.0, rbas=0.0, cosTheta=0.0, theta=0.0, x_in=0.0, x_out=0.0;
	
	x_in   = NozBas->Bsp_center->Coefs[0];
	x_out  = NozBas->Bsp_center->Coefs[NozBas->Bsp_center->NbrCoefs/2-1];
	
	if ( CrdOld[0] < x_in-1e-6 || CrdOld[0] > x_out+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	alp = (CrdOld[0]-x_in)/(x_out-x_in);
	
	//--- Evaluate baseline nozzle at crd
	
	x   = CrdOld[0];
	Evaluate_Nozzle ( NozBas, &x, &r1_bas, &r2_bas, &zbas );
	
	//---
	
	double theta1=0.0, sinTheta=0.0;
	
	theta1 = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	
	sinTheta = CrdOld[1]/r1_bas;
	sinTheta = min(1.0, sinTheta);
	sinTheta = max(-1.0, sinTheta);
	double theta2 = PI_NUMBER-asin(sinTheta);
	
	theta = alp*alp*theta2+(1-alp*alp)*theta1;
	
	//---	
	
	r1   = 1.0;
	r2   = 1.0;
	zcenter = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcenter+r2*cos(theta);
	
	return 1;
	
}

int ProjectNozzleWall_Down_DV_3 (double *CrdOld, double *CrdNew, CadNozzle * Noz, CadNozzle *NozBas)
{
	
	double thetaCut, zcut, alp, theta;
	double x, r1, r2, zcen, xin, xout;
	double x_bas, r1_bas, r2_bas, zbas, xin_bas, xout_bas;
		
	//--- Get alp
	
	x_bas     = CrdOld[0];
	xin_bas   = NozBas->Bsp_center->Coefs[0];
	xout_bas  = NozBas->Bsp_center->Coefs[NozBas->Bsp_center->NbrCoefs/2-1];
	
	if ( x_bas < xin_bas-1e-6 || x_bas > xout_bas+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	alp = (x_bas-xin_bas)/(xout_bas-xin_bas);
	
	//--- Evaluate baseline nozzle
	
	Evaluate_Nozzle ( NozBas, &x_bas, &r1_bas, &r2_bas, &zbas );
	
	//--- Evaluate nozzle at crd
	
	xin   = Noz->Bsp_center->Coefs[0];
	xout  = Noz->Bsp_center->Coefs[Noz->Bsp_center->NbrCoefs/2-1];
	x = xin + alp*(xout-xin);
		
	Evaluate_Nozzle ( Noz, &x, &r1, &r2, &zcen );
	
	thetaCut  = alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
	zcut      = zcen+r2*cos(thetaCut);	
	
	//--- Newton
	
	int    count;
	double dist, ldist;
	double eps=1e-30;
	double U1[3], V1[3], U2[3], V2[3], P0[3], dx[3], UV[3];
	double uv[2],uvs[2];
	double b0, b1, a00, a10, a11, det;
	
	// Initial guess
	
	uv[0] = 0.0;
	uv[1] = PI_NUMBER-asin(CrdOld[1]/r1_bas);

	ldist = 1e308;
	
	//--- Newton-Raphson
	
	double a, b, t, t1[3], t2[3], nrm;
	double tmin = thetaCut;
	double tmax = PI_NUMBER;
	
	// Initial guess
	t = PI_NUMBER-asin(CrdOld[1]/r1_bas);
	
	CrdNew[0] = x;
	CrdNew[1] = r1*sin(t);
	CrdNew[2] = (1.0-alp)*(zcen+r2*cos(t)) + alp*zcut;
	
	return 1;
	
	for (count = 0; count < 10; count++) {
		
		if ((t < tmin) || (t > tmax))	break;
		
		P0[1] = r1*sin(t);
		P0[2] = (1.0-alp)*(zcen+r2*cos(t)) + alp*zcut;
		
		dx[0] = 0.0; 
		dx[1] = P0[1]-CrdOld[1];
		dx[2] = (P0[2]-zcen)-(CrdOld[2]-zbas);
		
		t1[0] = 0.0;
		t1[1] = r1*cos(t);
		t1[2] = -(1.0-alp)*r2*sin(t);
		nrm = sqrt(t1[1]*t1[1]+t1[2]*t1[2]);
		t1[1] /= nrm;
		t1[2] /= nrm;
		
		t2[0] = 0.0;
		t2[1] = -r1*sin(t);
		t2[2] = -(1.0-alp)*r2*cos(t);
		nrm = sqrt(t2[1]*t2[1]+t2[2]*t2[2]);
		t2[1] /= nrm;
		t2[2] /= nrm;
		
		dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
		
		if ( count == 0 )
			ldist = dist;
		
		if (dist < eps) break;
		
		b     = -( dx[0]*t1[0] +  dx[1]*t1[1] +  dx[2]*t1[2]);
		a     =  ( t1[0]*t1[0] +  t1[1]*t1[1] + t1[2]*t1[2]) +
		         ( dx[0]*t2[0] +  dx[1]*t2[1] +  dx[2]*t2[2]);
		
		if (a == 0.0)	break;
		
		printf("   %d: t=%lf   d=%le \n", count, t, ldist);
		
		b /= a;
    if (fabs(b) < 1.e-10*(tmax-tmin)) { printf("BREAK DIST %.3le\n", dist); break;}
    t += b;
		
	}
	
  if (t < tmin) t = tmin;
  if (t > tmax) t = tmax;
	 
	CrdNew[0] = x;
	CrdNew[1] = r1*sin(t);
	CrdNew[2] = (1.0-alp)*(zcen+r2*cos(t)) + alp*zcut;
	
	return 1;
	
}



double NormalEllipseProjection (double r1, double r2, double *Crd, double t0)
{
	
	double y0=Crd[1], z0=Crd[2];
	double t[2], alp[2];
	int ite=0, iteMax=30, idx=-1, j, N=100;
	double crd_cen[3]={0.0,0.0,0.0};
	double theta = -1, dtheta = -1;
	
	double vec[2], nor[2];
	
	double eps = 1e-15;
	
	t[0] = t0;//0.5*PI_NUMBER;
	t[1] = PI_NUMBER;
		
	//--- Get initial theta bounds
	
	//dtheta = (t[1]-t[0])/100;
	//
	//t[0] = GetThetaProj(Crd[1], Crd[2]);
	//
	//t[1] = min(t[0]+dtheta,PI_NUMBER);
	//t[0] = max(t[0]-dtheta,t0);
	
	double tmin=t[0], tmax=t[1];
	
	
	GetEllipseNormal (r1, r2, t[0], nor);
	//GetVecCur_theta (crd_cen, r1, r2, t[0], nor);
	GetVecCur_theta (Crd, r1, r2, t[0], vec);
	
	alp[0] = acos(nor[0]*vec[0]+nor[1]*vec[1]);	
	
	for (ite=0; ite<N-1; ite++) {
		
		t[1] = t[0]+dtheta;
		
		GetEllipseNormal (r1, r2, t[1], nor);
		//GetVecCur_theta (crd_cen, r1, r2, t[1], nor);
		GetVecCur_theta (Crd, r1, r2, t[1], vec);
		alp[1] = acos(nor[0]*vec[0]+nor[1]*vec[1]);	
		
		
		printf("INI BOUND ite %d : t : [%lf %lf], alp : [%lf,%lf]\n", ite, t[0], t[1], alp[0], alp[1]);
		
		if ( alp[1] > alp[0] ) {
			break;
		}
		
		t[0] = t[1];
		alp[0] = alp[1];
		
	}
	
	for (j=0; j<2; j++) {
		GetEllipseNormal (r1, r2, t[j], nor);
		GetVecCur_theta (Crd, r1, r2, t[j], vec);
		alp[j] = acos(nor[0]*vec[0]+nor[1]*vec[1]);		
	}
	
	//for (ite=0; ite<iteMax; ite++) {
	//	
	//	if ( alp[0] < eps ) {
	//		printf("FOUND IT! ite %d theta %lf break\n",ite, t[0]);
	//		break;
	//	}
	//	
	//	if ( alp[0] < alp[1] ) {
	//		idx = 1;
	//	}
	//	else {
	//		idx = 0;
	//	}
	//	
	//	t[idx] = 0.5*(t[0]+t[1]);
	//	
	//	GetEllipseNormal (r1, r2, t[idx], nor);
	//	GetVecCur_theta (Crd, r1, r2, t[idx], vec);
	//	alp[idx] = acos(nor[0]*vec[0]+nor[1]*vec[1]);		
	//	
	//	printf("ite %d : t : [%lf %lf], alp : [%lf,%lf]\n", ite, t[0], t[1], alp[0], alp[1]);
	//	
	//}
	
	printf("\n----------\nr1=%lf\nr2=%lf\ny0=%lf\nz0=%lf \nt0=%lf\ntfin=%lf\ntmin=%lf\ntmax=%lf\n----------\n\n", r1, r2, Crd[1], Crd[2], t0, t[0], tmin, tmax);
	
	theta = t[0];
	
	
	
	return theta;
	
	
}

int ProjectNozzleWall_Down_DV_2 (double *CrdOld, double *CrdNew, CadNozzle * NozBas)
{
	
	double x=0.0, r1=0.0, r2=0.0, zcenter=0.0, alp=0.0;
	double r1_bas=0.0, r2_bas=0.0, zbas=0.0, rbas=0.0, cosTheta=0.0, theta=0.0, x_in=0.0, x_out=0.0;
	double crd[3];
	
	x_in   = NozBas->Bsp_center->Coefs[0];
	x_out  = NozBas->Bsp_center->Coefs[NozBas->Bsp_center->NbrCoefs/2-1];
	
	if ( CrdOld[0] < x_in-1e-6 || CrdOld[0] > x_out+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	alp = (CrdOld[0]-x_in)/(x_out-x_in);
	
	//--- Evaluate baseline nozzle at crd
	
	x   = CrdOld[0];
	Evaluate_Nozzle ( NozBas, &x, &r1_bas, &r2_bas, &zbas );
	
	crd[0] = 0.0;
	crd[1] = CrdOld[1];
	crd[2] = CrdOld[2]-zbas;
	
	double t0 = alp*NozBas->ThetaCutOut + (1.0-alp)*NozBas->ThetaCutIn;
	
	theta = NormalEllipseProjection (r1_bas, r2_bas, crd, t0);
	
	//printf("EXIT\n");
	//exit(1);
	
	//---
	
	//double theta1=0.0, sinTheta=0.0;
	//
	//theta1 = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	//
	//sinTheta = CrdOld[1]/r1_bas;
	//sinTheta = min(1.0, sinTheta);
	//sinTheta = max(-1.0, sinTheta);
	//double theta2 = PI_NUMBER-asin(sinTheta);
	//
	//theta = alp*alp*theta2+(1-alp*alp)*theta1;
	//
	//---	
	
	r1   = 1.0;
	r2   = 1.0;
	zcenter = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcenter+r2*cos(theta);
	
	return 1;
	
}




int ProjectNozzleWall_Up_DV (double *CrdOld, double *CrdNew, CadNozzle * NozBas)
{
	
	double x=0.0, r1=0.0, r2=0.0, zcenter=0.0, alp=0.0;
	double r1_bas=0.0, r2_bas=0.0, zbas=0.0, rbas=0.0, cosTheta=0.0, theta=0.0, x_in=0.0, x_out=0.0;
	
	x_in   = NozBas->Bsp_center->Coefs[0];
	x_out  = NozBas->Bsp_center->Coefs[NozBas->Bsp_center->NbrCoefs/2-1];
	
	
	int verbose = 1;
	
	if ( CrdOld[0] < x_in-1e-6 || CrdOld[0] > x_out+1e-6  ) {
		if ( verbose > 0 )
			printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	alp = (CrdOld[0]-x_in)/(x_out-x_in);
	
	//--- Evaluate baseline nozzle at crd
	
	x   = CrdOld[0];
	Evaluate_Nozzle ( NozBas, &x, &r1_bas, &r2_bas, &zbas );
	
	rbas  = CrdOld[1]*CrdOld[1]/(r1_bas*r1_bas);
	rbas += (CrdOld[2]-zbas)*(CrdOld[2]-zbas)/(r2_bas*r2_bas);
	rbas = sqrt(rbas);
	cosTheta = (CrdOld[2]-zbas)/r2_bas/rbas;
	theta = acos(max(-1.0,min(1.0, cosTheta)));
	
	//--- CrdNew
	
	r1   = 1.0;
	r2   = 1.0;
	zcenter = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcenter+r2*cos(theta);
	
	return 0;
	
}

int ProjectNozzleWall_Up (double *CrdOld, double *CrdNew, double *BasParam)
{
	double r_inlet      = BasParam[3];              // r inlet
	
	double r_outlet[2] = {BasParam[7], BasParam[8]};            // r1 r2 outlet
	
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]};
	
	double zbas = fz_bas (CrdOld[0],  BasParam);
	
	double cosTheta = 0.0, theta = 0.0, r1=0.0, r2=0.0, zcen=0.0;
	
	double r1_bas = fr1_bas (CrdOld[0],  BasParam);
	double r2_bas = fr2_bas (CrdOld[0],  BasParam);
	
	double rbas = 0.0;
	
	double alp = 0.0;
	
	double thetaMax = alp*BasParam[BasOutletTheta] + (1.0-alp)*BasParam[BasInletTheta];
	
	rbas += CrdOld[1]*CrdOld[1]/(r1_bas*r1_bas);
	rbas += (CrdOld[2]-zbas)*(CrdOld[2]-zbas)/(r2_bas*r2_bas);
	rbas = sqrt(rbas);
	
	cosTheta = (CrdOld[2]-zbas)/r2_bas/rbas;
	
	cosTheta = min(1.0, cosTheta);
	cosTheta = max(-1.0, cosTheta);
	
	theta = acos(cosTheta);
	
	r1 = 1.0;
	r2 = 1.0;
	zcen = 0.0;
	
	alp = (CrdOld[0]-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcen+r2*cos(theta);
	
	return 0;
	
}


double GetThetaProj(double y, double z)
{
	
	double  len, cosTheta, theta;
	
	len = sqrt(y*y+z*z);
	
	cosTheta = min(max(z/len, -1.0), 1.0);
	theta = acos(cosTheta);
	
	return theta;	
}

double GetThetaProj_2(double y, double z, double r1,  double r2)
{
	
	double  len, cosTheta, theta;
	
	
	z /= r2;
	y /= r1;
	
	len = sqrt(y*y+z*z);
	
	cosTheta = min(max(z/len, -1.0), 1.0);
	theta = acos(cosTheta);
	
	return theta;	
}

int ProjectNozzleWall_Down_Save (double *CrdOld, double *CrdNew, double *BasParam)
{
	double cen_inlet[3]  = {BasParam[0], BasParam[1], BasParam[2]}; 
	double cen_outlet[3] = {BasParam[4], BasParam[5], BasParam[6]}; 
	
	
	if ( CrdOld[0] < cen_inlet[0]-1e-6 || CrdOld[0] > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	double alp = (CrdOld[0]-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	double zbas = fz_bas (CrdOld[0],  BasParam);
	
	double cosTheta = 0.0, theta = 0.0, r1=0.0, r2=0.0, zcen=0.0, sinTheta=0.0;
	double y0, alp0, theta0, theta1;
	double r1_bas = fr1_bas (CrdOld[0],  BasParam);
	double r2_bas = fr2_bas (CrdOld[0],  BasParam);
	double zcut = fzcut (CrdOld[0],  BasParam);
	
	double z1, z2, z3;
	double nor[4];
	
	// ---- 
	
	theta1 = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	
	
	//theta = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	//theta = GetThetaProj_2(CrdOld[1], CrdOld[2]-zbas, r1_bas, r2_bas); 
	
	sinTheta = CrdOld[1]/r1_bas;
	sinTheta = min(1.0, sinTheta);
	sinTheta = max(-1.0, sinTheta);
	double theta2 = PI_NUMBER-asin(sinTheta);
	
	theta = alp*alp*theta2+(1-alp*alp)*theta1;
	
	r1 = 1.0;
	r2 = 1.0;
	zcen = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcen+r2*cos(theta);
	
	return 1;
	
}


double ComputeLength (double *crd0, double *crd1) 
{
	int j;
	double len=0.0;
	
	for (j=0; j<3; j++) {
		len += (crd1[j]-crd0[j])*(crd1[j]-crd0[j]);
	}
	
	return sqrt(len);
}

int GetVecCur (double *crd0, double *crd1, double *vec) 
{
	int j;
	double nrm=0.0;
	
	vec[0] = crd1[1]-crd0[1];
	vec[1] = crd1[2]-crd0[2];
			
	for (j=0; j<2; j++) {
		nrm += vec[j]*vec[j];
	}
	
	nrm = sqrt(nrm);
	for (j=0; j<2; j++)
		vec[j] /= nrm;	
	
	return 1;
}


int GetVecCur_theta (double *crd0, double r1, double r2, double theta, double *vec) 
{
	int j;
	double nrm=0.0;
	
	double crd[3];
	
	crd[1] = r1*sin(theta);
	crd[2] = r2*cos(theta);

	vec[0] = crd0[1]-crd[1];
	vec[1] = crd0[2]-crd[2];
			
	for (j=0; j<2; j++) {
		nrm += vec[j]*vec[j];
	}
	
	nrm = sqrt(nrm);
	for (j=0; j<2; j++)
		vec[j] /= nrm;	
	
	return 1;
}

int GetEllipseNormal (double r1, double r2, double theta, double *nor) 
{
  
	double nrm;
	
  nor[0] = -r1*sin(theta);
  nor[1] = -r2*cos(theta);
  nrm = sqrt(nor[0]*nor[0]+nor[1]*nor[1]);
  nor[0] /= nrm;
  nor[1] /= nrm;
  
  return nor;
}


int ProjectNozzleWall_Down (double *CrdOld, double *CrdNew, double *BasParam)
{
	
	
	double cen_inlet[3]  = { BasParam[BasInletx], \
													 BasParam[BasInlety], \
													 BasParam[BasInletz]}; 
												
	double cen_outlet[3]  = { BasParam[BasOutletx], \
	                     		  BasParam[BasOutlety], \
	                     		  BasParam[BasOutletz]}; 
	
	
	if ( CrdOld[0] < cen_inlet[0]-1e-6 || CrdOld[0] > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
	
	double alp = (CrdOld[0]-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	double zbas   = fz_bas  (CrdOld[0],  BasParam);
	double r1_bas = fr1_bas (CrdOld[0],  BasParam);
	double r2_bas = fr2_bas (CrdOld[0],  BasParam);
	double zcut   = fzcut   (CrdOld[0],  BasParam);
	
	//---
	//--- Perform dichotomy to find theta
	//---
	
	double theta1 = GetThetaProj(CrdOld[1], CrdOld[2]-zbas); 
	
	double sinTheta = CrdOld[1]/r1_bas;
	sinTheta = min(1.0, sinTheta);
	sinTheta = max(-1.0, sinTheta);
	double theta2 = PI_NUMBER-asin(sinTheta);
	
	double theta = alp*alp*theta2+(1-alp*alp)*theta1;
	
	double thetaMin = alp*BasParam[BasOutletTheta] + (1.0-alp)*BasParam[BasInletTheta];
	
	//--- Compute thetaMin based on zcut, not thetacut
	
	double cosTheta = max(-1.0, min(1.0,(zcut-zbas)/r2_bas));
	thetaMin = acos(cosTheta);
	
	//printf("thetaMin %lf New %lf\n", thetaMin, theta);
	
	double thetaMax = PI_NUMBER-1e-6;
	
	double dTheta = 0.03*PI_NUMBER;
	
	int ite=0, maxIte=30, j;
		
	double t0 = thetaMin, t1 = thetaMax, t01 = 0.0;
	double l0=0.0, l1=0.0, n0=0.0, n1=0.0, nrm=0.0;
	
	t0 = max(thetaMin,theta-dTheta);
	t1 = min(thetaMax,theta+dTheta);
	
	double crdCur[3] = {0.0,0.0,0.0};
	double crd[3] = {0.0,0.0,0.0};
	double vec[3] = {0.0,0.0,0.0};
	double nor[3]    = {0.0,0.0,0.0};
	
	crdCur[0] = CrdOld[0];
	crdCur[1] = CrdOld[1];
	crdCur[2] = CrdOld[2] - zbas;
	
	double t[2] = {0.0,0.0};
	double n[2] = {0.0,0.0};
	
	t[0] = max(thetaMin,theta-dTheta);
	t[1] = min(thetaMax,theta+dTheta);
	
	
	FILE *ellHdl = fopen("ellipse.dat", "wb");
	
	fprintf(ellHdl, "%le %le %le %le\n", r1_bas, r2_bas, crdCur[1], crdCur[2]);
	fclose(ellHdl);
	
	FILE *cvgHdl = fopen("cvgNor.dat", "wb");
	
	while( ite < maxIte ) {
		
		printf(" -- ITE %d : t0 %lf t1 %lf\n", ite, t[0], t[1]);
		
		for (j=0; j<2; j++) {
			
			//--- 
			
			crd[0] = CrdOld[0];
			crd[1] = r1_bas*sin(t[j]);
			crd[2] = r2_bas*cos(t[j]);
			
			//--- Normal to ellipse
			
			nor[0] = 0.0;
			nor[1] = -r1_bas*sin(t[j]);
			nor[2] = -r2_bas*cos(t[j]);
			nrm = sqrt(nor[1]*nor[1]+nor[2]*nor[2]);
			nor[1] /= nrm;
			nor[2] /= nrm;
			
			//--- Vector
			
			GetVecCur (crd, crdCur, vec);
			
			printf ("vec %d = %lf %lf\n", j, vec[1], vec[2]);
			printf ("nor %d = %lf %lf\n", j, nor[1], nor[2]);
			
			n[j] = (vec[1]-nor[1])*(vec[1]-nor[1])+(vec[2]-nor[2])*(vec[2]-nor[2]);
			fprintf(cvgHdl, "%le %le %le %le %le %le ", crd[1], crd[2], nor[1], nor[2], vec[1], vec[2]);
		}
		fprintf(cvgHdl, "\n");
		
		if ( fabs(n[0]-n[1]) < 1e-4 ) {
			theta = t[0];
			break;
		}
		
		if ( n[0] > n[1] ) {
			t[0] = 0.5*(t[0]+t[1]);
			theta = t[1];
		}
		else {
			t[1] = 0.5*(t[0]+t[1]);
			theta = t[0];
		}
		
		
		printf("n0 %lf n1 %lf\n", n[0], n[1]);
		
		ite++;
		
	}
	
	fclose(cvgHdl);
	
	
	//printf("thetamin %lf thetamax %lf theta %lf\n", thetaMin, thetaMax, theta);
	//
	//if ( theta < 0.21 )
	//	exit(1);
	

	if ( theta < thetaMin ) {
		printf("ERROR\n");
		exit(1);
	}
	
	//exit(1);
	
	// --- Project to cylinder
	
	double r1, r2, zcen;
	
	r1 = 1.0;
	r2 = 1.0;
	zcen = 0.0;
	
	CrdNew[0] = alp;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcen+r2*cos(theta);
	
	return 1;
	
}


int ProjectToDV(double *Crd, CadNozzle *Noz, double *BasParam, int Patch) 
{
	
	int j=0;
	double  CrdNew[3] = {0.0,0.0,0.0}, alp=0.0, theta=0.0, zcut=0.0, thetaBas=0.0, ycut=0.0;
	
	CadBspline * Bsp_center = Noz->Bsp_center;
	CadBspline * Bsp_r1     = Noz->Bsp_r1;
	CadBspline * Bsp_r2     = Noz->Bsp_r2;
		
	double len = Bsp_center->Coefs[Bsp_center->NbrCoefs/2-1];
	double x = Bsp_center->Coefs[0]+Crd[0]*len, zcenter=0.0, r1=0.0, r2=0.0, dydx=0.0;
	
	
	//--- Centerline
	bSplineGeo3 (Bsp_center->Knots, Bsp_center->Coefs, &x, &zcenter, &dydx, 1, Bsp_center->NbrKnots, Bsp_center->NbrCoefs/2);
	
	//--- R1
	bSplineGeo3 (Bsp_r1->Knots, Bsp_r1->Coefs, &x, &r1, &dydx, 1, Bsp_r1->NbrKnots, Bsp_r1->NbrCoefs/2);
	
	//--- R2
	bSplineGeo3 (Bsp_r2->Knots, Bsp_r2->Coefs, &x, &r2, &dydx, 1, Bsp_r2->NbrKnots, Bsp_r2->NbrCoefs/2);
	
	//printf("x %lf ycenter %lf r1 %lf r2 %lf\n", x, zcenter, r1, r2);
	
	//--- Theta correction
	double thetaMaxBas, thetaMaxDV;
	
	alp = Crd[0];
	
	//double xbas = Crd[0]*(BasParam[BasOutletx]-BasParam[BasInletx]);
	//double zbas   = fz_bas  (xbas,  BasParam);
	//double r2_bas = fr2_bas (xbas,  BasParam);
	//zcut   = fzcut   (xbas,  BasParam);
	//double cosTheta = max(-1.0, min(1.0,(zcut-zbas)/r2_bas));
	//thetaMaxBas = acos(cosTheta);
	
	thetaMaxBas = alp*BasParam[BasOutletTheta] + (1.0-alp)*BasParam[BasInletTheta];
	thetaMaxDV  = alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
	
	theta = acos(Crd[2]);
	theta *= thetaMaxDV/thetaMaxBas;
	
	//printf("CMP : acos %lf sin %lf sin(acos) %lf\n", acos(Crd[2]), Crd[1], sin(acos(Crd[2])));
	//printf(" sin(theta) %lf, cos(theta) %lf, Crd[1] %lf, Crd[2] %lf",  sin(theta), cos(theta), Crd[1], Crd[2]);
	//printf("thetaMaxDV %lf, thetaMaxBas %lf\n", thetaMaxDV, thetaMaxBas);
	
	//printf("thetaMaxBas %lf thetaMaxDV %lf\n", thetaMaxBas, thetaMaxDV);
	//CrdNew[0] = x;
	//CrdNew[1] = r1*Crd[1];
	//CrdNew[2] = zcenter+r2*Crd[2];
	
	CrdNew[0] = x;
	CrdNew[1] = r1*sin(theta);
	CrdNew[2] = zcenter+r2*cos(theta);
	
	if ( Patch == NOZZLEDOWN ) {
		
		//theta = thetaMaxBas;//alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
		//zcut  = zcenter+r2*cos(theta);		
		
		zcut = zcenter+r2*cos(thetaMaxDV);
		
		//printf("alp %lf in %lf out %lf cut %lf zcut %lf theta(crd) = %lf\n", alp, Noz->ThetaCutIn/PI_NUMBER, Noz->ThetaCutOut/PI_NUMBER, theta/PI_NUMBER, zcut, acos(Crd[2])/PI_NUMBER);
		alp = Crd[0];
		//CrdNew[2] = zcut;
		CrdNew[2] = alp*zcut + (1.0-alp)*CrdNew[2];
		CrdNew[2] = zcut;
	}
	
	//---
	
	Crd[0] = CrdNew[0];
	Crd[1] = CrdNew[1];
	Crd[2] = CrdNew[2];
	
	return 0;
	
}


int ProjectToDV_DV(double *Crd, CadNozzle *Noz, CadNozzle *Noz_bas, int Patch) 
{
	
	int j=0;
	double  CrdNew[3] = {0.0,0.0,0.0}, alp=0.0, theta=0.0, zcut=0.0, thetaBas=0.0, ycut=0.0;
	
	CadBspline * Bsp_center = Noz->Bsp_center;
	CadBspline * Bsp_r1     = Noz->Bsp_r1;
	CadBspline * Bsp_r2     = Noz->Bsp_r2;
		
	double len = Bsp_center->Coefs[Bsp_center->NbrCoefs/2-1];
	double x = Bsp_center->Coefs[0]+Crd[0]*len, zcenter=0.0, r1=0.0, r2=0.0, dydx=0.0;
	
	//--- Centerline
	bSplineGeo3 (Bsp_center->Knots, Bsp_center->Coefs, &x, &zcenter, &dydx, 1, Bsp_center->NbrKnots, Bsp_center->NbrCoefs/2);
	
	//--- R1
	bSplineGeo3 (Bsp_r1->Knots, Bsp_r1->Coefs, &x, &r1, &dydx, 1, Bsp_r1->NbrKnots, Bsp_r1->NbrCoefs/2);
	
	//--- R2
	bSplineGeo3 (Bsp_r2->Knots, Bsp_r2->Coefs, &x, &r2, &dydx, 1, Bsp_r2->NbrKnots, Bsp_r2->NbrCoefs/2);
	
	//printf("x %lf ycenter %lf r1 %lf r2 %lf\n", x, zcenter, r1, r2);
	
	//--- Theta correction
	double thetaMaxBas, thetaMaxDV;
	
	alp = Crd[0];
	
	//double xbas = Crd[0]*(BasParam[BasOutletx]-BasParam[BasInletx]);
	//double zbas   = fz_bas  (xbas,  BasParam);
	//double r2_bas = fr2_bas (xbas,  BasParam);
	//zcut   = fzcut   (xbas,  BasParam);
	//double cosTheta = max(-1.0, min(1.0,(zcut-zbas)/r2_bas));
	//thetaMaxBas = acos(cosTheta);
	
	thetaMaxBas = alp*Noz_bas->ThetaCutOut + (1.0-alp)*Noz_bas->ThetaCutIn;
	thetaMaxDV  = alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
	
	//printf("thetaMaxBas %lf thetaMaxDV %lf\n", thetaMaxBas, thetaMaxDV);
	
	CrdNew[0] = x;
	CrdNew[1] = r1*Crd[1];
	CrdNew[2] = zcenter+r2*Crd[2];
	
	if ( Patch == NOZZLEDOWN ) {
		
		theta = thetaMaxBas;//alp*Noz->ThetaCutOut + (1.0-alp)*Noz->ThetaCutIn;
		zcut  = zcenter+r2*cos(theta);		
		
		//printf("alp %lf in %lf out %lf cut %lf zcut %lf theta(crd) = %lf\n", alp, Noz->ThetaCutIn/PI_NUMBER, Noz->ThetaCutOut/PI_NUMBER, theta/PI_NUMBER, zcut, acos(Crd[2])/PI_NUMBER);
		alp = Crd[0];
		//CrdNew[2] = zcut;
		CrdNew[2] = alp*zcut + (1.0-alp)*CrdNew[2];
		
	}
	
	//---
	
	Crd[0] = CrdNew[0];
	Crd[1] = CrdNew[1];
	Crd[2] = CrdNew[2];
	
	return 0;
	
}



int ProjectToDV_test(double *Crd, CadNozzle *Noz, double *BasParam, int Patch) 
{
	
	//--- Baseline parameters
	
	double cen_inlet[3]  = { BasParam[BasInletx], \
													 BasParam[BasInlety], \
													 BasParam[BasInletz]}; 
												
	double cen_outlet[3]  = { BasParam[BasOutletx], \
	                     		  BasParam[BasOutlety], \
	                     		  BasParam[BasOutletz]}; 
	
	double cut_inlet[3]  = {0,  0.439,  0.0999};
	double cut_outlet[3] = {2.33702,  0.88272, 0.122373};
	
	if ( Crd[0] < cen_inlet[0]-1e-6 || Crd[0] > cen_outlet[0]+1e-6  ) {
		printf("  ## ERROR ProjectNozzleWall_Down : x out of range!\n");
		exit(1);
	}
		
	double zbas=0, r1_bas=0, r2_bas=0, alp=0;
	double ycut_bas=0, zcut_bas=0;
	double x=0, y=0, z=0;
	double r1=0, r2=0, dydx=0, zcenter=0;
	double ycut=0, zcut=0;
	double alpy=0, alpz=0;
	
	//--- Evaluate baseline at current x location
	
	x = Crd[0];
	y = Crd[1];
	z = Crd[2];
	
	alp = (x-cen_inlet[0])/(cen_outlet[0]-cen_inlet[0]);
	
	r1_bas = fr1_bas (x,  BasParam);
	r2_bas = fr2_bas (x,  BasParam);
	zbas   = fz_bas  (x,  BasParam);
	
	ycut_bas = alp*cut_outlet[1] + (1.0-alp)*cut_inlet[1];
	zcut_bas = alp*cut_outlet[2] + (1.0-alp)*cut_inlet[2];
	
	//--- Evaluate bsplines to project on
	
	CadBspline * Bsp_center = Noz->Bsp_center;
	CadBspline * Bsp_r1     = Noz->Bsp_r1;
	CadBspline * Bsp_r2     = Noz->Bsp_r2;
	
	//--- Centerline
	bSplineGeo3 (Bsp_center->Knots, Bsp_center->Coefs, &x, &zcenter, &dydx, 1, Bsp_center->NbrKnots, Bsp_center->NbrCoefs/2);
	
	//--- R1
	bSplineGeo3 (Bsp_r1->Knots, Bsp_r1->Coefs, &x, &r1, &dydx, 1, Bsp_r1->NbrKnots, Bsp_r1->NbrCoefs/2);
	
	//--- R2
	bSplineGeo3 (Bsp_r2->Knots, Bsp_r2->Coefs, &x, &r2, &dydx, 1, Bsp_r2->NbrKnots, Bsp_r2->NbrCoefs/2);
	
	Crd[1] = y/r1_bas*r1;
	Crd[2] = zcenter+(z-zbas)/r2_bas*r2;
	
	if ( Patch == NOZZLEDOWN ) {
					
		ycut = ycut_bas/r1_bas*r1;
		zcut = zcenter+(zcut_bas-zbas)/r2_bas*r2;
		
		alpy = y/ycut_bas;
		alpz = z/zcut_bas;
		
		Crd[1] = alp*alpy*ycut + (1.0-alp)*Crd[1];
		Crd[2] = alp*alpz*zcut + (1.0-alp)*Crd[2];
					
	}
	
	return 0;
}



int NozzleWallProjection_test (Options *mshopt, Mesh *Msh, CadNozzle * CadNoz, int refUp, int refDown, char *OutNam)
{
		
	//---
	
	int WrtMsh = 1, patch=-1, WrtFullMsh=0;
		
	int iTri, ref=0, i, j, iVer, vid=-1, is[3];
	
	double CrdNew[3];
	
	int NbrVerPrj=0 , NbrTriPrj=0;
	
	int *Tag = (int*) malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(Msh->NbrVer+1));
	
	//--- Parameters
	
	double BasParam[BasMaxKwd];
	memset(BasParam, 0.0, sizeof(double)*BasMaxKwd);
	
	BasParam[BasInletx    ]  = 0.0       ;
	BasParam[BasInlety    ]  = 0.0       ;
	BasParam[BasInletz    ]  = 0.099908 ;
	BasParam[BasInletr    ]  = 0.439461   ;
	BasParam[BasOutletx   ]  = 2.33702   ;
	BasParam[BasOutlety   ]  = 0.0       ;
	BasParam[BasOutletz   ]  = 0.19      ;
	BasParam[BasOutletr1  ]  = 0.92      ;
	BasParam[BasOutletr2  ]  = 0.24      ;
	BasParam[BasOutletzcut]  = 0.122638  ;
	BasParam[BasInletzcut ]  = 0.099     ;
	BasParam[BasOutletycut]  = 0.879257  ;
	BasParam[BasInletycut ]  = 0.438972  ;
	BasParam[BasInletTheta]  = 1.572865  ;
	BasParam[BasOutletTheta] = 1.855294  ;
	
	
	//--- Tag one ref after the other for consistent line re-projection
	
	NbrTriPrj = 0;
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refDown )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	} 
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refUp  )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	}	
	
	//--- Create mesh for visu
	
	NbrVerPrj = 0;
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] > 0 )
			NbrVerPrj++;
	}
	
	int SizMsh[GmfMaxKwd+1];
	
	for (i=0; i<GmfMaxKwd; i++)
		SizMsh[i] = 0;
	
	SizMsh[GmfDimension]  = 3;
	SizMsh[GmfVertices]   = NbrVerPrj;
	SizMsh[GmfTriangles]  = NbrTriPrj;
	
	//printf("NbrVerPrj %d NbrRTriPrj %d\n", NbrVerPrj, NbrTriPrj);
	
	Mesh *MshViz = NULL;
	
	if ( WrtMsh == 1 ) {
		MshViz = AllocMesh(SizMsh);
		MshViz->Dim = 3;		
	}	
		
	FILE *motionHdl = fopen(OutNam, "wb");
	
	if ( motionHdl ) {
		printf("%%%% %s OPENED (WRITE).\n", OutNam);
	}
	else
	{
		printf("  ## ERROR : Could not open %s.\n", OutNam);
		exit(1);
	}
	
	//--- Project points
		
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		if ( Tag[iVer] <= 0 )
			continue;
		
		CrdNew[0] = Msh->Ver[iVer][0];
		CrdNew[1] = Msh->Ver[iVer][1];
		CrdNew[2] = Msh->Ver[iVer][2];
		
		patch = -1;
		
		
		if ( Tag[iVer] == refUp ) {
			//ProjectNozzleWall_Up (Msh->Ver[iVer], CrdNew, BasParam);			
			patch = NOZZLEUP;
		}
		else if ( Tag[iVer] == refDown ) {
			
			//---
			//ProjectNozzleWall_Down_Save (Msh->Ver[iVer], CrdNew, BasParam);
			patch = NOZZLEDOWN;
		}
		
		if ( patch == -1 ) {
			printf("  ## ERROR : Wrong nozzle CAD patch.\n");
			exit(1);
		}
				
		ProjectToDV_test(CrdNew, CadNoz, BasParam, patch);
		
		if ( WrtMsh == 1 ) {
			MshViz->NbrVer++;
			AddVertex(MshViz, MshViz->NbrVer, CrdNew);
			Tag[iVer] = MshViz->NbrVer;
		}
		
		if (WrtFullMsh == 1) {
			Msh->Ver[iVer][0] = CrdNew[0];
			Msh->Ver[iVer][1] = CrdNew[1];
			Msh->Ver[iVer][2] = CrdNew[2];
		}
		
		fprintf(motionHdl, "%d %le %le %le\n", iVer-1, CrdNew[0], CrdNew[1], CrdNew[2]);
		
	}
	
	//--- Write visu mesh
	
	if ( WrtMsh == 1 ) {
		
		for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
			ref = Msh->Tri[iTri][3];
			if ( ref != refUp && ref != refDown )
				continue;
			
			for (j=0; j<3; j++) {
				vid = Msh->Tri[iTri][j];
				is[j] = Tag[vid];
			}
			
			MshViz->NbrTri++;
			AddTriangle(MshViz,MshViz->NbrTri,is,ref);
		}
		
		printf("%d ver, %d tri\n", MshViz->NbrVer, MshViz->NbrTri);
		WriteGMFMesh("visu", MshViz, 1);
	}
	
	if (WrtFullMsh == 1) {
		WriteGMFMesh("visu_full", Msh, 1);
	}
	
	if (motionHdl)
		fclose(motionHdl);
	
	if ( Tag )
		free(Tag);
	
	if ( MshViz )
		FreeMesh(MshViz);
	
	if ( CadNoz )
		FreeCadNozzle (CadNoz);
	
	return 0;
	
}




int NozzleWallProjection (Options *mshopt, Mesh *Msh, CadNozzle * CadNoz, int refUp, int refDown, char *OutNam)
{
		
	//---
	
	int WrtMsh = 1, patch=-1, WrtFullMsh=1;
		
	int iTri, ref=0, i, j, iVer, vid=-1, is[3];
	
	double CrdNew[3];
	
	int NbrVerPrj=0 , NbrTriPrj=0;
	
	int *Tag = (int*) malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(Msh->NbrVer+1));
	
	//--- Parameters
	
	double BasParam[BasMaxKwd];
	memset(BasParam, 0.0, sizeof(double)*BasMaxKwd);
	
	//BasParam[BasInletx    ]  = 0.0       ;
	//BasParam[BasInlety    ]  = 0.0       ;
	//BasParam[BasInletz    ]  = 0.0999079 ;
	//BasParam[BasInletr    ]  = 0.438972  ;
	//BasParam[BasOutletx   ]  = 2.33702   ;
	//BasParam[BasOutlety   ]  = 0.0       ;
	//BasParam[BasOutletz   ]  = 0.19      ;
	//BasParam[BasOutletr1  ]  = 0.92      ;
	//BasParam[BasOutletr2  ]  = 0.24      ;
	//BasParam[BasOutletzcut]  = 0.122638  ;
	//BasParam[BasInletzcut ]  = 0.099     ;
	//BasParam[BasOutletycut]  = 0.879257  ;
	//BasParam[BasInletycut ]  = 0.438972  ;
	//BasParam[BasInletTheta]  = 1.572865  ;
	//BasParam[BasOutletTheta] = 1.855294  ;
	
	
	BasParam[BasInletx    ]  = 0.0       ;
	BasParam[BasInlety    ]  = 0.0       ;
	BasParam[BasInletz    ]  = 0.099908 ;
	BasParam[BasInletr    ]  = 0.439461   ;
	BasParam[BasOutletx   ]  = 2.33702   ;
	BasParam[BasOutlety   ]  = 0.0       ;
	BasParam[BasOutletz   ]  = 0.19      ;
	BasParam[BasOutletr1  ]  = 0.92      ;
	BasParam[BasOutletr2  ]  = 0.24      ;
	BasParam[BasOutletzcut]  = 0.122638  ;
	BasParam[BasInletzcut ]  = 0.099     ;
	BasParam[BasOutletycut]  = 0.879257  ;
	BasParam[BasInletycut ]  = 0.438972  ;
	BasParam[BasInletTheta]  = 1.572865  ;
	BasParam[BasOutletTheta] = 1.855294  ;
	
	
	//--- Tag one ref after the other for consistent line re-projection
	
	NbrTriPrj = 0;
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refDown )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	} 
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refUp  )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	}	
	
	//--- Create mesh for visu
	
	NbrVerPrj = 0;
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] > 0 )
			NbrVerPrj++;
	}
	
	int SizMsh[GmfMaxKwd+1];
	
	for (i=0; i<GmfMaxKwd; i++)
		SizMsh[i] = 0;
	
	SizMsh[GmfDimension]  = 3;
	SizMsh[GmfVertices]   = NbrVerPrj;
	SizMsh[GmfTriangles]  = NbrTriPrj;
	
	//printf("NbrVerPrj %d NbrRTriPrj %d\n", NbrVerPrj, NbrTriPrj);
	
	Mesh *MshViz = NULL;
	
	if ( WrtMsh == 1 ) {
		MshViz = AllocMesh(SizMsh);
		MshViz->Dim = 3;		
	}	
		
	FILE *motionHdl = fopen(OutNam, "wb");
	
	if ( motionHdl ) {
		printf("%%%% %s OPENED (WRITE).\n", OutNam);
	}
	else
	{
		printf("  ## ERROR : Could not open %s.\n", OutNam);
		exit(1);
	}
	
	//--- Project points
		
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		if ( Tag[iVer] <= 0 )
			continue;
		
		CrdNew[0] = Msh->Ver[iVer][0];
		CrdNew[1] = Msh->Ver[iVer][1];
		CrdNew[2] = Msh->Ver[iVer][2];
		
		patch = -1;
		
		if ( Tag[iVer] == refUp ) {
			ProjectNozzleWall_Up (Msh->Ver[iVer], CrdNew, BasParam);			
			patch = NOZZLEUP;
		}
		else if ( Tag[iVer] == refDown ) {
			
			//---
			ProjectNozzleWall_Down_Save (Msh->Ver[iVer], CrdNew, BasParam);
			patch = NOZZLEDOWN;
		}
		
		if ( patch == -1 ) {
			printf("  ## ERROR : Wrong nozzle CAD patch.\n");
			exit(1);
		}
		
		//ProjectToDV(CrdNew, CadNoz, BasParam, patch);
		
		if ( WrtMsh == 1 ) {
			MshViz->NbrVer++;
			AddVertex(MshViz, MshViz->NbrVer, CrdNew);
			Tag[iVer] = MshViz->NbrVer;
		}
		
		if (WrtFullMsh == 1) {
			Msh->Ver[iVer][0] = CrdNew[0];
			Msh->Ver[iVer][1] = CrdNew[1];
			Msh->Ver[iVer][2] = CrdNew[2];
		}
		
		fprintf(motionHdl, "%d %le %le %le\n", iVer-1, CrdNew[0], CrdNew[1], CrdNew[2]);
		
	}
	
	//--- Write visu mesh
	
	if ( WrtMsh == 1 ) {
		
		for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
			ref = Msh->Tri[iTri][3];
			if ( ref != refUp && ref != refDown )
				continue;
			
			for (j=0; j<3; j++) {
				vid = Msh->Tri[iTri][j];
				is[j] = Tag[vid];
			}
			
			MshViz->NbrTri++;
			AddTriangle(MshViz,MshViz->NbrTri,is,ref);
		}
		
		printf("%d ver, %d tri\n", MshViz->NbrVer, MshViz->NbrTri);
		WriteGMFMesh("visu", MshViz, 1);
	}
	
	if (WrtFullMsh == 1) {
		WriteGMFMesh("visu_full", Msh, 1);
	}
	
	if (motionHdl)
		fclose(motionHdl);
	
	if ( Tag )
		free(Tag);
	
	if ( MshViz )
		FreeMesh(MshViz);
	
	if ( CadNoz )
		FreeCadNozzle (CadNoz);
	
	return 0;
	
}



//int GetThetaInOut(Mesh *Msh, CadNozzle *Noz_bas)
//{
//	
//	int iVer, iTri, j, ref;
//	double x_in=1e6, x_out=-1e6, ymin_in=1e6, ymin_out=-1e6;
//	double xmin=0, xmax=0;
//	
//	int *tag = (int*)malloc(sizeof(int)*(Msh->NbrVer+1));
//	memset(tag,0,sizeof(int)*(Msh->NbrVer+1));
//	
//	//x_in = Noz_bas->Bsp_center->Coefs[0];
//	//x_out = x_in + Noz_bas->Bsp_center->Coefs[Noz_bas->Bsp_center->NbrCoefs/2-1];
//	
//	xmin = 1e6;
//	xmax = -1e6;
//	
//	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
//		ref = Msh->Tri[iTri][3];
//		if ( ref == refUp  ) {
//			for (j=0; j<3; j++) {
//				iVer = Msh->Tri[iTri][j];
//				
//				xmin = min(xmin, Msh->Ver[iVer][0]);
//				
//				if ( fabs(Msh->Ver[iVer][0]-x_in) < 1e-6 ) {
//					zmin_in = min(zmin_in, Msh->Ver[iVer][2]);
//				}
//			}
//		}
//	}
//	
//	
//	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
//		
//		if ( )
//		
//		if ( fabs(Msh->Ver[iVer][0]-x)in) <  )
//	}
//	
//}



int NozzleWallProjection_DV (Options *mshopt, Mesh *Msh, CadNozzle * CadNoz, CadNozzle *CadNoz_bas,  int refUp, int refDown, char *OutNam, int verbose)
{
	
	//---
	
	int WrtMsh = 1, patch=-1;
		
	int iTri, ref=0, i, j, iVer, vid=-1, is[3];
	
	double CrdNew[3];
	
	int NbrVerPrj=0 , NbrTriPrj=0;
	
	int *Tag = (int*) malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag, 0, sizeof(int)*(Msh->NbrVer+1));
	
	//--- Tag one ref after the other for consistent line re-projection
	
	NbrTriPrj = 0;
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refUp  )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	}	
	
	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		ref = Msh->Tri[iTri][3];
		if ( ref != refDown )
			continue;
		for (j=0; j<3; j++) {
			vid = Msh->Tri[iTri][j];
			Tag[vid] = ref;
		}
		NbrTriPrj++;
	}
	
	//--- Create mesh for visu
	
	NbrVerPrj = 0;
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		if ( Tag[iVer] > 0 )
			NbrVerPrj++;
	}
	
	int SizMsh[GmfMaxKwd+1];
	
	for (i=0; i<GmfMaxKwd; i++)
		SizMsh[i] = 0;
	
	SizMsh[GmfDimension]  = 3;
	SizMsh[GmfVertices]   = NbrVerPrj;
	SizMsh[GmfTriangles]  = NbrTriPrj;
	
	//printf("NbrVerPrj %d NbrRTriPrj %d\n", NbrVerPrj, NbrTriPrj);
	
	Mesh *MshViz = NULL;
	
	if ( WrtMsh == 1 ) {
		MshViz = AllocMesh(SizMsh);
		MshViz->Dim = 3;		
	}	
		
	FILE *motionHdl = fopen(OutNam, "wb");
	
	if ( motionHdl ) {
		if ( verbose > 0 )
			printf("%%%% %s OPENED (WRITE).\n", OutNam);
	}
	else
	{
		if ( verbose > 0  )
			printf("  ## ERROR : Could not open %s.\n", OutNam);
		exit(1);
	}
	
	//--- Project points
		
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
		if ( Tag[iVer] <= 0 )
			continue;
		
		CrdNew[0] = Msh->Ver[iVer][0];
		CrdNew[1] = Msh->Ver[iVer][1];
		CrdNew[2] = Msh->Ver[iVer][2];
		
		patch = -1;
		
		//if ( Tag[iVer] == refUp ) {
		//	//---
		//	ProjectNozzleWall_Up_DV (Msh->Ver[iVer], CrdNew, CadNoz_bas);			
		//	patch = NOZZLEUP;
		//}
		//else if ( Tag[iVer] == refDown ) {			
		//	//---
		//	ProjectNozzleWall_Down_DV (Msh->Ver[iVer], CrdNew, CadNoz_bas);
		//	patch = NOZZLEDOWN;
		//}
		//
		//if ( patch == -1 ) {
		//	if ( verbose > 0 )
		//		printf("  ## ERROR : Wrong nozzle CAD patch.\n");
		//	exit(1);
		//}
		//
		//ProjectToDV_DV(CrdNew, CadNoz, CadNoz_bas, patch);

		if ( Tag[iVer] == refUp ) {
			//---
			ProjectNozzleWall_Up_DV (Msh->Ver[iVer], CrdNew, CadNoz_bas);	
			ProjectToDV_DV(CrdNew, CadNoz, CadNoz_bas, patch);		
			patch = NOZZLEUP;
		}
		else if ( Tag[iVer] == refDown ) {			
			//---
			ProjectNozzleWall_Down_DV_3 (Msh->Ver[iVer], CrdNew, CadNoz, CadNoz_bas);
			patch = NOZZLEDOWN;
		}
		
		if ( patch == -1 ) {
			if ( verbose > 0 )
				printf("  ## ERROR : Wrong nozzle CAD patch.\n");
			exit(1);
		}
		

		//CrdNew[0] = Msh->Ver[iVer][0];
		//CrdNew[1] = Msh->Ver[iVer][1];
		//CrdNew[2] = Msh->Ver[iVer][2];
		//	
		if ( WrtMsh == 1 ) {
			MshViz->NbrVer++;
			AddVertex(MshViz, MshViz->NbrVer, CrdNew);
			Tag[iVer] = MshViz->NbrVer;
		}
		
		fprintf(motionHdl, "%d %le %le %le\n", iVer-1, CrdNew[0], CrdNew[1], CrdNew[2]);
		
	}
	
	//--- Write visu mesh
	
	if ( WrtMsh == 1 ) {
		
		for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
			ref = Msh->Tri[iTri][3];
			if ( ref != refUp && ref != refDown )
				continue;
			
			for (j=0; j<3; j++) {
				vid = Msh->Tri[iTri][j];
				is[j] = Tag[vid];
			}
			
			MshViz->NbrTri++;
			AddTriangle(MshViz,MshViz->NbrTri,is,ref);
		}
		
		WriteGMFMesh("visu", MshViz, 1);
		
	}
	
	if (motionHdl)
		fclose(motionHdl);
	
	if ( Tag )
		free(Tag);
	
	if ( MshViz )
		FreeMesh(MshViz);
	
	if ( CadNoz )
		FreeCadNozzle (CadNoz);
	
	return 0;
	
}



