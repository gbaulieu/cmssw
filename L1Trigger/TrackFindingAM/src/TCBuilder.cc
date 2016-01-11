/**
   C++ implementation of the TC builder

   S.Viret, G.Galbit : 09/12/15
**/

#include "../interface/TCBuilder.h"

TCBuilder::TCBuilder():TrackFitter(0){

}

TCBuilder::TCBuilder(int nb):TrackFitter(nb)
{
  m_nMissingHits = 1;           //Maximum number of missing layers in a TC (from the number of layers in the pattern)
  m_nFractionnalPartWidth = 0;  //Binning
  updateThresholds();
}

TCBuilder::~TCBuilder()
{
}

void TCBuilder::initialize(){

}

void TCBuilder::setFractionnalPartWidth(int nbFloatingPoint){
  if(nbFloatingPoint>=0){
    m_nFractionnalPartWidth = nbFloatingPoint;
    updateThresholds();
  }
}

void TCBuilder::updateThresholds(){
  //If m_nFractionnalPartWidth != 10 bits, the floating point thresholds are used
  if (m_nFractionnalPartWidth == 10){
    /* 10 bits fractionnal part Thresholds */
  
    //Barrel
    addThresholds( 0,  1,  2, SEC_BARREL, 0.031250, 3.665039);
    addThresholds( 0,  1,  8, SEC_BARREL, 0.074219, 38.231445);
    addThresholds( 0,  1,  9, SEC_BARREL, 0.118164, 39.138672);
    addThresholds( 0,  1, 10, SEC_BARREL, 0.160156, 42.259766);
    addThresholds( 0,  2,  8, SEC_BARREL, 0.065430, 74.090820);
    addThresholds( 0,  2,  9, SEC_BARREL, 0.145508, 74.028320);
    addThresholds( 0,  2, 10, SEC_BARREL, 0.223633, 79.116211);
    addThresholds( 1,  2,  8, SEC_BARREL, 0.031250, 41.812500);
    addThresholds( 1,  2,  9, SEC_BARREL, 0.079102, 43.343750);
    addThresholds( 1,  2, 10, SEC_BARREL, 0.096680, 45.228516);

    //Hybrid
    addThresholds( 0,  1,  2, SEC_HYBRID, 0.038086, 4.089844);
    addThresholds( 0,  1,  8, SEC_HYBRID, 0.061523, 36.709961);
    addThresholds( 0,  1,  9, SEC_HYBRID, 0.125000, 38.314453);
    addThresholds( 0,  1, 10, SEC_HYBRID, 0.125000, 42.623047);
    addThresholds( 0,  1, 11, SEC_HYBRID, 0.134766, 62.708984);
    addThresholds( 0,  1, 12, SEC_HYBRID, 0.126953, 67.507812);
    addThresholds( 0,  1, 13, SEC_HYBRID, 0.172852, 70.373047);
    addThresholds( 0,  2,  8, SEC_HYBRID, 0.093750, 74.268555);
    addThresholds( 0,  2,  9, SEC_HYBRID, 0.129883, 76.775391);
    addThresholds( 0,  2, 10, SEC_HYBRID, 0.125977, 79.952148);
    addThresholds( 0,  2, 11, SEC_HYBRID, 0.180664, 144.421875);
    addThresholds( 0,  2, 12, SEC_HYBRID, 0.205078, 140.078125);
    addThresholds( 0,  2, 13, SEC_HYBRID, 0.201172, 146.367188);
    addThresholds( 1,  2,  8, SEC_HYBRID, 0.042969, 41.355469);
    addThresholds( 1,  2,  9, SEC_HYBRID, 0.084961, 39.357422);
    addThresholds( 1,  2, 10, SEC_HYBRID, 0.093750, 46.836914);
    addThresholds( 1,  2, 11, SEC_HYBRID, 0.106445, 87.046875);
    addThresholds( 1,  2, 12, SEC_HYBRID, 0.117188, 79.388672);
    addThresholds( 1,  2, 13, SEC_HYBRID, 0.097656, 81.233398);

    //Endcap
    addThresholds( 0,  1,  3, SEC_ENDCAP, 0.040039, 5.756836);
    addThresholds( 0,  1,  4, SEC_ENDCAP, 0.082031, 6.778320);
    addThresholds( 0,  1,  5, SEC_ENDCAP, 0.045898, 6.506836);
    addThresholds( 0,  1, 11, SEC_ENDCAP, 0.078125, 64.822266);
    addThresholds( 0,  1, 12, SEC_ENDCAP, 0.116211, 81.608398);
    addThresholds( 0,  1, 13, SEC_ENDCAP, 0.143555, 92.575195);
    addThresholds( 0,  1, 14, SEC_ENDCAP, 0.203125, 97.450195);
    addThresholds( 0,  1, 15, SEC_ENDCAP, 0.309570, 96.867188);
    addThresholds( 0,  3,  4, SEC_ENDCAP, 0.056641, 10.769531);
    addThresholds( 0,  3,  5, SEC_ENDCAP, 0.053711, 19.661133);
    addThresholds( 0,  3,  6, SEC_ENDCAP, 0.049805, 13.232422);
    addThresholds( 0,  3,  7, SEC_ENDCAP, 0.051758, 11.861328);
    addThresholds( 0,  3, 11, SEC_ENDCAP, 0.148438, 190.636719);
    addThresholds( 0,  3, 12, SEC_ENDCAP, 0.167969, 195.087891);
    addThresholds( 0,  3, 13, SEC_ENDCAP, 0.219727, 193.606445);
    addThresholds( 0,  3, 14, SEC_ENDCAP, 0.246094, 208.910156);
    addThresholds( 0,  3, 15, SEC_ENDCAP, 0.405273, 210.363281);
    addThresholds( 0,  4,  5, SEC_ENDCAP, 0.040039, 18.907227);
    addThresholds( 0,  4,  6, SEC_ENDCAP, 0.052734, 13.366211);
    addThresholds( 0,  4,  7, SEC_ENDCAP, 0.050781, 12.603516);
    addThresholds( 0,  4, 12, SEC_ENDCAP, 0.228516, 220.043945);
    addThresholds( 0,  4, 13, SEC_ENDCAP, 0.137695, 227.261719);
    addThresholds( 0,  4, 14, SEC_ENDCAP, 0.236328, 229.281250);
    addThresholds( 0,  4, 15, SEC_ENDCAP, 0.365234, 200.322266);
    addThresholds( 1,  3,  4, SEC_ENDCAP, 0.019531, 7.486328);
    addThresholds( 1,  3,  5, SEC_ENDCAP, 0.020508, 6.088867);
    addThresholds( 1,  3, 11, SEC_ENDCAP, 0.102539, 128.959961);
    addThresholds( 1,  3, 12, SEC_ENDCAP, 0.096680, 118.861328);
    addThresholds( 1,  3, 13, SEC_ENDCAP, 0.118164, 122.890625);
    addThresholds( 1,  3, 14, SEC_ENDCAP, 0.106445, 130.061523);
    addThresholds( 1,  3, 15, SEC_ENDCAP, 0.131836, 124.733398);
    addThresholds( 1,  4,  5, SEC_ENDCAP, 0.014648, 6.622070);
    addThresholds( 1,  4, 12, SEC_ENDCAP, 0.107422, 150.178711);
    addThresholds( 1,  4, 13, SEC_ENDCAP, 0.096680, 141.776367);
    addThresholds( 1,  4, 14, SEC_ENDCAP, 0.148438, 156.682617);
    addThresholds( 1,  4, 15, SEC_ENDCAP, 0.196289, 137.015625);
    addThresholds( 3,  4,  5, SEC_ENDCAP, 0.012695, 7.040039);
    addThresholds( 3,  4,  6, SEC_ENDCAP, 0.021484, 10.887695);
    addThresholds( 3,  4,  7, SEC_ENDCAP, 0.022461, 15.625977);
    addThresholds( 3,  4, 12, SEC_ENDCAP, 0.056641, 60.820312);
    addThresholds( 3,  4, 13, SEC_ENDCAP, 0.044922, 61.541016);
    addThresholds( 3,  4, 14, SEC_ENDCAP, 0.053711, 63.903320);
    addThresholds( 3,  4, 15, SEC_ENDCAP, 0.069336, 76.256836);

    /* End of 10 bits fractionnal part Thresholds */

  }
  else{
    /* Floating point Thresholds */

    //Barrel
    addThresholds( 0,  1,  2, SEC_BARREL, 0.025571, 3.656240);
    addThresholds( 0,  1,  8, SEC_BARREL, 0.075301, 38.258154);
    addThresholds( 0,  1,  9, SEC_BARREL, 0.124700, 39.175377);
    addThresholds( 0,  1, 10, SEC_BARREL, 0.155787, 42.231133);
    addThresholds( 0,  2,  8, SEC_BARREL, 0.056070, 74.100319);
    addThresholds( 0,  2,  9, SEC_BARREL, 0.132452, 74.038652);
    addThresholds( 0,  2, 10, SEC_BARREL, 0.210517, 79.126298);
    addThresholds( 1,  2,  8, SEC_BARREL, 0.017905, 41.821861);
    addThresholds( 1,  2,  9, SEC_BARREL, 0.059329, 43.340059);
    addThresholds( 1,  2, 10, SEC_BARREL, 0.089593, 45.167844);

    //Hybrid
    addThresholds( 0,  1,  2, SEC_HYBRID, 0.030079, 4.083949);
    addThresholds( 0,  1,  8, SEC_HYBRID, 0.059964, 36.661217);
    addThresholds( 0,  1,  9, SEC_HYBRID, 0.124758, 38.346169);
    addThresholds( 0,  1, 10, SEC_HYBRID, 0.129082, 42.516624);
    addThresholds( 0,  1, 11, SEC_HYBRID, 0.121676, 62.673522);
    addThresholds( 0,  1, 12, SEC_HYBRID, 0.122195, 67.448446);
    addThresholds( 0,  1, 13, SEC_HYBRID, 0.193056, 70.364759);
    addThresholds( 0,  2,  8, SEC_HYBRID, 0.069076, 74.257499);
    addThresholds( 0,  2,  9, SEC_HYBRID, 0.117036, 76.740289);
    addThresholds( 0,  2, 10, SEC_HYBRID, 0.134230, 79.913063);
    addThresholds( 0,  2, 11, SEC_HYBRID, 0.179327, 144.488991);
    addThresholds( 0,  2, 12, SEC_HYBRID, 0.194377, 140.078556);
    addThresholds( 0,  2, 13, SEC_HYBRID, 0.243792, 146.363807);
    addThresholds( 1,  2,  8, SEC_HYBRID, 0.036669, 41.356791);
    addThresholds( 1,  2,  9, SEC_HYBRID, 0.102266, 39.336151);
    addThresholds( 1,  2, 10, SEC_HYBRID, 0.046263, 46.853872);
    addThresholds( 1,  2, 11, SEC_HYBRID, 0.092028, 87.089641);
    addThresholds( 1,  2, 12, SEC_HYBRID, 0.096336, 79.450331);
    addThresholds( 1,  2, 13, SEC_HYBRID, 0.091826, 81.263897);

    //Endcap
    addThresholds( 0,  1,  3, SEC_ENDCAP, 0.049016, 5.789638);
    addThresholds( 0,  1,  4, SEC_ENDCAP, 0.098305, 6.791889);
    addThresholds( 0,  1,  5, SEC_ENDCAP, 0.025410, 6.441114);
    addThresholds( 0,  1, 11, SEC_ENDCAP, 0.075897, 64.845979);
    addThresholds( 0,  1, 12, SEC_ENDCAP, 0.094022, 81.627829);
    addThresholds( 0,  1, 13, SEC_ENDCAP, 0.149944, 92.453032);
    addThresholds( 0,  1, 14, SEC_ENDCAP, 0.213731, 97.360846);
    addThresholds( 0,  1, 15, SEC_ENDCAP, 0.334518, 96.698723);
    addThresholds( 0,  3,  4, SEC_ENDCAP, 0.056367, 10.747764);
    addThresholds( 0,  3,  5, SEC_ENDCAP, 0.048138, 19.651452);
    addThresholds( 0,  3,  6, SEC_ENDCAP, 0.032072, 13.317437);
    addThresholds( 0,  3,  7, SEC_ENDCAP, 0.037592, 11.832096);
    addThresholds( 0,  3, 11, SEC_ENDCAP, 0.157577, 190.592070);
    addThresholds( 0,  3, 12, SEC_ENDCAP, 0.163619, 195.157998);
    addThresholds( 0,  3, 13, SEC_ENDCAP, 0.196538, 193.688462);
    addThresholds( 0,  3, 14, SEC_ENDCAP, 0.226246, 208.811893);
    addThresholds( 0,  3, 15, SEC_ENDCAP, 0.368564, 210.417263);
    addThresholds( 0,  4,  5, SEC_ENDCAP, 0.033095, 18.907360);
    addThresholds( 0,  4,  6, SEC_ENDCAP, 0.030136, 13.398305);
    addThresholds( 0,  4,  7, SEC_ENDCAP, 0.044999, 12.633459);
    addThresholds( 0,  4, 12, SEC_ENDCAP, 0.195552, 220.013405);
    addThresholds( 0,  4, 13, SEC_ENDCAP, 0.134733, 227.209487);
    addThresholds( 0,  4, 14, SEC_ENDCAP, 0.231035, 229.285057);
    addThresholds( 0,  4, 15, SEC_ENDCAP, 0.342347, 200.509661);
    addThresholds( 1,  3,  4, SEC_ENDCAP, 0.011823, 7.458285);
    addThresholds( 1,  3,  5, SEC_ENDCAP, 0.007385, 6.078510);
    addThresholds( 1,  3, 11, SEC_ENDCAP, 0.109007, 128.931298);
    addThresholds( 1,  3, 12, SEC_ENDCAP, 0.103831, 118.932845);
    addThresholds( 1,  3, 13, SEC_ENDCAP, 0.112873, 122.867231);
    addThresholds( 1,  3, 14, SEC_ENDCAP, 0.124153, 130.008205);
    addThresholds( 1,  3, 15, SEC_ENDCAP, 0.101961, 124.808447);
    addThresholds( 1,  4,  5, SEC_ENDCAP, 0.007987, 6.613035);
    addThresholds( 1,  4, 12, SEC_ENDCAP, 0.107208, 150.096417);
    addThresholds( 1,  4, 13, SEC_ENDCAP, 0.093686, 141.718773);
    addThresholds( 1,  4, 14, SEC_ENDCAP, 0.158431, 156.732945);
    addThresholds( 1,  4, 15, SEC_ENDCAP, 0.151566, 137.073274);
    addThresholds( 3,  4,  5, SEC_ENDCAP, 0.007473, 7.082068);
    addThresholds( 3,  4,  6, SEC_ENDCAP, 0.010956, 10.968808);
    addThresholds( 3,  4,  7, SEC_ENDCAP, 0.021753, 15.601805);
    addThresholds( 3,  4, 12, SEC_ENDCAP, 0.048389, 60.788335);
    addThresholds( 3,  4, 13, SEC_ENDCAP, 0.036618, 61.534191);
    addThresholds( 3,  4, 14, SEC_ENDCAP, 0.045263, 63.905185);
    addThresholds( 3,  4, 15, SEC_ENDCAP, 0.050656, 76.211291);

    /* End of Floating point Thresholds */
  }
}

void TCBuilder::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}

// Tentative duplicate removal
// Not used for the moment

void TCBuilder::mergeTracks(){
}

/* From the hits of the best TC, process the track parameters, create and fill a Track object and return its pointer */
Track* TCBuilder::createFittedTrack(vector <Hit*> &bestTC)
{
  int lay;
  bool barrelmod;

  int size = bestTC.size();
    
  float r;
  float rmin,rmax;
  
  float phi_est,z_est,eta_est,pt_est;
  
  double parZX[2][2];
  double resZX[2];
  double invZX[2][2];
  double detZX = 0;
     
  double x,y,z;
  
  int npt=0;
  
  float x1,x2,y1,y2;
  float x0,y0;
  
  for (int i=0;i<2;++i)
    {
      resZX[i] = 0.;
      for (int jj=0;jj<2;++jj) parZX[i][jj] = 0.;
      for (int jj=0;jj<2;++jj) invZX[i][jj] = 0.;
    }
    
  rmin = 1000;
  rmax = 0;
  int imax=-1;
  int imin=-1;
  
  x0=0;
  y0=0;
  
  x2=0;
  y2=0;

  int n2Sb=0;
  int nPSb=0;
  int n2Se=0;
  int nPSe=0;

  // Loop over stubs in the TC
  // In order to determine the lowest/largest radius of the 
  // TC, and the TC composition
  //
  // The stub with the largest radius is our reference in 
  // the conformal space


  for (int kk=0;kk<size;++kk) 
    {
      lay  = bestTC.at(kk)->getLayer();
    
      if (lay<=7)
	{
	  (lay<=2)
	    ? nPSb++
	    : nPSe++;
	}
      else
	{
	  (lay<11)
	    ? n2Sb++
	    : n2Se++;
	}
            
      if (lay>10) continue; // 2S disks stubs not used
      ++npt;
        
      x = bestTC.at(kk)->getX();
      y = bestTC.at(kk)->getY();
      r = sqrt(x*x+y*y);
    
      if (r<rmin)
	{
	  rmin = r;
	  imin = kk;
	  x2   = x;
	  y2   = y;
	}
        
      if (r>rmax)
	{
	  rmax = r;
	  imax = kk;
	  x0   = x;
	  y0   = y;
	}
    }

  float rmax2 = 0;
    
  x1 = 0;
  y1 = 0;
  
  int nc=0;
    
  float xtemp1,ytemp1;
  float xtemp2,ytemp2;
  
  xtemp1=0.;
  ytemp1=0.;

  // Loop 2 over stubs in the TC
  // In order to determine the point with the second largest radius 

  for (int kk=0;kk<size;++kk) // Loop over stubs in the TC
    {
      if (kk==imax || kk==imin) continue; // Already used

      lay  = bestTC.at(kk)->getLayer();

      barrelmod=0;
      if (lay<=2 || (lay>=8 && lay<=10)) barrelmod=1;
      if (!barrelmod && (nPSb+n2Sb)>=3) continue; // Barrel modules have a priority
      if (lay>10 && (nPSb+nPSe)>=3) continue;     // Don't use 2S modules in the disks if possible
    
      x = bestTC.at(kk)->getX();
      y = bestTC.at(kk)->getY();
      r = sqrt(x*x+y*y);
    
      if (r>rmax2)
	{
	  rmax2  = r;
	  x1     = x;
	  y1     = y;
	  nc++;
	}
    }

  // Now get the coordinates in the conformal space.

  xtemp1 = (x1-x0)/((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  ytemp1 = (y1-y0)/((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  
  xtemp2 = (x2-x0)/((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
  ytemp2 = (y2-y0)/((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
  
  x1=xtemp1;
  y1=ytemp1;
  x2=xtemp2;
  y2=ytemp2;
  

  // Now we got everything for the r/phi plane   

  double a = x0-0.5*(y1-y2)/(y2*x1-y1*x2);
  double b = y0-0.5*(x1-x2)/(y1*x2-y2*x1);
    
  int charge =-b/fabs(b);
   
  phi_est = atan2(charge*a,-charge*b);

  //Apply the rotation
  phi_est += sec_phi;

  //Set the value betweend -Pi and Pi
  phi_est = fmod(phi_est + M_PI, 2 * M_PI) - M_PI;

  pt_est  = 0.003*3.833*sqrt(a*a+b*b);

  // Then we do the RZ fit (LS)

  float wght;
  int cnt=0;
  for (int kk=0;kk<size;++kk) // Loop over stubs in the TC
    {
      lay  = bestTC.at(kk)->getLayer();

      if (lay>7) continue; // Don't use 2S modules
      if (lay>2 && nPSb>=2) continue; // Don't use PS modules of the disks if number of good point in the barrel is sufficient        

      ++cnt;
      x = bestTC.at(kk)->getX();
      y = bestTC.at(kk)->getY();
      z = bestTC.at(kk)->getZ();        
      r = sqrt(x*x+y*y);
            
      wght=1;
      if (lay>2) wght=1/7.;

      parZX[0][0] += wght*r*r;
      parZX[1][1] += wght*1;
      parZX[1][0] += wght*r;
    
      resZX[0] += wght*r*z;
      resZX[1] += wght*z;       
    } // End of stub loop
    

  detZX = parZX[0][0]*parZX[1][1]-parZX[1][0]*parZX[1][0];

  invZX[0][0] =  parZX[1][1]/detZX;
  invZX[1][0] = -parZX[1][0]/detZX;
  invZX[1][1] =  parZX[0][0]/detZX;

  // Finally estimate track parameters in the R/Z plane 

  eta_est = std::asinh((invZX[0][0]*resZX[0] + invZX[1][0]*resZX[1]));
  z_est   = invZX[1][0]*resZX[0] + invZX[1][1]*resZX[1];

  Track* fit_track = new Track();
  fit_track->setCurve(pt_est);
  fit_track->setPhi0(phi_est);
  fit_track->setEta0(eta_est);
  fit_track->setZ0(z_est);
      
  for(unsigned int hitIndex=0;hitIndex < bestTC.size();hitIndex++)
    {
      fit_track->addStubIndex(bestTC[hitIndex]->getID());
    }

  return fit_track;
}

/* Function which apply a binning on a float number to simulate the Hardware fixed point representation */
double TCBuilder::binning(double fNumber, int nFractionnalPartWidth)
{
  double fDivRes;
  double fBinnedNumber;

  if (nFractionnalPartWidth == 0)
    {
      //If this parameter is set to zero, no binning is applied
      fBinnedNumber = fNumber;
    }
  else
    {
      //The number is divided by the power of 2 corresponding to the fractionnal part width
      fDivRes = fNumber / pow(2 , -nFractionnalPartWidth);

      //The result is rounded to the nearest integer and then multiplied by the power of 2
      fBinnedNumber = round(fDivRes) * pow(2 , -nFractionnalPartWidth);
    }

  return fBinnedNumber;
}

/* Process the alignment scores (on RPHI and on RZ plan) between the 2 seeds and an other stub */
void TCBuilder::alignScore(Hit& hSeed1, Hit& hSeed2, Hit& hTestStub, double tScores[])
{
  double fRPHI_Score , fRZ_Score;

  double X1, Y1, Z1, R1, PHI1;
  double X2, Y2, Z2, R2, PHI2;
  double X3, Y3, Z3, R3, PHI3;

  double RPHI_S1, RPHI_S2, RZ_S1, RZ_S2;

  //Coordinates X, Y, Z are already binned
  X1 = hSeed1.getX();
  Y1 = hSeed1.getY();
  Z1 = hSeed1.getZ();
	
  X2 = hSeed2.getX();
  Y2 = hSeed2.getY();
  Z2 = hSeed2.getZ();

  X3 = hTestStub.getX();
  Y3 = hTestStub.getY();
  Z3 = hTestStub.getZ();
	
  R1 = binning(sqrt(X1*X1 + Y1*Y1), m_nFractionnalPartWidth);
  R2 = binning(sqrt(X2*X2 + Y2*Y2), m_nFractionnalPartWidth);
  R3 = binning(sqrt(X3*X3 + Y3*Y3), m_nFractionnalPartWidth);

  //RPHI plan
  PHI1 = binning(atan(Y1/X1), m_nFractionnalPartWidth);
  PHI2 = binning(atan(Y2/X2), m_nFractionnalPartWidth);
  PHI3 = binning(atan(Y3/X3), m_nFractionnalPartWidth);

  RPHI_S1 = binning((PHI2 - PHI1) * (R3 - R2), m_nFractionnalPartWidth);
  RPHI_S2 = binning((PHI2 - PHI3) * (R2 - R1), m_nFractionnalPartWidth);

  fRPHI_Score = binning(RPHI_S1 + RPHI_S2, m_nFractionnalPartWidth);

  //RZ plan
  RZ_S1 = binning((Z2 - Z1) * (R3 - R2), m_nFractionnalPartWidth);
  RZ_S2 = binning((Z2 - Z3) * (R2 - R1), m_nFractionnalPartWidth);

  fRZ_Score = binning(RZ_S1 + RZ_S2, m_nFractionnalPartWidth);

  tScores[0] = fRPHI_Score;
  tScores[1] = fRZ_Score;
}


/* Fill thresholds data with new thresholds */
void TCBuilder::addThresholds(int nLaySeed1, int nLaySeed2, int nLayTestStub, SEC_TYPE secType, double fRPHI_Thresh, double fRZ_Thresh)
{

  //Fill the tab corresponding to the sector type with the pair of thresholds
  switch (secType)
    {
    case SEC_BARREL : m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0] = fRPHI_Thresh;
      m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1] = fRZ_Thresh;
      break;
    case SEC_HYBRID : m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0] = fRPHI_Thresh;
      m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1] = fRZ_Thresh;
      break;
    case SEC_ENDCAP : m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0] = fRPHI_Thresh;
      m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1] = fRZ_Thresh;
    }
}

/* Get the thresholds corresponding to the 3 layerID for a given sector type */
void TCBuilder::getThresholds(int nLaySeed1, int nLaySeed2, int nLayTestStub, SEC_TYPE secType, double tabThresh[])
{
  switch (secType)
    {
    case SEC_BARREL : tabThresh[0] = m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0];
      tabThresh[1] = m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1];
      break;
    case SEC_HYBRID : tabThresh[0] = m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0];
      tabThresh[1] = m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1];
      break;
    case SEC_ENDCAP : tabThresh[0] = m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0];
      tabThresh[1] = m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1];
    }
}


/* Operate the layerID transcoding (to get the layerID 4bits Hardawre representation) */
char TCBuilder::transcodeLayer(Hit * pHit)
{
  int nOrigLayer = pHit->getLayer();

  //layer transcoding of the disks is based on the radius, it can be optimized by using the ladder_ID
  double X = pHit->getX();
  double Y = pHit->getY();

  int nTransLayer;
		
  if (nOrigLayer <= 10)
    {
      //If the stub is in the barrel

      if (nOrigLayer <= 7)
	{    
	  //If layer 5, 6, 7
	  nTransLayer = nOrigLayer - 5;     //5->0, 6->1, ...
	}    
      else
	{
	  //If layer 8, 9, 10
	  nTransLayer = nOrigLayer;         //no change
	}
    }
  else if (sqrt(X*X + Y*Y) >= 62)
    {
      //If the stub is on an outer ring of a disk !!! (2S modules)

      if (nOrigLayer <= 15)
	{
	  //If layer 11, 12, 13, 14, 15
	  nTransLayer = nOrigLayer;         //no change
	}
      else
	{
	  //If layer 18, 19, 20, 21, 22
	  nTransLayer = nOrigLayer - 7;     //18->11, 19->12, ...
	}
    }
  else
    {
      //If the stub is on an inner ring of a disk !!! (PS modules)

      if (nOrigLayer <= 15)
	{
	  //If layer 11, 12, 13, 14, 15
	  nTransLayer = nOrigLayer - 8;     //11->3, 12->4, ...
	}
      else
	{
	  //If layer 18, 19, 20, 21, 22
	  nTransLayer = nOrigLayer - 15;    //18->3, 19->4, ...
	}
    }

  return nTransLayer;	
}

// TC builder module
/* Take as input the list of stubs contained in a matched road */

void TCBuilder::fit(vector<Hit*> originalHits)
{

  //cout<<"trying to fit "<<originalHits.size()<<" points"<<endl;

  int tow = sector_id; // The tower ID, necessary to get the phi shift

  SEC_TYPE currentSec;

  //Get the sec_type from the sector_id
  if (tow >= 16 && tow <= 31)
    currentSec = SEC_BARREL; 
  else if (tow >=8 && tow <=39)
    currentSec = SEC_HYBRID;
  else
    currentSec = SEC_ENDCAP;

  //Process the starting phi of the tower
  sec_phi = (tow%8) * M_PI / 4.0 - 0.4;

  //cos and sin values for a rotation of an angle -sec_phi
  double ci = cos(-sec_phi);
  double si = sin(-sec_phi);

  double rotatedX, rotatedY;

  //Create a new vector to store the custom hits parameters
  vector <Hit> hits;
  
  Hit* pOrigHit;
  
  //For each hit of the lists
  for (unsigned int origHitIndex = 0; origHitIndex<originalHits.size(); origHitIndex++)
    {
      pOrigHit = originalHits[origHitIndex];

      //Process the rotated coordinnates
      rotatedX = pOrigHit->getX() * ci - pOrigHit->getY() * si;
      rotatedY = pOrigHit->getX() * si + pOrigHit->getY() * ci;

      //Add the modified hit to the hits vector
      hits.push_back( Hit(transcodeLayer(pOrigHit),
			  pOrigHit->getLadder(),
			  pOrigHit->getModule(),
			  pOrigHit->getSegment(),
			  pOrigHit->getStripNumber(),
			  pOrigHit->getID(),
			  pOrigHit->getParticuleID(),
			  pOrigHit->getParticulePT(),
			  pOrigHit->getParticuleIP(),
			  pOrigHit->getParticuleETA(),
			  pOrigHit->getParticulePHI0(),
			  binning(rotatedX, m_nFractionnalPartWidth),
			  binning(rotatedY, m_nFractionnalPartWidth),
			  binning(pOrigHit->getZ(), m_nFractionnalPartWidth),
			  pOrigHit->getX0(),
			  pOrigHit->getY0(),
			  pOrigHit->getZ0(),
			  pOrigHit->getBend())
		      );

    }

  //Sort the hits by ascending order of layerID
  //(using a lambda definition of the sorting criteria which return a boolean)
  sort(hits.begin(), hits.end(), [ ]( const Hit& lhs, const Hit& rhs ) { return lhs.getLayer() < rhs.getLayer(); });


  int nLayersCurrentPattern = 0;
  int lastAddedLayer = -1;
  //Count the number of layers present in the pattern
  for (unsigned int hitIndex=0; hitIndex < hits.size(); hitIndex++)
    {
      if (lastAddedLayer != hits[hitIndex].getLayer())
	{
	  nLayersCurrentPattern++;
	  lastAddedLayer = hits[hitIndex].getLayer();
	}
    }


  vector <Hit*> vecCurrentCandidateHits;
  vector <double> vecCurrentCandidateScore;

  vector <Hit*> vecBestCandidateHits;
  double fBestCandidateScore = 0.0;

  for (unsigned int seed1Index=0; seed1Index<hits.size(); seed1Index++)
    {
      Hit& hSeed1   = hits[seed1Index];
      int nLaySeed1 = hSeed1.getLayer();

      if (nLaySeed1 == 2) continue; //layer 2 can't be the innermost seed stub
      if (nLaySeed1 > 3) break;     //no more possible combinations for this pattern

      //We have a correct Seed1


      for (unsigned int seed2Index = seed1Index+1; seed2Index<hits.size(); seed2Index++)
	{
	  Hit& hSeed2   = hits[seed2Index];
	  int nLaySeed2 = hSeed2.getLayer();

	  if (nLaySeed1 == nLaySeed2) continue; //The seed layers have to be differents
	  if (nLaySeed2 > 4) break;             //no more possible combinations for the current seed1


	  //We have a correct Seed1/Seed2 combination !!!

	  //Current candidate initialization (the 2 seeds)
	  vecCurrentCandidateHits.clear();
	  vecCurrentCandidateHits.push_back(&hSeed1);
	  vecCurrentCandidateHits.push_back(&hSeed2);

	  vecCurrentCandidateScore.clear();


	  for (unsigned int testStubIndex = seed2Index+1; testStubIndex<hits.size(); testStubIndex++)
	    {

	      Hit& hTestStub    = hits[testStubIndex];
	      int nLayTestStub  = hTestStub.getLayer();

	      if (nLayTestStub == nLaySeed2) continue; //The layers have to be differents

        
	      //Score processing of the Seed1/Seed2/testStub combination
	      double tabScore[2];
	      alignScore(hSeed1, hSeed2, hTestStub, tabScore);
        
	      //cout<< "Score RPHI = "<<tabScore[0]<<"   Score RZ = "<<tabScore[1]<<endl;

	      //Get the thresholds corresponding to the current layer combination
	      double tabThresh[2];
	      getThresholds(nLaySeed1, nLaySeed2, nLayTestStub, currentSec, tabThresh);

	      if (tabScore[0] <= tabThresh[0] && tabScore[1] <= tabThresh[1])
		{
		  //The stub is in the window defined by the seed projection (correct stub candidate !)
          
		  if (nLayTestStub != vecCurrentCandidateHits.back()->getLayer())
		    {
		      //The current testStub layer is not yet in the TC
		      vecCurrentCandidateHits.push_back(&hTestStub);
		      vecCurrentCandidateScore.push_back(tabScore[0]);
		    }
		  else if (tabScore[0] < vecCurrentCandidateScore.back())
		    {
		      //The layer is already in the TC but the Phi score of the current stub is better than the previous one
		      vecCurrentCandidateHits.back()   = &hTestStub;            
		      vecCurrentCandidateScore.back()  = tabScore[0];
		    }
		}
	    }

	  //If the current candidate own more than 6 stubs, the lasts (outtermost) are removed
	  while (vecCurrentCandidateHits.size() > 6)
	    {
	      vecCurrentCandidateHits.pop_back();
	      vecCurrentCandidateScore.pop_back();
	    }

	  //All the stubs have been tested for the current Seeds combination

	  if (int(vecCurrentCandidateHits.size()) >= nLayersCurrentPattern - m_nMissingHits)
	    {
	      //The current candidate has enough stubs to be a candidate
 
	      //Process the score of the track candidate
	      double fCurrentCandidateScore = 0.0;
	      while (vecCurrentCandidateScore.empty() == false)
		{
		  fCurrentCandidateScore += abs(vecCurrentCandidateScore.back());  //TODO norm?
		  vecCurrentCandidateScore.pop_back();
		}

	      if (vecCurrentCandidateHits.size() > vecBestCandidateHits.size() || (vecCurrentCandidateHits.size() == vecBestCandidateHits.size() && fCurrentCandidateScore < fBestCandidateScore))
		{
		  //The current candidate is better than the previous best one
		  vecBestCandidateHits = vecCurrentCandidateHits;
		  fBestCandidateScore = fCurrentCandidateScore;
		}
	    }
	}
    }

  //All the Seeds combinations have been tested

  if (vecBestCandidateHits.empty() == false)
    {
      //If there is a recorded best candidate

      //Fit the parameters and create the corresponding track object
      Track * fit_track;
      fit_track = createFittedTrack(vecBestCandidateHits);
    
      //cout<<"adding one track..."<<endl;
      tracks.push_back(fit_track);
    }
}

void TCBuilder::fit(){
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    fit(allHits);
  }
}

TrackFitter* TCBuilder::clone(){
  TCBuilder* fit = new TCBuilder(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

