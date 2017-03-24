/**
   C++ implementation of the TC builder

   S.Viret, G.Galbit : 09/12/15
**/

#include "../interface/TCBuilder.h"

TCBuilder::TCBuilder():TrackFitter(0){
  l2gConverter=NULL;
}

TCBuilder::TCBuilder(int nb):TrackFitter(nb)
{
  l2gConverter=NULL;
  m_nMissingHits = 1;
  m_minimum_number_for_TC = 0;                 //Minimum number of layers with stub to a create a TC (defined per pattern)
  maxseeds=-1; // By default all seeds are used.
  updateThresholds();
}

TCBuilder::~TCBuilder()
{
}

void TCBuilder::initialize(){

}

void TCBuilder::setLocalToGlobalConverter(LocalToGlobalConverter* l){
  l2gConverter = l;
}

void TCBuilder::setMaxSeeds(int nmax){
  maxseeds = nmax;
}

void TCBuilder::setPatternSize(int rs){
  if(rs>0)
    m_minimum_number_for_TC = rs-m_nMissingHits;
}

void TCBuilder::setHardwareEmulation(bool hardwareEmulation)
{
  CommonTools::hardwareSimulation = hardwareEmulation;
  updateThresholds();
}

void TCBuilder::updateThresholds(){

  if (CommonTools::hardwareSimulation){
    //// Hardware format LUT constants ////

    addThresholds( 0,  1,  3, SEC_ENDCAP, 0.00408554077148437500, 0.63134765625000000000);
    addThresholds( 0,  1,  4, SEC_ENDCAP, 0.00606918334960937500, 1.22338867187500000000);
    addThresholds( 0,  1,  5, SEC_ENDCAP, 0.00756835937500000000, 1.37548828125000000000);
    addThresholds( 0,  1, 11, SEC_ENDCAP, 0.00791549682617187500, 5.45458984375000000000);
    addThresholds( 0,  1, 12, SEC_ENDCAP, 0.01086425781250000000, 6.17846679687500000000);
    addThresholds( 0,  1, 13, SEC_ENDCAP, 0.01330947875976562500, 7.08740234375000000000);
    addThresholds( 0,  1, 14, SEC_ENDCAP, 0.01737594604492187500, 7.40600585937500000000);
    addThresholds( 0,  1, 15, SEC_ENDCAP, 0.01712799072265625000, 8.18334960937500000000);
    addThresholds( 0,  3,  4, SEC_ENDCAP, 0.00185775756835937500, 0.97607421875000000000);
    addThresholds( 0,  3,  5, SEC_ENDCAP, 0.00394439697265625000, 1.46850585937500000000);
    addThresholds( 0,  3,  6, SEC_ENDCAP, 0.00727081298828125000, 2.31469726562500000000);
    addThresholds( 0,  3,  7, SEC_ENDCAP, 0.01189422607421875000, 3.46313476562500000000);
    addThresholds( 0,  3, 11, SEC_ENDCAP, 0.00584793090820312500, 5.71972656250000000000);
    addThresholds( 0,  3, 12, SEC_ENDCAP, 0.00525283813476562500, 6.21728515625000000000);
    addThresholds( 0,  3, 13, SEC_ENDCAP, 0.00622558593750000000, 7.14184570312500000000);
    addThresholds( 0,  3, 14, SEC_ENDCAP, 0.00809860229492187500, 8.37329101562500000000);
    addThresholds( 0,  3, 15, SEC_ENDCAP, 0.01074600219726562500, 10.42944335937500000000);
    addThresholds( 0,  4,  5, SEC_ENDCAP, 0.00188446044921875000, 0.79711914062500000000);
    addThresholds( 0,  4,  6, SEC_ENDCAP, 0.00409317016601562500, 1.14843750000000000000);
    addThresholds( 0,  4,  7, SEC_ENDCAP, 0.00736618041992187500, 1.62402343750000000000);
    addThresholds( 0,  4, 12, SEC_ENDCAP, 0.00636672973632812500, 6.70434570312500000000);
    addThresholds( 0,  4, 13, SEC_ENDCAP, 0.00527191162109375000, 7.42749023437500000000);
    addThresholds( 0,  4, 14, SEC_ENDCAP, 0.00624465942382812500, 8.50512695312500000000);
    addThresholds( 0,  4, 15, SEC_ENDCAP, 0.00770568847656250000, 10.11230468750000000000);
    addThresholds( 1,  3,  4, SEC_ENDCAP, 0.00173950195312500000, 0.81323242187500000000);
    addThresholds( 1,  3,  5, SEC_ENDCAP, 0.00442504882812500000, 1.23291015625000000000);
    addThresholds( 1,  3, 11, SEC_ENDCAP, 0.00558471679687500000, 5.72119140625000000000);
    addThresholds( 1,  3, 12, SEC_ENDCAP, 0.00519180297851562500, 6.19433593750000000000);
    addThresholds( 1,  3, 13, SEC_ENDCAP, 0.00575637817382812500, 7.19091796875000000000);
    addThresholds( 1,  3, 14, SEC_ENDCAP, 0.00707244873046875000, 7.69384765625000000000);
    addThresholds( 1,  3, 15, SEC_ENDCAP, 0.01002883911132812500, 8.91235351562500000000);
    addThresholds( 1,  4,  5, SEC_ENDCAP, 0.00331115722656250000, 0.80493164062500000000);
    addThresholds( 1,  4, 12, SEC_ENDCAP, 0.00575256347656250000, 6.63671875000000000000);
    addThresholds( 1,  4, 13, SEC_ENDCAP, 0.00526046752929687500, 7.41308593750000000000);
    addThresholds( 1,  4, 14, SEC_ENDCAP, 0.00607299804687500000, 7.76953125000000000000);
    addThresholds( 1,  4, 15, SEC_ENDCAP, 0.00779342651367187500, 7.99926757812500000000);
    addThresholds( 3,  4,  5, SEC_ENDCAP, 0.00134277343750000000, 1.21142578125000000000);
    addThresholds( 3,  4,  6, SEC_ENDCAP, 0.00341415405273437500, 2.15625000000000000000);
    addThresholds( 3,  4,  7, SEC_ENDCAP, 0.00611495971679687500, 3.18896484375000000000);
    addThresholds( 3,  4, 12, SEC_ENDCAP, 0.00630187988281250000, 6.79858398437500000000);
    addThresholds( 3,  4, 13, SEC_ENDCAP, 0.00521850585937500000, 7.54370117187500000000);
    addThresholds( 3,  4, 14, SEC_ENDCAP, 0.00571441650390625000, 8.86938476562500000000);
    addThresholds( 3,  4, 15, SEC_ENDCAP, 0.00725173950195312500, 10.91772460937500000000);

    addThresholds( 0,  1,  2, SEC_BARREL, 0.00342941284179687500, 0.57543945312500000000);
    addThresholds( 0,  1,  8, SEC_BARREL, 0.00617980957031250000, 2.71679687500000000000);
    addThresholds( 0,  1,  9, SEC_BARREL, 0.01032638549804687500, 2.93017578125000000000);
    addThresholds( 0,  1, 10, SEC_BARREL, 0.01612091064453125000, 3.89453125000000000000);
    addThresholds( 0,  2,  8, SEC_BARREL, 0.00391769409179687500, 2.55566406250000000000);
    addThresholds( 0,  2,  9, SEC_BARREL, 0.00702285766601562500, 2.69750976562500000000);
    addThresholds( 0,  2, 10, SEC_BARREL, 0.01240158081054687500, 3.70800781250000000000);
    addThresholds( 1,  2,  8, SEC_BARREL, 0.00452804565429687500, 2.60107421875000000000);
    addThresholds( 1,  2,  9, SEC_BARREL, 0.00805282592773437500, 2.85424804687500000000);
    addThresholds( 1,  2, 10, SEC_BARREL, 0.01434326171875000000, 4.43383789062500000000);

    addThresholds( 0,  1,  2, SEC_HYBRID, 0.00379562377929687500, 0.51196289062500000000);
    addThresholds( 0,  1,  8, SEC_HYBRID, 0.00775527954101562500, 2.78686523437500000000);
    addThresholds( 0,  1,  9, SEC_HYBRID, 0.01086807250976562500, 3.05126953125000000000);
    addThresholds( 0,  1, 10, SEC_HYBRID, 0.01525878906250000000, 5.78320312500000000000);
    addThresholds( 0,  1, 11, SEC_HYBRID, 0.01112365722656250000, 5.12768554687500000000);
    addThresholds( 0,  1, 12, SEC_HYBRID, 0.01467132568359375000, 5.28295898437500000000);
    addThresholds( 0,  1, 13, SEC_HYBRID, 0.02006149291992187500, 5.51098632812500000000);
    addThresholds( 0,  2,  8, SEC_HYBRID, 0.00384140014648437500, 2.58593750000000000000);
    addThresholds( 0,  2,  9, SEC_HYBRID, 0.00765609741210937500, 2.80224609375000000000);
    addThresholds( 0,  2, 10, SEC_HYBRID, 0.01053619384765625000, 4.22973632812500000000);
    addThresholds( 0,  2, 11, SEC_HYBRID, 0.00677490234375000000, 4.82226562500000000000);
    addThresholds( 0,  2, 12, SEC_HYBRID, 0.00902557373046875000, 5.01293945312500000000);
    addThresholds( 0,  2, 13, SEC_HYBRID, 0.01129150390625000000, 5.31420898437500000000);
    addThresholds( 1,  2,  8, SEC_HYBRID, 0.00382232666015625000, 2.62939453125000000000);
    addThresholds( 1,  2,  9, SEC_HYBRID, 0.00746154785156250000, 2.88012695312500000000);
    addThresholds( 1,  2, 10, SEC_HYBRID, 0.01178359985351562500, 4.67285156250000000000);
    addThresholds( 1,  2, 11, SEC_HYBRID, 0.00666427612304687500, 4.91992187500000000000);
    addThresholds( 1,  2, 12, SEC_HYBRID, 0.00848007202148437500, 5.12963867187500000000);
    addThresholds( 1,  2, 13, SEC_HYBRID, 0.01080703735351562500, 5.34667968750000000000);

    //// Hardware format LUT constants ////

  }
  else{
    ///// Floating point Thresholds ////

    addThresholds( 0,  1,  3, SEC_ENDCAP, 0.00408554077148437500, 0.63134765625000000000);
    addThresholds( 0,  1,  4, SEC_ENDCAP, 0.00606918334960937500, 1.22338867187500000000);
    addThresholds( 0,  1,  5, SEC_ENDCAP, 0.00756835937500000000, 1.37548828125000000000);
    addThresholds( 0,  1, 11, SEC_ENDCAP, 0.00791549682617187500, 5.45458984375000000000);
    addThresholds( 0,  1, 12, SEC_ENDCAP, 0.01086425781250000000, 6.17846679687500000000);
    addThresholds( 0,  1, 13, SEC_ENDCAP, 0.01330947875976562500, 7.08740234375000000000);
    addThresholds( 0,  1, 14, SEC_ENDCAP, 0.01737594604492187500, 7.40600585937500000000);
    addThresholds( 0,  1, 15, SEC_ENDCAP, 0.01712799072265625000, 8.18334960937500000000);
    addThresholds( 0,  3,  4, SEC_ENDCAP, 0.00185775756835937500, 0.97607421875000000000);
    addThresholds( 0,  3,  5, SEC_ENDCAP, 0.00394439697265625000, 1.46850585937500000000);
    addThresholds( 0,  3,  6, SEC_ENDCAP, 0.00727081298828125000, 2.31469726562500000000);
    addThresholds( 0,  3,  7, SEC_ENDCAP, 0.01189422607421875000, 3.46313476562500000000);
    addThresholds( 0,  3, 11, SEC_ENDCAP, 0.00584793090820312500, 5.71972656250000000000);
    addThresholds( 0,  3, 12, SEC_ENDCAP, 0.00525283813476562500, 6.21728515625000000000);
    addThresholds( 0,  3, 13, SEC_ENDCAP, 0.00622558593750000000, 7.14184570312500000000);
    addThresholds( 0,  3, 14, SEC_ENDCAP, 0.00809860229492187500, 8.37329101562500000000);
    addThresholds( 0,  3, 15, SEC_ENDCAP, 0.01074600219726562500, 10.42944335937500000000);
    addThresholds( 0,  4,  5, SEC_ENDCAP, 0.00188446044921875000, 0.79711914062500000000);
    addThresholds( 0,  4,  6, SEC_ENDCAP, 0.00409317016601562500, 1.14843750000000000000);
    addThresholds( 0,  4,  7, SEC_ENDCAP, 0.00736618041992187500, 1.62402343750000000000);
    addThresholds( 0,  4, 12, SEC_ENDCAP, 0.00636672973632812500, 6.70434570312500000000);
    addThresholds( 0,  4, 13, SEC_ENDCAP, 0.00527191162109375000, 7.42749023437500000000);
    addThresholds( 0,  4, 14, SEC_ENDCAP, 0.00624465942382812500, 8.50512695312500000000);
    addThresholds( 0,  4, 15, SEC_ENDCAP, 0.00770568847656250000, 10.11230468750000000000);
    addThresholds( 1,  3,  4, SEC_ENDCAP, 0.00173950195312500000, 0.81323242187500000000);
    addThresholds( 1,  3,  5, SEC_ENDCAP, 0.00442504882812500000, 1.23291015625000000000);
    addThresholds( 1,  3, 11, SEC_ENDCAP, 0.00558471679687500000, 5.72119140625000000000);
    addThresholds( 1,  3, 12, SEC_ENDCAP, 0.00519180297851562500, 6.19433593750000000000);
    addThresholds( 1,  3, 13, SEC_ENDCAP, 0.00575637817382812500, 7.19091796875000000000);
    addThresholds( 1,  3, 14, SEC_ENDCAP, 0.00707244873046875000, 7.69384765625000000000);
    addThresholds( 1,  3, 15, SEC_ENDCAP, 0.01002883911132812500, 8.91235351562500000000);
    addThresholds( 1,  4,  5, SEC_ENDCAP, 0.00331115722656250000, 0.80493164062500000000);
    addThresholds( 1,  4, 12, SEC_ENDCAP, 0.00575256347656250000, 6.63671875000000000000);
    addThresholds( 1,  4, 13, SEC_ENDCAP, 0.00526046752929687500, 7.41308593750000000000);
    addThresholds( 1,  4, 14, SEC_ENDCAP, 0.00607299804687500000, 7.76953125000000000000);
    addThresholds( 1,  4, 15, SEC_ENDCAP, 0.00779342651367187500, 7.99926757812500000000);
    addThresholds( 3,  4,  5, SEC_ENDCAP, 0.00134277343750000000, 1.21142578125000000000);
    addThresholds( 3,  4,  6, SEC_ENDCAP, 0.00341415405273437500, 2.15625000000000000000);
    addThresholds( 3,  4,  7, SEC_ENDCAP, 0.00611495971679687500, 3.18896484375000000000);
    addThresholds( 3,  4, 12, SEC_ENDCAP, 0.00630187988281250000, 6.79858398437500000000);
    addThresholds( 3,  4, 13, SEC_ENDCAP, 0.00521850585937500000, 7.54370117187500000000);
    addThresholds( 3,  4, 14, SEC_ENDCAP, 0.00571441650390625000, 8.86938476562500000000);
    addThresholds( 3,  4, 15, SEC_ENDCAP, 0.00725173950195312500, 10.91772460937500000000);

    addThresholds( 0,  1,  2, SEC_BARREL, 0.00342941284179687500, 0.57543945312500000000);
    addThresholds( 0,  1,  8, SEC_BARREL, 0.00617980957031250000, 2.71679687500000000000);
    addThresholds( 0,  1,  9, SEC_BARREL, 0.01032638549804687500, 2.93017578125000000000);
    addThresholds( 0,  1, 10, SEC_BARREL, 0.01612091064453125000, 3.89453125000000000000);
    addThresholds( 0,  2,  8, SEC_BARREL, 0.00391769409179687500, 2.55566406250000000000);
    addThresholds( 0,  2,  9, SEC_BARREL, 0.00702285766601562500, 2.69750976562500000000);
    addThresholds( 0,  2, 10, SEC_BARREL, 0.01240158081054687500, 3.70800781250000000000);
    addThresholds( 1,  2,  8, SEC_BARREL, 0.00452804565429687500, 2.60107421875000000000);
    addThresholds( 1,  2,  9, SEC_BARREL, 0.00805282592773437500, 2.85424804687500000000);
    addThresholds( 1,  2, 10, SEC_BARREL, 0.01434326171875000000, 4.43383789062500000000);

    addThresholds( 0,  1,  2, SEC_HYBRID, 0.00379562377929687500, 0.51196289062500000000);
    addThresholds( 0,  1,  8, SEC_HYBRID, 0.00775527954101562500, 2.78686523437500000000);
    addThresholds( 0,  1,  9, SEC_HYBRID, 0.01086807250976562500, 3.05126953125000000000);
    addThresholds( 0,  1, 10, SEC_HYBRID, 0.01525878906250000000, 5.78320312500000000000);
    addThresholds( 0,  1, 11, SEC_HYBRID, 0.01112365722656250000, 5.12768554687500000000);
    addThresholds( 0,  1, 12, SEC_HYBRID, 0.01467132568359375000, 5.28295898437500000000);
    addThresholds( 0,  1, 13, SEC_HYBRID, 0.02006149291992187500, 5.51098632812500000000);
    addThresholds( 0,  2,  8, SEC_HYBRID, 0.00384140014648437500, 2.58593750000000000000);
    addThresholds( 0,  2,  9, SEC_HYBRID, 0.00765609741210937500, 2.80224609375000000000);
    addThresholds( 0,  2, 10, SEC_HYBRID, 0.01053619384765625000, 4.22973632812500000000);
    addThresholds( 0,  2, 11, SEC_HYBRID, 0.00677490234375000000, 4.82226562500000000000);
    addThresholds( 0,  2, 12, SEC_HYBRID, 0.00902557373046875000, 5.01293945312500000000);
    addThresholds( 0,  2, 13, SEC_HYBRID, 0.01129150390625000000, 5.31420898437500000000);
    addThresholds( 1,  2,  8, SEC_HYBRID, 0.00382232666015625000, 2.62939453125000000000);
    addThresholds( 1,  2,  9, SEC_HYBRID, 0.00746154785156250000, 2.88012695312500000000);
    addThresholds( 1,  2, 10, SEC_HYBRID, 0.01178359985351562500, 4.67285156250000000000);
    addThresholds( 1,  2, 11, SEC_HYBRID, 0.00666427612304687500, 4.91992187500000000000);
    addThresholds( 1,  2, 12, SEC_HYBRID, 0.00848007202148437500, 5.12963867187500000000);
    addThresholds( 1,  2, 13, SEC_HYBRID, 0.01080703735351562500, 5.34667968750000000000);
   
    //// End of Floating point Thresholds ////
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
      for (int j=0;j<2;++j) parZX[i][j] = 0.;
      for (int j=0;j<2;++j) invZX[i][j] = 0.;
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


  for (int i=0;i<size;++i) 
    {
      lay  = bestTC.at(i)->getLayer();

      if (lay<=7) //If PS
        {
          (lay<=2) ? nPSb++ : nPSe++;
        }
      else  //else 2S
        {
          (lay<11) ? n2Sb++ : n2Se++;
        }

      if (lay>10) continue; // 2S disks stubs not used
      ++npt;

      x = bestTC.at(i)->getX();
      y = bestTC.at(i)->getY();
      r = sqrt(x*x+y*y);

      if (r<rmin)
        {
          rmin = r;
          imin = i;
          x2   = x;
          y2   = y;
        }

      if (r>rmax)
        {
          rmax = r;
          imax = i;
          x0   = x;
          y0   = y;
        }
    }

  float rmax2 = 0;

  x1 = 0;
  y1 = 0;

  int nc=0;

  // Loop 2 over stubs in the TC
  // In order to determine the point with the second largest radius 

  for (int i=0;i<size;++i) // Loop over stubs in the TC
    {
      if (i==imax || i==imin) continue; // Already used

      lay  = bestTC.at(i)->getLayer();

      barrelmod=0;
      if (lay<=2 || (lay>=8 && lay<=10)) barrelmod=1;
      if (!barrelmod && (nPSb+n2Sb)>=3) continue; // Barrel modules have a priority
      if (lay>10 && (nPSb+nPSe)>=3) continue;     // Don't use 2S modules in the disks if possible

      x = bestTC.at(i)->getX();
      y = bestTC.at(i)->getY();
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

  double sqR1  = CommonTools::binning((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0), 13, 18, UNSIGNED);
  double sqR2  = CommonTools::binning((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0), 13, 18, UNSIGNED);

  x1 = CommonTools::binning((x1-x0)/sqR1, -3, 18, SIGNED);
  y1 = CommonTools::binning((y1-y0)/sqR1, -3, 18, SIGNED);

  x2 = CommonTools::binning((x2-x0)/sqR2, -3, 18, SIGNED);
  y2 = CommonTools::binning((y2-y0)/sqR2, -3, 18, SIGNED);


  double mult_x1_y2 = CommonTools::binning(x1*y2, -8, 18, SIGNED);
  double mult_x2_y1 = CommonTools::binning(x2*y1, -8, 18, SIGNED);

  // Now we got everything for the r/phi plane   

  double sub_of_mult = CommonTools::binning(mult_x1_y2 - mult_x2_y1, -11, 18, SIGNED);

  double divA = CommonTools::binning((y1-y2)/sub_of_mult, 17, 18, SIGNED);
  double divB = CommonTools::binning((x2-x1)/sub_of_mult, 17, 18, SIGNED);

  double a = CommonTools::binning(x0-0.5*divA, 17, 18, SIGNED);
  double b = CommonTools::binning(y0-0.5*divB, 17, 18, SIGNED);

  double a2_b2 = CommonTools::binning(a*a + b*b, 29, 18, UNSIGNED);

  pt_est  = 0.003*3.833*sqrt(a2_b2);

  int charge =-b/fabs(b);

  phi_est = atan2(charge*a,-charge*b);

  //Apply the rotation
  phi_est += sec_phi;

  //Set the value betweend -Pi and Pi
  phi_est = fmod(phi_est + M_PI, 2 * M_PI) - M_PI;

  // Then we do the RZ fit (LS)

  float wght;
  int cnt=0;
  for (int i=0;i<size;++i) // Loop over stubs in the TC
    {
      lay  = bestTC.at(i)->getLayer();

      if (lay>7) continue; // Don't use 2S modules
      if (lay>2 && nPSb>=2) continue; // Don't use PS modules of the disks if number of good point in the barrel is sufficient        

      ++cnt;
      x = bestTC.at(i)->getX();
      y = bestTC.at(i)->getY();
      z = bestTC.at(i)->getZ();        
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
  fit_track->setCharge(charge);
      
  for(unsigned int hitIndex=0;hitIndex < bestTC.size();hitIndex++)
    {
      fit_track->addStubIndex(bestTC[hitIndex]->getID());
    }

  return fit_track;
}

/* Process the alignment scores (on RPHI and on RZ plan) between the 2 seeds and an other stub */
void TCBuilder::alignScore(Hit& hSeed1, Hit& hSeed2, Hit& hTestStub, double tScores[])
{
  double fRPHI_Score , fRZ_Score;

  double X1, Y1, Z1, R1, PHI1;
  double X2, Y2, Z2, R2, PHI2;
  double X3, Y3, Z3, R3, PHI3;

  double result_R, result_PHI;

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
	
  /*

  //Old way

  R1 = CommonTools::binning(sqrt(X1*X1 + Y1*Y1), 6, 18, SIGNED);
  R2 = CommonTools::binning(sqrt(X2*X2 + Y2*Y2), 6, 18, SIGNED);
  R3 = CommonTools::binning(sqrt(X3*X3 + Y3*Y3), 6, 18, SIGNED);

  //RPHI plan
  PHI1 = CommonTools::binning(atan(Y1/X1), 0, 17, SIGNED);
  PHI2 = CommonTools::binning(atan(Y2/X2), 0, 17, SIGNED);
  PHI3 = CommonTools::binning(atan(Y3/X3), 0, 17, SIGNED);
  */

  //New bitwise way to get the polar coordinates
  CommonTools::binCordic(X1, Y1, result_R, result_PHI);
  R1 = result_R;
  PHI1 = result_PHI;
  
  CommonTools::binCordic(X2, Y2, result_R, result_PHI);
  R2 = result_R;
  PHI2 = result_PHI;

  CommonTools::binCordic(X3, Y3, result_R, result_PHI);
  R3 = result_R;
  PHI3 = result_PHI;

  //cout<<"Polar "<<hSeed1.getID()<<" : "<<R1<<"/"<<PHI1<<"/"<<Z1<<endl;
  //cout<<"Polar "<<hSeed2.getID()<<" : "<<R2<<"/"<<PHI2<<"/"<<Z2<<endl;
  //cout<<"Polar "<<hTestStub.getID()<<" : "<<R3<<"/"<<PHI3<<"/"<<Z3<<endl;

  RPHI_S1 = CommonTools::binning((PHI2 - PHI1) * (R3 - R2), 8, 20, SIGNED);
  RPHI_S2 = CommonTools::binning((PHI2 - PHI3) * (R2 - R1), 8, 20, SIGNED);

  fRPHI_Score = CommonTools::binning(fabs(RPHI_S1 + RPHI_S2), 7, 18, UNSIGNED);

  //RZ plan
  RZ_S1 = CommonTools::binning((Z2 - Z1) * (R3 - R2), 12, 20, SIGNED);
  RZ_S2 = CommonTools::binning((Z2 - Z3) * (R2 - R1), 12, 20, SIGNED);

  fRZ_Score = CommonTools::binning(fabs(RZ_S1 + RZ_S2), 10, 18, UNSIGNED);

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

void TCBuilder::fit(vector<Hit*> originalHits, int pattern_id)
{

  //cout<<"trying to fit "<<originalHits.size()<<" points"<<endl;

  int tow = sector_id; // The tower ID, necessary to get the phi shift

  int nseeds=0;

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
      
      /**************** FROM LOCAL TO GLOBAL COORDINATES ****************/
      vector<float> coords;
      if(l2gConverter!=NULL){
	coords = l2gConverter->toGlobal(pOrigHit);
	rotatedX = coords[0];
	rotatedY = coords[1];
      }
      else{
	// If we do not have a converter, use the coordinates from CMSSW
	coords.push_back(pOrigHit->getX());
	coords.push_back(pOrigHit->getY());
	coords.push_back(pOrigHit->getZ());
	//Process the rotated coordinnates
	rotatedX = coords[0] * ci - coords[1] * si;
	rotatedY = coords[0] * si + coords[1] * ci;
      }
      /*****************************************************************/

      //cout<<"Cartesian "<<pOrigHit->getID()<<" : "<<CommonTools::binning(rotatedX, 6, 18, SIGNED)<<"/"<<CommonTools::binning(rotatedY, 6, 18, SIGNED)<<"/"<<CommonTools::binning((double)coords[2], 8, 18, SIGNED)<<endl;

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
			  CommonTools::binning(rotatedX, 6, 18, SIGNED),
			  CommonTools::binning(rotatedY, 6, 18, SIGNED),
			  CommonTools::binning((double)coords[2], 8, 18, SIGNED),
			  pOrigHit->getX0(),
			  pOrigHit->getY0(),
			  pOrigHit->getZ0(),
			  pOrigHit->getBend())
		      );
    }

  //Sort the hits by ascending order of layerID
  //(using a lambda definition of the sorting criteria which return a boolean)
  sort(hits.begin(), hits.end(), [ ]( const Hit& lhs, const Hit& rhs ) { return lhs.getLayer() < rhs.getLayer(); });


  int lastAddedLayer = -1;
  //Count the number of layers present in the pattern
  for (unsigned int hitIndex=0; hitIndex < hits.size(); hitIndex++)
    {
      if (lastAddedLayer != hits[hitIndex].getLayer())
	{
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

      //Get the radius of the seed1
      double fRseed1 = CommonTools::binning(sqrt(hSeed1.getX()*hSeed1.getX() + hSeed1.getY()*hSeed1.getY()), 6, 18, SIGNED);

      for (unsigned int seed2Index = seed1Index+1; seed2Index<hits.size(); seed2Index++)
	{
	  Hit& hSeed2   = hits[seed2Index];
	  int nLaySeed2 = hSeed2.getLayer();

	  if (nLaySeed1 == nLaySeed2) continue; //The seed layers have to be differents
	  if (nLaySeed2 > 4) break;             //no more possible combinations for the current seed1

	  ++nseeds;
	  if (nseeds>maxseeds && maxseeds>0) break;

	  //We have a correct Seed1/Seed2 combination !!!

	  //Get the radius of the seed2
	  double fRseed2 = CommonTools::binning(sqrt(hSeed2.getX()*hSeed2.getX() + hSeed2.getY()*hSeed2.getY()), 6, 18, SIGNED);

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
        

	      //Get the thresholds corresponding to the current layer combination
	      double tabNormThresh[2];
	      getThresholds(nLaySeed1, nLaySeed2, nLayTestStub, currentSec, tabNormThresh);


	      //Process the real thresholds from the normalized one stored in the LUT
	      double fThreshRPHI = CommonTools::binning(fabs(tabNormThresh[0] * (fRseed2 - fRseed1)), 4, 18, SIGNED);
	      double fThreshRZ = CommonTools::binning(fabs(tabNormThresh[1] * (fRseed2 - fRseed1)), 8, 18, SIGNED);

	      if (tabScore[0] <= fThreshRPHI && tabScore[1] <= fThreshRZ)
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

	  //All the stubs have been tested for the current Seeds combination

	  if (vecCurrentCandidateHits.size() >= m_minimum_number_for_TC)
	    {
	      //The current candidate has enough stubs to be a candidate
 
	      //Process the score of the track candidate
	      double fCurrentCandidateScore = 0.0;
	      while (vecCurrentCandidateScore.empty() == false)
		{
		  fCurrentCandidateScore += vecCurrentCandidateScore.back();
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

  if ( vecBestCandidateHits.size() >= m_minimum_number_for_TC )
    {
      //If there is a recorded best candidate

      //If the candidate owns more than 6 stubs, the outtermost ones are removed
      while (vecBestCandidateHits.size() > 6){
	vecBestCandidateHits.pop_back();
      }

      // If we have a 6 stubs TC with a stub on the last endcap disk -> we remove it
      // The 5 stubs TC will be easier to handle for the PCA
      if(vecBestCandidateHits.size()==6){
	if(vecBestCandidateHits[vecBestCandidateHits.size()-1]->getLayer()==7 || vecBestCandidateHits[vecBestCandidateHits.size()-1]->getLayer()==15){
	  vecBestCandidateHits.pop_back();
	}
      }

      //Fit the parameters and create the corresponding track object
      Track * fit_track;
      fit_track = createFittedTrack(vecBestCandidateHits);
      fit_track->setOriginPatternID(pattern_id);
    
      //cout<<"adding one track..."<<endl;
      tracks.push_back(fit_track);
    }
}

void TCBuilder::fit(){
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    setPatternSize(patterns[i]->getNbLayers()-patterns[i]->getNbFakeSuperstrips());
    fit(allHits, patterns[i]->getOrderInChip());
  }
}

TrackFitter* TCBuilder::clone(){
  TCBuilder* fit = new TCBuilder(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

