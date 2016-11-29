#ifndef _TCBUILDER_H_
#define _TCBUILDER_H_

#include "TrackFitter.h"
#include "LocalToGlobalConverter.h"
#include <math.h>

#include <iomanip>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>


enum SEC_TYPE {SEC_BARREL, SEC_HYBRID, SEC_ENDCAP};

/**
   \brief Implementation of a Seed Clustering  fitter
**/
class TCBuilder:public TrackFitter{

 private:

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
    {
      ar << boost::serialization::base_object<TrackFitter>(*this);
      ar << sec_phi;
    }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
    {
      ar >> boost::serialization::base_object<TrackFitter>(*this);
      ar >> sec_phi;
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  double m_tabBarrelThresholds[2][3][10][2];  //Barrel LayMaxSeed1 = 1, LayMaxSeed2 = 2, LayMaxTestStub = 10
  double m_tabHybridThresholds[2][3][14][2];  //Hybrid LayMaxSeed1 = 1, LayMaxSeed2 = 2, LayMaxTestStub = 13
  double m_tabEndcapThresholds[4][5][16][2];  //Endcap LayMaxSeed1 = 3, LayMaxSeed2 = 4, LayMaxTestStub = 15

  int m_nMissingHits;
  int maxseeds;
  int rsize;
  LocalToGlobalConverter* l2gConverter;

  /**
     \brief Updates the Thresholds with respect to the fractionnal part width
   **/
  void updateThresholds();
  void addThresholds(int, int, int, SEC_TYPE, double, double);
  void getThresholds(int, int, int, SEC_TYPE, double[]);
  char transcodeLayer(Hit *);
  void alignScore(Hit& , Hit& , Hit& , double []);

 public:

  /**
     \brief Default constructor   
  **/
  TCBuilder();
  /**
     \brief Constructor   
     \param nb Layers number
  **/  
  TCBuilder(int nb);
  ~TCBuilder();

  void initialize();
  void mergePatterns();
  void mergeTracks();

  Track* createFittedTrack(vector <Hit*>&);

  void setLocalToGlobalConverter(LocalToGlobalConverter* l);
  void setMaxSeeds(int nmax);
  void setSize(int rs);

  /**
     \brief Configure the way the computing is done
     \param hardwareEmulation If true the computing will be done using the hardware precision. If false, float computing is used.
   **/
  void setHardwareEmulation(bool hardwareEmulation);

  void fit();
  void fit(vector<Hit*> hits, int pattern_id=-1);
  TrackFitter* clone();
};
#endif
