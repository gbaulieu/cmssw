#ifndef _GRADEDPATTERN_H_
#define _GRADEDPATTERN_H_

#include "Pattern.h"

#include <boost/serialization/base_object.hpp>

using namespace std;

/**
   \brief A Pattern wih a grade (the number of tracks corresponding to this pattern).
   Can also contain (optionnaly) the average Pt of the tracks
**/

class GradedPattern : public Pattern{
 public:
  /**
     \brief Constructor
  **/
  GradedPattern();
  /**
     \brief Constructor
  **/
  GradedPattern(int nbLayers);
  /**
     \brief Copy Constructor
  **/
  GradedPattern(const Pattern& p);
  GradedPattern(const GradedPattern& p);

  virtual ~GradedPattern(){};
  /**
     \brief Get the grade of the Pattern
     \return The number of tracks having generated the pattern
  **/
  int getGrade() const;
  /**
     \brief Get the average Pt of the tracks having generated the pattern (if used)
     \return The average Pt
  **/
  float getAveragePt() const;
  /**
     \brief Get the minimal Pt of the tracks having generated the pattern
     \return The minimal Pt
  **/
  float getMinPt() const;
  /**
     \brief Get the maximal Pt of the tracks having generated the pattern
     \return The maximal Pt
  **/
  float getMaxPt() const;
  /**
     \brief Get the charge of the generating particles
     \return -1 for PDG<0, 1 for PDG>0, 0 if different charge were used
  **/
  int getCharge() const;
  /**
     Increment the grade (tracks occurences + 1)
  **/
  void increment();
  /**
     Increment the grade (tracks occurences + 1), add a Pt value to the average Pt, add the charge according to the PDG
     @param pt The Pt value of the last track
  **/
  void increment(float pt, int pdg);
  /**
     \brief Allows to compare 2 patterns on their grade
     \param gp The second pattern
     \return -1 if the pattern has a lower grade
  **/
  int operator<(const GradedPattern& gp);
  /**
     \brief Update the data of the pattern according to the data of the argument
   **/
  void mergeData(const GradedPattern& gp);

 private:
  int grade;
  float averagePt;
  float minPT;
  float maxPT;
  int charge;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const{
    ar << boost::serialization::base_object<Pattern>(*this);
    ar << grade;
    ar << averagePt;
    ar << minPT;
    ar << maxPT;
    ar << charge;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version){
    ar >> boost::serialization::base_object<Pattern>(*this);
    ar >> grade;
    ar >> averagePt;
    if(version>0){
      ar >> minPT;
      ar >> maxPT;
      ar >> charge;
    }
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

};
BOOST_CLASS_VERSION(GradedPattern, 1)
#endif
