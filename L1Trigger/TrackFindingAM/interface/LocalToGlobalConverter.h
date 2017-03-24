#ifndef _LOCALTOGLOBALCONVERTER_H_
#define _LOCALTOGLOBALCONVERTER_H_

#include <vector>
#include <stdexcept>
#include "CommonTools.h"
#include "Hit.h"

using namespace std;

/**
   \brief Converts local coordinates (layer, ladder, module, segment, hit local to the sector) to global cartesian coordinates
 **/
class LocalToGlobalConverter{

 protected:
  /**
     module_pos[prbf2_layer][prbf2_ladder][prbf2_module] contains 
     [
     module_center_X, 
     module_center_Y, 
     module_center_Y, 
     strip_pitch_X, 
     strip_pitch_Y, 
     strip_pitch_Z,
     seg_pitch_X, 
     seg_pitch_Y, 
     seg_pitch_Z
     ]
  **/
  //TEC+ or TEC-?
  bool tracker_side;

 public:
  /**
     \brief Default constructor of the abstract class
  **/
  LocalToGlobalConverter();

  /**
     \brief Destructor
  **/
  virtual ~LocalToGlobalConverter();

  /**
     \brief Convertion of local coordinates (layer/ladder/module/segment/strip) to cartesian coordinates global to the tracker. Whether you need PRBF2 or CMSSW local values relies upon the underlying implementation.
     \param layer The layer information
     \param ladder The ladder information
     \param module The module information
     \param segment The segment information
     \param strip The strip information
     \return A vector of 3 float containing the X/Y/Z cartesian coordinates
   **/
  virtual vector<float> toGlobal(int layer, int ladder, int module, int segment, float strip) const throw (std::runtime_error) =0;

  /**
     \brief Convertion of local coordinates contained in a Hit object to cartesian coordinates global to the tracker
     \param h A Hit object
     \return A vector of 3 float containing the X/Y/Z cartesian coordinates
   **/
  virtual vector<float> toGlobal(const Hit* h) const throw (std::runtime_error) =0;

};
#endif
