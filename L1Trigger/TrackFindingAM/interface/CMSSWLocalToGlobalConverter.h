#ifndef _CMSSWLOCALTOGLOBALCONVERTER_H_
#define _CMSSWLOCALTOGLOBALCONVERTER_H_

#include "LocalToGlobalConverter.h"
#include "Sector.h"
#include "CommonTools.h"

using namespace std;

/**
   \brief Converts local coordinates (layer, ladder, module, segment, hit local to the sector) to global cartesian coordinates
 **/
class CMSSWLocalToGlobalConverter:public LocalToGlobalConverter{

 private:
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
  map<int, map<int, map<int, vector<float>>>> module_pos;

 public:
  /**
     \brief Constructor
     \param sectorID The sector concerned by the convertion
     \param geometryFile Text file containing the global coordinates of all modules
  **/
  CMSSWLocalToGlobalConverter(int sectorID, string geometryFile);

  /**
     \brief Destructor
  **/
  ~CMSSWLocalToGlobalConverter();

  /**
     \brief Convertion of geometric coordinates (CMSSW layer/ladder/module/segment/strip format) to cartesian coordinates global to the tracker
     \param layer The CMSSW value of the layer (5 to 22)
     \param ladder The CMSSW value of the ladder
     \param module The CMSSW value of the module
     \param segment The CMSSW value of the segment
     \param strip The CMSSW strip value
     \return A vector of 3 float containing the X/Y/Z cartesian coordinates
   **/
  vector<float> toGlobal(int layer, int ladder, int module, int segment, float strip) const throw (std::runtime_error);

  /**
     \brief Convertion of local coordinates contained in a Hit object to cartesian coordinates global to the tracker
     \param h A Hit object
     \return A vector of 3 float containing the X/Y/Z cartesian coordinates
   **/
  vector<float> toGlobal(const Hit* h) const throw (std::runtime_error);

};
#endif
