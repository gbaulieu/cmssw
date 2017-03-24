#include "../interface/PRBF2LocalToGlobalConverter.h"

PRBF2LocalToGlobalConverter::PRBF2LocalToGlobalConverter(const Sector* sectorDefinition, string geometryFile):LocalToGlobalConverter(){
  string line;
  ifstream myfile (geometryFile.c_str());
  int layer = -1;
  int ladder = -1;
  int local_ladder = -1;
  int module = -1;
  int local_module = -1;
  float coef_value = -1;
  stringstream val;

  module_pos = NULL;
  sector = sectorDefinition;

  //Which side of the tracker are we on?
  if(sectorDefinition->getOfficialID()<24)
    tracker_side=true;
  else
    tracker_side=false;

  //Parse the geometry file containing the cartesian coordinates of all modules
  if (myfile.is_open()){
        
    //Memory allocation of the structure
    module_pos = new vector<float>**[16];
    for(int i=0;i<16;i++){//we can have up to 16 layers
      module_pos[i] = new vector<float>*[16];
      for(int j=0;j<16;j++){//16 ladders
	module_pos[i][j] = new vector<float>[32];
	for(int k=0;k<32;k++){//32 modules
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	  module_pos[i][j][k].push_back(-1);
	}
      }
    }
    
    while ( myfile.good() ){
      getline (myfile,line);
      if(line.length()>0 && line.find("#")!=0){
	stringstream ss(line);
	std::string item;
	vector<string> items;
	while (getline(ss, item, '/')) {
	  std::string::iterator end_pos = std::remove(item.begin(), item.end(), ' ');
	  item.erase(end_pos, item.end());
	  items.push_back(item);
	}
	if(items.size()==11){
	  val.clear();
	  val.str(items[0]);
	  val >> layer;
	  val.clear();
	  val.str(items[1]);
	  val >> ladder;
	  ladder = ladder-1;//numbering in file starts with 1 -> corrects to 0
	  val.clear();
	  val.str(items[2]);
	  val >> module;
	  local_ladder = sectorDefinition->getLadderCode(layer,ladder);
	  if(local_ladder==-1)//not in the sector
	    continue;
	  local_module = sectorDefinition->getModuleCode(layer,ladder,module);
	  if(local_module==-1)//not in the sector
	    continue;
	  bool isPSModule = false;
	  if((layer>=5 && layer<=7) || (layer>10 && ladder<=8)){
	    isPSModule=true;
	  }
	  int prbf2_layer = CMSPatternLayer::cmssw_layer_to_prbf2_layer(layer,isPSModule);

	  val.clear();
	  val.str(items[8]);
	  val >> coef_value;
	  module_pos[prbf2_layer][local_ladder][local_module][0]=coef_value;
	  val.clear();
	  val.str(items[9]);
	  val >> coef_value;
	  module_pos[prbf2_layer][local_ladder][local_module][1]=coef_value;
	  val.clear();
	  val.str(items[10]);
	  val >> coef_value;
	  module_pos[prbf2_layer][local_ladder][local_module][2]=coef_value;

	  //Process the starting phi of the tower
	  double sec_phi = (sectorDefinition->getOfficialID()%8) * M_PI / 4.0 - 0.4;
	  
	  //cos and sin values for a rotation of an angle -sec_phi
	  double ci = cos(-sec_phi);
	  double si = sin(-sec_phi);

	  //Rotates the module center to the first sector
	  double rotatedX = module_pos[prbf2_layer][local_ladder][local_module][0] * ci - module_pos[prbf2_layer][local_ladder][local_module][1] * si;
	  double rotatedY = module_pos[prbf2_layer][local_ladder][local_module][0] * si + module_pos[prbf2_layer][local_ladder][local_module][1] * ci;
	  module_pos[prbf2_layer][local_ladder][local_module][0] = rotatedX;
	  module_pos[prbf2_layer][local_ladder][local_module][1] = rotatedY;

	  //Computes the angle between X axis and the module center
	  float PhiMod = atan2(module_pos[prbf2_layer][local_ladder][local_module][1],module_pos[prbf2_layer][local_ladder][local_module][0]);//atan2(y,x)
	  float fModuleWidth,fModuleHeight,nModuleStrips,nModuleSegments;
	  if (isPSModule){
	    //PS
	    fModuleWidth 	= 4.48144;
	    fModuleHeight 	= 9.59;
	    nModuleStrips	= 959.0;
	    nModuleSegments	= 31.0;
	  }
	  else{
	    //2S
	    fModuleWidth 	= 5.025;
	    fModuleHeight 	= 9.135;
	    nModuleStrips	= 1015.0;
	    nModuleSegments	= 1.0;
	  }
	  float fStripPitch = fModuleHeight / nModuleStrips;
	  float fSegmentPitch = fModuleWidth / nModuleSegments;

	  if(layer<11){
	    module_pos[prbf2_layer][local_ladder][local_module][3]=fStripPitch * sin(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][4]=-fStripPitch * cos(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][5]=0.0;
	    module_pos[prbf2_layer][local_ladder][local_module][6]=0.0;
	    module_pos[prbf2_layer][local_ladder][local_module][7]=0.0;
	    module_pos[prbf2_layer][local_ladder][local_module][8]=-fSegmentPitch;
	  }
	  else{
	    if(tracker_side){
	      module_pos[prbf2_layer][local_ladder][local_module][3]=fStripPitch * sin(PhiMod);
	      module_pos[prbf2_layer][local_ladder][local_module][4]=-fStripPitch * cos(PhiMod);
	    }
	    else{
	      module_pos[prbf2_layer][local_ladder][local_module][3]=-fStripPitch * sin(PhiMod);
	      module_pos[prbf2_layer][local_ladder][local_module][4]=fStripPitch * cos(PhiMod);
	    }
	    module_pos[prbf2_layer][local_ladder][local_module][5]=0.0;
	    module_pos[prbf2_layer][local_ladder][local_module][6]=-fSegmentPitch * cos(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][7]=-fSegmentPitch * sin(PhiMod);
	    module_pos[prbf2_layer][local_ladder][local_module][8]=0.0;
	  }

	  //Apply the HW binning to the coefficients 
	  module_pos[prbf2_layer][local_ladder][local_module][0] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][0], 6, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][1] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][1], 6, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][2] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][2], 8, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][3] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][3], -7, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][4] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][4], -7, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][5] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][5], -7, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][6] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][6], 2, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][7] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][7], 2, 18, SIGNED);
	  module_pos[prbf2_layer][local_ladder][local_module][8] = CommonTools::binning(module_pos[prbf2_layer][local_ladder][local_module][8], 2, 18, SIGNED);
	  
	}
      }
    }
    myfile.close();
  }
  else{
    cout << "Can not find file "<<geometryFile<<" to load the modules position lookup table!"<<endl;
  }
}

PRBF2LocalToGlobalConverter::~PRBF2LocalToGlobalConverter(){ 
  if(module_pos!=NULL){
    for(int i=0;i<16;i++){//we can have up to 16 layers
      for(int j=0;j<16;j++){//16 ladders
	for(int k=0;k<32;k++){//32 modules
	  module_pos[i][j][k].clear();
	}
	delete [] module_pos[i][j];
      }
      delete [] module_pos[i];
    }
    delete [] module_pos;
  }
}

vector<float> PRBF2LocalToGlobalConverter::toGlobal(const Hit* h) const throw (std::runtime_error){
  int hit_layer = h->getLayer();
  int hit_ladder = h->getLadder();
  int hit_module = h->getModule();	    
  bool isPSModule = false;
  if((hit_layer>=5 && hit_layer<=7) || (hit_layer>10 && hit_ladder<=8)){
    isPSModule=true;
  }
  int prbf2_layer = CMSPatternLayer::cmssw_layer_to_prbf2_layer(hit_layer,isPSModule);
  int prbf2_ladder = sector->getLadderCode(hit_layer, hit_ladder);
  int prbf2_module = sector->getModuleCode(hit_layer, hit_ladder, hit_module);

  return toGlobal(prbf2_layer,prbf2_ladder, prbf2_module, h->getSegment(), h->getHDStripNumber());
}

vector<float> PRBF2LocalToGlobalConverter::toGlobal(int layer, int ladder, int module, int segment, float strip) const throw (std::runtime_error){

  if(module_pos==NULL){
    throw std::runtime_error("Modules position lookup table not found");
  }

  //Correction of a strip number bug (when strip decimal is not zero, it has to be 0.5)
  if (strip != floor(strip)){
    strip = floor(strip)+0.5;
  }

  bool isPS = (layer<8);
  float relatStrip=0.0;
  float relatSeg=0.0;
  float X;
  float Y;
  float Z;
  
  if(isPS){
    relatStrip = strip-(959.0/2.0);
    relatSeg   = segment-(31.0/2.0);
  }
  else{
    relatStrip = strip-(1015.0/2.0);
    relatSeg   = segment-(1.0/2.0);
  }

  vector<float> res;
  vector<float> positions = module_pos[layer][ladder][module];


  X = positions[0];
  Y = positions[1];
  Z = positions[2];

  X += CommonTools::binning(relatStrip*positions[3], 6, 18, SIGNED);
  Y += CommonTools::binning(relatStrip*positions[4], 6, 18, SIGNED);
  Z += CommonTools::binning(relatStrip*positions[5], 8, 18, SIGNED);

  X += CommonTools::binning(relatSeg*positions[6], 6, 18, SIGNED);
  Y += CommonTools::binning(relatSeg*positions[7], 6, 18, SIGNED);
  Z += CommonTools::binning(relatSeg*positions[8], 8, 18, SIGNED);

  res.push_back(CommonTools::binning(X, 6, 18, SIGNED));
  res.push_back(CommonTools::binning(Y, 6, 18, SIGNED));
  res.push_back(CommonTools::binning(Z, 8, 18, SIGNED));

  return res;
}
