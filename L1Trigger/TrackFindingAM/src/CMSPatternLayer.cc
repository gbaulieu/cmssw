#include "../interface/CMSPatternLayer.h"

short CMSPatternLayer::MOD_START_BIT = 10;
short CMSPatternLayer::PHI_START_BIT = 15;
short CMSPatternLayer::STRIP_START_BIT = 0;
short CMSPatternLayer::SEG_START_BIT = 0;
short CMSPatternLayer::MOD_MASK = 0x1F;
short CMSPatternLayer::PHI_MASK = 0x1;
short CMSPatternLayer::STRIP_MASK = 0x3FF;
short CMSPatternLayer::SEG_MASK = 0x0;
short CMSPatternLayer::OUTER_LAYER_SEG_DIVIDE = 2;
short CMSPatternLayer::INNER_LAYER_SEG_DIVIDE = 2;

map<string, int> CMSPatternLayer::phi_lut = loadPhiLUT("/scratch/lut.txt");
map<string, int> CMSPatternLayer::z_lut = loadZLUT("/scratch/lut.txt");

CMSPatternLayer::CMSPatternLayer():PatternLayer(){

}

CMSPatternLayer* CMSPatternLayer::clone(){
  CMSPatternLayer* p = new CMSPatternLayer();
  p->bits=this->bits;
  memcpy(p->dc_bits,this->dc_bits, DC_BITS*sizeof(char));
  return p;
}


bool CMSPatternLayer::isFake(){
  return (getPhi()==1);
}

vector<SuperStrip*> CMSPatternLayer::getSuperStrip(int l, const vector<int>& ladd, const map<int, vector<int> >& modules, Detector& d){
  int nb_dc = getDCBitsNumber();
  vector<SuperStrip*> v;

  if(isFake()){ // this is a fake superstrip! We link it to the dump superstrip
    vector<short> positions = getPositionsFromDC();
    if(positions.size()==0){
      SuperStrip* patternStrip = d.getDump();
      v.push_back(patternStrip);
    }
    else{
      for(unsigned int i=0;i<positions.size();i++){
	SuperStrip* patternStrip = d.getDump();
	v.push_back(patternStrip);
      }
    }
    return v;
  }
  else{
    Layer* la = d.getLayerFromAbsolutePosition(l);
    if(la!=NULL){
      //int ladderID = ladd[getPhi()];//getPhi() is the position in the sector;ladd[getPhi()] gives the ID of the ladder
      Ladder* patternLadder = la->getLadder(0);
      if(patternLadder!=NULL){
	//map<int, vector<int> >::const_iterator iterator = modules.find(ladderID); // get the vector of module IDs for this ladder
	//int moduleID = iterator->second[getModule()];// get the module ID from its position
	Module* patternModule = patternLadder->getModule(0);
	if(patternModule!=NULL){
	  Segment* patternSegment = patternModule->getSegment(getModule());
	  if(patternSegment!=NULL){
	    int base_index = getStripCode()<<nb_dc;
	    if(nb_dc>0){
	      vector<short> positions=getPositionsFromDC();
	      for(unsigned int i=0;i<positions.size();i++){
		int index = base_index | positions[i];
		SuperStrip* patternStrip = patternSegment->getSuperStripFromIndex(grayToBinary(index));
		v.push_back(patternStrip);
	      }
	    }
	    else{
	      SuperStrip* patternStrip = patternSegment->getSuperStripFromIndex(grayToBinary(base_index));
	      v.push_back(patternStrip);
	    }
	    return v;
	  }
	  else{
	    cout<<"cannot find segment "<<getModule()<<endl;
	  }
	}
	else{
	  cout<<"cannot find module 0"<<endl;
	}
      }
      else{
	cout<<"cannot find ladder 0"<<endl;
      }
    }
    cout<<"Error : can not link layer "<<l<<" ladder "<<ladd[0]<<" module "<<getModule()<<" segment "<<getSegment()<<" strip "<<getStrip()<<endl;
  }
  return v;
}

#ifdef IPNL_USE_CUDA
void CMSPatternLayer::getSuperStripCuda(int l, const vector<int>& ladd, const map<int, vector<int> >& modules, int layerID, unsigned int* v){
  int nb_dc = getDCBitsNumber();

  if(isFake()){ // this is a fake superstrip! -> No index
    return;
  }
  else{
    int layer_index = cuda_layer_index[layerID];
    if(layer_index!=-1){
      int ladderID = ladd[getPhi()];//getPhi() is the position in the sector;ladd[getPhi()] gives the ID of the ladder
      map<int, vector<int> >::const_iterator iterator = modules.find(ladderID); // get the vector of module IDs for this ladder
      int moduleID = iterator->second[getModule()];// get the module ID from its position
      int segment = getSegment();
      int base_index = getStripCode()<<nb_dc;
      vector<short> positions=getPositionsFromDC();
      for(unsigned int i=0;i<positions.size();i++){
	int strip_index = base_index | positions[i];
	int index = layer_index*SIZE_LAYER+ladderID*SIZE_LADDER+moduleID*SIZE_MODULE+segment*SIZE_SEGMENT+grayToBinary(strip_index);
	v[i]=index;
      }
      return;
    }
    cout<<"Error : can not link layer "<<l<<" ladder "<<ladd[getPhi()]<<" module "<<getModule()<<" segment "<<getSegment()<<" strip "<<getStrip()<<endl;
  }
}
#endif

short CMSPatternLayer::binaryToGray(short num)
{
  return (num >> 1) ^ num;
}

short CMSPatternLayer::grayToBinary(short gray)
{
  gray ^= (gray >> 8);
  gray ^= (gray >> 4);
  gray ^= (gray >> 2);
  gray ^= (gray >> 1);
  return(gray);
}

void CMSPatternLayer::setValues(short m, short phi, short strip, short seg){
  strip=binaryToGray(strip);
  bits |= (m&MOD_MASK)<<MOD_START_BIT |
    (phi&PHI_MASK)<<PHI_START_BIT |
    (strip&STRIP_MASK)<<STRIP_START_BIT |
    (seg&SEG_MASK)<<SEG_START_BIT;
}

void CMSPatternLayer::computeSuperstrip(short layerID, short module, short phi, short strip, short seg, int sstripSize, bool fake){
    if(fake){
      setValues(0, 1, 0, 0);
      return;
    }

    ostringstream oss;
    oss<<std::setfill('0');    
    oss<<layerID<<setw(2)<<phi<<module;
    int delta_z = z_lut.at(oss.str());
    //cout<<"z_lut["<<oss.str()<<") = "<<delta_z<<endl;
    //if(module%2!=0)
    //  delta_z++;
    
    int z = delta_z-seg;
 
    //if(layerID>7 && layerID<11) 
    //  z = z-seg;
    
    //cout<<"-> module = "<<z<<endl;
    
    if(layerID>=5 && layerID<=7)
      z = z/(32*INNER_LAYER_SEG_DIVIDE);
    else if (layerID>10 && phi<=8)
      z = z/(32*INNER_LAYER_SEG_DIVIDE);
    else
      z = z/OUTER_LAYER_SEG_DIVIDE;

    /****************/
    if(module%2!=0){
      module--;
    
      switch(layerID){
      case 5:{
	strip = (1.026*strip)-13;
      }break;
      case 6:{
	strip = (1.017*strip)-8;
      }break;
      case 7:{
	strip = (0.988*strip)+7;
      }break;
      case 8:{
	strip = (1.009*strip)-5;
      }break;
      case 9:{
	strip = (1.006*strip)-4;
      }break;
      case 10:{
	strip = (1.006*strip)-3;
      }break;
      }
      oss.clear();
      oss.str("");
      oss<<layerID<<setw(2)<<phi<<module;
      
      try{
	phi_lut.at(oss.str());
      }
      catch (const std::out_of_range& oor) {
	module+=2;
	oss.clear();
	oss.str("");
	oss<<layerID<<setw(2)<<phi<<module;
      }
    }
    /***************/

    int delta_phi = phi_lut.at(oss.str());
    //cout<<"phi_lut["<<oss.str()<<") = "<<delta_phi<<endl;

    int superStrip = delta_phi-strip;
    //cout<<"-> Superstrip = "<<superStrip<<endl;
    //cout<<"utilisation de "<<sstripSize<<" strips par superstrip"<<endl;
    superStrip = superStrip/sstripSize;
    //cout<<"  -> superstrip apres resolution = "<<superStrip<<endl;
    
    //z=0;
    //  if(layerID>7 && layerID<=10)
    //cout<<"  -> module apres resolution = "<<z<<endl;
    
    setValues(z, 0, superStrip, 0);
}

short CMSPatternLayer::getModule(){
  int val = bits.to_ulong();
  short r = (val>>MOD_START_BIT)&MOD_MASK;
  return r;
}

short CMSPatternLayer::getPhi(){
  int val = bits.to_ulong();
  short r = (val>>PHI_START_BIT)&PHI_MASK;
  return r;
}

short CMSPatternLayer::getStrip(){
  int val = bits.to_ulong();
  short r = (val>>STRIP_START_BIT)&STRIP_MASK;
  r=grayToBinary(r);
  return r;
}

short CMSPatternLayer::getStripCode(){
  int val = bits.to_ulong();
  short r = (val>>STRIP_START_BIT)&STRIP_MASK;
  return r;
}

short CMSPatternLayer::getSegment(){
  int val = bits.to_ulong();
  short r = (val>>SEG_START_BIT)&SEG_MASK;
  return r;
}

string CMSPatternLayer::toString(){
  ostringstream oss;
  oss<<"Ladder "<<getPhi()<<" Module "<<getModule()<<" Segment "<<getSegment()<<" strip "<<getStrip();
  if(dc_bits[0]!=3){
    oss<<" (";
    for(int i=0;i<DC_BITS;i++){
      if(dc_bits[i]==2)
	oss<<"X";
      else if(dc_bits[i]!=3)
	oss<<(int)dc_bits[i];
    }
    oss<<")";
  }
  return oss.str();
}

string CMSPatternLayer::toStringSuperstripBinary(){
  short seg = getSegment();
  short sstrip = getStripCode();
  int initialValue = getIntValue();
  int moduleLadder = (initialValue>>PHI_START_BIT);//remove the superstrip and segment informations
  int newValue = 0;
  newValue |= (moduleLadder&0x1FF)<<PHI_START_BIT |
    (seg&SEG_MASK)<<(PHI_START_BIT-1) |
    (sstrip&STRIP_MASK)<<SEG_START_BIT;
  ostringstream oss;
  oss<<newValue;
  if(dc_bits[0]!=3){
    oss<<" (";
    for(int i=0;i<DC_BITS;i++){
      if(dc_bits[i]==2)
	oss<<"X";
      else if(dc_bits[i]!=3)
	oss<<(int)dc_bits[i];
    }
    oss<<")";
  }
  return oss.str();
}

string CMSPatternLayer::toStringBinary(){
  ostringstream oss;
  oss<<getIntValue();
  if(dc_bits[0]!=3){
    oss<<" (";
    for(int i=0;i<DC_BITS;i++){
      if(dc_bits[i]==2)
	oss<<"X";
      else if(dc_bits[i]!=3)
	oss<<(int)dc_bits[i];
    }
    oss<<")";
  }
  return oss.str();
}

string CMSPatternLayer::toAM05Format(){

  /**
     The input superstrip is 16 bits long
     The stored superstrip in the AM05 is 18 bits long (room for DC bits)
     At least 2 DC bits will be used to fit in the 18 bits : 16 bits are becoming 14 bits + 2DC bits -> 18 bits
  **/
  int nb_dc_bits = 0;
  int used_dc_bits = getDCBitsNumber();
  if(used_dc_bits<3)
    nb_dc_bits=2;
  else
    nb_dc_bits=used_dc_bits;

  //Default organization of bits (0 to 2 DC bits used)
  short AM05_MOD_START_BIT = 13;
  short AM05_PHI_START_BIT = 9;
  short AM05_SEG_START_BIT = 8;
  short AM05_STRIP_START_BIT = 4;
  short AM05_STRIP_DC0_BIT = 2;
  short AM05_STRIP_DC1_BIT = 0;
  short AM05_STRIP_DC2_BIT = 0;

  short AM05_MOD_MASK = 0x1F;
  short AM05_PHI_MASK = 0xF;
  short AM05_SEG_MASK = 0x1;
  short AM05_STRIP_MASK = 0xF;//4 bits + 2 DC bits
  short AM05_STRIP_DC0_MASK = 0x3;
  short AM05_STRIP_DC1_MASK = 0x3;
  short AM05_STRIP_DC2_MASK = 0;

  short z = getModule();
  short ladder = getPhi();
  short seg = getSegment();
  short sstrip = getStripCode();

  int am_format=0;//18 bits value for the AM05 chip  

  if(used_dc_bits==0){

    /*
      we need to encode the 2 last bits of the sstrip as DC bits
      First we retrieve the values and remove them from the sstrip position      
    */
    short dcbit0_val = (sstrip>>1)&0x1;
    short dcbit1_val = sstrip&0x1;
    sstrip = sstrip>>2;

    // Now we encode the values as DC bit values
    if(dcbit0_val==0)
      dcbit0_val=1;//01
    else
      dcbit0_val=2;//10

    if(dcbit1_val==0)
      dcbit1_val=1;//01
    else
      dcbit1_val=2;//10

    // in case this is a fake superstrip, it must not be activable : we use the 11 value of the DC bits
    if(isFake()){
      dcbit0_val=3;//11
      dcbit1_val=3;//11
    }

    //5 bits for Z + 4 bits for ladder + 1 bit for seg + 4 bits for sstrip + 2 bits for sstrips DC bit 0 + 2 bits for sstrips DC bit 1 = 18 bits
    am_format |= (z&AM05_MOD_MASK)<<AM05_MOD_START_BIT |
      (ladder&AM05_PHI_MASK)<<AM05_PHI_START_BIT |
      (seg&AM05_SEG_MASK)<<AM05_SEG_START_BIT |
      (sstrip&AM05_STRIP_MASK)<<AM05_STRIP_START_BIT |
      (dcbit0_val&AM05_STRIP_DC0_MASK)<<AM05_STRIP_DC0_BIT |
      (dcbit1_val&AM05_STRIP_DC1_MASK)<<AM05_STRIP_DC1_BIT;
  }
  else if(used_dc_bits==1){

    /*
      we need to encode the last bit of the sstrip as a DC bit
      First we retrieve the value and remove it from the sstrip position      
    */
    short dcbit0_val = sstrip&0x1;
    short dcbit1_val = dc_bits[0];
    sstrip = sstrip>>1;

    // Now we encode the values as DC bit values
    if(dcbit0_val==0)
      dcbit0_val=1;//01
    else
      dcbit0_val=2;//10

    switch(dcbit1_val){
    case 2 : dcbit1_val=0;//00:X
      break;
    case 0 : dcbit1_val=1;//01
      break;
    case 1 : dcbit1_val=2;//10
      break;
    }
    
    // in case this is a fake superstrip, it must not be activable : we use the 11 value of the DC bits
    if(isFake()){
      dcbit0_val=3;//11
      dcbit1_val=3;//11
    }

    //5 bits for Z + 4 bits for ladder + 1 bit for seg + 4 bits for sstrip + 2 bits for sstrips DC bit 0 + 2 bits for sstrips DC bit 1 = 18 bits
    am_format |= (z&AM05_MOD_MASK)<<AM05_MOD_START_BIT |
      (ladder&AM05_PHI_MASK)<<AM05_PHI_START_BIT |
      (seg&AM05_SEG_MASK)<<AM05_SEG_START_BIT |
      (sstrip&AM05_STRIP_MASK)<<AM05_STRIP_START_BIT |
      (dcbit0_val&AM05_STRIP_DC0_MASK)<<AM05_STRIP_DC0_BIT |
      (dcbit1_val&AM05_STRIP_DC1_MASK)<<AM05_STRIP_DC1_BIT;
  }
  else if(used_dc_bits==2){

    short dcbit0_val = dc_bits[0];
    short dcbit1_val = dc_bits[1];

    // we encode the values as DC bit values
    switch(dcbit0_val){
    case 2 : dcbit0_val=0;//00:X
      break;
    case 0 : dcbit0_val=1;//01
      break;
    case 1 : dcbit0_val=2;//10
      break;
    }

    switch(dcbit1_val){
    case 2 : dcbit1_val=0;//00:X
      break;
    case 0 : dcbit1_val=1;//01
      break;
    case 1 : dcbit1_val=2;//10
      break;
    }

    // in case this is a fake superstrip, it must not be activable : we use the 11 value of the DC bits
    if(isFake()){
      dcbit0_val=3;//11
      dcbit1_val=3;//11
    }    

    //5 bits for Z + 4 bits for ladder + 1 bit for seg + 4 bits for sstrip + 2 bits for sstrips DC bit 0 + 2 bits for sstrips DC bit 1 = 18 bits
    am_format |= (z&AM05_MOD_MASK)<<AM05_MOD_START_BIT |
      (ladder&AM05_PHI_MASK)<<AM05_PHI_START_BIT |
      (seg&AM05_SEG_MASK)<<AM05_SEG_START_BIT |
      (sstrip&AM05_STRIP_MASK)<<AM05_STRIP_START_BIT |
      (dcbit0_val&AM05_STRIP_DC0_MASK)<<AM05_STRIP_DC0_BIT |
      (dcbit1_val&AM05_STRIP_DC1_MASK)<<AM05_STRIP_DC1_BIT;
  }
  else if(used_dc_bits==3){

    AM05_MOD_START_BIT = 14;
    AM05_PHI_START_BIT = 10;
    AM05_SEG_START_BIT = 9;
    AM05_STRIP_START_BIT = 6;
    AM05_STRIP_DC0_BIT = 4;
    AM05_STRIP_DC1_BIT = 2;
    AM05_STRIP_DC2_BIT = 0;
    
    AM05_MOD_MASK = 0xF;
    AM05_PHI_MASK = 0xF;
    AM05_SEG_MASK = 0x1;
    AM05_STRIP_MASK = 0x7;//3 bits + 3 DC bits
    AM05_STRIP_DC0_MASK = 0x3;
    AM05_STRIP_DC1_MASK = 0x3;
    AM05_STRIP_DC2_MASK = 0x3;

    short dcbit0_val = dc_bits[0];
    short dcbit1_val = dc_bits[1];
    short dcbit2_val = dc_bits[2];

    // we encode the values as DC bit values
    switch(dcbit0_val){
    case 2 : dcbit0_val=0;//00:X
      break;
    case 0 : dcbit0_val=1;//01
      break;
    case 1 : dcbit0_val=2;//10
      break;
    }

    switch(dcbit1_val){
    case 2 : dcbit1_val=0;//00:X
      break;
    case 0 : dcbit1_val=1;//01
      break;
    case 1 : dcbit1_val=2;//10
      break;
    }

    switch(dcbit2_val){
    case 2 : dcbit2_val=0;//00:X
      break;
    case 0 : dcbit2_val=1;//01
      break;
    case 1 : dcbit2_val=2;//10
      break;
    }
    
    //we are using 4 bits for the Z value so it must be below 16 (should be ok with official trigger towers)
    if(z>15){
      cout<<"The module value is too high ("<<z<<">15) : pattern can not be stored in an AM05 chip"<<endl;
      exit(-1);
    }

    // in case this is a fake superstrip, it must not be activable : we use the 11 value of the DC bits
    if(isFake()){
      dcbit0_val=3;//11
      dcbit1_val=3;//11
      dcbit2_val=3;//11
    }

    //4 bits for Z + 4 bits for ladder + 1 bit for seg + 3 bits for sstrip + 2 bits for sstrips DC bit 0 + 2 bits for sstrips DC bit 1 + 2 bits for sstrips DC bit 2 = 18 bits
    am_format |= (z&AM05_MOD_MASK)<<AM05_MOD_START_BIT |
      (ladder&AM05_PHI_MASK)<<AM05_PHI_START_BIT |
      (seg&AM05_SEG_MASK)<<AM05_SEG_START_BIT |
      (sstrip&AM05_STRIP_MASK)<<AM05_STRIP_START_BIT |
      (dcbit0_val&AM05_STRIP_DC0_MASK)<<AM05_STRIP_DC0_BIT |
      (dcbit1_val&AM05_STRIP_DC1_MASK)<<AM05_STRIP_DC1_BIT |
      (dcbit2_val&AM05_STRIP_DC2_MASK)<<AM05_STRIP_DC2_BIT;
  }
  ostringstream oss;
  oss<<am_format<<" "<<nb_dc_bits;
  return oss.str();
}

vector<int> CMSPatternLayer::getLayerIDs(){
  vector<int> layers;

  //BARREL
  layers.push_back(5);
  layers.push_back(6);
  layers.push_back(7);
  layers.push_back(8);
  layers.push_back(9);
  layers.push_back(10);

  //ENDCAP 1
  layers.push_back(11);
  layers.push_back(12);
  layers.push_back(13);
  layers.push_back(14);
  layers.push_back(15);

  //ENDCAP 2
  layers.push_back(18);
  layers.push_back(19);
  layers.push_back(20);
  layers.push_back(21);
  layers.push_back(22);

  return layers;
}

int CMSPatternLayer::getNbStripsInSegment(){
  return 1024;
}

int CMSPatternLayer::getSegmentCode(int layerID, int ladderID, int segmentID){
  return segmentID;
  /*
  if(layerID>7 && layerID<11)
    return segmentID;
  if(layerID>=5 && layerID<=7)
    return segmentID/16;
  if(ladderID<=8)
    return segmentID/16;
  return segmentID;
  */
}


int CMSPatternLayer::getModuleCode(int layerID, int moduleID){
  switch(layerID){
  case 5 : return (moduleID);
  case 6 : return (moduleID);
  case 7 : return (moduleID);
  case 8 : return moduleID;
  case 9 : return moduleID;
  case 10 : return moduleID;
  default : return moduleID;
  }
}

int CMSPatternLayer::getLadderCode(int layerID, int ladderID){
  return ladderID;
}

 int CMSPatternLayer::getNbLadders(int layerID){
   if(layerID<5 || layerID>24)
     return -1;
   switch(layerID){
   case 5 : return 16;
   case 6 : return 24;
   case 7 : return 34;
   case 8 : return 48;
   case 9 : return 62;
   case 10 : return 76;
   default : return 15;
   }
 }

int CMSPatternLayer::getNbModules(int layerID, int ladderID){
  if(layerID==5)
    return 64;
  if(layerID==6)
    return 56;
  if(layerID==7)
    return 54;
  if(layerID>=8 && layerID<=10)
    return 24;
  if(layerID>=11 && layerID<=24){
    switch(ladderID){
    case 0:return 20;
    case 1:return 24;
    case 2:return 28;
    case 3:return 28;
    case 4:return 32;
    case 5:return 36;
    case 6:return 36;
    case 7:return 40;
    case 8:return 40;
    case 9:return 52;
    case 10:return 56;
    case 11:return 64;
    case 12:return 68;
    case 13:return 76;
    case 14:return 80;
    default:return 80;
    }
  }
  return -1;
}

map<int, pair<float,float> > CMSPatternLayer::getLayerDefInEta(){
  map<int,pair<float,float> > eta;
  eta[5]=pair<float,float>(-2.2,2.2);
  eta[6]=pair<float,float>(-1.72,1.72);
  eta[7]=pair<float,float>(-1.4,1.4);
  eta[8]=pair<float,float>(-1.2,1.2);
  eta[9]=pair<float,float>(-1.1,1.1);
  eta[10]=pair<float,float>(-0.9,0.9);
  eta[11]=pair<float,float>(1.08,2.24);
  eta[12]=pair<float,float>(1.21,2.45);
  eta[13]=pair<float,float>(1.36,2.5);
  eta[14]=pair<float,float>(1.49,2.5);
  eta[15]=pair<float,float>(1.65,2.5);
  eta[18]=pair<float,float>(-2.24,-1.08);
  eta[19]=pair<float,float>(-2.45,-1.21);
  eta[20]=pair<float,float>(-2.5,-1.36);
  eta[21]=pair<float,float>(-2.5,-1.49);
  eta[22]=pair<float,float>(-2.5,-1.65);
  return eta;
}

map< string, int > CMSPatternLayer::loadPhiLUT(string name){
  string line;
  ifstream myfile (name.c_str());
  map< string, int > phi_lut;
  if (myfile.is_open()){
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
	if(items.size()==6){
	  istringstream buffer(items[2]);
	  istringstream buffer2(items[3]);
	  int startStrip;
	  int superStripIndex;
	  buffer >> startStrip;
	  buffer2 >> superStripIndex;
	  phi_lut[items[1]]=startStrip+superStripIndex;
	}
      }
    }
    myfile.close();
  }
  else{
    cout << "Can not find file "<<name<<" to load the lookup table!"<<endl;
    exit(-1);
  }
  return phi_lut;
}

map< string, int > CMSPatternLayer::loadZLUT(string name){
  string line;
  ifstream myfile (name.c_str());
  map< string, int > phi_lut;
  if (myfile.is_open()){
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
	if(items.size()==6){
	  istringstream buffer(items[4]);
	  istringstream buffer2(items[5]);
	  int startStrip;
	  int superStripIndex;
	  buffer >> startStrip;
	  buffer2 >> superStripIndex;
	  phi_lut[items[1]]=startStrip+superStripIndex;
	}
      }
    }
    myfile.close();
  }
  else{
    cout << "Can not find file "<<name<<" to load the lookup table!"<<endl;
    exit(-1);
  }
  return phi_lut;
}

vector<int> CMSPatternLayer::getHDSuperstrips(){
  vector<int> array;
  int nb_dc = getDCBitsNumber();
  int base_index = getStripCode()<<nb_dc;
  if(nb_dc>0){
    vector<short> positions=getPositionsFromDC();
    for(unsigned int i=0;i<positions.size();i++){
      int index = base_index | positions[i];
      array.push_back(grayToBinary(index));
    }
  }
  else{
    array.push_back(grayToBinary(base_index));
  }
  return array;
}
