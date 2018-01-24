#include "../interface/PatternFinder.h"

#ifndef IPNL_USE_CUDA
PatternFinder::PatternFinder(int at, SectorTree* st, string& f, string& of):
  active_threshold(at), max_nb_missing_hit(0), useMissingHits(false),
  max_road_number(1000000), sectors(st), eventsFileName(f),
  outputFileName(of), tracker(), converter(NULL), hardware_limitations(false)
{

  //we don't need the map of patterns, a vector will be enough and uses less memory
  sectors->getAllSectors()[0]->getPatternTree()->switchToVector();

  map< int, vector<int> > detector_config = Sector::readConfig("detector.cfg");

  //We add the layers corresponding to the sector structure
  vector<Sector*> sector_list = sectors->getAllSectors();
  if(sector_list.size()>0){
    for(int i=0;i<sector_list[0]->getNbLayers();i++){
      int layerID = sector_list[0]->getLayerID(i);
      if(!detector_config.empty()){
	if(detector_config.find(layerID)!=detector_config.end())
	  if(layerID<11)//barrel : 1 module with 2*nb_modules segments
	    tracker.addLayer(detector_config[layerID][0],detector_config[layerID][1],1, detector_config[layerID][2]*2, detector_config[layerID][3], SectorTree::getSuperstripSize(layerID), true);
	  else // endcap
	    tracker.addLayer(detector_config[layerID][0],detector_config[layerID][1], detector_config[layerID][2],2, detector_config[layerID][3], 8, false);
	else
	  cout<<"WARNING : Layer "<<layerID<<" is used in the sector definition of the bank but is missing in the configuration of the virtual detector"<<endl;
      }
    }
  }

  tracker.setSectorMaps(sector_list[0]->getLadderCodeMap(),sector_list[0]->getModuleCodeMap());

  //Link the patterns with the tracker representation
  cout<<"linking..."<<endl;
  sectors->link(tracker);
  cout<<"done."<<endl;
}
#endif

PatternFinder::~PatternFinder(){
  if(converter!=NULL)
    delete converter;
}

#ifdef IPNL_USE_CUDA
PatternFinder::PatternFinder(int at, SectorTree* st, string& f, string& of, patternBank* p, deviceDetector* d, deviceParameters* dp):
  active_threshold(at), max_nb_missing_hit(0), useMissingHits(false), max_road_number(1000000),
  sectors(st), eventsFileName(f), outputFileName(of),converter(NULL), hardware_limitations(false),d_detector(d),
  d_p_bank(p), d_parameters(dp)
{

  //we don't need the map of patterns, a vector will be enough and uses less memory
  sectors->getAllSectors()[0]->getPatternTree()->switchToVector();

  //we use 1024 threads and want to compute the number of blocks
  int nb_patterns = sectors->getAllSectors()[0]->getLDPatternNumber();
  nb_threads = 1024;
  int nbIter = nb_patterns/nb_threads/1024+1;//we want nb_blocks<1025
  nb_blocks = nb_patterns/nbIter/nb_threads+1;

  if(cudaSuccess !=cudaMemcpy(d_parameters->threshold,&active_threshold, sizeof(int),cudaMemcpyHostToDevice))
    cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  
  if(nb_patterns>0){
    int nbIter = nb_patterns/(nb_blocks*nb_threads)+1;
    if(nbIter==0)
      nbIter = 1;

    if(cudaSuccess !=cudaMemcpy(d_parameters->iter,&nbIter, sizeof(int), cudaMemcpyHostToDevice))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;

    if(cudaSuccess !=cudaMemcpy(d_parameters->nbPatterns,&nb_patterns, sizeof(int), cudaMemcpyHostToDevice))
      cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  }

}
#endif

void PatternFinder::setSectorTree(SectorTree* s){
  sectors = s;
}

void PatternFinder::setMaxRoadNumber(unsigned int m){
  max_road_number=m;
}

void PatternFinder::setHardwareLimitations(bool b){
  hardware_limitations = b;
  if(hardware_limitations)
    Detector::setHWPatternLimitations(HW_LIMIT_PATTERN_LAYER_STUBS);
  else
    Detector::setHWPatternLimitations(0);
}

void PatternFinder::setEventsFile(string& f){
  eventsFileName = f;
}

#ifdef IPNL_USE_CUDA
void PatternFinder::findCuda(int start, int& stop, deviceStubs* d_stubs){

  /**************** OUTPUT FILE ****************/

  TTree* Out = new TTree("Patterns", "Active patterns");
  TTree* SectorOut = new TTree("Sectors", "Used Sectors");
  TFile *t = new TFile(outputFileName.c_str(),"recreate");

  const int MAX_NB_PATTERNS = 1500;
  const int MAX_NB_HITS = 100;
  const int MAX_NB_LADDERS_PER_LAYER = 16;
  const int MAX_NB_LAYERS = 9;
  
  int nb_layers;
  int nb_patterns=0;
  int ori_nb_stubs=0;
  int sel_nb_stubs=0;
  int nb_tc=0;
  int event_id;
  int superStrip_layer_0[MAX_NB_PATTERNS];
  int superStrip_layer_1[MAX_NB_PATTERNS];
  int superStrip_layer_2[MAX_NB_PATTERNS];
  int superStrip_layer_3[MAX_NB_PATTERNS];
  int superStrip_layer_4[MAX_NB_PATTERNS];
  int superStrip_layer_5[MAX_NB_PATTERNS];
  int superStrip_layer_6[MAX_NB_PATTERNS];
  int superStrip_layer_7[MAX_NB_PATTERNS];
  int superStrip_layer_8[MAX_NB_PATTERNS];
  float pattern_pt[MAX_NB_PATTERNS];
  int pattern_sector_id[MAX_NB_PATTERNS];
  
  //Array containing the strips arrays
  int* superStrips[MAX_NB_LAYERS];
  superStrips[0]=superStrip_layer_0;
  superStrips[1]=superStrip_layer_1;
  superStrips[2]=superStrip_layer_2;
  superStrips[3]=superStrip_layer_3;
  superStrips[4]=superStrip_layer_4;
  superStrips[5]=superStrip_layer_5;
  superStrips[6]=superStrip_layer_6;
  superStrips[7]=superStrip_layer_7;
  superStrips[8]=superStrip_layer_8;

  int sector_id=0;
  int sector_layers=0;
  int nb_ladders_layer[MAX_NB_LAYERS];
  int sector_layer_list[MAX_NB_LAYERS];
  int sector_layer_0[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_1[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_2[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_3[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_4[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_5[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_6[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_7[MAX_NB_LADDERS_PER_LAYER];
  int sector_layer_8[MAX_NB_LADDERS_PER_LAYER];

  int* sector_layers_detail[MAX_NB_LAYERS];
  sector_layers_detail[0]=sector_layer_0;
  sector_layers_detail[1]=sector_layer_1;
  sector_layers_detail[2]=sector_layer_2;
  sector_layers_detail[3]=sector_layer_3;
  sector_layers_detail[4]=sector_layer_4;
  sector_layers_detail[5]=sector_layer_5;
  sector_layers_detail[6]=sector_layer_6;
  sector_layers_detail[7]=sector_layer_7;
  sector_layers_detail[8]=sector_layer_8;

  int nbHitPerPattern[MAX_NB_PATTERNS];

  float* track_pt = new float[MAX_NB_PATTERNS];
  float* track_phi = new float[MAX_NB_PATTERNS];
  float* track_d0 = new float[MAX_NB_PATTERNS];
  float* track_eta = new float[MAX_NB_PATTERNS];
  float* track_z0 = new float[MAX_NB_PATTERNS];

  int totalNbHits=0;
  short* hit_layer = new short[MAX_NB_PATTERNS*MAX_NB_HITS];
  short* hit_ladder = new short[MAX_NB_PATTERNS*MAX_NB_HITS];
  short* hit_zPos = new short[MAX_NB_PATTERNS*MAX_NB_HITS];
  short* hit_segment = new short[MAX_NB_PATTERNS*MAX_NB_HITS];
  short* hit_strip = new short[MAX_NB_PATTERNS*MAX_NB_HITS];
  int* hit_tp = new int[MAX_NB_PATTERNS*MAX_NB_HITS];
  int* hit_idx = new int[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_ptGEN = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_etaGEN = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_phi0GEN = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_ip = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_x = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_y = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_z = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_x0 = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_y0 = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_z0 = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  float* hit_bend = new float[MAX_NB_PATTERNS*MAX_NB_HITS];
  
  SectorOut->Branch("sectorID",            &sector_id);
  SectorOut->Branch("nbLayers",            &sector_layers);
  SectorOut->Branch("layer",                  sector_layer_list, "layer[nbLayers]/I");
  SectorOut->Branch("nb_ladders_layer",       nb_ladders_layer, "nb_ladders_layer[nbLayers]/I");
  SectorOut->Branch("sectorLadders_layer_0",  sector_layer_0, "sectorLadders_layer_0[16]/I");
  SectorOut->Branch("sectorLadders_layer_1",  sector_layer_1, "sectorLadders_layer_1[16]/I");
  SectorOut->Branch("sectorLadders_layer_2",  sector_layer_2, "sectorLadders_layer_2[16]/I");
  SectorOut->Branch("sectorLadders_layer_3",  sector_layer_3, "sectorLadders_layer_3[16]/I");
  SectorOut->Branch("sectorLadders_layer_4",  sector_layer_4, "sectorLadders_layer_4[16]/I");
  SectorOut->Branch("sectorLadders_layer_5",  sector_layer_5, "sectorLadders_layer_5[16]/I");
  SectorOut->Branch("sectorLadders_layer_6",  sector_layer_6, "sectorLadders_layer_6[16]/I");
  SectorOut->Branch("sectorLadders_layer_7",  sector_layer_7, "sectorLadders_layer_7[16]/I");
  SectorOut->Branch("sectorLadders_layer_8",  sector_layer_8, "sectorLadders_layer_8[16]/I");

  Out->Branch("nbLayers",            &nb_layers);
  Out->Branch("nbPatterns",          &nb_patterns);
  Out->Branch("nbStubsInEvt",        &ori_nb_stubs);
  Out->Branch("nbStubsInPat",        &sel_nb_stubs);

  Out->Branch("nbTracks",            &nb_tc);

  Out->Branch("eventID",             &event_id);
  Out->Branch("sectorID",            pattern_sector_id, "sectorID[nbPatterns]/I");
  Out->Branch("patternPT",           pattern_pt, "patternPT[nbPatterns]/F");
  Out->Branch("superStrip0",         superStrip_layer_0, "superStrip0[nbPatterns]/I");
  Out->Branch("superStrip1",         superStrip_layer_1, "superStrip1[nbPatterns]/I");
  Out->Branch("superStrip2",         superStrip_layer_2, "superStrip2[nbPatterns]/I");
  Out->Branch("superStrip3",         superStrip_layer_3, "superStrip3[nbPatterns]/I");
  Out->Branch("superStrip4",         superStrip_layer_4, "superStrip4[nbPatterns]/I");
  Out->Branch("superStrip5",         superStrip_layer_5, "superStrip5[nbPatterns]/I");
  Out->Branch("superStrip6",         superStrip_layer_6, "superStrip6[nbPatterns]/I");
  Out->Branch("superStrip7",         superStrip_layer_7, "superStrip7[nbPatterns]/I");
  Out->Branch("superStrip8",         superStrip_layer_8, "superStrip8[nbPatterns]/I");
  Out->Branch("nbStubs",             nbHitPerPattern, "nbStubs[nbPatterns]/I");
  
  Out->Branch("track_pt",             track_pt, "track_pt[nbTracks]/F");
  Out->Branch("track_phi",            track_phi, "track_phi[nbTracks]/F");
  Out->Branch("track_d0",             track_d0, "track_d0[nbTracks]/F");
  Out->Branch("track_eta",            track_eta, "track_eta[nbTracks]/F");
  Out->Branch("track_z0",             track_z0, "track_z0[nbTracks]/F");
  
  Out->Branch("total_nb_stubs",      &totalNbHits, "total_nb_stubs/I");
  Out->Branch("stub_layers",         hit_layer,"stub_layers[total_nb_stubs]/S");
  Out->Branch("stub_ladders",        hit_ladder, "stub_ladders[total_nb_stubs]/S");
  Out->Branch("stub_module",         hit_zPos, "stub_module[total_nb_stubs]/S");
  Out->Branch("stub_segment",        hit_segment, "stub_segment[total_nb_stubs]/S");
  Out->Branch("stub_strip",          hit_strip, "stub_strip[total_nb_stubs]/S");
  Out->Branch("stub_tp",             hit_tp,    "stub_tp[total_nb_stubs]/I");
  Out->Branch("stub_idx",            hit_idx,    "stub_idx[total_nb_stubs]/I");
  Out->Branch("stub_ptGEN",          hit_ptGEN, "stub_ptGEN[total_nb_stubs]/F");
  Out->Branch("stub_etaGEN",         hit_etaGEN, "stub_etaGEN[total_nb_stubs]/F");
  Out->Branch("stub_phi0GEN",        hit_phi0GEN, "stub_phi0GEN[total_nb_stubs]/F");
  Out->Branch("stub_IP",             hit_ip, "stub_IP[total_nb_stubs]/F");
  Out->Branch("stub_x",              hit_x, "stub_x[total_nb_stubs]/F");
  Out->Branch("stub_y",              hit_y, "stub_y[total_nb_stubs]/F");
  Out->Branch("stub_z",              hit_z, "stub_z[total_nb_stubs]/F");
  Out->Branch("stub_X0",             hit_x0, "stub_X0[total_nb_stubs]/F");
  Out->Branch("stub_Y0",             hit_y0, "stub_Y0[total_nb_stubs]/F");
  Out->Branch("stub_Z0",             hit_z0, "stub_Z0[total_nb_stubs]/F");
  Out->Branch("stub_bend",           hit_bend, "stub_bend[total_nb_stubs]/F");

  /*********************************************/

  /******************* SAVING SECTORS **************/

  vector<Sector*> all_sectors = sectors->getAllSectors();
  for(unsigned int i=0;i<all_sectors.size();i++){
    Sector* tmpSec = all_sectors[i];
    sector_id=tmpSec->getOfficialID();
    if(sector_id==-1)
      sector_id=tmpSec->getKey();
    sector_layers = tmpSec->getNbLayers();
    for(int j=0;j<sector_layers;j++){
      vector<int> sec_l = tmpSec->getLadders(j);
      sector_layer_list[j] = tmpSec->getLayerID(j);
      nb_ladders_layer[j] = sec_l.size();
      for(unsigned int l=0;l<sec_l.size();l++){
	sector_layers_detail[j][l]=sec_l[l];
      }

    }
    cout<<tmpSec->getIDString()<<" -> "<<sector_id<<endl;
    SectorOut->Fill();
  }
  SectorOut->Write();
  delete SectorOut;

  /*************************************************/

  /***************** INPUT FILE ****************/

  TChain* TT = new TChain("TkStubs");
  TT->Add(eventsFileName.c_str());
  
  int               n_evt;

  int m_stub;

  vector<int>           m_stub_layer;  // Layer du stub (5 a 10 pour les 6 layers qui nous interessent)
  vector<int>           m_stub_module; // Position en Z du module contenant le stub
  vector<int>           m_stub_ladder; // Position en PHI du module contenant le stub
  vector<int>           m_stub_seg;    // Segment du module contenant le stub
  vector<float>         m_stub_strip;  // Strip du cluster interne du stub
  vector<int>           m_stub_tp;     // particule du stub
  vector<float>         m_stub_px_gen; // pt initial de la particule ayant genere le stub
  vector<float>         m_stub_py_gen; // pt initial de la particule ayant genere le stub
  vector<float>         m_stub_x0;     // utilise pour calculer la distance au point d'interaction
  vector<float>         m_stub_y0;     // utilise pour calculer la distance au point d'interaction
  vector<float>         m_stub_z0;
  vector<float>         m_stub_phi0;
  vector<float>         m_stub_eta_gen;
  vector<float>         m_stub_x;      // x coordinate of the hit
  vector<float>         m_stub_y;      // y coordinate of the hit
  vector<float>         m_stub_z;      // z coordinate of the hit

  vector<int>           *p_m_stub_layer =  &m_stub_layer;
  vector<int>           *p_m_stub_module = &m_stub_module;
  vector<int>           *p_m_stub_ladder = &m_stub_ladder;
  vector<int>           *p_m_stub_seg =    &m_stub_seg;
  vector<float>         *p_m_stub_strip =  &m_stub_strip;
  vector<int>           *p_m_stub_tp =     &m_stub_tp;
  vector<float>         *p_m_stub_pxGEN = &m_stub_px_gen;  
  vector<float>         *p_m_stub_pyGEN = &m_stub_py_gen;  
  vector<float>         *p_m_stub_x0 =     &m_stub_x0;
  vector<float>         *p_m_stub_y0 =     &m_stub_y0;
  vector<float>         *p_m_stub_z0 =     &m_stub_z0;
  vector<float>         *p_m_stub_phi0 =   &m_stub_phi0;
  vector<float>         *p_m_stub_etaGEN = &m_stub_eta_gen;
  vector<float>         *p_m_stub_x =      &m_stub_x;
  vector<float>         *p_m_stub_y =      &m_stub_y;
  vector<float>         *p_m_stub_z =      &m_stub_z;


  TT->SetBranchAddress("L1Tkevt",            &n_evt);
  TT->SetBranchAddress("L1TkSTUB_n",         &m_stub);
  TT->SetBranchAddress("L1TkSTUB_layer",     &p_m_stub_layer);
  TT->SetBranchAddress("L1TkSTUB_module",    &p_m_stub_module);
  TT->SetBranchAddress("L1TkSTUB_ladder",    &p_m_stub_ladder);
  TT->SetBranchAddress("L1TkSTUB_seg",       &p_m_stub_seg);
  TT->SetBranchAddress("L1TkSTUB_strip",     &p_m_stub_strip);
  TT->SetBranchAddress("L1TkSTUB_tp",        &p_m_stub_tp);
  TT->SetBranchAddress("L1TkSTUB_X0",        &p_m_stub_x0);
  TT->SetBranchAddress("L1TkSTUB_Y0",        &p_m_stub_y0);
  TT->SetBranchAddress("L1TkSTUB_Z0",        &p_m_stub_z0);
  TT->SetBranchAddress("L1TkSTUB_PHI0",      &p_m_stub_phi0);
  TT->SetBranchAddress("L1TkSTUB_etaGEN",    &p_m_stub_etaGEN);
  TT->SetBranchAddress("L1TkSTUB_pxGEN",     &p_m_stub_pxGEN);
  TT->SetBranchAddress("L1TkSTUB_pyGEN",     &p_m_stub_pyGEN);
  TT->SetBranchAddress("L1TkSTUB_x",         &p_m_stub_x);
  TT->SetBranchAddress("L1TkSTUB_y",         &p_m_stub_y);
  TT->SetBranchAddress("L1TkSTUB_z",         &p_m_stub_z);

  /*******************************************************/

  int n_entries_TT = TT->GetEntries();
  int num_evt = start;
  if(stop>n_entries_TT){
    stop=n_entries_TT-1;
    cout<<"Last event index too high : reset to "<<stop<<endl;
  }

  char* cuda_hits = new char[5000];
  int cuda_nb_hits = 0;


  while(num_evt<n_entries_TT && num_evt<=stop){
    TT->GetEntry(num_evt);

    cout<<"Event "<<n_evt<<endl;
    cuda_nb_hits = 0;
    vector<Hit*> hits;

    for(int i=0;i<m_stub;i++){
      int layer = m_stub_layer[i];
      int module = CMSPatternLayer::getModuleCode(layer, m_stub_module[i]);
      if(module<0)
	continue;
      int ladder = CMSPatternLayer::getLadderCode(layer, m_stub_ladder[i]);
      int segment =  CMSPatternLayer::getSegmentCode(layer, ladder, m_stub_seg[i]);
      if(segment<0 || segment>1){
	cout<<"Invalid segment on event "<<n_evt<<endl;
	continue;
      }
      int strip = m_stub_strip[i];
      int tp = m_stub_tp[i];
      float eta = m_stub_eta_gen[i];
      float phi0 = m_stub_phi0[i];
      float spt = sqrt(m_stub_px_gen[i]*m_stub_px_gen[i]+m_stub_py_gen[i]*m_stub_py_gen[i]);
      float x = m_stub_x[i];
      float y = m_stub_y[i];
      float z = m_stub_z[i];
      float x0 = m_stub_x0[i];
      float y0 = m_stub_y0[i];
      float z0 = m_stub_z0[i];
      
      float ip = sqrt(m_stub_x0[i]*m_stub_x0[i]+m_stub_y0[i]*m_stub_y0[i]);

      Hit* h = new Hit(layer,ladder, module, segment, strip, i, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0);

      if(sectors->getSector(*h)!=NULL){
	hits.push_back(h);

	int cuda_idx = cuda_nb_hits*CUDA_STUB_SIZE;
	cuda_hits[cuda_idx]=cuda_layer_index[layer];
	cuda_hits[cuda_idx+1]=ladder;
	cuda_hits[cuda_idx+2]=module;
	cuda_hits[cuda_idx+3]=segment;
	//cuda_hits[cuda_idx+4]=(char)(strip/sectors->getSuperStripSize(layer));
	cuda_nb_hits++;

      }
      else
	delete(h);
    }

    cudaCopyStubs(cuda_hits,d_stubs,cuda_nb_hits); 
    nb_patterns = findCuda(cuda_nb_hits, d_stubs);
    
    bool* active_stubs = new bool[cuda_nb_hits];
    cudaGetActiveStubs(active_stubs,d_stubs,&cuda_nb_hits);
    
    cout<<"Selected stubs : "<<endl;
    int total_nb_stubs = 0;
    for(int i=0;i<cuda_nb_hits;i++){
      if(active_stubs[i]){
	cout<<*(hits[i])<<endl;
	total_nb_stubs++;
      }
    }
    cout<<"Total nb stubs : "<<total_nb_stubs<<endl;

    //Traitement des patterns actif : enregistrement, affichage...
    nb_layers = tracker.getNbLayers();
    event_id=num_evt;//we use the index in the file as event_id (most of our input files do not have a valid event_id)
    ori_nb_stubs = (int)hits.size();
    
    nb_tc = 0;
    for(int i=0;i<nb_layers;i++){
      memset(superStrips[i],0,MAX_NB_PATTERNS*sizeof(int));
    }
    memset(hit_layer,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(short));
    memset(hit_ladder,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(short));
    memset(hit_tp,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(int));
    memset(hit_idx,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(int));
    memset(hit_ptGEN,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_etaGEN,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_phi0GEN,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_ip,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_x,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_y,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_z,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_x0,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_y0,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_z0,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));
    memset(hit_bend,0,MAX_NB_PATTERNS*MAX_NB_HITS*sizeof(float));

    int stubIndex = 0;
    totalNbHits = 0;

    for(unsigned int i=0;i<hits.size();i++){
      if(active_stubs[i]){
	//stubs of the patterns
	  
	hit_layer[stubIndex]=hits[i]->getLayer();
	hit_ladder[stubIndex]=hits[i]->getLadder();
	hit_zPos[stubIndex]=hits[i]->getModule();
	hit_segment[stubIndex]=hits[i]->getSegment();
	hit_strip[stubIndex]=hits[i]->getStripNumber();
	hit_tp[stubIndex]=hits[i]->getParticuleID();
	hit_idx[stubIndex]=hits[i]->getID();
	hit_ptGEN[stubIndex]=hits[i]->getParticulePT();
	hit_etaGEN[stubIndex]=hits[i]->getParticuleETA();
	hit_phi0GEN[stubIndex]=hits[i]->getParticulePHI0();
	hit_ip[stubIndex]=hits[i]->getParticuleIP();
	hit_x[stubIndex]=hits[i]->getX();
	hit_y[stubIndex]=hits[i]->getY();
	hit_z[stubIndex]=hits[i]->getZ();
	hit_x0[stubIndex]=hits[i]->getX0();
	hit_y0[stubIndex]=hits[i]->getY0();
	hit_z0[stubIndex]=hits[i]->getZ0();
	hit_bend[stubIndex]=hits[i]->getBend();
	
	stubIndex++;
      }
    }
    totalNbHits=stubIndex;
    sel_nb_stubs = stubIndex;

    Out->Fill();

    //////////////////////////////////
    for(unsigned int i=0;i<hits.size();i++){
      delete hits[i];
    }

    delete [] active_stubs;

    num_evt++;

  }

  delete[] cuda_hits;

  //Out->Print();
  Out->Write();
  t->Close();
  delete Out;
  delete t;

  delete TT;

  delete[] track_pt;
  delete[] track_phi;
  delete[] track_d0;
  delete[] track_eta;
  delete[] track_z0;
  delete[] hit_layer;
  delete[] hit_ladder;
  delete[] hit_zPos;
  delete[] hit_segment;
  delete[] hit_strip;
  delete[] hit_tp;
  delete[]  hit_idx;
  delete[]  hit_ptGEN;
  delete[]  hit_etaGEN;
  delete[]  hit_phi0GEN;
  delete[]  hit_ip;
  delete[]  hit_x;
  delete[]  hit_y;
  delete[]  hit_z;
  delete[]  hit_x0;
  delete[]  hit_y0;
  delete[] hit_z0;
  delete[] hit_bend;
}
#endif

vector<Sector*> PatternFinder::find(vector<Hit*> hits){
  tracker.clear();
  for(unsigned int i=0;i<hits.size();i++){
    tracker.receiveHit(*hits[i]);
  }
  if(useMissingHits){
    return sectors->getActivePatternsPerSectorUsingMissingHit(max_nb_missing_hit, active_threshold, max_road_number);
  }
  else{
    return sectors->getActivePatternsPerSector(active_threshold, max_road_number);
  }
}

#ifdef IPNL_USE_CUDA
int PatternFinder::findCuda(int nb, deviceStubs* d_stubs, cudaStream_t* stream){
  cudaSetHitsWrapper(d_stubs,nb,d_detector, stream);
  int res = cudaGetActivePatternsWrapper(d_detector, d_p_bank, d_stubs, d_parameters, active_threshold, nb_blocks, nb_threads, stream);
  resetDetector(d_detector, stream);
  return res;
}
#endif

void PatternFinder::useMissingHitThreshold(int max_nb_missing_hit){
  useMissingHits=true;
  this->max_nb_missing_hit = max_nb_missing_hit;
}

void PatternFinder::setVerboseMode(bool m){
  tracker.setVerboseMode(m);
}
