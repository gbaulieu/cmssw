#include "../interface/PatternTree.h"

PatternTree::PatternTree(){
}

PatternTree::~PatternTree(){
  clear();
}

void PatternTree::clear(){
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      delete (itr->second);
    }
    patterns.clear();
  }
  if(v_patterns.size()!=0){
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      delete (*itr);
    }
    v_patterns.clear();
  }
}

void PatternTree::addPattern(Pattern* ldp, Pattern* fdp){
  if(patterns.size()==0)
    switchToMap();
  string key = ldp->getKey();
  map<string, PatternTrunk*>::iterator it = patterns.find(key);
  if(it==patterns.end()){//not found
    PatternTrunk* pt = new PatternTrunk(ldp);
    pt->addFDPattern(fdp);
    patterns[key]=pt;
  }
  else{
    (it->second)->addFDPattern(fdp);
  }
}

void PatternTree::addPattern(Pattern* ldp, Pattern* fdp, float new_pt, int new_pdg){
  if(patterns.size()==0)
    switchToMap();
  string key = ldp->getKey();
  map<string, PatternTrunk*>::iterator it = patterns.find(key);
  if(it==patterns.end()){//not found
    PatternTrunk* pt = new PatternTrunk(ldp);
    pt->addFDPattern(fdp, new_pt, new_pdg);
    patterns[key]=pt;
  }
  else{
    (it->second)->addFDPattern(fdp, new_pt, new_pdg);
  }
}

void PatternTree::switchToVector(){
  for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
    v_patterns.push_back(itr->second);
  }
  patterns.clear();
}

void PatternTree::switchToMap(){
  for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
    GradedPattern* gp = (*itr)->getLDPattern();
    string key = gp->getKey();
    patterns[key]=*itr;
    delete gp;
  }
  v_patterns.clear();
}

vector<GradedPattern*> PatternTree::getFDPatterns(){
  vector<GradedPattern*> res;
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      vector<GradedPattern*> fdp = itr->second->getFDPatterns();
      res.insert(res.end(), fdp.begin(), fdp.end());
    }
  }
  if(v_patterns.size()!=0){
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      vector<GradedPattern*> fdp = (*itr)->getFDPatterns();
      res.insert(res.end(),fdp.begin(),fdp.end());
    }
  }
  return res;
}

vector<GradedPattern*> PatternTree::getLDPatterns(){
  vector<GradedPattern*> res;
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      GradedPattern* ldp = itr->second->getLDPattern();
      res.push_back(ldp);
    }
  }
  if(v_patterns.size()!=0){
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      GradedPattern* ldp = (*itr)->getLDPattern();
      res.push_back(ldp);
    }
  }
  return res;
}

vector<int> PatternTree::getPTHisto(){
  vector<int> h;
  for(int i=0;i<150;i++){
    h.push_back(0);
  }
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      float pt = itr->second->getLDPatternPT();
      if(pt>149)
	pt=149;
      if(pt<0)
	pt=0;
      h[(int)pt]=h[(int)pt]+1;
    }
  }
  if(v_patterns.size()!=0){
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      float pt = (*itr)->getLDPatternPT();
      if(pt>149)
	pt=149;
      if(pt<0)
	pt=0;
      h[(int)pt]=h[(int)pt]+1;
    }
  }
  return h;
}

int PatternTree::getFDPatternNumber(){
  vector<GradedPattern*> res;
  int num=0;
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      num += itr->second->getFDPatternNumber();
    }
  }
  if(v_patterns.size()!=0){
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      num += (*itr)->getFDPatternNumber();
    }
  }
  return num;
}

int PatternTree::getLDPatternNumber(){
  if(patterns.size()!=0)
    return patterns.size();
  else
    return v_patterns.size();
}

void PatternTree::computeAdaptativePatterns(vector<int> r){
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      itr->second->computeAdaptativePattern(r);
    }
  }
  else{
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      (*itr)->computeAdaptativePattern(r);
    }
  }
}

void PatternTree::link(Detector& d){
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      itr->second->link(d);
    }
  }
  else{
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      (*itr)->link(d);
    }
  }
}

#ifdef IPNL_USE_CUDA
void PatternTree::linkCuda(patternBank* p, deviceDetector* d, const vector< vector<int> >& sec, const vector<map<int, vector<int> > >& modules, vector<int> layers){
  int counter=0;
  unsigned int* cache = new unsigned int[PATTERN_LAYERS*PATTERN_SSTRIPS];
  if(patterns.size()!=0){
    cudaSetNbPatterns(p, patterns.size());
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      itr->second->linkCuda(p, d, counter, sec, modules, layers, cache);
      counter++;
      if(counter%10000==0){
	cout<<counter*100/patterns.size()<<"%\r";
	cout.flush();
      }
    }
  }
  else{
    cudaSetNbPatterns(p, v_patterns.size());
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      (*itr)->linkCuda(p, d, counter, sec, modules, layers, cache);
      counter++;
      if(counter%10000==0){
	cout<<counter*100/v_patterns.size()<<"%\r";
	cout.flush();
      }
    }
  }
  delete[] cache;
}
#endif

void PatternTree::getActivePatterns(int active_threshold, vector<GradedPattern*>& active_patterns){
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      GradedPattern* p = itr->second->getActivePattern(active_threshold);
      if(p!=NULL)
	active_patterns.push_back(p);
    }
  }
  else{
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      GradedPattern* p = (*itr)->getActivePattern(active_threshold);
      if(p!=NULL)
	active_patterns.push_back(p);
    }
  }
}

void PatternTree::getActivePatternsUsingMissingHit(int max_nb_missing_hit, int active_threshold, vector<GradedPattern*>& active_patterns){
  if(patterns.size()!=0){
    for(map<string, PatternTrunk*>::iterator itr = patterns.begin(); itr != patterns.end(); ++itr){
      GradedPattern* p = itr->second->getActivePatternUsingMissingHit(max_nb_missing_hit, active_threshold);
      if(p!=NULL)
	active_patterns.push_back(p);
    }
  }
  else{
    for(vector<PatternTrunk*>::iterator itr = v_patterns.begin(); itr != v_patterns.end(); ++itr){
      GradedPattern* p = (*itr)->getActivePatternUsingMissingHit(max_nb_missing_hit, active_threshold);
      if(p!=NULL)
	active_patterns.push_back(p);
    }
  }
}

void PatternTree::addPatternsFromTree(PatternTree* p){
  if(patterns.size()==0)
    switchToMap();
  vector<GradedPattern*> ld = p->getLDPatterns();
  for(unsigned int i=0;i<ld.size();i++){
    GradedPattern* patt = ld[i];

    addPatternForMerging(patt);
    
    delete patt;
  }
}

void PatternTree::addPatternForMerging(GradedPattern* ldp){
  if(patterns.size()==0)
    switchToMap();
  string key = ldp->getKey();
  map<string, PatternTrunk*>::iterator it = patterns.find(key);
  if(it==patterns.end()){//not found
    PatternTrunk* pt = new PatternTrunk(ldp);
    patterns[key]=pt;
  }
  else{
    (it->second)->updateWithPattern(ldp);
  }
}

bool PatternTree::checkPattern(Pattern* lp, Pattern* hp){
  if(patterns.size()==0)
    switchToMap();
  if(lp==NULL || hp==NULL)
    return false;
  string key = lp->getKey();
  map<string, PatternTrunk*>::iterator it = patterns.find(key);
  if(it==patterns.end()){//not found
    return false;
  }
  else{
    return (it->second)->checkPattern(hp);
  }
}

bool comparePatterns(PatternTrunk* p1, PatternTrunk* p2){
  if(p1->getLDPatternGrade()==p2->getLDPatternGrade())
    return p1->getLDPatternPT()>p2->getLDPatternPT();
  else
    return p1->getLDPatternGrade()>p2->getLDPatternGrade();
}

bool comparePatternsbyPT(PatternTrunk* p1, PatternTrunk* p2){
  if(p1->getLDPatternPT()==p2->getLDPatternPT())
    return p1->getLDPatternGrade()>p2->getLDPatternGrade();
  else
    return p1->getLDPatternPT()>p2->getLDPatternPT();
}

struct ScoreComparer {
  ScoreComparer(float minPT, float pt_span, int minGrade, float grade_span) { this->min_pt = minPT; this->min_grade = minGrade; this->pt_span = pt_span; this->grade_span=grade_span;}
  bool operator () (PatternTrunk* p1, PatternTrunk* p2) {
    float p1_score = (((p1->getLDPatternPT()-min_pt)/pt_span)*100+((p1->getLDPatternGrade()-min_grade)/grade_span)*100);
    float p2_score = (((p2->getLDPatternPT()-min_pt)/pt_span)*100+((p2->getLDPatternGrade()-min_grade)/grade_span)*100);
    return p1_score>p2_score;;
  }

  float min_pt;
  int min_grade;
  float pt_span;
  float grade_span;
};

void PatternTree::desactivateModules(PatternTree* ref_pt, set<unsigned int> defective_modules){
  clear();
  ref_pt->switchToVector();
  CMSPatternLayer fakeSuperstrip;
  fakeSuperstrip.computeSuperstrip(0,0,0,0,0,0,false,true);

  for(unsigned int i=0;i<ref_pt->v_patterns.size();i++){
    //If we have dead modules : we replace the corresponding superstrips by fake ones
    if(defective_modules.size()>0){
      int nbLayers = ref_pt->v_patterns[i]->getLDPattern()->getNbLayers();
      GradedPattern* nPattern = ref_pt->v_patterns[i]->getLDPattern();
      // Loop over the layers
      for(int k=0;k<nbLayers;k++){
	CMSPatternLayer* pl = (CMSPatternLayer*)nPattern->getLayerStrip(k);
	// If the superstrip is already a fake one : nothing to do
	if(!pl->isFake()){
	  unsigned int module = k*100000+pl->getPhi()*1000+pl->getModule()*10+pl->getSegment();
	  if(defective_modules.find(module) != defective_modules.end()){
	    nPattern->setLayerStrip(k,&fakeSuperstrip);
	  }
	}
      }
      addPatternForMerging(nPattern);
      delete nPattern;
    }
  }
}

void PatternTree::truncate(int nbPatterns, int sorting_algo, vector<unsigned int> defective_addresses){
  switchToVector();
  switch(sorting_algo){
  case 0:
    sort(v_patterns.begin(),v_patterns.end(), comparePatterns);//sort by popularity then PT
    break;
  case 1:
    sort(v_patterns.begin(),v_patterns.end(), comparePatternsbyPT);//sort by PT then popularity
    break;
  case 2:
    // Sort by score
    float min_pt=200;
    int min_grade=5000;
    float max_pt=0;
    int max_grade=0;
    float pt_span=0;
    int grade_span=0;
    for(unsigned int i=0;i<v_patterns.size();i++){
      int current_grade = v_patterns[i]->getLDPatternGrade();
      float current_pt = v_patterns[i]->getLDPatternPT();
      if(current_pt<min_pt)
	min_pt=current_pt;
      if(current_pt>max_pt)
	max_pt=current_pt;
      if(current_grade<min_grade)
	min_grade=current_grade;
      if(current_grade>max_grade)
	max_grade=current_grade;
    }
    pt_span=max_pt-min_pt;
    grade_span=max_grade-min_grade;
    ScoreComparer sp(min_pt, pt_span, min_grade, grade_span);
    sort(v_patterns.begin(),v_patterns.end(), ScoreComparer(min_pt, pt_span, min_grade, grade_span));//sort by a score computed with PT and Popularity
    break;
  }

  if(nbPatterns>0){
    switch(sorting_algo){
    case 0:
      cout<<"Popularity ranging from  : "<<v_patterns[0]->getLDPatternGrade()<<" to "<<v_patterns[v_patterns.size()-1]->getLDPatternGrade()<<endl;
      break;
    case 1:
      cout<<"PT ranging from  : "<<v_patterns[0]->getLDPatternPT()<<" to "<<v_patterns[v_patterns.size()-1]->getLDPatternPT()<<endl;
      break;
    case 2:
      cout<<"Popularity/PT ranging from  : "<<v_patterns[0]->getLDPatternGrade()<<"/"<<v_patterns[0]->getLDPatternPT()<<" to "<<v_patterns[v_patterns.size()-1]->getLDPatternGrade()<<"/"<<v_patterns[v_patterns.size()-1]->getLDPatternPT()<<endl;
      break;
    }
    if(defective_addresses.size()>2048){
      cout<<"Too many defective addresses in the chip : can not handle more than 2048!"<<endl;
      cout<<"DEFECTIVE ADDRESSES WILL BE IGNORED!"<<endl;
      defective_addresses.clear();
    }

    for(unsigned int k=0;k<defective_addresses.size();k++){
      //Creating an invalid pattern to fill the defective addresses of the chip
      int nbLayers = v_patterns[0]->getLDPattern()->getNbLayers();
      Pattern* pat = new GradedPattern(nbLayers);
      for(int i=0;i<nbLayers;i++){
	CMSPatternLayer pl;
	pl.computeSuperstrip(5, 31, k/128, k%128, 0, 1, false, false);
	int nb_dc = v_patterns[0]->getLDPattern()->getLayerStrip(i)->getDCBitsNumber();
	for(int j=0;j<nb_dc;j++){
	  pl.setDC(j,4);
	}
	pat->setLayerStrip(i,&pl);
      }
      PatternTrunk* invalid_pattern = new PatternTrunk(pat);
      delete pat;
      
      vector<PatternTrunk*>::iterator it;
      it = v_patterns.begin()+defective_addresses[k];
      v_patterns.insert(it,invalid_pattern);
    }

    int nbToDelete = v_patterns.size()-nbPatterns;
    for(int i=0;i<nbToDelete;i++){
      v_patterns.pop_back();
    }
    cout<<"Keep "<<v_patterns.size()<<" patterns with ";
    switch(sorting_algo){
    case 0:
      cout<<"popularity ranging from  : "<<v_patterns[0]->getLDPatternGrade()<<" to "<<v_patterns[v_patterns.size()-1]->getLDPatternGrade()<<endl;
      break;
    case 1:
      cout<<"PT ranging from  : "<<v_patterns[0]->getLDPatternPT()<<" to "<<v_patterns[v_patterns.size()-1]->getLDPatternPT()<<endl;
      break;
    case 2:
      cout<<"popularity/PT ranging from  : "<<v_patterns[0]->getLDPatternGrade()<<"/"<<v_patterns[0]->getLDPatternPT()<<" to "<<v_patterns[v_patterns.size()-1]->getLDPatternGrade()<<"/"<<v_patterns[v_patterns.size()-1]->getLDPatternPT()<<endl;
      break;
    }
  }

  for(unsigned int i=0;i<v_patterns.size();i++){
    v_patterns[i]->setOrderInChip(i);
  }

  switchToMap();
}

void PatternTree::removePatterns(int minFS, int maxFS){
  switchToVector();

  vector<PatternTrunk*>::iterator itr = v_patterns.begin();

  int min = minFS;
  int max = maxFS;

  auto it = std::remove_if (v_patterns.begin(), v_patterns.end(), [min,max](PatternTrunk* p){
      GradedPattern* ldp = p->getLDPattern();
      int nbFS = ldp->getNbFakeSuperstrips();
      delete ldp;
      if (nbFS>max || nbFS<min) {
        return true;
      }
      return false;
    });
  v_patterns.erase(it, v_patterns.end());

  switchToMap();
}
