#include "../interface/GradedPattern.h"

GradedPattern::GradedPattern():Pattern(){
  grade=0;
  averagePt=0;
  minPT=200;
  maxPT=0;
  charge = -2;
}

GradedPattern::GradedPattern(int nbLayers):Pattern(nbLayers){
  grade=0;
  averagePt=0;
  minPT=200;
  maxPT=0;
  charge = -2;
}

GradedPattern::GradedPattern(const Pattern& p):Pattern(p){
  grade=0;
  averagePt=0;
  minPT=200;
  maxPT=0;
  charge = -2;
}

GradedPattern::GradedPattern(const GradedPattern& p):Pattern(p){
  grade=p.getGrade();
  averagePt=p.getAveragePt();
  minPT=p.getMinPt();
  maxPT=p.getMaxPt();
  charge = p.getCharge();
}

int GradedPattern::getGrade() const{
  return grade;
}

float GradedPattern::getAveragePt() const{
  return averagePt;
}

float GradedPattern::getMinPt() const{
  return minPT;
}

float GradedPattern::getMaxPt() const{
  return maxPT;
}

int GradedPattern::getCharge() const{
  return charge;
}

void GradedPattern::increment(){
  grade++;
}

void GradedPattern::increment(float pt, int pdg){
  averagePt=(averagePt*grade+pt)/(grade+1);
  if(pt>maxPT)
    maxPT=pt;
  if(pt<minPT)
    minPT=pt;
  int partCharge = (pdg>0)?1:-1;
  if(charge==-2)//this is a new pattern
    charge=partCharge;
  else if(charge!=partCharge)
    charge = 0;
  grade++;
}

int GradedPattern::operator<(const GradedPattern& gp){
  return grade<gp.getGrade();
}

void GradedPattern::mergeData(const GradedPattern& gp){
  averagePt=(grade*averagePt+gp.grade*gp.averagePt)/(grade+gp.grade);
  if(gp.minPT<minPT)
    minPT=gp.minPT;
  if(gp.maxPT>maxPT)
    maxPT=gp.maxPT;
  grade+=gp.grade;
  if(gp.charge!=-2){
    if(charge==-2)
      charge=gp.charge;
    if(charge!=gp.charge)
      charge=0;
  }
}
