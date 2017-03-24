#include "../interface/Track.h"

Track::Track(){
  charge= 0.;
  curve = 0.;
  d0    = 0.;
  phi0  = 0.;
  eta0  = 0.;
  z0    = 0.;
  w_xy  = 0.;
  w_rz  = 0.;
  chi2  = -1;
  pattern_id=-1;
}

Track::Track(double c, double d, double p, double p_a, double p_b, double cg, double Wxy, double Wrz, double chi){
  charge= cg;
  curve = c;
  d0    = d;
  phi0  = p;
  eta0  = p_a;
  z0    = p_b;
  w_xy  = Wxy;
  w_rz  = Wrz;
  chi2 = chi;
  pattern_id = -1;
}

Track::Track(const Track& ref){
  charge= ref.charge;
  curve = ref.curve;
  d0    = ref.d0;
  phi0  = ref.phi0;
  eta0  = ref.eta0;
  z0    = ref.z0;
  w_xy  = ref.w_xy;
  w_rz  = ref.w_rz;
  for(unsigned int i=0;i<ref.stub_ids.size();i++){
    stub_ids.push_back(ref.stub_ids[i]);
  }
  pattern_id = ref.pattern_id;
}

void Track::setCharge(double cg){
  charge=cg;
}

void Track::setCurve(double c){
  curve=c;
}

void Track::setD0(double d){
  d0=d;
}

void Track::setPhi0(double p){
  phi0=p;
}

void Track::setEta0(double p_a){
  eta0=p_a;
}

void Track::setZ0(double p_b){
  z0=p_b;
}

void Track::setWxy(double Wxy){
  w_xy=Wxy;
}

void Track::setWrz(double Wrz){
  w_rz=Wrz;
}

void Track::setChi2(double c){
  chi2=c;
}

void Track::addStubIndex(int s){
  if(s>=0)
    stub_ids.push_back(s);
}

vector<int> Track::getStubs() const{
  return stub_ids;
}

void Track::clearStubList(){
  stub_ids.clear();
}

double Track::getCharge() const{
  return charge;
}

double Track::getCurve() const{
  return curve;
}

double Track::getD0() const{
  return d0;
}

double Track::getPhi0() const{
  return phi0;
}

double Track::getEta0() const{
  return eta0;
}

double Track::getZ0() const{
  return z0;
}

double Track::getWxy() const{
  return w_xy;
}

double Track::getWrz() const{
  return w_rz;
}

double Track::getChi2() const{
  return chi2;
}

void Track::setOriginPatternID(int id){
  if(id>=0)
    pattern_id = id;
}

int Track::getOriginPatternID() const{
  return pattern_id;
}
