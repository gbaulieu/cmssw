#ifndef _TRACK_H_
#define _TRACK_H_

#include <vector>

using namespace std;

/**
   \brief A track is defined with its curve (PT), D0, Phi0, Eta0 and Z0 parameters
**/
class Track{

 private:
  double charge;
  double curve;
  double d0;
  double phi0;
  double eta0;
  double z0;
  double w_xy;
  double w_rz;
  double chi2;
  vector<int> stub_ids;
  int pattern_id;//optional pattern_id leading to this track

 public:
  /**
     \brief Default Constructor : all values set to 0
  **/
  Track();
  /**
     \brief Constructor
     \param c The PT of the track
     \param d The D0 of the track
     \param p The PHI0 of the track
     \param p_a The Eta0 of the track
     \param p_b The Z0 of the track
     \param cg The charge of the track
     \param Wxy The weight of the XY-retina maximum
     \param Wrz The weight of the RZ-retina maximum
     \param chi The chi2 parameter of the track
  **/
  Track(double c, double d, double p, double p_a, double p_b, double cg=0., double Wxy=-1., double Wrz=-1., double chi=-1);
  /**
     \brief Copy Constructor
  **/
  Track(const Track&);

  /**
     \brief Set the charge of the track
     \param cg The charge of the track
  **/
  void setCharge(double cg);
  /**
     \brief Get the charge of the track estimated from Fit Module
     \return The charge of the track
  **/
  double getCharge() const;
  /**
     \brief Set the PT of the track
     \param p The PT of the track
  **/
  void setCurve(double p);
  /**
     \brief Set the D0 of the track
     \param d The D0 of the track
  **/
  void setD0(double d);
  /**
     \brief Set the Phi of the track
     \param p The Phi of the track
  **/
  void setPhi0(double p);
  /**
     \brief Set the Eta of the track
     \param e The Eta of the track
  **/
  void setEta0(double e);
  /**
     \brief Set the Z0 of the track
     \param z The Z0 of the track
  **/
  void setZ0(double z);
  /**
     \brief Set the weight of the XY-retina maximum
     \param Wxy The weight of the XY-retina maximum
  **/
  void setWxy(double Wxy);
  /**
     \brief Set the weight of the RZ-retina maximum
     \param Wrz The weight of the RZ-retina maximum
  **/
  void setWrz(double Wrz);
 /**
     \brief Set the Chi2 of the track
     \param c The Chi2 value
  **/
  void setChi2(double c);

  /**
     \brief Add a stub to the list of stubs used to create the track
     \param s The ID of the stub
  **/
  void addStubIndex(int s);
  /**
     \brief Get the list of the index of stubs used to compute the track
     \return A vector with the list of index
  **/
  vector<int> getStubs() const;

  /**
     \brief CLear the list of stubs used to create the track
  **/
  void clearStubList();

  /**
     \brief Get the PT of the track
     \return The PT of the track
  **/
  double getCurve() const;
  /**
     \brief Get the D0 of the track
     \return The D0 of the track
  **/
  double getD0() const;
  /**
     \brief Get the Phi of the track
     \return The Phi of the track
  **/
  double getPhi0() const;
  /**
     \brief Get the Eta of the track
     \return The Eta of the track
  **/
  double getEta0() const;
  /**
     \brief Get the Z0 of the track
     \return The Z0 of the track
  **/
  double getZ0() const;
  /**
     \brief Get the weight of the XY-retina maximum
     \return The weight of the XY-retina maximum
  **/
  double getWxy() const;
  /**
     \brief Get the weight of the RZ-retina maximum
     \return The weight of the RZ-retina maximum
  **/
  double getWrz() const;
  /**
     \brief Get the Chi2 value
     \return The Chi2 value of the track
  **/
  double getChi2() const;

  /**
     \brief Set the ID of the pattern which led to this track.
     \param id The address of the pattern
  **/
  void setOriginPatternID(int id);

  /**
     \brief Get the ID of the pattern which led to this track.
     \return The ID of the pattern, -1 if unknown
  **/
  int getOriginPatternID() const;

};
#endif
