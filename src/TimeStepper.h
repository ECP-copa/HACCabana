#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <cmath>

#define TS_FLOAT double

class TimeStepper {
 public:

  TimeStepper(TS_FLOAT alpha_, 
	      TS_FLOAT ain_, 
	      TS_FLOAT afin_,
	      int nsteps_, 
	      TS_FLOAT omega_matter_, 
	      TS_FLOAT omega_cdm_,
	      TS_FLOAT omega_baryon_,
	      TS_FLOAT omega_cb_,
	      TS_FLOAT omega_nu_,
	      TS_FLOAT omega_radiation_,
	      TS_FLOAT f_nu_massless_,
	      TS_FLOAT f_nu_massive_,
	      TS_FLOAT w_,
              TS_FLOAT wa_);
  ~TimeStepper();

  void advanceHalfStep();
  void advanceFullStep();
  void reverseHalfStep();
  void reverseFullStep();
  void reduceTimestep(int factor);
  void increaseTimestep(int factor);

  TS_FLOAT aa() { return m_aa; }
  TS_FLOAT pp() { return m_pp; }
  TS_FLOAT zz() { return m_zz; }
  TS_FLOAT alpha() { return m_alpha; }
  TS_FLOAT tau() { return m_tau; }
  TS_FLOAT tau2() { return m_tau2; }
  TS_FLOAT adot() { return m_adot; }
  TS_FLOAT H_ratio(); // H/H_0
  TS_FLOAT omega_matter() { return m_omega_matter; }
  TS_FLOAT omega_cdm() { return m_omega_cdm; }
  TS_FLOAT omega_baryon() { return m_omega_baryon; }
  TS_FLOAT omega_cb() { return m_omega_cb; }
  TS_FLOAT omega_nu() { return m_omega_nu; }
  TS_FLOAT omega_radiation() { return m_omega_radiation; }
  TS_FLOAT f_nu_massless() { return m_f_nu_massless; }
  TS_FLOAT f_nu_massive() { return m_f_nu_massive; }
  TS_FLOAT ain() { return m_ain; }
  TS_FLOAT afin() { return m_afin; }
  TS_FLOAT pin() { return m_pin; }
  TS_FLOAT pfin() { return m_pfin; }
  TS_FLOAT zin() { return m_zin; }
  TS_FLOAT zfin() { return m_zfin; }
  int nsteps() { return m_nsteps; }
  TS_FLOAT phiscal() { return m_phiscal; }
  TS_FLOAT fscal() { return m_fscal; }
  TS_FLOAT w() { return m_w; }
  TS_FLOAT wa() { return m_wa; }
  int currentHalfStep() { return m_currentHalfStep; }

 private:

  TimeStepper();
  TimeStepper( const TimeStepper& );
  TimeStepper& operator = (const TimeStepper& );

  TS_FLOAT omega_nu_massive(TS_FLOAT a);
  void set_adot();
  void set_scal();

  TS_FLOAT m_aa;
  TS_FLOAT m_pp;
  TS_FLOAT m_zz;
  TS_FLOAT m_alpha;
  TS_FLOAT m_tau;
  TS_FLOAT m_tau2;
  TS_FLOAT m_adot;

  TS_FLOAT m_omega_matter;
  TS_FLOAT m_omega_cdm;
  TS_FLOAT m_omega_baryon;
  TS_FLOAT m_omega_cb;
  TS_FLOAT m_omega_nu;
  TS_FLOAT m_omega_radiation;
  TS_FLOAT m_f_nu_massless;
  TS_FLOAT m_f_nu_massive;

  TS_FLOAT m_ain;
  TS_FLOAT m_afin;
  TS_FLOAT m_pin;
  TS_FLOAT m_pfin;
  TS_FLOAT m_zin;
  TS_FLOAT m_zfin;
  int m_nsteps;

  TS_FLOAT m_phiscal;
  TS_FLOAT m_fscal;
  TS_FLOAT m_w;
  TS_FLOAT m_wa;

  int m_currentHalfStep;
};

#endif
