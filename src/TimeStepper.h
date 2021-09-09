/*
 *                 Copyright (C) 2017, UChicago Argonne, LLC
 *              as Operator of Argonne National Laboratory under
 *  Prime Contract No. DE-AC02-06CH11357 with the U.S. Department of Energy
 *                            All Rights Reserved
 *
 *       Hardware/Hybrid Accelerated Cosmology Code (HACC), Version 1.0
 *
 * Salman Habib, Adrian Pope, Hal Finkel, Nicholas Frontiere, Katrin Heitmann,
 *     Vitali Morozov, Jeffrey Emberson, Thomas Uram, Esteban Rangel
 *                       (Argonne National Laboratory)
 *
 * David Daniel, Patricia Fasel, Chung-Hsing Hsu, Zarija Lukic, James Ahrens
 *                     (Los Alamos National Laboratory)
 *
 *                              George Zagaris
 *                                (Kitware)
 *
 *                           OPEN SOURCE LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * ****************************************************************************
 *
 *                               DISCLAIMER
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * ****************************************************************************
 */

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
