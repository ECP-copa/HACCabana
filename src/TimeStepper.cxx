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

#include "TimeStepper.h"
#include "assert.h"

TimeStepper::TimeStepper(TS_FLOAT alpha_, 
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
                         TS_FLOAT wa_) :
  m_aa(-1.0),
  m_pp(-1.0),
  m_zz(-1.0),
  m_alpha(-1.0),
  m_tau(-1.0),
  m_tau2(-1.0),
  m_adot(-1.0),
  m_omega_matter(-1.0),
  m_omega_cdm(-1.0),
  m_omega_baryon(-1.0),
  m_omega_cb(-1.0),
  m_omega_nu(-1.0),
  m_omega_radiation(-1.0),
  m_f_nu_massless(-1.0),
  m_f_nu_massive(-1.0),
  m_ain(-1.0),
  m_afin(-1.0),
  m_pin(-1.0),
  m_pfin(-1.0),
  m_zin(-1.0),
  m_zfin(-1.0),
  m_nsteps(-1),
  m_phiscal(-1.0),
  m_fscal(-1.0),
  m_w(-1.0),
  m_wa(-1.0),
  m_currentHalfStep(0)
{
  m_alpha = alpha_;
  m_ain = ain_;
  m_afin = afin_;
  m_nsteps = nsteps_;
  m_omega_matter = omega_matter_;
  m_omega_cdm = omega_cdm_;
  m_omega_baryon = omega_baryon_;
  m_omega_cb = omega_cb_;
  m_omega_nu = omega_nu_;
  m_omega_radiation = omega_radiation_;
  m_f_nu_massless = f_nu_massless_;
  m_f_nu_massive = f_nu_massive_;
  m_w = w_;
  m_wa = wa_;

  m_pin = pow(m_ain, m_alpha);
  m_pfin = pow(m_afin, m_alpha);
  m_zin = 1.0/m_ain - 1.0;
  m_zfin = 1.0/m_afin - 1.0;
  m_tau = (m_pfin - m_pin)/(1.0*m_nsteps);
  m_tau2 = 0.5*m_tau;

  m_pp = m_pin;
  m_aa = m_ain;
  m_zz = m_zin;

  set_adot();
  set_scal();
}

TimeStepper::~TimeStepper() {}

TS_FLOAT TimeStepper::omega_nu_massive(TS_FLOAT a) {
  TS_FLOAT mat = m_omega_nu/pow(a,3.0);
  TS_FLOAT rad = m_f_nu_massive*m_omega_radiation/pow(a,4.0);
  return (mat>=rad)*mat + (rad>mat)*rad;
}

TS_FLOAT TimeStepper::H_ratio() {
  return sqrt(m_omega_cb/pow(m_aa,3.0f)
           + (1.0 + m_f_nu_massless)*m_omega_radiation/pow(m_aa,4.0f)
           + omega_nu_massive(m_aa)
           + (1.0 - m_omega_matter - (1.0+m_f_nu_massless)*m_omega_radiation)
           *pow(m_aa,(float)(-3.0*(1.0+m_w+m_wa)))*exp(-3.0*m_wa*(1.0-m_aa))
         );
}

void TimeStepper::set_adot() {//NOTE! Adot is in units of H0^-1, meaning adot = H0 adot_code. 
  // compared to normal equation for H(a) terms are multiplied by a^3
  TS_FLOAT pp1 = pow(m_aa, -3.0*(m_w+m_wa))*exp(-3.0*m_wa*(1.0-m_aa));
  TS_FLOAT tmp = m_omega_cb
    + (1.0+m_f_nu_massless)*m_omega_radiation/m_aa
    + omega_nu_massive(m_aa)*pow(m_aa,3.0)
    + (1.0-m_omega_matter-(1.0+m_f_nu_massless)*m_omega_radiation)*pp1;
  tmp /= m_aa;
  m_adot = sqrt(tmp);
  return;
}

// DOES THIS NEED TO BE CHANGED FOR W != -1?
void TimeStepper::set_scal() {
  set_adot();
  float dtdy = m_aa/(m_alpha*m_adot*m_pp);
  m_phiscal = 1.5*m_omega_cb;//Poisson equation is grad^2 phi = 3/2 omega_m (rho-1)
  m_fscal = m_phiscal*dtdy*(1.0/m_aa);//1/a from dp/dy = -gradphi/a, and dt/dy to time transform from t to y units. 
  return;
}

void TimeStepper::advanceHalfStep() {
  m_pp += m_tau2;
  m_aa = pow(m_pp, 1.0/m_alpha);
  m_zz = 1.0/m_aa - 1.0;
  set_adot();
  set_scal();
  m_currentHalfStep++;
  return;
}

void TimeStepper::reverseHalfStep() {
  m_pp -= m_tau2;
  m_aa = pow(m_pp, 1.0/m_alpha);
  m_zz = 1.0/m_aa - 1.0;
  set_adot();
  set_scal();
  m_currentHalfStep--;
  return;
}

void TimeStepper::advanceFullStep() {
  advanceHalfStep();
  advanceHalfStep();
  return;
}

void TimeStepper::reverseFullStep() {
  reverseHalfStep();
  reverseHalfStep();
  return;
}

void TimeStepper::reduceTimestep(int factor){
	assert(factor > 0);
	m_tau /= factor; //reduce dt by factor
	m_tau2 = 0.5*m_tau;//recalc half dt
	m_currentHalfStep *= factor;//increase current half step count by factor (it is what the half steps would be if the stepper was always at this reduced level)
}

void TimeStepper::increaseTimestep(int factor){
	assert(factor > 0);
        m_tau *= factor; //increase dt by factor
        m_tau2 = 0.5*m_tau;//recalc half dt
        m_currentHalfStep /= factor;//decrease current half step count by factor (it is what the half steps would be if the stepper was always at this increased level)
}
