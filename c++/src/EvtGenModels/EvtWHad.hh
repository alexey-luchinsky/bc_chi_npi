//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtWHad.hh
//
// Description: Routine to calculate W -> (n pi) + (m K) current
//			according to [Kuhn, Was, Acta.Phys.Polon B39 (2008) 147]
//
// Modification history:
//	A V Luchinsky  	20 Jan, 2013	Module created
//
//------------------------------------------------------------------------

#ifndef EvtWHad_HH
#define EvtWHad_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"

#include <vector>

class EvtWHad {

public:

    EvtWHad();

    EvtVector4C WCurrent(const EvtVector4R& q1) const;

    EvtVector4C WCurrent(const EvtVector4R& q1, const EvtVector4R& q2) const;

    EvtVector4C WCurrent(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3) const;

    EvtVector4C WCurrent(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3,
			 const EvtVector4R& q4, const EvtVector4R& q5) const;

    EvtVector4C WCurrent_5pi(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3,
			 const EvtVector4R& q4, const EvtVector4R& q5) const;


    EvtVector4C WCurrent_KKP(const EvtVector4R& pKplus, const EvtVector4R& pKminus,
			     const EvtVector4R& pPiPlus) const;

    EvtVector4C WCurrent_KPP(const EvtVector4R& pKplus, const EvtVector4R& pPiPlus,
			     const EvtVector4R& pPiMinus) const;

    EvtVector4C WCurrent_KSK(const EvtVector4R& pKS, const EvtVector4R& pKplus) const;

    EvtVector4C WCurrent_KKPPP(const EvtVector4R& pKplus, const EvtVector4R& pKminus,
                               const EvtVector4R& pPi1Plus, const EvtVector4R& pPi2Plus, const EvtVector4R& pPiMinus) const;

  // 1=pi+ 2=pi+ 3=pi+ 4=pi+ 5=pi- 6=pi- 7=pi- with symmetrization of the identical particles
    EvtVector4C WCurrent_7pi(const EvtVector4R& p1, const EvtVector4R& p2,
                             const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5,
                             const EvtVector4R& p6, const EvtVector4R& p7) const;

  // 1=K+ 2 = pi+ 3 = pi+ 4 = pi- 5 = pi- with symmetrizatiom
   EvtVector4C WCurrent_K4pi(const EvtVector4R& p1, const EvtVector4R& p2,
                                    const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5) const;

    // a1 -> p+ p- pi+
    EvtVector4C WCurrent_ppPi(const EvtVector4R& p1, const EvtDiracSpinor& sp1, const EvtVector4R& p2, const EvtDiracSpinor& sp2, const EvtVector4R& k) const;
  
  

protected:
  

    EvtVector4C JB(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3, 
		   const EvtVector4R& q4, const EvtVector4R& q5) const;
    EvtComplex Den(double q,  double mR, double gammaR, double gR) const;

    EvtComplex BWa(const EvtVector4R& q) const;

    EvtComplex BWf(const EvtVector4R& q) const;

    EvtComplex BWr(const EvtVector4R& q) const;

    EvtComplex BWKK(double s, int i) const;

    double pi3G(double Q2) const;

    EvtComplex pcm(double s) const;

    EvtComplex BW( double s, double m, double gamma, double xm1, double xm2 ) const;

    EvtVector4C WCurrent_KKPPP_nosym(const EvtVector4R& pKplus, const EvtVector4R& pKminus,
                               const EvtVector4R& pPi1Plus, const EvtVector4R& pPi2Plus, const EvtVector4R& pPiMinus) const;

  // a1 -> a1(1=pi+ 2=pi+ 3=pi+ 5=pi- 6=pi-) f0(4=pi+ 7=pi-) without symmetrization of the identical particles
   EvtVector4C WCurrent_7pi_nosymm(const EvtVector4R& p1, const EvtVector4R& p2,
                             const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5,
                             const EvtVector4R& p6, const EvtVector4R& p7) const;

  // a1 -> K*0 (1=K+ 4=pi-) a1(2=pi+ 3=pi+ 5=pi-)
   EvtVector4C WCurrent_K4pi_nosymm(const EvtVector4R& p1, const EvtVector4R& p2,
                                    const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5) const;

private:

    std::vector<double> mRho_, gamma0_, cK_, mK_, gammaK_, gKRho_, gKPi_, mPi_;

};

#endif
