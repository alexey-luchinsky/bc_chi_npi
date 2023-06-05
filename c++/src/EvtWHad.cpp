//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtWHad.cc
//
// Description: Routine to calculate W -> (n pi) + (m K) current
//			according to [Kuhn, Was, Acta.Phys.Polon B39 (2008) 147]
//
// Modification history:
//	A V Luchinsky	 20 Jan, 2013	Module created
//	A V Luchinsky	 09 May, 2019	Ks K mode added
//      
//------------------------------------------------------------------------

#include "EvtGenModels/EvtWHad.hh"
#include "EvtGenBase/EvtTensor4C.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtTensor4C.hh"

EvtWHad::EvtWHad() :
    mRho_(),
    gamma0_(),
    cK_(0),
    mK_(),
    gammaK_(),
    gKRho_(),
    gKPi_(),
    mPi_()
{
    // cK coefficients from Eur. Phys. J. C39, 41 (2005), arXiv:hep-ph/0409080 [hep-ph]  

    // rho(770)
    mRho_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("rho0")));
    gamma0_.push_back(EvtPDL::getWidth( EvtPDL::getId("rho0" )));
    cK_.push_back(1.195);

    // rho(1450)
    mRho_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("rho(2S)0")));
    gamma0_.push_back(EvtPDL::getWidth( EvtPDL::getId("rho(2S)0" )));
    cK_.push_back(-0.112);

    // rho(1700)
    mRho_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("rho(3S)0")));
    gamma0_.push_back(EvtPDL::getWidth( EvtPDL::getId("rho(3S)0" )));
    cK_.push_back(-0.083);

    // rho(2150), PRD 76 092005
    mRho_.push_back(2.150);
    gamma0_.push_back(0.310);
    cK_.push_back(0.0);

    // Storing K - resonance  information

   // K(892)
    
     mK_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("K*0")));
     gammaK_.push_back(EvtPDL::getWidth( EvtPDL::getId( "K*0" )));
     gKRho_.push_back(0.0); 
     gKPi_.push_back(3.26);
   

    // K1(1270)

     mK_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("K_10")));
     gammaK_.push_back(EvtPDL::getWidth( EvtPDL::getId( "K_10" )));
     gKRho_.push_back(2.71); 
     gKPi_.push_back(0.792);


     // K1(1400)
     mK_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("K'_10")));
     gammaK_.push_back(EvtPDL::getWidth( EvtPDL::getId( "K'_10" )));
     gKRho_.push_back(0.254); 
     gKPi_.push_back(2.509);

     //m Pi
     mPi_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("pi+")));
     mPi_.push_back(EvtPDL::getMeanMass(EvtPDL::getId("pi-")));

}


EvtComplex EvtWHad::BWKK(double s, int i) const {
    double m2 = pow(mRho_[i], 2);
    EvtComplex qs = pcm(s);
    EvtComplex qm = pcm(m2);
    if (abs(qm) < 1e-10) {return 0;}

    EvtComplex rat = qs/qm;
    EvtComplex rat3 = rat*rat*rat;
    if (abs(s) < 1e-10) {return 0;}

    EvtComplex gamma = m2*rat3*gamma0_[i]/s;
    EvtComplex I(0.0, 1.0);

    EvtComplex denBW = m2 - s - I*sqrt(s)*gamma;
    if (abs(denBW) < 1e-10) {return 0;}

    return cK_[i]*m2/denBW;
}

EvtVector4C EvtWHad::WCurrent_KSK(const EvtVector4R& pKS, const EvtVector4R& pKplus) const {
    double s = (pKS + pKplus).mass2();
    EvtComplex f = BWKK(s, 0) + BWKK(s, 1) + BWKK(s, 2);

    return f*(pKS - pKplus);

}

EvtComplex EvtWHad::pcm(double s) const {
    double mpi2 = pow(0.140, 2);
    if (abs(s) < 1e-10) return 0;

    double pcm2 = 1.0 - 4.0*mpi2/s;
    EvtComplex result;

    if (pcm2 >= 0.0) {
	result = EvtComplex(sqrt(pcm2), 0.0);
    } else {
	result = EvtComplex(0.0, sqrt(-pcm2));
    }

    return result;
}

// =================== W+ -> pi_ current ========================================

EvtVector4C EvtWHad::WCurrent(const EvtVector4R& q1) const {
    return q1;
}

//====================== W+ -> pi+ pi0 current =========================================

EvtVector4C EvtWHad::WCurrent(const EvtVector4R& q1, const EvtVector4R& q2) const {
    return BWr(q1 + q2)*(q1 - q2);
}

//========================= W+ -> pi+ pi+ pi- current ==============================================

EvtVector4C EvtWHad::WCurrent(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3) const {
    EvtVector4R Q = q1 + q2 + q3;
    EvtVector4R q13 = q1 - q3, q23 = q2 - q3;
    double Q2 = Q.mass2();
    return BWa(Q)*(q13 - (Q*(Q*q13)/Q2) * BWr(q2 + q3) + q23 - (Q*(Q*q23)/Q2) * BWr(q1 + q3));
}

// ================= W+ -> pi+ pi+ pi- pi- pi+ current with symmetrization ================================

EvtVector4C EvtWHad::WCurrent(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3, 
			      const EvtVector4R& q4, const EvtVector4R& q5) const {

    EvtVector4C term1 = JB(q1, q2, q3, q4, q5);
    EvtVector4C term2 = JB(q5, q2, q3, q4, q1);
    EvtVector4C term3 = JB(q1, q5, q3, q4, q2);
    EvtVector4C term4 = JB(q1, q2, q4, q3, q5);
    EvtVector4C term5 = JB(q5, q2, q4, q3, q1);
    EvtVector4C term6 = JB(q1, q5, q4, q3, q2);

    EvtVector4C V = term1 + term2 + term3 + term4 + term5 + term6;
    return V;

}

// W+ -> pi+ pi+ pi+ pi- pi-
EvtVector4C EvtWHad::WCurrent_5pi(const EvtVector4R& q1, const EvtVector4R& q2, const EvtVector4R& q3,
                                  const EvtVector4R& q4, const EvtVector4R& q5) const 
{
  return EvtWHad::WCurrent(q1, q2, q4, q5, q3); // WCurrent(++--+)
}



// =========================W+ -> K+ K- pi+ current ==================================================

EvtVector4C EvtWHad::WCurrent_KKP(const EvtVector4R& pKplus, const EvtVector4R& pKminus, 
				  const EvtVector4R& pPiPlus) const {


    const double mA1(1.239), gammaA1(0.600);

    EvtVector4R q = pKplus + pKminus + pPiPlus;
    double q2 = q.mass2();
    EvtVector4R pK = pKminus + pPiPlus;
    double pK2 = pK.mass2();

    EvtComplex I(0.0, 1.0), den1, den2;

    den1 = 1.0/(q2 - mA1*mA1 + I*mA1*gammaA1);
    den2 = 1.0/(pK2 - mK_[0]*mK_[0] + I*mK_[0]*gammaK_[0]);   //K(892)

    EvtTensor4C ten = EvtTensor4C::g() - (1.0/q2)*EvtGenFunctions::directProd(q, q);

    EvtVector4C vec = den1*den2*(pKminus - pPiPlus);
    vec = ten.cont2(vec);

    return vec;
}


// hadronic hurrent W -> K+ K- pi+ pi+ pi- with identical pi+ symmetry
EvtVector4C EvtWHad::WCurrent_KKPPP(const EvtVector4R& pKplus, const EvtVector4R& pKminus,
                                    const EvtVector4R& pPi1Plus, const EvtVector4R& pPi2Plus, const EvtVector4R& pPiMinus) const 
{
  return EvtWHad::WCurrent_KKPPP_nosym(pKplus, pKminus, pPi1Plus, pPi2Plus, pPiMinus) + 
    EvtWHad::WCurrent_KKPPP_nosym(pKplus, pKminus, pPi2Plus, pPi1Plus, pPiMinus);
}


// hadronic hurrent W -> a1(K+ K- pi1+) f0(pi2+ pi-) without identical pi+ symmetry
EvtVector4C EvtWHad::WCurrent_KKPPP_nosym(const EvtVector4R& pKplus, const EvtVector4R& pKminus,
                                    const EvtVector4R& pPi1Plus, const EvtVector4R& pPi2Plus, const EvtVector4R& pPiMinus) const 
{
  //EvtVector4C EvtWHad::WCurrent_KKP(const EvtVector4R& pKplus, const EvtVector4R& pKminus, const EvtVector4R& pPiPlus) const 
  EvtVector4R pf0 = pPi2Plus + pPiMinus;
  EvtVector4C epsA1 = EvtWHad::WCurrent_KKP(pKplus, pKminus, pPi1Plus);
  EvtVector4R q = pKplus + pKminus + pPi1Plus + pPi2Plus + pPiMinus;
  return BWa(q)*epsA1*BWf(pf0);
}

  
 // 1=pi+ 2=pi+ 3=pi+ 4=pi+ 5=pi- 6=pi- 7=pi- with symmetrization of the identical particles
EvtVector4C EvtWHad::WCurrent_7pi(const EvtVector4R& p1, const EvtVector4R& p2,
                             const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5,
                                  const EvtVector4R& p6, const EvtVector4R& p7) const 
{
// a1 -> a1(1=pi+ 2=pi+ 3=pi+ 5=pi- 6=pi-) f0(4=pi+ 7=pi-) without symmetrization of the identical particles
// making p4 symmetric with p1, p2, p3
//        p7                p5, p6
    EvtVector4C eps;
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p2, p3, p4, p5, p6, p7);
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p2, p4, p3, p5, p6, p7);
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p4, p3, p2, p5, p6, p7);
    eps += EvtWHad::WCurrent_7pi_nosymm(p4, p2, p3, p1, p5, p6, p7);
    //
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p2, p3, p4, p5, p7, p6);
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p2, p4, p3, p5, p7, p6);
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p4, p3, p2, p5, p7, p6);
    eps += EvtWHad::WCurrent_7pi_nosymm(p4, p2, p3, p1, p5, p7, p6);
    //
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p2, p3, p4, p7, p6, p5);
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p2, p4, p3, p7, p6, p5);
    eps += EvtWHad::WCurrent_7pi_nosymm(p1, p4, p3, p2, p7, p6, p5);
    eps += EvtWHad::WCurrent_7pi_nosymm(p4, p2, p3, p1, p7, p6, p5);

   return eps;
}



// a1 -> a1(1=pi+ 2=pi+ 3=pi+ 5=pi- 6=pi-) f0(4=pi+ 7=pi-) without symmetrization of the identical particles
EvtVector4C EvtWHad::WCurrent_7pi_nosymm(const EvtVector4R& p1, const EvtVector4R& p2,
                             const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5,
                                    const EvtVector4R& p6, const EvtVector4R& p7) const {
  EvtVector4R qTot = p1+p2+p3+p4+p5+p6+p7;
  //EvtVector4R qa = p1+p2+p3+p5+p6;
  EvtVector4C eps1 = EvtWHad::WCurrent_5pi(p1, p2, p3, p5, p6); // pi+ pi+ pi+ pi- p-
  EvtVector4R pf0 = p4 + p7;
  return eps1*BWa(qTot)*BWf(pf0);
}

// hadronic current W+ -> K+ pi+ pi-

EvtVector4C EvtWHad::WCurrent_KPP(const EvtVector4R& pKplus, const EvtVector4R& pPiPlus, 
				  const EvtVector4R& pPiMinus) const {

    const double cK1p = 0.210709, cK1r = -0.0152997, cK2p = 0.0945309, cK2r = 0.504315;
    const double gRho_PiPi = 6.02;

    EvtVector4R q = pKplus + pPiPlus + pPiMinus;
    double q2 = q.mass2(), pp2(0.0);

    EvtVector4C curr(0, 0, 0, 0), curr1;

    // W+ -> K1+(1270) -> K+ rho0 -> K+ pi+ pi-


    pp2 = (pPiPlus + pPiMinus).mass2();
    curr1 = (pPiPlus - pPiMinus)*Den(q2,mK_[1],gammaK_[1],gKRho_[1])*Den(pp2,mRho_[0],gamma0_[0], gRho_PiPi); //K1(1270) and rho(770) 
    curr = curr + cK1r*curr1;

    // W+ -> K1+(1270) -> K*(892)0 pi+ -> K+ pi- pi-

  
    pp2 = (pKplus + pPiMinus).mass2();
    curr1 = (pKplus - pPiMinus)*Den(q2,mK_[1],gammaK_[1],gKPi_[1])*Den(pp2,mK_[0],gammaK_[0],gKPi_[0]); //K1(1270) and K(892) 
    curr = curr + cK1p*curr1;

    // W+ -> K1+(1400) -> K+ rho0 -> K+ pi+ pi-

    pp2 = (pPiMinus + pPiPlus).mass2();
    curr1 = (pPiPlus - pPiMinus) *Den(q2,mK_[2],gammaK_[2],gKRho_[2])*Den(pp2,mRho_[0],gamma0_[0],gRho_PiPi);//K1(1400) and rho(770) 
    curr = curr + cK2r*curr1;

    // W+ -> K1+(1400) -> K*(892)0 pi+ -> K+ pi- pi+

	
    pp2 = (pKplus + pPiMinus).mass2();
    curr1 = (pKplus - pPiPlus)*Den(q2,mK_[2],gammaK_[2],gKPi_[2])*Den(pp2,mK_[0],gammaK_[0],gKPi_[0]);//K1(1400) and K(892) 
    curr = curr + cK2p*curr1;

    EvtTensor4C ten = EvtTensor4C::g() - (1.0/q2)*EvtGenFunctions::directProd(q, q);
    curr = ten.cont2(curr);

    return curr;
}

EvtComplex EvtWHad::Den(double qSq, const double mR, const double gammaR, const double gR) const {
		EvtComplex I(0.0, 1.0), tmp;
	  tmp = qSq - mR*mR + I*mR*gammaR;

		if(abs(tmp)<1e-10) return 0;
 		return gR / tmp;
}

// a1 -> pi+ pi+ pi- BW
EvtComplex EvtWHad::BWa(const EvtVector4R& q) const {

    double _mA1(1.26), _GA1(0.4);
    double _mA1Sq = _mA1*_mA1;
    EvtComplex I(0.0, 1.0);
    double Q2 = q.mass2();
    double GA1 = _GA1*pi3G(Q2)/pi3G(_mA1Sq);

    EvtComplex denBA1(_mA1Sq - Q2, -1.0*_mA1*GA1);
    if(abs(denBA1)<1e-10) return 0;
    return _mA1Sq/denBA1;
}

EvtComplex EvtWHad::BWf(const EvtVector4R& q) const {
    double mf(0.8), Gf(0.6);
    double mfSq = mf*mf;
    EvtComplex I(0.0, 1.0);
    double Q2 = q.mass2();
    return mfSq/(mfSq - Q2 - I*mf*Gf);
}

EvtComplex EvtWHad::BWr(const EvtVector4R& q) const {
    double beta(-0.108);
  
    double s = q.mass2();
    EvtComplex BW_rho = BW(s, mRho_[0], gamma0_[0], mPi_[0], mPi_[1]);
    EvtComplex BW_rhopr = BW(s, mRho_[1],gamma0_[1], mPi_[0], mPi_[1]);
    return (BW_rho + beta*BW_rhopr)/(1.0 + beta);
}

double EvtWHad::pi3G(double Q2) const {

    double mpi2 = pow( mPi_[0],2.);
    if (Q2 < pow(mRho_[0] + mPi_[0], 2.)) {
        double arg = Q2 - 9. * mpi2;
        return 4.1 * pow(arg, 3.) * (1. - 3.3 * arg + 5.8 * pow(arg, 2.));
    } else
        return Q2 * (1.623 + 10.38 / Q2 - 9.32 / pow(Q2, 2.) + 0.65 / pow(Q2, 3.));

}

EvtVector4C EvtWHad::JB(const EvtVector4R& p1, const EvtVector4R& p2, const EvtVector4R& p3, 
			const EvtVector4R& p4, const EvtVector4R& p5) const {

    EvtVector4R Qtot = p1 + p2 + p3 + p4 + p5, Qa = p1 + p2 + p3;
    EvtTensor4C T = (1.0/Qtot.mass2()) * EvtGenFunctions::directProd(Qtot, Qtot) - EvtTensor4C::g();

    EvtVector4R p13 = p1 - p3, p23 = p2 - p3;
    EvtVector4R V13 = Qa*(p2*p13)/Qa.mass2() - p13;
    EvtVector4R V23 = Qa*(p1*p23)/Qa.mass2() - p23;

    return BWa(Qtot)*BWa(Qa)*BWf(p4 + p5)*(T.cont1(V13)*BWr(p1 + p3) + T.cont1(V23)*BWr(p2 + p3));

}

EvtComplex EvtWHad::BW( double s, double m, double gamma, double xm1, double xm2 ) const {
  double m2 = pow( m, 2.);
  if ( s > pow( xm1+xm2, 2.) ) {
    double qs = sqrt( fabs( (s-pow(xm1+xm2,2.)) * (s-pow(xm1-xm2,2.)) ) ) / sqrt(s); 
    double qm = sqrt( fabs( (m2-pow(xm1+xm2,2.)) * (m2-pow(xm1-xm2,2.)) ) ) / m;
    
    gamma *= m2/s * pow( qs/qm, 3.); 
  }
  else
    gamma = 0.;
  EvtComplex denBW( m2 - s, -1.* sqrt(s) * gamma );
  return m2 / denBW;
}



EvtVector4C EvtWHad::WCurrent_K4pi(const EvtVector4R& p1, const EvtVector4R& p2,
                                   const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5) const
{
  return EvtWHad::WCurrent_K4pi_nosymm(p1, p2, p3, p4, p5) +EvtWHad::WCurrent_K4pi_nosymm(p1, p2, p3, p5, p4);
}


// a1 -> K*0 (1=K+ 4=pi-) a1(2=pi+ 3=pi+ 5=pi-)
EvtVector4C EvtWHad::WCurrent_K4pi_nosymm(const EvtVector4R& p1, const EvtVector4R& p2,
                                          const EvtVector4R& p3, const EvtVector4R& p4, const EvtVector4R& p5) const
{
  const EvtComplex I(0, 1);
  EvtVector4R pKstar = p1 + p4, pa1 = p2+p3+p5;//, pTot = pKstar + pa1;
  EvtComplex denKstar = pKstar*pKstar -  mK_[0]*mK_[0] + I*mK_[0]*gammaK_[0]; //K(892)
  if(abs(denKstar)<1e-10) {
    denKstar = 1e10;
  }
  EvtVector4C epsKstar = 1/denKstar*(p1 - p4);
  EvtVector4C epsA1 = WCurrent(p2, p3, p5);
  EvtVector4C eps = dual(EvtGenFunctions::directProd(epsKstar, epsA1)).cont2(pKstar-pa1);
  //  return BWa(pTot)*eps;
  return eps;
}

EvtVector4C EvtWHad::WCurrent_ppPi(const EvtVector4R& p1, const EvtDiracSpinor& sp1, const EvtVector4R& p2,
				   const EvtDiracSpinor& sp2, const EvtVector4R& k) const {
  const EvtVector4R q = p1 + p2 + k;
  const double q2 = q.mass2();
  const double mp = p1.mass(), mp2 = mp*mp, mpi = k.mass(), mpi2 = mpi*mpi;
  const double mn = EvtPDL::getMeanMass( EvtPDL::getId( "n0" ) );
  const double mn2 = mn*mn;
  const EvtComplex II(0, 1);

  const double kp2 = k*p2;
  const double p1p2 = p1*p2;

  const double f1 = 1.0;
  const double f2 = 3.7/(2.0 + mp);
  const double g1 = 1.25;
  const double g3 = 2.0*mp*g1/(p1p2 + mp2);

  const EvtComplex curS = EvtLeptonSCurrent(sp1, sp2);
  const EvtComplex curP = EvtLeptonPCurrent(sp1, sp2);
  const EvtVector4C curV = EvtLeptonVCurrent(sp1, sp2);
  const EvtVector4C curA = EvtLeptonACurrent(sp1, sp2);
  const EvtTensor4C curT = EvtLeptonTCurrent(sp1, sp2);

  const double D1 = 1/(2*kp2 + mp2 - mn2 + mpi2);  //(k+p2)^2-mn^2

  // Amplitude: ~U(p1).GA5.Vpp(alpha).prop(-p2-k).V(p2) + permutations
  // U() and V() are proton and antiproton spinors, prop() is the proton's propagator,
  // and Vpp(alpha) is the W->pp vertex.
  // Expand terms to use basic spinor currents (Scalar, Vector, Axial, etc)

  EvtVector4C current;
  current += curA*(-(f2*II));
  current += (D1*g1*II)*curT.cont2(k);
  current += -(-0.5*(D1*f1))*dual(curT).cont2(k);
  current += -(D1*f2*II*mp)*dual(curT).cont2(k);
  current += k*(-(D1*g1))*curS;
  current += k*(-(D1*f1))*curP;
  current += k*(2*D1*f2*II*mp)*curP;
  current += k*(curV*k)*(D1*g3);
  current += k*(curA*k)*(D1*f2*II);
  current += p1*(curV*k)*(D1*g3);
  current += p1*(curA*k)*(-(D1*f2*II));
  current += p2*(curV*k)*(D1*g3);
  current += p2*(curA*k)*(D1*f2*II);

  const EvtTensor4C JT = (1.0/q2)*EvtGenFunctions::directProd(q, q) - EvtTensor4C::g();
  current = JT.cont2(current);

  return BWa(q)*current;
}

