//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtBcVHad.cc
//
// Description: Module to implement Bc -> psi + (n pi) + (m K) decays
//
// Modification history:
//
//    A V Luchinsky     Jan 29, 2013        Module created
//    A V Luchinsky     Apr 30, 2019        psi K_S K mode added
//.   A V Luchinsky.    Nov 4, 2021         K+4pi, KK+3pi, 7pi modes added
//
//------------------------------------------------------------------------


// The following decays are supported:
// out_code=1: B_c+ -> V pi+
// out_code=2: B_c+ -> V pi+ pi0              (auto sym)
// out_code=3: B_c+ -> V + 2pi+ pi-           (auto sym)
// out_code=5: B_c+ -> V + 3pi+ + 2 pi-       (auto sym)
// out_code=6: B_c+ -> V K+ K- pi+            (auto sym)
// out_code=7: B_c+ -> V K+ pi+ pi-           (auto sym)
// out_code=8: B_c+ -> V K_S0 K+
// out_code=9: B_c+ -> V K+ K- pi+ pi+ pi-    (auto sym)
// out_code=10: B_c+ -> V + 4 pi+ + 3 pi-     (auto sym)
// out_code=11: B_c+ -> V K+ pi+ pi+ pi- pi-  (auto sym)


#include "EvtGenModels/EvtBcVHad.hh"

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtBCVFF2.hh"
#include "EvtGenModels/EvtWHad.hh"

#include <iostream>

EvtBcVHad::EvtBcVHad() :  
    whichfit(0), 
    idVector(0), 
    out_code(0),
    ffmodel(0), 
    wcurr(0) 
{
}

EvtBcVHad::~EvtBcVHad()
{

    if (ffmodel) {delete ffmodel;}
    if (wcurr) {delete wcurr;}

}

std::string EvtBcVHad::getName() {
    return "BC_VHAD";
}

EvtDecayBase* EvtBcVHad::clone() {
    return new EvtBcVHad;
}

void EvtBcVHad::parseDecay() 
{
    EvtIdSet BcPlusID("B_c+"), BcMinusID("B_c-");
    EvtIdSet theK("K+", "K-", "K_S0");
    EvtId parentId = getParentId();
    int cMode = BcMinusID.contains(parentId);
    
    EvtIdSet PiPlusID(cMode == 0 ? "pi+" : "pi-");
    EvtIdSet PiMinusID(cMode == 0 ? "pi-" : "pi+");
    EvtIdSet PiZeroID("pi0");
    EvtIdSet KPlusID(cMode == 0 ? "K+" : "K-");
    EvtIdSet KMinusID(cMode == 0 ? "K-" : "K+");
    EvtGenReport(EVTGEN_INFO, "EvtBcVHad") << "parentId = " << getParentId() << std::endl;
        
    int PiPlusFound = 0, PiMinusFound = 0, PiZeroFound = 0, KPlusFound = 0, KMinusFound = 0;
    for (int iDaughter = 0; iDaughter < getNDaug(); ++iDaughter) 
    {
      const EvtId daugId = getDaug(iDaughter);
      EvtGenReport(EVTGEN_INFO, "EvtBcVHad") << "iDaughter = " << iDaughter
					     << " id = " <<getDaug(iDaughter).getName() <<std::endl;
      if (PiPlusID.contains(daugId) && PiPlusFound < 4)
      {
        iPiPlus[PiPlusFound] = iDaughter;
        PiPlusFound++;
      }
      else if (PiMinusID.contains(daugId) && PiMinusFound < 4)
      {
        iPiMinus[PiMinusFound] = iDaughter;
        PiMinusFound++;
      }
      else if (PiZeroID.contains(daugId) && PiZeroFound < 4)
      {
        iPiZero[PiZeroFound] = iDaughter;
        PiZeroFound++;
      }
      else if (KPlusID.contains(daugId) && KPlusFound < 4)
      {
        iKPlus[KPlusFound] = iDaughter;
        KPlusFound++;
      }
      else if (KMinusID.contains(daugId) && KMinusFound < 4)
      {
        iKMinus[KMinusFound] = iDaughter;
        KMinusFound++;
      }
    }    
    
    if (getNDaug() == 2 && PiPlusFound == 1)
    {
      out_code = 1; // pi+
    }
    else if (getNDaug() == 3 && PiPlusFound == 1 && PiZeroFound == 1)
    {     
      out_code = 2; // pi+ pi0
    }
    else if (getNDaug() == 4 && PiPlusFound == 2 && PiMinusFound == 1)
    {
       out_code = 3; // pi+ pi+ pi-
    } 
    else if (getNDaug() == 5 && PiPlusFound == 2 && PiMinusFound == 1 && PiZeroFound == 1)
    {
      out_code = 4; // pi+ pi+ pi- pi0
    }
    else if (getNDaug() == 6 && PiPlusFound == 3 && PiMinusFound == 2)
    {
      out_code = 5; // 5pi
    }
    else if (getNDaug() == 4 && KPlusFound == 1 && KMinusFound == 1 && PiPlusFound == 1)
    {
      out_code = 6; // KKpi
    }
    else if (getNDaug() == 4 && KPlusFound == 1 && PiPlusFound == 1 && PiMinusFound == 1)
    {
      out_code = 7; // K+ pi+ pi- 
    }
    else if (getNDaug() == 3 && theK.contains(getDaug(1)) && theK.contains(getDaug(2)))
    {
      out_code = 8; // KK
    }
    else if (getNDaug() == 6 && KPlusFound == 1 && KMinusFound == 1 && PiPlusFound == 2 && PiMinusFound == 1)
    {
      out_code = 9; // K+ K- pi+ pi+ pi-
    }
    else if (getNDaug() == 8 && PiPlusFound == 4 && PiMinusFound == 3)
    {
      out_code = 10; // 7pi
    }
    else if (getNDaug() == 6 && KPlusFound == 1 && PiPlusFound == 2 && PiMinusFound == 2)
    {
      out_code = 11; // K+ pi+ pi+ pi- pi-
    }
    else
    {
      EvtGenReport(EVTGEN_ERROR, "EvtBcVHad")<<"Init: unknown decay"<<std::endl;
    }
    
    EvtGenReport(EVTGEN_INFO, "EvtBcVHad") << "out_code = " <<out_code << ", whichfit = " << whichfit << std::endl;
    for (int i = 0; i < 4; ++i) 
    {
      EvtGenReport(EVTGEN_INFO, "EvtBcVHad") << " i = "<< i
					     << ", iPiPlus = " << iPiPlus[i]
					     << ", iPiMinus = " << iPiMinus[i]
					     << ", iPiZero = " << iPiZero[i]
					     << ", iKPlus = " << iKPlus[i]
					     << ", iKMinus = " << iKMinus[i] << std::endl;
    }
    EvtGenReport(EVTGEN_INFO, "EvtBcVHad") << "PiPlusFound = " << PiPlusFound
					   << ", PiMinusFound = " << PiMinusFound
					   << ", PiZeroFound = " << PiZeroFound
					   << ", KPlusFound = " << KPlusFound
					   << ", KMinusFound = " << KMinusFound << std::endl;
}

//======================================================
void EvtBcVHad::init() {

    checkNArg(1);
    checkSpinParent(EvtSpinType::SCALAR);
    checkSpinDaughter(0, EvtSpinType::VECTOR);
    for (int i = 1; i <= (getNDaug() - 1); i++) {
        checkSpinDaughter(i, EvtSpinType::SCALAR);
    }

    idVector = getDaug(0).getId();
    whichfit = int(getArg(0) + 0.1);
    ffmodel = new EvtBCVFF2(idVector, whichfit);

    wcurr = new EvtWHad();
    // Determine the code of final hadronic state    
    parseDecay();
    
}

//======================================================

void EvtBcVHad::initProbMax() {
    if (out_code == 1) {
    if (idVector == EvtPDL::getId("J/psi").getId()        && whichfit == 1 && getNDaug() == 2) setProbMax(500.);
    else if (idVector == EvtPDL::getId("J/psi").getId()   && whichfit == 2 && getNDaug() == 2) setProbMax(300.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1 && getNDaug() == 2) setProbMax(17.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2 && getNDaug() == 2) setProbMax(40.);
  } else if (out_code == 2) {
    if (idVector == EvtPDL::getId("J/psi").getId()        && whichfit == 1 && getNDaug() == 3) setProbMax(10950.);
    else if (idVector == EvtPDL::getId("J/psi").getId()   && whichfit == 2 && getNDaug() == 3) setProbMax(4200.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1 && getNDaug() == 3) setProbMax(500.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2 && getNDaug() == 3) setProbMax(700.);
  } else if (out_code == 3) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1 && getNDaug() == 4) setProbMax(42000.);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2 && getNDaug() == 4) setProbMax(90000.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1 && getNDaug() == 4) setProbMax(1660.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2 && getNDaug() == 4) setProbMax(2600.);
  } else if (out_code == 5) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1 && getNDaug() == 6) setProbMax(720000.);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2 && getNDaug() == 6) setProbMax(519753.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1 && getNDaug() == 6) setProbMax(40000.);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2 && getNDaug() == 6) setProbMax(30000.);
  } else if (out_code == 6) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1) setProbMax(50000.);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2) setProbMax(22000.0);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1) setProbMax(2300.0);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2) setProbMax(1700.00);
  } else if (out_code == 7) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1) setProbMax(2.2e+06);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2) setProbMax(930000);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1) setProbMax(92000.0);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2) setProbMax(93000.0);
  } else if (out_code == 8) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1) setProbMax(2e2);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2) setProbMax(80);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1) setProbMax(10);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2) setProbMax(10);
  } else if(out_code == 9) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1) setProbMax(3e4);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2) setProbMax(18540);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1) setProbMax(0.15*1e4);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2) setProbMax(2*500);
  } else if(out_code == 10) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1) setProbMax(2e6);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2) setProbMax(5e6);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1) setProbMax(1.5e5);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2) setProbMax(1e5);
  } else if(out_code == 11) {
    if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 1) setProbMax(2.5e8);
    else if (idVector == EvtPDL::getId("J/psi").getId() && whichfit == 2) setProbMax(1.4e7);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 1) setProbMax(2e6);
    else if (idVector == EvtPDL::getId("psi(2S)").getId() && whichfit == 2) setProbMax(8e4);
  } else {
    EvtGenReport(EVTGEN_ERROR, "EvtBcVHad") << "probax: Have not yet implemented this final state in BC_VHAD model, out_code = " << out_code << std::endl;
    for (int id = 0; id < (getNDaug() - 1); id++) {
      EvtGenReport(EVTGEN_ERROR, "EvtBcVHad") << "Daug " << id << " " << EvtPDL::name(getDaug(id)).c_str() << std::endl;
    }
    ::abort();
  }

}

//======================================================

EvtVector4C EvtBcVHad::hardCurr(EvtParticle *root_particle) const {

    EvtVector4C hardCur;

    if (out_code == 1) {

	// pi+
        hardCur = wcurr->WCurrent(root_particle->getDaug(iPiPlus[0])->getP4());

    } else if (out_code == 2) {

	// pi+ pi0
        hardCur = wcurr->WCurrent(
                                  root_particle->getDaug(iPiPlus[0])->getP4(),
                                  root_particle->getDaug(iPiZero[0])->getP4());

    } else if (out_code == 3) {

	// pi+ pi+ pi-
        hardCur = wcurr->WCurrent(
          root_particle->getDaug(iPiPlus[0])->getP4(),
          root_particle->getDaug(iPiPlus[1])->getP4(),
          root_particle->getDaug(iPiMinus[0])->getP4()
        );
    } else if (out_code == 5) {
	// Bc -> psi pi+ pi+ pi- pi- pi+ from Kuhn, Was, hep-ph/0602162
        hardCur = wcurr->WCurrent_5pi(
                                      root_particle->getDaug(iPiPlus[0])->getP4(),
                                      root_particle->getDaug(iPiPlus[1])->getP4(),
                                      root_particle->getDaug(iPiPlus[2])->getP4(),
                                      root_particle->getDaug(iPiMinus[0])->getP4(),
                                      root_particle->getDaug(iPiMinus[1])->getP4()
                                      );

    } else if (out_code == 6) {
	// K+ K- pi+
        hardCur = wcurr->WCurrent_KKP(
                                      root_particle->getDaug(iKPlus[0])->getP4(),
                                      root_particle->getDaug(iKMinus[0])->getP4(),
                                      root_particle->getDaug(iPiPlus[0])->getP4());

    } else if (out_code == 7) {
	// K+ pi+ pi-
        hardCur = wcurr->WCurrent_KPP(
                                      root_particle->getDaug(iKPlus[0])->getP4(),
                                      root_particle->getDaug(iPiPlus[0])->getP4(),
                                      root_particle->getDaug(iPiMinus[0])->getP4());

    } else if (out_code == 8) {
	
	// K_S0 K+
        hardCur = wcurr->WCurrent_KSK(root_particle->getDaug(1)->getP4(),
				      root_particle->getDaug(2)->getP4());

    } 
    else if (out_code == 9) {
      // K+ K- pi+ pi+ pi-
      hardCur = wcurr -> WCurrent_KKPPP(
                                        root_particle->getDaug(iKPlus[0])->getP4(), // K+
                                        root_particle->getDaug(iKMinus[0])->getP4(), // K-
                                        root_particle->getDaug(iPiPlus[0])->getP4(), // pi+
                                        root_particle->getDaug(iPiPlus[1])->getP4(), // pi+
                                        root_particle->getDaug(iPiMinus[0])->getP4() // pi-
                                        );
    }
    else if (out_code == 10) {
 // 1=pi+ 2=pi+ 3=pi+ 4=pi+ 5=pi- 6=pi- 7=pi- with symmetrization of the identical particles
      hardCur = wcurr -> WCurrent_7pi(
                                        root_particle->getDaug(iPiPlus[0])->getP4(), // pi+
                                        root_particle->getDaug(iPiPlus[1])->getP4(), // pi+
                                        root_particle->getDaug(iPiPlus[2])->getP4(), // pi+
                                        root_particle->getDaug(iPiPlus[3])->getP4(), // pi+
                                        root_particle->getDaug(iPiMinus[0])->getP4(), // pi-
                                        root_particle->getDaug(iPiMinus[1])->getP4(), // pi-
                                        root_particle->getDaug(iPiMinus[2])->getP4() // pi-
                                        );
    }
    else if (out_code == 11) {
 // 1=K+ 2 = pi+ 3 = pi+ 4 = pi- 5 = pi- with symmetrizatiom
      hardCur = wcurr -> WCurrent_K4pi(
                                        root_particle->getDaug(iKPlus[0])->getP4(), // pi+
                                        root_particle->getDaug(iPiPlus[0])->getP4(), // pi+
                                        root_particle->getDaug(iPiPlus[1])->getP4(), // pi+
                                        root_particle->getDaug(iPiMinus[0])->getP4(), // pi-
                                        root_particle->getDaug(iPiMinus[1])->getP4() // pi-
                                        );
    } else {
      EvtGenReport(EVTGEN_ERROR, "EvtBcVHad") << "hardCurr: Have not yet implemented this final state in BC_VHAD model" << std::endl;
      for (int id = 0; id < (getNDaug() - 1); id++) {
	EvtGenReport(EVTGEN_ERROR, "EvtBcVHad") << "Daug " << id << " " << EvtPDL::name(getDaug(id)).c_str() << std::endl;
      }
      ::abort();

    }

    return hardCur;
}

//======================================================

void EvtBcVHad::decay(EvtParticle *root_particle) {

    root_particle->initializePhaseSpace(getNDaug(), getDaugs());

    // Calculate hadronic current
    EvtVector4C hardCur = hardCurr(root_particle);

    EvtParticle* Jpsi = root_particle->getDaug(0);
    
    EvtVector4R
	p4b(root_particle->mass(), 0., 0., 0.), // Bc momentum
	p4meson = Jpsi->getP4(), // J/psi momenta
	Q = p4b - p4meson, 
	p4Sum = p4meson + p4b;
    double Q2 = Q.mass2();

    // Calculate Bc -> V W form-factors
    double a1f(0.0), a2f(0.0), vf(0.0), a0f(0.0);

    double m_meson = Jpsi->mass();
    double m_b = root_particle->mass();
    double mVar = m_b + m_meson;

    ffmodel->getvectorff(root_particle->getId(),
			 Jpsi->getId(), 
			 Q2, m_meson, &a1f, &a2f, &vf, &a0f);

    double a3f = (mVar/(2.0*m_meson))*a1f - ((m_b - m_meson)/(2.0*m_meson))*a2f;

    // Calculate Bc -> V W current
    EvtTensor4C H = a1f*mVar*EvtTensor4C::g();
    H.addDirProd((-a2f/mVar)*p4b, p4Sum);
    H += EvtComplex(0.0, vf/mVar)*dual(EvtGenFunctions::directProd(p4Sum, Q));
    H.addDirProd((a0f - a3f)*2.0*(m_meson/Q2)*p4b, Q);
    EvtVector4C Heps = H.cont2(hardCur);

    for (int i = 0; i < 4; i++) {
        EvtVector4C eps = Jpsi->epsParent(i).conj(); // psi-meson polarization vector
        EvtComplex amp = eps*Heps;
        vertex(i, amp);
    }

}
