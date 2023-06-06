//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtBcVHad.hh
//
// Description: Module to implement Bc -> psi + (n pi) + (m K) decays
//
// Modification history:
//
//    A V Luchinsky     Jan 29, 2013        Module created
//    A V Luchinsky     Apr 30, 2019        psi K_S K node added
//
//------------------------------------------------------------------------

#ifndef EvtBcVHad_HH
#define EvtBcVHad_HH

#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

#include <array>
#include <string>

class EvtBCVFF2;
class EvtParticle;
class EvtWHad;

class EvtBcVHad : public EvtDecayAmp {
public:

    EvtBcVHad();
    virtual ~EvtBcVHad();

    std::string getName();
    EvtDecayBase* clone();
    void initProbMax();
    void init();
    void decay(EvtParticle *p);

protected:

    // Hadronic current function
    EvtVector4C hardCurr(EvtParticle *root_particle) const;


    void parseDecay(void);


    /// @brief amplitude of bc -> psi W decay
    /// @param root_particle 
    /// @param hadCur effective polrization of virtual W
    /// @param iPol index of psi's polarization (0,1,2)
    /// @return amplitude of the decay
    EvtComplex amp_psi(EvtParticle *root_particle, EvtVector4C hadCur, int iPol);
    void decay_Psi(EvtParticle *p);
    void decay_Psi_mn(EvtParticle *p);

    /// @brief amplitude of bc -> chi_c0 W decay
    /// @param root_particle 
    /// @param hadCur effective polrization of virtual W
    /// @return amplitude of the decay
    EvtComplex amp_chiC0(EvtParticle *root_particle, EvtVector4C hadCur);
    void decay_chiC0(EvtParticle *p);
    void decay_chiC0_mn(EvtParticle *p);

    /// @brief amplitude of bc -> chi_c1 W decay
    /// @param root_particle 
    /// @param hadCur effective polrization of virtual W
    /// @param iPol index of chi_c1 polarization (0,1,2)
    /// @return amplitude of the decay
    EvtComplex amp_chiC1(EvtParticle *root_particle, EvtVector4C hadCur, int iPol);
    void decay_chiC1(EvtParticle *p);
    void decay_chiC1_mn(EvtParticle *p);

    /// @brief amplitude of bc -> chi_c2 W decay
    /// @param root_particle 
    /// @param hadCur effective polrization of virtual W
    /// @param iPol index of chi polarization (0,1,2,3,4)
    /// @return amplitude of the decay
    EvtComplex amp_chiC2(EvtParticle *root_particle, EvtVector4C harCur, int iPol);
    void decay_chiC2(EvtParticle *root_particle);  
    void decay_chiC2_mn(EvtParticle *p);

private:
    // whichfit --- code of the Bc -> VW formfactor set:
    //   1 - SR
    //   2 - PM
    int whichfit;

    // idVector --- final vector particle code
    int idVector;

    // out_code: code of the hadronic final state
    //   1 - pi+
    //   2 - pi+ pi0
    //   3 - pi+ pi+ pi-
    //   4 - 4pi
    //   5 - pi+ pi+ pi- pi- pi+
    //   6 - K+ K- pi+
    //   7 - K+ pi+ pi-
    //   8 - K_S0 K+
    int out_code;  

    EvtBCVFF2 *ffmodel;
    EvtWHad *wcurr;

    std::array<int, 4> iPiPlus{-1, -1, -1, -1};
    std::array<int, 4> iPiMinus{-1, -1, -1, -1};
    std::array<int, 4> iPiZero{-1, -1, -1, -1};
    std::array<int, 4> iKPlus{-1, -1, -1, -1};
    std::array<int, 4> iKMinus{-1, -1, -1, -1};
    std::array<int, 4> iLepton{-1, -1, -1, -1};
    std::array<int, 4> iNeutrino{-1, -1, -1, -1};

};

#endif
