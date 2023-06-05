/*******************************************************************************
* Copyright 2019 Alexey Luchinsky                                              *
*   email: alexey.luchinsky@gmail.com                                          *
*                                                                              *
* This program is free software: you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License as published by         *
* the Free Software Foundation, either version 3 of the License, or            *
* (at your option) any later version.                                          *
*                                                                              *
* It is distributed in the hope that it will be useful,                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see   <https://www.gnu.org/licenses/>.     *
*******************************************************************************/


#include "segR.h"
using namespace std;

double EvtRootRandomEngine::random() {
    static TRandom3 randengine(seed);
    return randengine.Rndm();
}

float rlu_() {
    return EvtRandom::Flat();
}

float begran_(int *) {
    return EvtRandom::Flat();
}

void segR(int argc, char *argv[]) {

    EvtStdHep evtstdhep;
    EvtParticle *parent;
    int seed = -1;
    EvtGen *myGenerator;
    EvtRandomEngine* myRandomEngine = 0;
    EvtId PART;
    int number = 10;
    double val = -3.0969;
    double P = 0.0;
    double E = 0.0;
    double mp = 0.938272;
    double md = 1.875613;
    double mtarg;

    Int_t ev = 0;
    Int_t nLine[MAX], pdgID[MAX], nDau[MAX], nM1[MAX], nM2[MAX], nDF[MAX], nDL[MAX];
    Double_t fPx[MAX], fPy[MAX], fPz[MAX], fE[MAX];
    Double_t fVx[MAX], fVy[MAX], fVz[MAX], fT[MAX];
    Double_t fTht[MAX], fM[MAX], fP[MAX], fPt[MAX];
    Double_t prob;

    Int_t nTrk = 0;
    TTree *ntp;

    TStopwatch timer;
    timer.Start();

    if (argc < 3) {
        cout << "\nUSAGE: simpleEvtGenRO <particle> <dec-file> <# events> <pbar-mom/cms-energy> <rand seed>\n" << endl;
        cout << "  <particle> = particle type to decay, e.g. 'eta_c', 'pbarpSystem' etc." << endl; //argv[1]
        cout << "  <dec-file> = EvtGen decay file (.DEC) to use; see directory 'test' for examples" << endl; //argv[2]
        cout << "  <# events> = number of events to produce; default value = 10" << endl; //argv[3]
        cout << "  <pbar-mom> = (>0) momentum of the pbar beam; (<0) negativ cms energy; default value = mass of <particle>" << endl; //argv[4]
        cout << "               mandatory, when <particle> = pbarpSystem" << endl;
        cout << "  <rand seed> = random seed for TRandom3. Value < 0 = use default random gen.; default = -1\n" << endl; //argv[3]
        cout << "Output is stored in file 'evtOutput.root'.\n\n" << endl;
        return;
    }

    if (argc > 5) seed = atoi(argv[5]);
    if (seed >= 0)
        myRandomEngine = new EvtRootRandomEngine(seed);

    myGenerator = new EvtGen(argv[2], "../src/evt.pdl", myRandomEngine);


    PART = EvtPDL::getId(std::string(argv[1]));

    if (argc >= 4) number = atoi(argv[3]);

    if (std::string(argv[1]) == "pbarpSystem" && argc < 5) {
        cout << "\n******  FATAL EVT_ERROR: <particle> is 'pbarpSystem'; MUST give pbar momentum or cms energy!\n\n" << endl;
        return;
    } else if (std::string(argv[1]) != "pbarpSystem" && argc >= 5) {
        cout << "\n****** WARNING: overriding given momentum, setting cms energy to mass of " << argv[1] << ".\n" << endl;
    }

    if (std::string(argv[1]) == "pbardSystem" && argc < 5) {
        cout << "\n******  FATAL EVT_ERROR: <particle> is 'pbardSystem'; MUST give pbar momentum!\n\n" << endl;
        return;
    }

    if (std::string(argv[1]) == "pbardSystem" && atof(argv[4]) < 0) {
        cout << "\n******  FATAL EVT_ERROR: <particle> is 'pbardSystem'; cms energy doesn't make sense, give the pbar momentum!\n\n" << endl;
        return;
    }


    if (std::string(argv[1]) == "pbarpSystem") mtarg = mp;
    if (std::string(argv[1]) == "pbardSystem") mtarg = md;



    if (argc >= 5)
        val = atof(argv[4]);
    else
        val = -EvtPDL::getMass(PART);

    // val is the momentum of the pbar beam
    if (val > 0) {
        P = val;
        E = mtarg + sqrt(P * P + mp * mp);
    } else //val is -E_cm (for p target only)
    {
        val = -val;
        E = val * val / (2 * mp);
        P = sqrt(E * E - val * val);
    }

    P = 0;
    E = EvtPDL::getMass(PART);

    cout << "\n\n############# Generating with following conditions:\n\n";
    cout << "particle       : '" << argv[1] << "'" << endl;
    cout << "decay file     : " << argv[2] << endl;
    cout << "incident 4-mom : (" << E << ", 0, 0, " << P << "), m = " << sqrt(E * E - P * P) << endl;
    cout << "# events       : " << number << "\n\n######################\n\n" << endl;


    ntp = new TTree("ntp", "ntp");

    ntp->Branch("ev", &ev, "ev/I");
    ntp->Branch("prob", &prob, " prod/D");
    ntp->Branch("nTrk", &nTrk, "nTrk/I");

    ntp->Branch("N", nLine, "N[nTrk]/I");
    ntp->Branch("Id", pdgID, "Id[nTrk]/I");

    ntp->Branch("M1", nM1, "M1[nTrk]/I");
    ntp->Branch("M2", nM2, "M2[nTrk]/I");
    ntp->Branch("DF", nDF, "DF[nTrk]/I");
    ntp->Branch("DL", nDL, "DL[nTrk]/I");
    ntp->Branch("nDau", nDau, "nDau[nTrk]/I");

    ntp->Branch("px", fPx, "px[nTrk]/D");
    ntp->Branch("py", fPy, "py[nTrk]/D");
    ntp->Branch("pz", fPz, "pz[nTrk]/D");
    ntp->Branch("E", fE, "E[nTrk]/D");

    ntp->Branch("t", fT, "t[nTrk]/D");
    ntp->Branch("x", fVx, "x[nTrk]/D");
    ntp->Branch("y", fVy, "y[nTrk]/D");
    ntp->Branch("z", fVz, "z[nTrk]/D");

    ntp->Branch("m", fM, "m[nTrk]/D");
    ntp->Branch("p", fP, "p[nTrk]/D");
    ntp->Branch("pt", fPt, "pt[nTrk]/D");
    ntp->Branch("tht", fTht, "tht[nTrk]/D");

    // Loop to create nEvents, starting from an Upsilon(4S)
    int i, j;
    for (i = 0; i < number; i++) {
        // Set up the parent particle

        ev = i;

        EvtVector4R pInit(E, 0.0000, -0.0000, P);
        parent = EvtParticleFactory::particleFactory(PART, pInit);
        parent->setDiagonalSpinDensity();

        // Generate the event
        myGenerator->generateDecay(parent);

        // Write out the results
        evtstdhep.init();
        parent->makeStdHep(evtstdhep);

        //print out some status info
        //    if (i<10) cout << "event Number\t"<< i << evtstdhep << endl;
        if (!((i + 1) % (number / 10))) cout << "event Number\t" << i + 1 << endl;

        // Write the output file
        //cout << i << "\t" << evtstdhep.getNPart();
        //cout <<evtstdhep<<endl;
        nTrk = evtstdhep.getNPart();

        if (nTrk >= MAX) {
            nTrk = MAX;
            cout << "-W- Event #" << i << ": exeeding maximum particle number " << MAX << endl;
        }

        for (j = 0; j < nTrk; j++) {
            nLine[j] = j;
            pdgID[j] = evtstdhep.getStdHepID(j); // PDG code
            nDF[j] = evtstdhep.getFirstDaughter(j); // Index first daughter (-1 if stable)
            nDL[j] = evtstdhep.getLastDaughter(j); // Index last daughter (-1 if stable)
            nDau[j] = nDL[j] > 0 ? (nDL[j] - nDF[j]) + 1 : 0; // Number of daughters 
            nM1[j] = evtstdhep.getFirstMother(j); // Index of first mother 
            nM2[j] = evtstdhep.getLastMother(j); // Index of second mother

            EvtVector4R p4 = evtstdhep.getP4(j);
            EvtVector4R x4 = evtstdhep.getX4(j);

            fE[j] = p4.get(0); // Energy
            fPx[j] = p4.get(1); // Px component
            fPy[j] = p4.get(2); // Py component
            fPz[j] = p4.get(3); // Pz component

            fT[j] = x4.get(0); // Time
            fVx[j] = x4.get(1); // Vertex x component
            fVy[j] = x4.get(2); // Vertex y component
            fVz[j] = x4.get(3); // Vertex z component

            double pt2 = p4.get(1) * p4.get(1) + p4.get(2) * p4.get(2);
            double p2 = pt2 + p4.get(3) * p4.get(3);

            fPt[j] = sqrt(pt2); // Transverse Momentum
            fP[j] = sqrt(p2); // Modulus of Momentum

            fM[j] = sqrt(p4.get(0) * p4.get(0) - p2); // Invariant Mass of particle
            fTht[j] = atan2(fPt[j], p4.get(3)); // Polar angle
        }
        prob = *parent->decayProb();
        ntp->Fill();

        parent->deleteTree();
    }

    TFile *f = new TFile("evtOutput.root", "RECREATE");
    ntp->Write();
    f->Close();

}

