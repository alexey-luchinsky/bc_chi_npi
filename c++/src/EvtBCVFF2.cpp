//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtBCVFF2.cc
//
// Description: form factors for B->Vlnu 
//
// Modification history:
//
//    AVL Jan 29, 2013, Module created
//
//------------------------------------------------------------------------
// 

#include "EvtGenModels/EvtBCVFF2.hh"


using namespace std;

EvtBCVFF2::EvtBCVFF2(int idV, int fit) {

  idVector = idV;
  whichfit = fit;
  //cout<<"==== EvtBCVFF2:: idVector="<<idVector<<" whichfit="<<whichfit<<endl;
  return;
}

void EvtBCVFF2::getvectorff(EvtId,EvtId,
			   double t, double, double *a1f,
			   double *a2f, double *vf, double *a0f ){
  double q2=t;

  if(whichfit == 0) {
    *vf = 0;
    *a0f = 0;
    *a1f = 1;
    *a2f=0;
    
    return;
  };

  if( idVector == EvtPDL::getId("J/psi").getId() ) { // Bc -> J/psi
      if(whichfit == 1) { // SR form factor set from [Kiselev, hep-ph/0211021]
	double Mbc = 6.277, Mpsi=3.0967; // Experimental values
	double Mpole2 = 4.5*4.5, den = 1./(1.-q2/Mpole2);
	double FV = 0.11*den,
	FAp = -0.071*den,
	FA0 = 5.9*den,
	FAm = 0.12*den;
	*vf = (Mbc + Mpsi)*FV;
	*a2f = -(Mbc+Mpsi)*FAp;
	*a1f = FA0/(Mbc+Mpsi);
	*a0f = (q2*FAm + (Mbc+Mpsi)*(*a1f)-(Mbc-Mpsi)*(*a2f))/(2*Mpsi);    
	return;
      }
      else if(whichfit == 2) {  // form factor set from  [Ebert, hep-ph/0306306] 
	*vf = (0.49077824756158533 - 0.0012925655191347828*q2)/(1 - 0.06292520325875656*q2);
	*a0f = (0.4160345034630221 - 0.0024720095310225023*q2)/(1 - 0.061603451915567785*q2);
	*a1f = (0.4970212860605933 - 0.0067519730024654745*q2)/(1 - 0.050487026667172176*q2);
	*a2f = (0.7315284919705497 + 0.0014263826220727142*q2 -  0.0006946090066269195*q2*q2)/(1 - 0.04885587273651653*q2);
	return;
      };
  }
  else if(idVector == EvtPDL::getId("psi(2S)").getId()) { // Bc -> psi((2S)
      if(whichfit == 1) {
	////cout<<"BC2:: psi2S, Kiselev, q2="<<q2<<endl;
	double Mbc = 6.277, Mpsi=3.0967, Mpsi2S = 3.686, kappa = Mpsi/Mpsi2S; // Experimental values
	double Mpole2 = 4.5*4.5, den=1./(1.-q2/Mpole2);
	double FV = 0.11*den*kappa/3.1,
	FAp = -0.071*den*kappa/4.9,
	FA0 = 5.9*den*kappa/3.5,
	FAm = 0.12*den*kappa/2.3;
	*vf = (Mbc + Mpsi2S)*FV;
	*a2f = -(Mbc+Mpsi2S)*FAp;
	*a1f = FA0/(Mbc+Mpsi2S);
	*a0f = (q2*FAm + (Mbc+Mpsi2S)*(*a1f)-(Mbc-Mpsi2S)*(*a2f))/(2*Mpsi2S);  
	return;
      }
      else if(whichfit == 2) {
	////cout<<"BC2:: psi2S, Ebert, q2="<<q2<<endl;
	*vf  =  (0.24177223968739653 - 0.053589051007278135*q2)/(1 - 0.0977848994260899*q2);
	*a0f = (0.23996026570086615 - 0.03530198514007337*q2)/(1 - 0.09371162519983989*q2);
	*a1f = (0.17418379258849329 - 0.004129699022085851*q2*q2)/(1 + 0.06607665248402918*q2);
	*a2f = (0.1352376939112041 - 0.040361722565209444*q2 + 0.003343515369431853*q2*q2)/(1 - 0.1463698128333418*q2);
	return;
      };
  }
  else {
//   report(ERROR,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF2.\n";  
    ::abort();
  };
}

void EvtBCVFF2::getscalarff(EvtId, EvtId, double q2, 
	double, 
	double* fPlus, double* fMinus)
{
	if(idVector == EvtPDL::getId("chi_c0").getId() && whichfit==1) {
		*fPlus = 0.2735114574549966 + 0.02698177783795792*q2 + 0.0004438591230092838*pow(q2,2) + 0.00023394089715202517*pow(q2,3);
		*fMinus = -0.7888608969785541 - 0.10052789894587946*q2 + 0.0028122053015600472*pow(q2,2) - 0.0015941397716725178*pow(q2,3);
	}
	else if(idVector == EvtPDL::getId("chi_c0").getId() && whichfit==2) {
		*fPlus=-0.33991380592027914 - 0.040024442563952885*q2 + 0.0025230010943415723*pow(q2,2) - 0.0004656159966881966*pow(q2,3);
		*fMinus=0.7701897692262628 + 0.0881680421897945*q2 - 0.00879204111309994*pow(q2,2) + 0.0010886713551476884*pow(q2,3);
	}
//   report(ERROR,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF2.\n";  
}


void EvtBCVFF2::getaxialff( EvtId parent, EvtId daught, double q2, double mass, 
					double *hV1,  double *hV2, double *hV3, double *hA ) {
	if(idVector == EvtPDL::getId("chi_c1").getId() && whichfit==1) {
		*hV1=-0.07249517884406438 + 0.006386114470446252*q2 - 0.0016598461472340164*pow(q2,2) + 0.00043125401922510393*pow(q2,3);
		*hV2=-0.42708574288736423 - 0.08150607461417403*q2 + 0.00655250512925111*pow(q2,2) - 0.002089171872485879*pow(q2,3);
		*hV3=0.9721190093276546 + 0.14785518940757297*q2 - 0.010128808927653912*pow(q2,2) + 0.0033475878397617523*pow(q2,3);
		*hA=-1.0807476041824315 - 0.1250995787006434*q2 + 0.00032468802756664515*pow(q2,2) - 0.0016479191668232192*pow(q2,3);
	}
	if(idVector == EvtPDL::getId("chi_c1").getId() && whichfit==2) {
			*hV1=-0.24435448046190159 + 0.010331648277938579*q2 - 0.0009134031378323107*pow(q2,2) + 0.00022159216585685204*pow(q2,3);
			*hV2=-0.7103876821583478 - 0.0602774820822653*q2 + 0.004266553207625208*pow(q2,2) - 0.000565962924451258*pow(q2,3);
			*hV3=1.4934249265375008 + 0.08998682725105182*q2 - 0.006857963573522601*pow(q2,2) + 0.000711550359979238*pow(q2,3);
			*hA=1.2720231435712308 + 0.09003897439415781*q2 -  0.005137318246822382*pow(q2,2) + 0.0006148691245052408*pow(q2,3);
	}
}

void EvtBCVFF2::gettensorff(EvtId, EvtId, double q2, double, 
	double* tV,  double* tA1, double* tA2, double* tA3) {
	if(idVector == EvtPDL::getId("chi_c2").getId() && whichfit==1) {
		*tV=-0.7593279619789294 - 0.14373638874539618*q2 + 0.015428531108854306*pow(q2,2) - 0.0037404399820982963*pow(q2,3);
		*tA1=-0.5845892322932432 - 0.05829150603138564*q2 + 0.0003828588805250468*pow(q2,2) - 0.0007578111067025861*pow(q2,3);
		*tA2=0.0939174867271461 - 0.029819354099933058*q2 + 0.01373975840177015*pow(q2,2) - 0.0019524446633583286*pow(q2,3);
		*tA3=-0.0020369655557486233 + 0.0033852618872049524*q2 - 0.0008299552762151167*pow(q2,2) + 0.00006460810519319177*pow(q2,3);  
	}
	if(idVector == EvtPDL::getId("chi_c2").getId() && whichfit==2) {
		*tV=2.0192706161363225 + 0.10153860806900702*q2 - 0.0052357919072080665*pow(q2,2) + 0.0005919186313734745*pow(q2,3);
		*tA1=-0.9888948657318037 - 0.0683950179746637*q2 + 0.005431661340191046*pow(q2,2) - 0.0006880049067587609*pow(q2,3);
		*tA2=-0.2921648565036161 - 0.04757032020969248*q2 + 0.009924447477444942*pow(q2,2) - 0.0012187229555837636*pow(q2,3);
		*tA3=1.3421413598569518 + 0.10110866108376225*q2 - 0.0010856339717409947*pow(q2,2) + 0.0003398028563210277*pow(q2,3);
	}
}



void EvtBCVFF2::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
//   report(ERROR,"EvtGen") << "Not implemented :getbaryonff in EvtBCVFF2.\n";  
  ::abort();

}

void EvtBCVFF2::getdiracff(EvtId, EvtId, double, double, double*, double*,
			       double*, double*, double*, double*) {
  
//   report(ERROR,"EvtGen") << "Not implemented :getdiracff in EvtBCVFF2.\n";
  ::abort();

}

void EvtBCVFF2::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				double*, double*, double*, double*, double*, double*) {
  
//   report(ERROR,"EvtGen") << "Not implemented :getraritaff in EvtBCVFF2.\n";
  ::abort();

}

