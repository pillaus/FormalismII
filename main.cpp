#include "PentaquarkMatrix.h"
#include "Isobar.h"
#include "TLorentzVector.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"


void shuffle(double p[4])
{
  double x = p[3] * .001;
  for (int i =3; i > 0; i-- ) p[i] = p[i-1] * .001;
  p[0] = x;
  return;
}

int main (int argc, char **argv)
{


	int nn = 100000, seed = 0;

	for (int ii=0;ii<argc;ii++)
	{
		if (strcmp(argv[ii],"-n")==0) nn = atoi(argv[ii+1]);
		if (strcmp(argv[ii],"-s")==0) seed = atoi(argv[ii+1]);
	}

	if (argc == 1)
	{
		printf("USAGE: main -n [EVENTS] -s [SEED]\n");
	}
	printf("Generating %d events with seed %d... \n", nn, seed);



	TFile *f = new TFile("mc.root", "RECREATE");
	TTree *t= new TTree("ntuple", "ntuple");


  double mb = 5.61958, mpsi = 3.0969, mp = .938272081, mK = .493677;
  double pmu_p[4], pmu_m[4], pK[4], pp[4], ss, tt;

  double weight;
  for (int mu = 0; mu < 4; mu++)
  {
    char suff[10];
    switch (mu)
    {
      case 0: sprintf(suff, "E"); break;
      case 1: sprintf(suff, "Px"); break;
      case 2: sprintf(suff, "Py"); break;
      case 3: sprintf(suff, "Pz"); break;
    }
    t->Branch(Form("mup_%s", suff), &pmu_p[mu]);
    t->Branch(Form("mum_%s", suff), &pmu_m[mu]);
    t->Branch(Form("K_%s", suff), &pK[mu]);
    t->Branch(Form("p_%s", suff), &pp[mu]);


  }
  t->Branch("s", &ss); t->Branch("t", &tt);
  t->Branch("weight", &weight);

  PentaquarkMatrix *p = new PentaquarkMatrix(mb, mpsi, mp, mK);

//   Isobar(double mass, double width, cd coupling unsigned int j, unsigned int spin, unsigned int L, bool eta, char channel)
  std::vector<Isobar> *Isobars = new std::vector<Isobar>();

 //Isobar iso; 1.519
  Isobars->push_back(Isobar(1.519, 0.156, 1., 3, 3, 2, true, true, 's', Isobar::BreitWigner));
//  Isobars->push_back(Isobar(1.,1., 1., 3, 3, 2, true, true, 's', Isobar::BreitWigner));
 //std::cout << Isobars->at(0).J() << std::endl;

  p->SetIsobars(Isobars);


  TLorentzVector Lambda_b(0.0, 0.0, 0.0, mb);
  //(Momentum, Energy units are Gev/C, GeV)
  double masses[3] = { mpsi, mp, mK } ;
  double masses2[2] = { .105, .105 };
  TGenPhaseSpace event, event2;
  event.SetDecay(Lambda_b, 3, masses);
  for (int n=0;n<nn ;n++)
  {
     weight = event.Generate();
     TLorentzVector *tpPsi = event.GetDecay(0);
     TLorentzVector *tpP    = event.GetDecay(1);
     TLorentzVector *tpK    = event.GetDecay(2);

     TLorentzVector tpLambdastar = *tpP + *tpK;
     TLorentzVector tpPc = *tpP + *tpPsi;


     // std::cout << tpPsi->E() << std::endl;
     event2.SetDecay(*tpPsi, 2, masses2);
     weight *= event2.Generate();
     TLorentzVector *tpmup = event2.GetDecay(0);
     TLorentzVector *tpmum = event2.GetDecay(1);

     pmu_p[0] = tpmup->E(); pmu_p[1] = tpmup->Px(); pmu_p[2] = tpmup->Py(); pmu_p[3] = tpmup->Pz();
     pmu_m[0] = tpmum->E(); pmu_m[1] = tpmum->Px(); pmu_m[2] = tpmum->Py(); pmu_m[3] = tpmum->Pz();
     pp[0] = tpP->E(); pp[1] = tpP->Px(); pp[2] = tpP->Py(); pp[3] = tpP->Pz();
     pK[0] = tpK->E(); pK[1] = tpK->Px(); pK[2] = tpK->Py(); pK[3] = tpK->Pz();

     ss = tpLambdastar.M2();
     tt = tpPc.M2();


     // printf("%lf %lf %lf %lf\n", pmu_p[0], pmu_p[1], pmu_p[2], pmu_p[3]);
     // break;

  // p->SetKinematics(4., 20.);
  // p->Print();

  // double pmu_p[4] = {-2291.164 , -1528.241 , 55768.683 , 55836.745 };
  // double pmu_m[4] = { -236.8123 , -2036.502 , 15932.630 , 16064.348 };
  // double pK[4] = { 423.32345 , -880.7687 , 20447.629 , 20476.919 };
  // double pp[4] = { -203.5984 , -820.3133 , 36541.246 , 36563.061 };
  // shuffle(pmu_p);shuffle(pmu_m);shuffle(pK);shuffle(pp);

  // double mpsi = 0., mb = 0., mp = 0., mK= 0.;
  // for (int mu = 0; mu < 4; mu++)
  // {
  //   mpsi += pow(pmu_p[mu] + pmu_m[mu], 2) * ( mu == 0 ? 1 : -1);
  //   mb += pow(pmu_p[mu] + pmu_m[mu] + pK[mu] + pp[mu], 2) * ( mu == 0 ? 1 : -1);
  //   mK += pow(pK[mu], 2) * ( mu == 0 ? 1 : -1);
  //   mp += pow(pp[mu], 2) * ( mu == 0 ? 1 : -1);
  // }
  // mpsi = sqrt(mpsi); mb = sqrt(mb); mK = sqrt(mK), mp = sqrt(mp);


    p->SetKinematics(pmu_p, pmu_m, pK, pp);
    // std::cout << weight << " " << p->Evaluate() << std::endl;

    std::complex<double> iso = Isobars->at(0).Evaluate(ss);
    // std::cout << abs(iso)*abs(iso) << std::endl;
    weight *= p->Evaluate();

    // std::cout << weight << std::endl;


    t->Fill();
  }

  t->Write();
  f->Close();
// p->Print();

   // std::cout << p->Evaluate() << std::endl;

 return 0;


}
