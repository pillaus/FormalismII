#include "PentaquarkMatrix.h"
#include "Isobar.h"



void shuffle(double p[4])
{
  double x = p[3] * .001;
  for (int i =3; i > 0; i-- ) p[i] = p[i-1] * .001;
  p[0] = x;
  return;
}

int main ()
{



  // p->SetKinematics(4., 20.);
  // p->Print();

  double pmu_p[4] = {-2291.164 , -1528.241 , 55768.683 , 55836.745 };
  double pmu_m[4] = { -236.8123 , -2036.502 , 15932.630 , 16064.348 };
  double pK[4] = { 423.32345 , -880.7687 , 20447.629 , 20476.919 };
  double pp[4] = { -203.5984 , -820.3133 , 36541.246 , 36563.061 };
  shuffle(pmu_p);shuffle(pmu_m);shuffle(pK);shuffle(pp);

  double mpsi = 0., mb = 0., mp = 0., mK= 0.;
  for (int mu = 0; mu < 4; mu++)
  {
    mpsi += pow(pmu_p[mu] + pmu_m[mu], 2) * ( mu == 0 ? 1 : -1);
    mb += pow(pmu_p[mu] + pmu_m[mu] + pK[mu] + pp[mu], 2) * ( mu == 0 ? 1 : -1);
    mK += pow(pK[mu], 2) * ( mu == 0 ? 1 : -1);
    mp += pow(pp[mu], 2) * ( mu == 0 ? 1 : -1);
  }
  mpsi = sqrt(mpsi); mb = sqrt(mb); mK = sqrt(mK), mp = sqrt(mp);

   PentaquarkMatrix *p = new PentaquarkMatrix(mb, mpsi, mp, mK);

//   Isobar(double mass, double width, cd coupling unsigned int j, unsigned int spin, unsigned int L, bool eta, char channel)
  std::vector<Isobar> *Isobars = new std::vector<Isobar>();

  //Isobar iso;
  Isobars->push_back(Isobar(1.519, 0.156, 1., 3, 3, 2, true, true, 's', Isobar::BreitWigner));
//  Isobars->push_back(Isobar(1.,1., 1., 3, 3, 2, true, true, 's', Isobar::BreitWigner));
  //std::cout << Isobars->at(0).J() << std::endl;

  p->SetIsobars(Isobars);

  p->SetKinematics(pmu_p, pmu_m, pK, pp);

//p->Print();

   std::cout << p->Evaluate() << std::endl;

 return 0;
  cd uno = 1.;
  cd due = 2.i;
  std::cout << "cococo " << uno/due << std::endl;

  double _mass = 1., x= .5, _width=1., _coupling = 1.;
  cd den = _mass * _mass - x - 1i * _mass * _width;
  std::cout << "cococo " << _coupling / den  << std::endl;

  std::cout << "cococo " << Isobars->at(1).Evaluate(.5)  << std::endl;

  // std::cout << "cg " <<SpecialFunc::ClebschGordan(1,-1,1,1,1,0) << std::endl;

  return 0;

}
