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


   PentaquarkMatrix *p = new PentaquarkMatrix();
  // p->SetKinematics(4., 20.);
  // p->Print();

  double pmu_p[4] = {-2291.164 , -1528.241 , 55768.683 , 55836.745 };
  double pmu_m[4] = { -236.8123 , -2036.502 , 15932.630 , 16064.348 };
  double pK[4] = { 423.32345 , -880.7687 , 20447.629 , 20476.919 };
  double pp[4] = { -203.5984 , -820.3133 , 36541.246 , 36563.061 };
  shuffle(pmu_p);shuffle(pmu_m);shuffle(pK);shuffle(pp);


//   Isobar(double mass, double width, cd coupling unsigned int j, unsigned int spin, unsigned int L, bool eta, char channel)
  std::vector<Isobar> *Isobars = new std::vector<Isobar>();
  Isobar iso(1.519, 0.156, 1., 3, 3, 2, true, true, 's');
  //Isobar iso;
  Isobars->push_back(Isobar(1.519, 0.156, 1., 3, 3, 2, true, true, 's'));
  std::cout << Isobars->at(0).J() << std::endl;

  p->SetIsobars(Isobars);

  p->SetKinematics(pmu_p, pmu_m, pK, pp);

  std::cout << p->Evaluate() << std::endl;



  return 0;

}
