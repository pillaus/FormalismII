/**
    PentaquarkMatrix.cpp
    Calculate the matrices in Appendix G of 1805.02113

    @author Alessandro Pilloni
    @version 1.0 2018-05-13
*/
#include <math.h>
#include <iostream>
#include <vector>
#include "Isobar.h"


#ifndef PENTA_H
#define PENTA_H


class PentaquarkMatrix {

private:
  double s, t, u, rs, ru, ps, qs, pu, qu, ps2, qs2, pu2, qu2, Eb_s, Ep_s, Epsi_s, num_s, num_u, Eb_u, Ep_u, Epsi_u, u2;
  const double sp, sm, up, um, mp2, mpsi2, mb2, mK2;
  const double mb, mpsi, mK, mp;

  cd Structures[6][4][4][4];
  double Matrices[4][6][6];
  double lepton[4][4];
  static const double metric[4];
  static const cd gamma[6][4][4];
  static const double over4PI;
  int g5(int i) { if (i > 1) return i - 2; else return i+2; };
  std::vector<Isobar> *Isobars;


  void slash(double p[4], cd ret[4][4]);
  void slashm(double p[4], double m, cd ret[4][4]);


public:

  PentaquarkMatrix(double _mb = 5.61958, double _mpsi = 3.0969, double _mp = .938272081, double _mK = .493677);

  void SetKinematics(double pmup[4], double pmum[4], double pK[4], double pp[4]);

  void SetKinematics(double _s, double _u, bool check = true);

  void SetIsobars(std::vector<Isobar> *_Isobars) { Isobars = _Isobars; }

  double Evaluate();

  double GetPCschannel(int i, int j);
  double GetPVschannel(int i, int j);
  double GetPCuchannel(int i, int j);
  double GetPVuchannel(int i, int j);

  double GetMatrix(int k, int i, int j)
  {
    switch (k)
    {
      case 0: return GetPCschannel(i,j);
      case 1: return GetPVschannel(i,j);
      case 2: return GetPCuchannel(i,j);
      case 3: return GetPVuchannel(i,j);
    }
  }


  ~PentaquarkMatrix() {};
  void Print();

};

#endif
