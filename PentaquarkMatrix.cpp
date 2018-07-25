/**
    PentaquarkMatrix.cpp
    Calculate the matrices in Appendix G of 1805.02113

    @author Alessandro Pilloni
    @version 1.0 2018-05-13
*/
#include "PentaquarkMatrix.h"
  /**
      Constructor. It gets the masses as optional input. Default from LivePDG, 2018-05-13

      @param _mb Lambda_b mass
      @param _mpsi J/psi mass
      @param _mp proton mass
      @param _mK Kaon mass

  */
const double PentaquarkMatrix::metric[4] = { 1., -1., -1., -1. };

const cd PentaquarkMatrix::gamma[6][4][4] = {
//gamma0
        { { 1., 0., 0., 0. }, { 0., 1., 0., 0. }, { 0., 0., -1., 0. }, { 0., 0., 0., -1. } },
//gamma1
        { { 0., 0., 0., 1. }, { 0., 0., 1., 0. }, { 0., -1., 0., 0. }, { -1., 0., 0., 0. } },
//gamma2
        { { 0., 0., 0., -1.i }, { 0., 0., 1.i, 0. }, { 0.,  1.i, 0., 0. }, { -1.i, 0., 0., 0. } },
//gamma3
        { { 0., 0., 1., 0. }, { 0., 0., 0., -1. }, { -1., 0., 0., 0. }, { 0., 1., 0., 0. } },
//gamma4 = gamma0
        { { 1., 0., 0., 0. }, { 0., 1., 0., 0. }, { 0., 0., -1., 0. }, { 0., 0., 0., -1. } },
//gamma5
        { { 0., 0., 1., 0. }, { 0., 0., 0., 1. }, { 1., 0., 0., 0. }, { 0., 1., 0., 0. } }
};

const double PentaquarkMatrix::over4PI = 1./std::atan(1.)/16.;

  PentaquarkMatrix::PentaquarkMatrix(double _mb, double _mpsi, double _mp, double _mK) :
    mb(_mb), mpsi(_mpsi), mp(_mp), mK(_mK), mb2(_mb*_mb),
    mpsi2(_mpsi*_mpsi), mp2(_mp*_mp), mK2(_mK*_mK),
    sp(pow(_mb + _mpsi, 2)), sm(pow(_mb - _mpsi, 2)),
    up(pow(_mpsi + _mp, 2)), um(pow(mpsi - mp, 2)) {  };

void PentaquarkMatrix::slash(double p[4], cd ret[4][4])
{
  for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
  {
    ret[i][j] = 0.;
    for (int mu = 0; mu < 4; mu++) ret[i][j]  += metric[mu]*p[mu]*gamma[mu][i][j];
   }
 }
void PentaquarkMatrix::slash(double p[4], double m, cd ret[4][4])
{
  for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
  {
    ret[i][j] = (i == j) ? m : 0.;
    for (int mu = 0; mu < 4; mu++) ret[i][j]  += metric[mu]*p[mu]*gamma[mu][i][j];
  }
}

    /**
        Sets the kinematics)

        @param pmu_p[4] mu+ four-momentum
        @param pmu_m[4] mu- four-momentum
        @param pK[4] kaon four-momentum
        @param pp[4] proton four-momentum
    */
  void PentaquarkMatrix::SetKinematics(double pmu_p[4], double pmu_m[4], double pK[4], double pp[4])
{
  double ppsi[4], pb[4];
  for (int mu =0; mu < 4; mu++) { ppsi[mu] = pmu_p[mu] + pmu_m[mu]; pb[mu] = ppsi[mu] + pp[mu] + pK[mu]; }
  // for (int mu =0; mu < 4; mu++) for (int nu =0; nu < 4; nu++)
  // {
  //   lepton[mu][nu] = pmu_p[mu]*pmu_m[nu] + pmu_p[nu]*pmu_m[mu];
  //   lepton[mu][nu] -= (mu == nu) ? metric[mu]*mpsi2/2. : 0.;
  // }
  for (int mu =0; mu < 4; mu++) for (int nu =0; nu < 4; nu++)
  {
    lepton[mu][nu]  = ppsi[mu]*ppsi[nu]/mpsi2;
    lepton[mu][nu] -= (mu == nu) ? metric[mu] : 0.;
  }
          // To add muons

  cd ppsislash[4][4];
  slash(ppsi, ppsislash);
  slash(pb, mb, pbm);
  slash(pp, mp, ppm);

  double _s = 0., _u = 0.;
  for (int mu = 0; mu < 4; mu++) { _s += metric[mu]*pow(pp[mu] + pK[mu], 2); _u += metric[mu]*pow(pp[mu] + ppsi[mu], 2);  }
  PentaquarkMatrix::SetKinematics(_s, _u);

  // Lorentz-Dirac structures. PC
  for (int i =0; i < 4; i++) for (int j =0; j < 4; j++) for (int mu =0; mu < 4; mu++)
  {
    Structures[0][mu][i][j] = gamma[5][i][j] * pb[mu];
    Structures[1][mu][i][j] = gamma[5][i][j] * pp[mu];
    Structures[2][mu][i][j] = ppsislash[PentaquarkMatrix::g5(i)][j] * pb[mu]; // left multiplication of gamma5 is the swapping of the rows
    Structures[3][mu][i][j] = ppsislash[PentaquarkMatrix::g5(i)][j] * pp[mu];
    Structures[4][mu][i][j] = gamma[mu][PentaquarkMatrix::g5(i)][j];
    Structures[5][mu][i][j] = 0.;
    for (int k = 0; k < 4; k++) Structures[5][mu][i][j] += gamma[mu][PentaquarkMatrix::g5(i)][k]*ppsislash[k][j];
  }

  for (int k =0; k < 4; k++) for (int i =0; i < 6; i++) for (int j =0; j < 6; j++)  Matrices[k][i][j] = PentaquarkMatrix::GetMatrix(k,i,j);
  return;
}

cd PentaquarkMatrix::Evaluate()
{
  cd FF[4][6] = {}; // First index is channel and PC/PV. second index the 6 guys




  unsigned int size = Isobars->size();
  for (unsigned int n = 0; n < size; n++)
  {
    int j = Isobars->at(n).J(); // remember j is twice the spin
    int S = Isobars->at(n).Spin();
    int L = Isobars->at(n).L();
    int PC = Isobars->at(n).PC()  ?1:-1;
    int eta  = Isobars->at(n).Eta() ?1:-1;
    double valP = 1., valM = 1., fact;
    // fact are the common values for both etabars. varP is for etabar +, varM for etabar -
    int lambdamin = (j == 1) ? 0 : -1;
    for (int lambda = 1; lambda >= lambdamin; lambda--)
    {
      int kl = 1 - lambda;
      if (Isobars->at(n).Channel() == 's' )
      {

        fact = over4PI * (j + 1) * pow(ps*qs, (lambda == -1) ? (j - 3)/2 : (j - 1)/2 )
                    * SpecialFunc::ClebschGordan(1,1,2, -2*lambda,S, 1 - 2*lambda) * SpecialFunc::ClebschGordan(S, 1 - 2*lambda, L, 0, j, 1 - 2*lambda);

        if (lambda == -1)
        {
          fact *=  mpsi*mp/rs*eta;
          valM *= -1.;
        }
        else if (lambda == 0 && (j != 1 || eta != PC))
        {
          fact *= Epsi_s / mpsi;
        }
        if (eta > 0)
        {
          // eta = +, etabar = +
          valP *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, true, costheta_s);
          if (j == 1 && PC > 0) valP *= ps2*s/mpsi2;
          // eta = +, etabar = -
          valM *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, false, costheta_s);
          valM /= (PC > 0) ? ps*mp/qs/mpsi :  mpsi*mp/s/ps/qs;
        }
        else
        {
          // eta = -, etabar = -
          valM *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, true, costheta_s);
          if (j == 1 && PC < 0) valM *= ps2*s/mpsi2;
          // eta = -, etabar = +
          valP *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, false, costheta_s);
          valP *= (PC > 0) ? ps*mp/qs/mpsi :  mpsi*mp/s/ps/qs;
        }
        if (2*L - j + 2 + eta*PC > 0) fact *= ps2;

//        fact *= sqrt((Eb_s + mb)/(Ep_s + mp));
        cd cfact = Isobars->at(n).Evaluate(s) * fact;
        FF[PC ? 0 : 1][kl  ] += cfact*valP;
        FF[PC ? 0 : 1][kl+3] += cfact*valM;

      }
      else if (Isobars->at(n).Channel() == 'u' )
      {
          fact = over4PI * (j + 1) * pow(pu*qu, (lambda == -1) ? (j - 3)/2 : (j - 1)/2 )
                      * SpecialFunc::ClebschGordan(1,1,2, -2*lambda,S, 1 - 2*lambda) * SpecialFunc::ClebschGordan(S, 1 - 2*lambda, L, 0, j, 1 - 2*lambda);

          if (lambda == -1)
          {
            fact *=  mb*mp/ru*eta;
            valM *= -1.;
          }
          if (eta > 0)
          {
            fact *= (Ep_u + mp)/2./mp;
            // eta = +, etabar = +
            valP *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, true, costheta_u);
            if (j == 1) valP *= qu2*u/mp2;
            // eta = +, etabar = -
            valM *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, false, costheta_u);
            valM /= (PC > 0) ? qu*mb/pu/mp/(u - um) :  qu*mb*pu*mp/(u - um);
          }
          else
          {
            // eta = -, etabar = -
            valM *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, true, costheta_u);
            // eta = -, etabar = +
            valP *= SpecialFunc::WignerDhat(j, 1, 1 - 2*lambda, false, costheta_u);
            valP *= (PC > 0) ? qu*mb/pu/mp/(u - um) :  qu*mb*pu*mp/(u - um);
          }
          if (2*L - j + 2 + eta > 0) fact *= qu2;

//          fact *= sqrt((Eb_s + mb));
          cd cfact = Isobars->at(n).Evaluate(u) * fact;
          FF[PC ? 2 : 3][kl  ] += cfact*valP;
          FF[PC ? 2 : 3][kl+3] += cfact*valM;

        }
      }
    }
    cd CC[4][6] = {};
    for (int i =0; i < 4; i++) for (int j =0; j < 6; j++) for (int k =0; k < 6; k++) CC[i][j] += Matrices[i][j][k]*FF[i][k];

    cd CM[4][4][4] = {};
    for (int i =0; i < 4; i++) for (int j =0; j < 4; j++) for (int mu =0; mu < 4; mu++) for (int k = 0; k < 6 ; k++)
    {
      // s- and u-channel, PC
      CM[mu][i][j] += Structures[k][mu][i][j] * (CC[0][k] + CC[2][k]);
      // s- and u-channel, PV
      CM[mu][i][j] += Structures[k][mu][g5(i)][j] * (CC[1][k] + CC[3][k]);
    }

    cd ret = 0.;

    for (int i =0; i < 4; i++) for (int j =0; j < 4; j++) for (int k =0; k < 4; k++) for (int l =0; l < 4; l++) for (int mu =0; mu < 4; mu++) for (int nu =0; nu < 4; nu++)
    {
      ret += pbm[i][j] * CM[mu][j][k] * pbm[k][l] * std::conj(CM[nu][i][l]) * g0fy(l,i) *lepton[mu][nu];
    }






    // else if (!Isobars[n].PC() && Isobars[n].Channel() == 's' ) sel = 1;
    // else if (Isobars[n].PC() && Isobars[n].Channel() == 'u' ) sel = 2;
    // else if (!Isobars[n].PC() && Isobars[n].Channel() == 'u' ) sel = 3;




  return ret;
}
/**
    Sets the kinematics and check whether it falls in the Dalitz (optional)

    @param _s pK invariant mass squared
    @param _u J/psi p invariant mass squared
    @param check flag to check the kinematics

*/
  void PentaquarkMatrix::SetKinematics(double _s, double _u, bool check)
  {
    s = _s;
    u = _u;
    t = mb2 + mK2 + mp2 + mpsi2 - s - u;
    rs = sqrt(s);
    ru = sqrt(u);
    u2 = u*u;
    ps2 = (s - sp)*(s - sm)/4./s;
    ps = sqrt(ps2);
    qs2 = (s - pow(mp + mK, 2))*(s - pow(mp - mK, 2))/4./s;
    qs = sqrt(qs2);
    num_s = (s*(t - u) + (mb2 - mpsi2)*(mp2 - mK2))/4./s;
    Epsi_s = (s + mpsi2 - mb2)/2./rs;
    Eb_s = (s - mpsi2 + mb2)/2./rs;
    Ep_s = (s + mp2 - mK2)/2./rs;

    pu2 = (u - pow(mb + mK, 2))*(u - pow(mb - mK, 2))/4./u;
    pu = sqrt(pu2);
    qu2 = (u - up)*(u - um)/4./u;
    qu = sqrt(qu2);
    num_u = (u*(t - s) + (mb2 - mK2)*(mp2 - mpsi2))/4./u;
    Epsi_u= (u + mpsi2 - mp2)/2./ru;
    Eb_u = (u - mK2 + mb2)/2./ru;
    Ep_u = (u + mp2 - mpsi2)/2./ru;

    costheta_s = (s*(t - u) + (mb2 - mpsi2)*(mp2 - mK2))/(4.*s*ps*qs);
    costheta_u = (u*(t - s) + (mb2 - mK2)*(mp2 - mpsi2))/(4.*u*pu*qu);

    if (check && !(s <= sm && s >= pow(mp + mK, 2) && u >= up && u <= pow(mb - mK, 2) && num_s <= ps*qs && num_s >= -ps*qs && num_u <= pu*qu && num_u >= -pu*qu ))
    {
      std::cout << "Warning: Kinematics out of the Dalitz region" << std::endl;
    }

  }

  /**
      Calculate the matrix for the parity-conserving s-channel (Lambda*)

      @param i row of the matrix, 0->5
      @param j column of the matrix, 0->5
      @return output of the matrix

  */
  double PentaquarkMatrix::GetPCschannel(int i, int j)
  {
    double ret;
    if (i == 0 && j == 0) ret = mpsi2*mpsi/(4.*mp*ps2*s);
    else if (i == 0 && j == 1) ret = -(mpsi2*(mb + rs))/(2.*sqrt(2.)*mp*ps2*s);
    else if (i == 0 && j == 2) ret = -(-((Eb_s + mb)*(mp*mpsi2 + Ep_s*(2*Epsi_s*mb + mb2 - s))) + num_s*(2*Epsi_s*mb - mb2 + s))/(4.*mp2*ps2*rs);
    else if (i == 0 && j == 3) ret = -((Eb_s + mb)*mpsi2)/(4.*(Ep_s + mp)*ps2*s);
    else if (i == 0 && j == 4) ret =     ((Eb_s + mb)*mpsi*(-mb + rs))/(2.*sqrt(2.)*(Ep_s + mp)*ps2*s);
    else if (i == 0 && j == 5) ret = ((Eb_s + mb)*num_s*(-(mb*(2*Epsi_s + mb)) + s) + ps2*(mp*mpsi2 + Ep_s*(2*Epsi_s*mb - mb2 + s)))/(4.*mp*(Ep_s + mp)*mpsi*ps2*rs);
    else if (i == 1 && j == 2) ret = (mb + rs)/(2.*mp2);
    else if (i == 1 && j == 5) ret = ((Eb_s + mb)*(mb - rs))/(2.*mp*(Ep_s + mp)*mpsi);
    else if (i == 2 && j == 0) ret = (mpsi*(mb + rs))/(4.*mp*ps2*s);
    else if (i == 2 && j == 1) ret = -mpsi2/(2.*sqrt(2.)*mp*ps2*s);
    else if (i == 2 && j == 2) ret =     ((-Eb_s - mb)*(Epsi_s*(-Ep_s + mp) + (Eb_s - mb)*(Ep_s + mp)) + (Eb_s - Epsi_s + mb)*num_s)/(4.*mp2*ps2*rs);
    else if (i == 2 && j == 3) ret = ((Eb_s + mb)*(-mb + rs))/(4.*(Ep_s + mp)*ps2*s);
    else if (i == 2 && j == 4) ret = -((Eb_s + mb)*mpsi)/(2.*sqrt(2.)*(Ep_s + mp)*ps2*s);
    else if (i == 2 && j == 5) ret =     -((Eb_s + mb)*(-Eb_s + Epsi_s + mb)*num_s + ps2*(Ep_s*(Eb_s - Epsi_s + mb) - mp*(mb + rs)))/(4.*mp*(Ep_s + mp)*mpsi*ps2*rs);
    else if (i == 3 && j == 2) ret = 1./(2.*mp2);
    else if (i == 3 && j == 5) ret = (Eb_s + mb)/(2*Ep_s*mp*mpsi + 2*mp2*mpsi);
    else if (i == 4 && j == 0) ret =    ((Eb_s + mb)*mpsi*(mb - rs))/(4.*mp*ps2*rs);
    else if (i == 4 && j == 2) ret = -((Eb_s + mb)*num_s*(-mb + rs) + (Ep_s - mp)*ps2*(mb + rs))/(4.*mp2*ps2);
    else if (i == 4 && j == 3) ret = -((mb + rs)/(4*Ep_s*rs + 4*mp*rs));
    else if (i == 4 && j == 5) ret = -((Eb_s + mb)*(Ep_s + mp)*(-mb + rs) + num_s*(mb + rs))/(4.*mp*(Ep_s + mp)*mpsi);
    else if (i == 5 && j == 0) ret =    -((Eb_s + mb)*mpsi)/(4.*mp*ps2*rs);
    else if (i == 5 && j == 2) ret = -(Eb_s*num_s + mb*num_s - Ep_s*ps2 + mp*ps2)/(4.*mp2*ps2);
    else if (i == 5 && j == 3) ret = 1./(4*Ep_s*rs + 4*mp*rs);
    else if (i == 5 && j == 5) ret = (-((Eb_s + mb)*(Ep_s + mp)) + num_s)/(4.*mp*(Ep_s + mp)*mpsi);
    else return 0.;

    return -ret*sqrt((Ep_s + mp)/(Eb_s + mb));


  }

  /**
      Calculate the matrix for the parity-violating s-channel (Lambda*)

      @param i row of the matrix, 0->5
      @param j column of the matrix, 0->5
      @return output of the matrix

  */
  double PentaquarkMatrix::GetPVschannel(int i, int j)
  {
    double ret;
    if (i == 0 && j == 0) ret = ((Eb_s + mb)*mpsi2)/(4.*mp*ps2*rs);
    else if (i == 0 && j == 1) ret = ((Eb_s + mb)*mpsi*(mb - rs))/(2.*sqrt(2.)*mp*ps2*rs);
    else if (i == 0 && j == 2) ret = -((Eb_s + mb)*num_s*(-(mb*(2*Epsi_s + mb)) + s) + ps2*(-(mp*mpsi2) + Ep_s*(2*Epsi_s*mb - mb2 + s)))/(4.*mp2*mpsi*ps2);
    else if (i == 0 && j == 3) ret =     -(mpsi2*mpsi/(4*Ep_s*ps2*rs*s + 4*mp*ps2*rs*s));
    else if (i == 0 && j == 4) ret = (mpsi2*(mb + rs))/(2.*sqrt(2.)*(Ep_s + mp)*ps2*rs*s);
    else if (i == 0 && j == 5) ret = -((Eb_s + mb)*(-(mp*mpsi2) + Ep_s*(2*Epsi_s*mb + mb2 - s)) + num_s*(-2*Epsi_s*mb + mb2 - s))/(4.*mp*(Ep_s + mp)*ps2*s);
    else if (i == 1 && j == 2) ret = ((Eb_s + mb)*rs*(-mb + rs))/(2.*mp2*mpsi);
    else if (i == 1 && j == 5) ret = -((mb + rs)/(2*Ep_s*mp*rs + 2*mp2*rs));
    else if (i == 2 && j == 0) ret =    ((Eb_s + mb)*(mb - rs))/(4.*mp*ps2*rs);
    else if (i == 2 && j == 1) ret = ((Eb_s + mb)*mpsi)/(2.*sqrt(2.)*mp*ps2*rs);
    else if (i == 2 && j == 2) ret = ((Eb_s + mb)*(-Eb_s + Epsi_s + mb)*num_s + ps2*(Ep_s*(Eb_s - Epsi_s + mb) + mp*(mb + rs)))/(4.*mp2*mpsi*ps2);
    else if (i == 2 && j == 3) ret =     -((mb*mpsi + mpsi*rs)/(4*Ep_s*ps2*rs*s + 4*mp*ps2*rs*s));
    else if (i == 2 && j == 4) ret = mpsi2/(sqrt(2.)*(2*Ep_s*ps2*rs*s + 2*mp*ps2*rs*s));
    else if (i == 2 && j == 5) ret = -((Eb_s - Epsi_s + mb)*num_s + (Eb_s + mb)*(Ep_s*(-Eb_s + Epsi_s + mb) + mp*(-mb + rs)))/(4.*mp*(Ep_s + mp)*ps2*s);
    else if (i == 3 && j == 2) ret = -((Eb_s + mb)*rs)/(2.*mp2*mpsi);
    else if (i == 3 && j == 5) ret = -1./(2.*(Eb_s + Epsi_s)*mp*(Ep_s + mp));
    else if (i == 4 && j == 0) ret = (mb + rs)/(4.*mp);
    else if (i == 4 && j == 2) ret = (rs*((Eb_s + mb)*(Ep_s - mp)*(-mb + rs) + num_s*(mb + rs)))/(4.*mp2*mpsi);
    else if (i == 4 && j == 3) ret = ((Eb_s + mb)*mpsi*(-mb + rs))/(4.*(Ep_s + mp)*ps2*s);
    else if (i == 4 && j == 5) ret =     -((Eb_s + mb)*num_s*(mb - rs) - (Ep_s + mp)*ps2*(mb + rs))/(4.*mp*(Ep_s + mp)*ps2*rs);
    else if (i == 5 && j == 0) ret = -1./(4.*mp);
    else if (i == 5 && j == 2) ret = -((-(Ep_s*mb) + mb*mp + Eb_s*(-Ep_s + mp) + num_s)*rs)/(4.*mp2*mpsi);
    else if (i == 5 && j == 3) ret = (Eb_s*mpsi + mb*mpsi)/(4*Ep_s*ps2*s + 4*mp*ps2*s);
    else if (i == 5 && j == 5) ret =     (Eb_s*num_s + mb*num_s - Ep_s*ps2 - mp*ps2)/(4*Ep_s*mp*ps2*rs + 4*mp2*ps2*rs);
    else return 0.;

    return -ret*sqrt((Ep_s + mp)/(Eb_s + mb));


  }

  /**
      Calculate the matrix for the parity-conserving u-channel (P_c)

      @param i row of the matrix, 0->5
      @param j column of the matrix, 0->5
      @return output of the matrix

  */
  double PentaquarkMatrix::GetPCuchannel(int i, int j)
  {
    double ret;
         if (i == 0 && j == 2) ret = -(mp + ru)/(2.*mb2);
    else if (i == 0 && j == 5) ret = ((Ep_u + mp)*(-mp + ru))/(4.*mb2*mp*(u - um));
    else if (i == 1 && j == 0) ret = -(mp*mpsi2)/(4.*mb*qu2*u);
    else if (i == 1 && j == 1) ret = -(mp*mpsi*(mp + ru))/(2.*sqrt(2.)*mb*qu2*u);
    else if (i == 1 && j == 2) ret =     (2*Ep_u*mb*num_u*(-mp + ru) + 2*mb2*qu2*(mp + ru) + mb*num_u*(-2*mp2 + mpsi2 + 2*mp*ru))/(4.*mb*mb2*qu2*ru);
    else if (i == 1 && j == 3) ret = (Ep_u*mpsi2 + mp*mpsi2)/(8*mb*qu2*u2 - 8*mb*qu2*u*um);
    else if (i == 1 && j == 4) ret = ((Ep_u + mp)*mpsi*(-mp + ru))/(4.*sqrt(2.)*mb*qu2*u*(u - um));
    else if (i == 1 && j == 5) ret =     -(mp2*mpsi2*num_u + 2*mb*mp*mp2*qu2 + 2*mb*mp*mpsi2*qu2 + 2*mp2*num_u*qu2 + 3*mpsi2*num_u*qu2 + 4*mb*mp*qu2*qu2 + 4*num_u*qu2*qu2 +         Epsi_u*(4*Ep_u*Ep_u*mb*qu2 + mb*mpsi2*qu2 + mp*num_u*(mpsi2 + 2*qu2) + Ep_u*(mpsi2*num_u + 4*(mb*mp + num_u)*qu2)) - 2*mb*mp2*qu2*ru + mb*mpsi2*qu2*ru + Ep_u*(mb*qu2*(2*mp2 + 3*mpsi2 + 4*qu2) + mp*(mpsi2*num_u + 2*qu2*(num_u - mb*ru))))/     (8.*mb2*mp*qu2*u*(u - um));
    else if (i == 2 && j == 2) ret = 1./(2.*mb2);
    else if (i == 2 && j == 5) ret = (Ep_u + mp)/(4*mb2*mp*u - 4*mb2*mp*um);
    else if (i == 3 && j == 0) ret = (mp*       (2*mp2*mp2 + mp2*mpsi2 + Ep_u*Ep_u*(2*mp2 + mpsi2) + 5*mp2*qu2 + 2*qu2*qu2 + Epsi_u*(Ep_u + mp)*(3*Ep_u*mp + 3*mp2 + 2*qu2) + Ep_u*mp*(4*mp2 + 2*mpsi2 + 5*qu2)))/(4.*mb*Ep_u + mp2*qu2*ru*u);
    else if (i == 3 && j == 1) ret =     (mp*mpsi)/(2.*sqrt(2.)*mb*qu2*u);
    else if (i == 3 && j == 2) ret = (Epsi_u*mp*num_u + Ep_u*mp*num_u + mp2*num_u - mpsi2*num_u - 2*mb*qu2*ru)/(4.*mb2*qu2*u);
    else if (i == 3 && j == 3) ret =     (Epsi_u*Ep_u*mp + Epsi_u*mp2 + Ep_u*mpsi2 + mp*mpsi2 + mp*qu2 + 2*qu2*ru)/(8*mb*qu2*ru*u2 - 8*mb*qu2*ru*u*um);
    else if (i == 3 && j == 4) ret = ((Ep_u + mp)*mpsi)/(4.*sqrt(2.)*mb*qu2*u*(u - um));
    else if (i == 3 && j == 5) ret =     -((Epsi_u*Ep_u*mp*num_u + Epsi_u*mp2*num_u + Ep_u*mpsi2*num_u + mp*mpsi2*num_u - 2*Epsi_u*Ep_u*mb*qu2 - 2*Ep_u*Ep_u*mb*qu2 + 2*mb*mp2*qu2 - 2*mb*mpsi2*qu2 + mp*num_u*qu2)/(8*mb2*mp*qu2*u2 - 8*mb2*mp*qu2*u*um));
    else if (i == 4 && j == 0) ret =    (-(Epsi_u*mp2*(Ep_u + mp)) + Ep_u*mp*mpsi2 + mp2*mpsi2 + mp2*qu2)/(4.*mb*(mp2 - mpsi2)*qu2);
    else if (i == 4 && j == 2) ret = ((Ep_u + mp)*num_u*(mp - ru))/(4.*mb2*qu2);
    else if (i == 4 && j == 3) ret = -((mp + ru)/(8*mb*ru*u - 8*mb*ru*um));
    else if (i == 4 && j == 5) ret =     (2*Ep_u*mb2*mp + 2*mb2*mp2 - mb*mp*num_u - 2*Ep_u*mb2*ru - 2*mb2*mp*ru - mb*num_u*ru)/(8*mb*mb2*mp*u - 8*mb*mb2*mp*um);
    else if (i == 5 && j == 0) ret =    -(mp*(Ep_u + mp))/(4.*mb*qu2*ru);
    else if (i == 5 && j == 2) ret = -((Ep_u + mp)*num_u)/(4.*mb2*qu2);
    else if (i == 5 && j == 3) ret = 1./(8*mb*ru*u - 8*mb*ru*um);
    else if (i == 5 && j == 5) ret = -((2*Ep_u*mb2 + 2*mb2*mp - mb*num_u)/(8*mb*mb2*mp*u - 8*mb*mb2*mp*um));


    else return 0.;

    return ret*sqrt((u - um)/(Eb_u + mb)/(Ep_u + mp));

  }

  /**
      Calculate the matrix for the parity-violating u-channel (P_c)

      @param i row of the matrix, 0->5
      @param j column of the matrix, 0->5
      @return output of the matrix

  */
  double PentaquarkMatrix::GetPVuchannel(int i, int j)
  {
    double ret;
         if (i == 0 && j == 2) ret = -((mp + ru)/(2*Eb_u*mb*ru + 2*mb2*ru));
    else if (i == 0 && j == 5) ret = ((Ep_u + mp)*(-(mp*ru) + u))/(2.*mb2*mp*(u - um));
    else if (i == 1 && j == 0) ret =    -((mp*mpsi2)/(4*Eb_u*qu2*ru*u + 4*mb*qu2*ru*u));
    else if (i == 1 && j == 1) ret = -(mp*mpsi*(mp + ru))/(2.*sqrt(2.)*(Eb_u + mb)*qu2*ru*u);
    else if (i == 1 && j == 2) ret =     (Eb_u*Ep_u*mpsi2 + Ep_u*mb*mpsi2 + Eb_u*mp*mpsi2 + mb*mp*mpsi2 + 2*Epsi_u*Ep_u*num_u + 2*Epsi_u*mp*num_u + mpsi2*num_u + 2*Eb_u*mp*qu2 + 2*num_u*qu2 + 2*Eb_u*qu2*ru)/(4*Eb_u*mb*qu2*u + 4*mb2*qu2*u);
    else if (i == 1 && j == 3) ret =     (Ep_u*mpsi2 + mp*mpsi2)/(4*mb*qu2*ru*u - 4*mb*qu2*ru*um);
    else if (i == 1 && j == 4) ret = ((Ep_u + mp)*mpsi*(-mp + ru))/(2.*sqrt(2.)*mb*qu2*ru*(u - um));
    else if (i == 1 && j == 5) ret =     (2*Epsi_u*Ep_u*mp*num_u + 2*Epsi_u*mp2*num_u + Ep_u*mp2*num_u + mp*mp2*num_u - 2*Eb_u*Epsi_u*mp*qu2 + Eb_u*mp2*qu2 + mb*mpsi2*qu2 - Ep_u*num_u*u - mp*num_u*u - Eb_u*qu2*u)/(4*mb2*mp*qu2*u - 4*mb2*mp*qu2*um);
    else if (i == 2 && j == 2) ret = -(1./(2*Eb_u*mb*ru + 2*mb2*ru));
    else if (i == 2 && j == 5) ret = -((Ep_u*ru + mp*ru)/(2*mb2*mp*u - 2*mb2*mp*um));
    else if (i == 3 && j == 0) ret =    -((mp2 + mp*ru)/(4*Eb_u*qu2*ru*u + 4*mb*qu2*ru*u));
    else if (i == 3 && j == 1) ret = -((mp*mpsi)/(sqrt(2.)*(2*Eb_u*qu2*ru*u + 2*mb*qu2*ru*u)));
    else if (i == 3 && j == 2) ret =     -(Epsi_u*mb*mp - mb*mp2 + Eb_u*(Epsi_u - Ep_u + mp)*(Ep_u + mp) - Epsi_u*num_u + Ep_u*num_u + mp*num_u + Ep_u*mb*ru)/(4.*mb*(Eb_u + mb)*qu2*u);
    else if (i == 3 && j == 3) ret = ((Ep_u + mp)*(mp - ru))/(4.*mb*qu2*ru*(u - um));
    else if (i == 3 && j == 4) ret = -((Ep_u + mp)*mpsi)/(2.*sqrt(2.)*mb*qu2*ru*(u - um));
    else if (i == 3 && j == 5) ret =     (-(Epsi_u*(Ep_u*mb*mp + mb*mp2 - Ep_u*num_u - mp*num_u + Eb_u*qu2)) + (Ep_u + mp)*(-(mb*mp2) - Ep_u*num_u + mp*num_u + Eb_u*qu2 + Ep_u*mb*ru))/(4.*mb2*mp*qu2*(u - um));
    else if (i == 4 && j == 0) ret =    (mp*(Ep_u + mp)*(-mp + ru))/(4.*(Eb_u + mb)*qu2*u);
    else if (i == 4 && j == 2) ret = (-(mp2*num_u) + (Eb_u + mb)*qu2*ru + Ep_u*num_u*(-mp + ru) + mp*(Eb_u*qu2 + mb*qu2 + num_u*ru))/(4.*mb*(Eb_u + mb)*qu2*ru);
    else if (i == 4 && j == 3) ret = (mp + ru)/(4*mb*u - 4*mb*um);
    else if (i == 4 && j == 5) ret =     (ru*((Eb_u - mb)*(Epsi_u + Ep_u - mp)*(Ep_u + mp) + num_u*(mp + ru)))/(4.*mb2*mp*(u - um));
    else if (i == 5 && j == 0) ret = -(mp*(Ep_u + mp))/(4.*(Eb_u + mb)*qu2*u);
    else if (i == 5 && j == 2) ret = -((Ep_u*num_u + mp*num_u - Eb_u*qu2 - mb*qu2)/(4*Eb_u*mb*qu2*ru + 4*mb2*qu2*ru));
    else if (i == 5 && j == 3) ret = 1./(4*mb*u - 4*mb*um);
    else if (i == 5 && j == 5) ret =     (((-Eb_u + mb)*(Ep_u + mp) + num_u)*ru)/(4.*mb2*mp*(u - um));



    else return 0.;

    return ret*sqrt((u - um)/(Eb_u + mb)/(Ep_u + mp));

  }

  /**
      Prints the whole matrices for the given kinematics

  */
  void PentaquarkMatrix::Print()
  {
    std::cout << "Check kinematics (s channel)" << std::endl;
    std::cout << "p: " << ps << ", q: " << qs << ", Epsi_s: " << - Epsi_s << ", Ep_s: " << Ep_s << ", Eb_s: " << Eb_s << "cos theta_s: " << num_s/ps/qs << std::endl;

    std::cout << "s-channel PC" << std::endl;
    std::cout << "{ ";
    for (int i=0; i<6; i++)
    {
      std::cout << "{ ";
      for (int j=0; j<6; j++) std::cout << GetPCschannel(i,j) << (j == 5 ? "" : ", ");
      std::cout << (i == 5 ? " } }" : " }," )<< std::endl;
    }
    std::cout << std::endl << std::endl;
    std::cout << "s-channel PV" << std::endl;
    std::cout << "{ ";
    for (int i=0; i<6; i++)
    {
     std::cout << "{ ";
     for (int j=0; j<6; j++) std::cout << GetPVschannel(i,j) << (j == 5 ? "" : ", ");
     std::cout << (i == 5 ? " } }" : " }," )<< std::endl;
    }
    std::cout << std::endl << std::endl;

    std::cout << "Check kinematics (u channel)" << std::endl;
    std::cout << "p: " << pu << ", q: " << qu << ", Epsi_u: " << Epsi_u << ", Ep_u: " << Ep_u << ", Eb_u: " << Eb_u << "cos theta_u: " << num_u/pu/qu << std::endl;
    std::cout << "u-channel PC" << std::endl;
    std::cout << "{ ";
    for (int i=0; i<6; i++)
    {
     std::cout << "{ ";
     for (int j=0; j<6; j++) std::cout << GetPCuchannel(i,j) << (j == 5 ? "" : ", ");
     std::cout << (i == 5 ? " } }" : " }," )<< std::endl;
    }
    std::cout << std::endl << std::endl;
    std::cout << "u-channel PV" << std::endl;
    std::cout << "{ ";
    for (int i=0; i<6; i++)
    {
     std::cout << "{ ";
     for (int j=0; j<6; j++) std::cout << GetPVuchannel(i,j) << (j == 5 ? "" : ", ");
     std::cout << (i == 5 ? " } }" : " }," )<< std::endl;
    }

  }
