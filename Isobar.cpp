/**
    Isobar.cpp
    Calculate the matrices in Appendix G of 1805.02113

    @author Alessandro Pilloni
    @version 1.0 2018-07-22
*/
#define DEBUG_H
#include "Isobar.h"

const int Isobar::BreitWigner = 0;

  /**
      Constructor. It gets the masses as optional input. Default from LivePDG, 2018-05-13

      @param _mb Lambda_b mass
      @param _mpsi J/psi mass
      @param _mp proton mass
      @param _mK Kaon mass

  */
Isobar::Isobar(double mass, double width, cd coupling, unsigned int j, unsigned int spin, unsigned int L, bool eta, bool PC, char channel, int amp)
{
  this->Mass(mass);
  this->Width(width);
  this->Coupling(coupling);
  this->J(j);
  this->Spin(spin);
  this->L(L);
  this->Eta(eta);
  this->PC(PC);
  this->Channel(channel);
  this->Amplitude(amp);

  CheckConsistency();

  return;
}

double Isobar::Mass(double x)
{
  if (x < 0.) { DEBUG("Negative mass, set mass = -mass"); _mass = -_mass; }
  return _mass = x;
}

double Isobar::Width(double x)
{
  if (x < 0.) { DEBUG("Negative width, set width = -width"); _width = -_width; }
  return _width = x;
}
cd Isobar::Coupling(cd x)
{
  return _coupling = x;
}
unsigned int Isobar::J(unsigned x)
{
  if (x < 0 || x % 2 != 1) { DEBUG("Nonvalid J, set default 3/2"); x = 3; }
  return _j = x;
}
unsigned int Isobar::Spin(unsigned x)
{
  if (x != 1 && x!= 3) { DEBUG("Nonvalid Spin, set default 3/2"); x = 3; }
  return _spin = x;
}
unsigned int Isobar::L(unsigned x)
{
  if (x < 0) { DEBUG("Nonvalid L, set default 0"); x = 0; }
  return _L = x;
}
bool Isobar::Eta(bool x)
{
  return _eta = x;
}
bool Isobar::PC(bool x)
{
  return _PC = x;
}
char Isobar::Channel(char x)
{
  if (x != 's' && x != 'u') { DEBUG("Nonvalid channel, set default s"); x = 's'; }
  return _channel = x;
}
int Isobar::Amplitude(int x)
{
  if (x!= 0) { DEBUG("Unknown amplitude, setting BreitWigner"); x = 0; }
  return _amplitude = x;
}

cd Isobar::Evaluate(double x)
{
  switch (_amplitude)
  {
    case 0:
    cd den = _mass * _mass - x - 1.i * _mass * _width;
    return _coupling / den;
    return 1.i ; //_coupling/(_mass * _mass - x - 1.i * _mass * _width);
  }

}
void Isobar::CheckConsistency()
{
  if ( _j <  2*_L - _spin || _j >  2*_L + _spin) { DEBUG("L failing consistency, set default 0"); _L = 0; }
  bool x = (_j + 1 + 2*_L) % 4 == 0;
  if (_channel == 's' && !_PC) x = !x;
  if (x != _eta)  { DEBUG("Eta failing consistency, set Eta = not(Eta)"); _eta = !_eta; }
  return;
}
