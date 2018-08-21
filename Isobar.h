#include <math.h>
#include <iostream>
#include <complex>
#include <vector>



#ifdef DEBUG_H
#define DEBUG(x) std::cout << x << std::endl
#else
#define DEBUG(x)
#endif

#ifndef ISOBAR_H
#define ISOBAR_H

typedef std::complex<double> cd;

class Isobar
{
public:
  Isobar(double mass, double width, cd coupling, unsigned int j, unsigned int spin, unsigned int L, bool eta, bool PC, char channel, int amp);
  ~Isobar () {};
  double Mass (double x);
  double Width (double x);
  cd Coupling(cd x);

  double Mass() { return _mass; };
  double Width() { return _width; };
  cd Coupling() { return _coupling; };
  unsigned int Spin() { return _spin; };
  unsigned int L() { return _L; };
  unsigned int J() { return _j; };
  bool Eta() { return _eta; };
  bool PC() { return _PC; };
  char Channel () { return _channel; };
  cd Evaluate(double x);

  static const int BreitWigner = 0;

private:
  double _mass, _width;
  cd _coupling;
  unsigned int _spin, _j, _L;
  bool _eta, _PC;
  char _channel;
  int _amplitude;
  void CheckConsistency();
  unsigned int Spin(unsigned int x);
  unsigned int L(unsigned int x);
  unsigned int J(unsigned int x);
  int Amplitude(int x);
  bool Eta(bool x);
  bool PC(bool x);
  char Channel (char x);







};

#endif
