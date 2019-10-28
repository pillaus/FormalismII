{
	// s+ 2.05047817057674431e+00
	// s- 1.97664786049396607e-01
	// s- 6.36391438239999907e+00
	// s+ 7.59770235904000089e+01

        double mb = 5.61958, mpsi = 3.0969, mp = .938272081, mK = .493677;
	printf("(s*(2*t + s - %lf) + %lf) / sqrt((s - %lf)*(s - %lf)*(s - %lf)*(s - %lf))\n", mb*mb + mpsi*mpsi + mp*mp + mK*mK, (mb*mb - mpsi*mpsi)*(mp*mp - mK*mK),
		pow(mb + mpsi,2), pow(mb - mpsi,2), pow(mp + mK,2), pow(mp - mK,2));

	ntuple->Draw("(s*(-2*t - s + 42.294540) + 13.998952) / sqrt((s - 75.977024)*(s - 6.363914)*(s - 2.050478)*(s - 0.197665))","weight","COLZ");
}
