#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
 template<typename T> T GC0(T x, T w) {return 1.0*pow(2, 1.0L/4.0L)*pow(w, 1.0L/4.0L)*exp(-w*pow(x, 2))/pow(M_PI, 1.0L/4.0L);}
 template<typename T> T GC1(T x, T w) {return 1.4142135623731*pow(2, 3.0L/4.0L)*pow(w, 3.0L/4.0L)*x*exp(-w*pow(x, 2))/pow(M_PI, 1.0L/4.0L);}
 template<typename T> T GC2(T x, T w) {return pow(2, 3.0L/4.0L)*pow(w, 1.0L/4.0L)*(1.15470053837925*sqrt(3)*w*pow(x, 2) - 0.5)*exp(-w*pow(x, 2))/(pow(M_PI, 1.0L/4.0L)*sqrt(-0.408248290463863*sqrt(6) + 2.0));}
 template<typename T> T GC3(T x, T w) {return pow(2, 1.0L/4.0L)*pow(w, 3.0L/4.0L)*x*(1.46059348668044*sqrt(5)*w*pow(x, 2) - 1.41421356237309*sqrt(3))*exp(-w*pow(x, 2))/(pow(M_PI, 1.0L/4.0L)*sqrt(-0.547722557505166*sqrt(30) + 4.0));}
 template<typename T> T GC4(T x, T w) {return pow(2, 3.0L/4.0L)*pow(w, 1.0L/4.0L)*(0.390360029179413*sqrt(35)*pow(w, 2)*pow(x, 4) - 3.46410161513775*w*pow(x, 2) + 0.249999999999999*sqrt(3))*exp(-w*pow(x, 2))/(pow(M_PI, 1.0L/4.0L)*sqrt(-0.896421457000792*sqrt(70) - 0.612372435695792*sqrt(6) + 0.0517549169506764*sqrt(210) + 9.24999999999997));}
 template<typename T> T GC5(T x, T w) {return pow(2, 1.0L/4.0L)*pow(w, 1.0L/4.0L)*(sqrt(M_PI)*sqrt(w)*x*(1.10410489494777*sqrt(7)*pow(w, 2)*pow(x, 4) - 7.30296743340223*w*pow(x, 2) + 0.707106781186548*sqrt(15)) + 4.81725822971645e-16*sqrt(7)*pow(w, 2)*pow(x, 4) - 1.09189729375224e-15*sqrt(5)*w*pow(x, 2) + 6.81434636196684e-17*sqrt(15))*exp(-w*pow(x, 2))/sqrt(sqrt(M_PI)*(-1.74344399036082e-31*sqrt(70) - 1.31531873072206e-31*sqrt(6) + 8.70443129110198e-33*sqrt(210) + 1.85364304437592e-30) + M_PI*(-1.86571608980725e-15*sqrt(14) - 6.39349589175072e-16*sqrt(70) - 1.3648716171903e-15*sqrt(6) - 1.75945666501245e-16*sqrt(30) + 1.47031940189871e-16*sqrt(210) + 4.46005260969297e-15*sqrt(2) + 2.81926668964759e-15*sqrt(10)) + pow(M_PI, 3.0L/2.0L)*(-4.67707173346744*sqrt(14) - 1.36930639376292*sqrt(30) + 0.258774584753383*sqrt(210) + 22.2500000000001));}
 template<typename T> T GC6(T x, T w) {return pow(2, 1.0L/4.0L)*pow(w, 1.0L/4.0L)*(sqrt(w)*x*(-2.30821208100296e-16*sqrt(21)*pow(w, 2)*pow(x, 4) + 3.00703615104905e-15*sqrt(3)*w*pow(x, 2) - 1.30345259857917e-15*sqrt(5)) + sqrt(M_PI)*(0.156930636356401*sqrt(231)*pow(w, 3)*pow(x, 6) - 1.95180014589708*sqrt(21)*pow(w, 2)*pow(x, 4) + 1.73205080756889*sqrt(15)*w*pow(x, 2) - 0.250000000000001*sqrt(5)))*exp(-w*pow(x, 2))/sqrt(M_PI*(-9.3377343759136e-15*sqrt(14) - 3.99099221302972e-15*sqrt(6) - 2.017060897162e-15*sqrt(22) - 6.36068454365811e-16*sqrt(70) - 1.08479891862402e-16*sqrt(2310) - 2.65786956709986e-16*sqrt(30) + 9.1986965170096e-16*sqrt(210) + 1.11868514770392e-14*sqrt(2) + 5.52427815612967e-15*sqrt(10) + 1.50156418506225e-15*sqrt(154)) + sqrt(M_PI)*(-1.03932209817624e-30*sqrt(30) - 1.2078136831379e-30*sqrt(14) + 9.97234362408154e-32*sqrt(210) + 9.51410915001654e-30) + pow(M_PI, 3.0L/2.0L)*(-8.39477820650226*sqrt(22) - 3.36158046375302*sqrt(70) - 0.76546554461975*sqrt(6) - 0.0130039121652576*sqrt(2310) + 0.129387292376693*sqrt(210) + 0.472992167915174*sqrt(770) + 56.0000000000006));}
#pragma GCC diagnostic pop
template<typename T> T GC(T x, T w, int n) {
   if (n > 6) {
       return -1;
   }
   switch(n) {
       case 0: return GC0(x, w);
       case 1: return GC1(x, w);
       case 2: return GC2(x, w);
       case 3: return GC3(x, w);
       case 4: return GC4(x, w);
       case 5: return GC5(x, w);
       case 6: return GC6(x, w);
       default: return 0;
   }
}
