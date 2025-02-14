OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.028539874) q[0];
sx q[0];
rz(-2.324335) q[0];
sx q[0];
rz(-0.28491268) q[0];
rz(-2.3218396) q[1];
sx q[1];
rz(-2.2567891) q[1];
sx q[1];
rz(1.1362145) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2819913) q[0];
sx q[0];
rz(-1.4930269) q[0];
sx q[0];
rz(0.95232173) q[0];
rz(-2.4769463) q[2];
sx q[2];
rz(-2.342939) q[2];
sx q[2];
rz(-2.2448774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6416671) q[1];
sx q[1];
rz(-1.2878152) q[1];
sx q[1];
rz(0.44264117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3435279) q[3];
sx q[3];
rz(-1.1385001) q[3];
sx q[3];
rz(0.55124084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39584407) q[2];
sx q[2];
rz(-0.48600799) q[2];
sx q[2];
rz(0.29224545) q[2];
rz(0.62119836) q[3];
sx q[3];
rz(-2.0665593) q[3];
sx q[3];
rz(-2.9228389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.578823) q[0];
sx q[0];
rz(-1.1172453) q[0];
sx q[0];
rz(-0.086409464) q[0];
rz(-2.7482765) q[1];
sx q[1];
rz(-1.6875279) q[1];
sx q[1];
rz(-1.4837861) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4230712) q[0];
sx q[0];
rz(-0.72038652) q[0];
sx q[0];
rz(-1.9299401) q[0];
x q[1];
rz(-1.9772524) q[2];
sx q[2];
rz(-0.33824846) q[2];
sx q[2];
rz(1.8601314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68951974) q[1];
sx q[1];
rz(-1.4175002) q[1];
sx q[1];
rz(-0.91233715) q[1];
rz(-pi) q[2];
rz(1.8350164) q[3];
sx q[3];
rz(-2.1253028) q[3];
sx q[3];
rz(-2.4305024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94574404) q[2];
sx q[2];
rz(-1.5429147) q[2];
sx q[2];
rz(3.1405385) q[2];
rz(2.4347608) q[3];
sx q[3];
rz(-2.2471434) q[3];
sx q[3];
rz(-2.8640174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9029249) q[0];
sx q[0];
rz(-2.213573) q[0];
sx q[0];
rz(0.4784041) q[0];
rz(2.2967285) q[1];
sx q[1];
rz(-1.076315) q[1];
sx q[1];
rz(2.1873651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0811998) q[0];
sx q[0];
rz(-2.0234479) q[0];
sx q[0];
rz(-1.1552108) q[0];
rz(-pi) q[1];
rz(-2.235844) q[2];
sx q[2];
rz(-0.85341382) q[2];
sx q[2];
rz(2.5274072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0597998) q[1];
sx q[1];
rz(-2.2471618) q[1];
sx q[1];
rz(-2.6790819) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8370127) q[3];
sx q[3];
rz(-0.93627083) q[3];
sx q[3];
rz(-0.75283829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93762952) q[2];
sx q[2];
rz(-2.8014247) q[2];
sx q[2];
rz(-0.071268737) q[2];
rz(1.010262) q[3];
sx q[3];
rz(-1.9805757) q[3];
sx q[3];
rz(-2.8205813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45692745) q[0];
sx q[0];
rz(-1.848897) q[0];
sx q[0];
rz(-2.7780823) q[0];
rz(0.89206308) q[1];
sx q[1];
rz(-1.0328247) q[1];
sx q[1];
rz(1.2537289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9867986) q[0];
sx q[0];
rz(-2.5240233) q[0];
sx q[0];
rz(2.2876431) q[0];
rz(3.0405696) q[2];
sx q[2];
rz(-1.2083078) q[2];
sx q[2];
rz(2.6355138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4895136) q[1];
sx q[1];
rz(-1.8207912) q[1];
sx q[1];
rz(1.5279624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8785973) q[3];
sx q[3];
rz(-1.9501163) q[3];
sx q[3];
rz(2.5130481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38972798) q[2];
sx q[2];
rz(-2.3675297) q[2];
sx q[2];
rz(-1.5285726) q[2];
rz(-2.6914237) q[3];
sx q[3];
rz(-1.7894141) q[3];
sx q[3];
rz(1.2990797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12440974) q[0];
sx q[0];
rz(-0.15394177) q[0];
sx q[0];
rz(0.85087878) q[0];
rz(-1.7393913) q[1];
sx q[1];
rz(-1.5910999) q[1];
sx q[1];
rz(-0.15241399) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440118) q[0];
sx q[0];
rz(-0.23939238) q[0];
sx q[0];
rz(1.5967861) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28976243) q[2];
sx q[2];
rz(-2.1524977) q[2];
sx q[2];
rz(-1.5982472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9375657) q[1];
sx q[1];
rz(-1.4051249) q[1];
sx q[1];
rz(1.8574287) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7431454) q[3];
sx q[3];
rz(-2.6227847) q[3];
sx q[3];
rz(-2.0386396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57205814) q[2];
sx q[2];
rz(-2.2172838) q[2];
sx q[2];
rz(2.7657236) q[2];
rz(2.87319) q[3];
sx q[3];
rz(-1.0207876) q[3];
sx q[3];
rz(1.0387897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.038045) q[0];
sx q[0];
rz(-2.2232942) q[0];
sx q[0];
rz(0.12575664) q[0];
rz(3.0962931) q[1];
sx q[1];
rz(-0.89729512) q[1];
sx q[1];
rz(-0.50517857) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99085) q[0];
sx q[0];
rz(-1.3162386) q[0];
sx q[0];
rz(-1.6360511) q[0];
rz(-pi) q[1];
rz(-2.061902) q[2];
sx q[2];
rz(-1.668108) q[2];
sx q[2];
rz(2.8081913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0263848) q[1];
sx q[1];
rz(-1.7167101) q[1];
sx q[1];
rz(2.3664631) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0374378) q[3];
sx q[3];
rz(-1.2705876) q[3];
sx q[3];
rz(-1.2128346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.81856412) q[2];
sx q[2];
rz(-2.8105141) q[2];
sx q[2];
rz(2.8884812) q[2];
rz(-0.74462914) q[3];
sx q[3];
rz(-1.5537477) q[3];
sx q[3];
rz(2.871992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.700106) q[0];
sx q[0];
rz(-1.8550669) q[0];
sx q[0];
rz(3.098068) q[0];
rz(-1.9684017) q[1];
sx q[1];
rz(-1.2930608) q[1];
sx q[1];
rz(-0.80165577) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4960413) q[0];
sx q[0];
rz(-1.6694549) q[0];
sx q[0];
rz(-2.990574) q[0];
rz(-pi) q[1];
rz(-0.86362401) q[2];
sx q[2];
rz(-0.94085521) q[2];
sx q[2];
rz(0.0090321909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55085631) q[1];
sx q[1];
rz(-1.4155861) q[1];
sx q[1];
rz(-0.89977818) q[1];
rz(-pi) q[2];
rz(-2.4658937) q[3];
sx q[3];
rz(-2.7609918) q[3];
sx q[3];
rz(2.8322647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7511071) q[2];
sx q[2];
rz(-2.0483053) q[2];
sx q[2];
rz(2.2170179) q[2];
rz(-3.0604002) q[3];
sx q[3];
rz(-1.743318) q[3];
sx q[3];
rz(-2.413522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8586332) q[0];
sx q[0];
rz(-1.8959683) q[0];
sx q[0];
rz(2.0661195) q[0];
rz(1.3772427) q[1];
sx q[1];
rz(-1.3721162) q[1];
sx q[1];
rz(-0.65518728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0690883) q[0];
sx q[0];
rz(-3.0720815) q[0];
sx q[0];
rz(-0.68060912) q[0];
rz(-pi) q[1];
rz(-2.8497208) q[2];
sx q[2];
rz(-2.3322189) q[2];
sx q[2];
rz(2.3152817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6020169) q[1];
sx q[1];
rz(-1.5028483) q[1];
sx q[1];
rz(-0.69700239) q[1];
rz(-0.029441802) q[3];
sx q[3];
rz(-1.9092776) q[3];
sx q[3];
rz(-1.728831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7738771) q[2];
sx q[2];
rz(-0.39153063) q[2];
sx q[2];
rz(-1.1565304) q[2];
rz(1.997442) q[3];
sx q[3];
rz(-0.63026989) q[3];
sx q[3];
rz(1.3327848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82041204) q[0];
sx q[0];
rz(-2.2207566) q[0];
sx q[0];
rz(1.9658827) q[0];
rz(1.8020449) q[1];
sx q[1];
rz(-0.99937719) q[1];
sx q[1];
rz(0.40750009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20243199) q[0];
sx q[0];
rz(-2.6337886) q[0];
sx q[0];
rz(-0.61384551) q[0];
rz(-pi) q[1];
rz(-0.48353075) q[2];
sx q[2];
rz(-1.3322416) q[2];
sx q[2];
rz(1.8339744) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.54267) q[1];
sx q[1];
rz(-2.5537468) q[1];
sx q[1];
rz(-1.6076615) q[1];
rz(-0.27176933) q[3];
sx q[3];
rz(-1.6394233) q[3];
sx q[3];
rz(-2.8420699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0743865) q[2];
sx q[2];
rz(-1.9830474) q[2];
sx q[2];
rz(2.1161533) q[2];
rz(-0.49753672) q[3];
sx q[3];
rz(-2.189744) q[3];
sx q[3];
rz(0.40646762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24791726) q[0];
sx q[0];
rz(-1.1578639) q[0];
sx q[0];
rz(1.7561308) q[0];
rz(-1.0476184) q[1];
sx q[1];
rz(-1.7966929) q[1];
sx q[1];
rz(0.97283831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097869748) q[0];
sx q[0];
rz(-2.1132054) q[0];
sx q[0];
rz(-3.0412251) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49373105) q[2];
sx q[2];
rz(-1.0876473) q[2];
sx q[2];
rz(-1.8918623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9094338) q[1];
sx q[1];
rz(-1.2175908) q[1];
sx q[1];
rz(-1.722419) q[1];
rz(-pi) q[2];
rz(-1.1249969) q[3];
sx q[3];
rz(-0.74927038) q[3];
sx q[3];
rz(2.7188735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4086548) q[2];
sx q[2];
rz(-0.32821566) q[2];
sx q[2];
rz(-1.1019863) q[2];
rz(-1.7105626) q[3];
sx q[3];
rz(-0.82508636) q[3];
sx q[3];
rz(1.903532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86193209) q[0];
sx q[0];
rz(-1.4623549) q[0];
sx q[0];
rz(-2.097492) q[0];
rz(-0.73061371) q[1];
sx q[1];
rz(-0.63332557) q[1];
sx q[1];
rz(0.81061737) q[1];
rz(-0.43888406) q[2];
sx q[2];
rz(-2.4861797) q[2];
sx q[2];
rz(-0.7143153) q[2];
rz(0.35351462) q[3];
sx q[3];
rz(-1.8302038) q[3];
sx q[3];
rz(-1.3448546) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
