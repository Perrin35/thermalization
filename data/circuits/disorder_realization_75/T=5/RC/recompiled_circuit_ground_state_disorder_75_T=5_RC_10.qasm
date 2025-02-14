OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(6.5919696) q[0];
sx q[0];
rz(6.0433521) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(1.5129369) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37555239) q[0];
sx q[0];
rz(-2.7514725) q[0];
sx q[0];
rz(1.2287771) q[0];
rz(-1.0983019) q[2];
sx q[2];
rz(-1.3368946) q[2];
sx q[2];
rz(-2.5462674) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1354243) q[1];
sx q[1];
rz(-1.0841771) q[1];
sx q[1];
rz(1.7371337) q[1];
x q[2];
rz(-1.4024847) q[3];
sx q[3];
rz(-1.5136079) q[3];
sx q[3];
rz(-2.7528499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1538887) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(1.055701) q[2];
rz(2.740247) q[3];
sx q[3];
rz(-1.515712) q[3];
sx q[3];
rz(1.6373985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317247) q[0];
sx q[0];
rz(-2.8724176) q[0];
sx q[0];
rz(1.4058231) q[0];
rz(2.3954605) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(3.0768118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9262652) q[0];
sx q[0];
rz(-0.61867061) q[0];
sx q[0];
rz(-0.72665213) q[0];
rz(-0.6798052) q[2];
sx q[2];
rz(-1.4379939) q[2];
sx q[2];
rz(-0.50253579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3105615) q[1];
sx q[1];
rz(-1.4906724) q[1];
sx q[1];
rz(-1.7145654) q[1];
rz(-pi) q[2];
rz(-1.719789) q[3];
sx q[3];
rz(-2.9284366) q[3];
sx q[3];
rz(0.085863559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88227162) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-3.1003013) q[2];
rz(-1.1193554) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(-2.0004499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(0.69931716) q[0];
rz(0.62272561) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(0.25845382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39614933) q[0];
sx q[0];
rz(-0.93702261) q[0];
sx q[0];
rz(1.5859912) q[0];
rz(-pi) q[1];
rz(3.0270789) q[2];
sx q[2];
rz(-2.2199759) q[2];
sx q[2];
rz(2.442292) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.041641673) q[1];
sx q[1];
rz(-1.1219684) q[1];
sx q[1];
rz(1.5484518) q[1];
rz(0.036470099) q[3];
sx q[3];
rz(-1.1410645) q[3];
sx q[3];
rz(3.1231073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93566018) q[2];
sx q[2];
rz(-1.460133) q[2];
sx q[2];
rz(2.4350731) q[2];
rz(-0.62890729) q[3];
sx q[3];
rz(-0.8258515) q[3];
sx q[3];
rz(-2.0188873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(2.1674147) q[0];
rz(1.6575419) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(-1.0008224) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324686) q[0];
sx q[0];
rz(-1.4147458) q[0];
sx q[0];
rz(-2.1740995) q[0];
x q[1];
rz(-2.328726) q[2];
sx q[2];
rz(-0.95693254) q[2];
sx q[2];
rz(2.6809106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.448424) q[1];
sx q[1];
rz(-2.1285372) q[1];
sx q[1];
rz(-1.0156469) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9934512) q[3];
sx q[3];
rz(-2.8175948) q[3];
sx q[3];
rz(1.2202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68690825) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(-0.81721133) q[2];
rz(-0.59598437) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(1.7433085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3572094) q[0];
sx q[0];
rz(-1.4056118) q[0];
sx q[0];
rz(1.0719517) q[0];
rz(-1.1373854) q[1];
sx q[1];
rz(-1.5147361) q[1];
sx q[1];
rz(-0.17328182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.539285) q[0];
sx q[0];
rz(-2.0208911) q[0];
sx q[0];
rz(1.1316687) q[0];
rz(1.0801804) q[2];
sx q[2];
rz(-1.2254493) q[2];
sx q[2];
rz(2.5290348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76643244) q[1];
sx q[1];
rz(-2.3262437) q[1];
sx q[1];
rz(-0.098618193) q[1];
x q[2];
rz(-0.51637465) q[3];
sx q[3];
rz(-0.18680113) q[3];
sx q[3];
rz(-1.9641701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6413573) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(-0.87453169) q[2];
rz(-2.4747961) q[3];
sx q[3];
rz(-1.6801445) q[3];
sx q[3];
rz(-0.82505208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85103971) q[0];
sx q[0];
rz(-0.14851004) q[0];
sx q[0];
rz(2.9669951) q[0];
rz(2.5945276) q[1];
sx q[1];
rz(-0.51536307) q[1];
sx q[1];
rz(-1.863716) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1924612) q[0];
sx q[0];
rz(-1.8572079) q[0];
sx q[0];
rz(0.064563036) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75210394) q[2];
sx q[2];
rz(-1.9255203) q[2];
sx q[2];
rz(-2.6238476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6656832) q[1];
sx q[1];
rz(-1.2194677) q[1];
sx q[1];
rz(0.76900469) q[1];
rz(-pi) q[2];
rz(2.2383591) q[3];
sx q[3];
rz(-2.0732911) q[3];
sx q[3];
rz(2.5046405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7586907) q[2];
sx q[2];
rz(-0.26472696) q[2];
sx q[2];
rz(0.36941377) q[2];
rz(-2.3139125) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5740042) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(-2.9045203) q[0];
rz(1.392662) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(1.0038092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061351731) q[0];
sx q[0];
rz(-2.6727242) q[0];
sx q[0];
rz(-2.7659225) q[0];
rz(-pi) q[1];
rz(-2.4280274) q[2];
sx q[2];
rz(-1.9305561) q[2];
sx q[2];
rz(0.52116115) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3442093) q[1];
sx q[1];
rz(-2.4996669) q[1];
sx q[1];
rz(2.5688237) q[1];
rz(-2.4455261) q[3];
sx q[3];
rz(-2.4579774) q[3];
sx q[3];
rz(-3.0974577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5281333) q[2];
sx q[2];
rz(-1.5789092) q[2];
sx q[2];
rz(2.6252739) q[2];
rz(0.83827072) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(-0.18394884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8539921) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(0.91424346) q[0];
rz(-0.5842579) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(0.90726888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0153025) q[0];
sx q[0];
rz(-1.5889822) q[0];
sx q[0];
rz(-0.92211266) q[0];
rz(-1.3970988) q[2];
sx q[2];
rz(-2.1790884) q[2];
sx q[2];
rz(2.4243958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0338194) q[1];
sx q[1];
rz(-2.2108452) q[1];
sx q[1];
rz(2.1442851) q[1];
rz(-0.26772883) q[3];
sx q[3];
rz(-0.98891034) q[3];
sx q[3];
rz(-1.8393593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82627901) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(-2.2798174) q[2];
rz(1.9050542) q[3];
sx q[3];
rz(-0.084241353) q[3];
sx q[3];
rz(-1.9305362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(-2.1114517) q[0];
rz(-2.8252699) q[1];
sx q[1];
rz(-1.373469) q[1];
sx q[1];
rz(-2.8533459) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0228393) q[0];
sx q[0];
rz(-0.42834474) q[0];
sx q[0];
rz(-2.785925) q[0];
rz(-pi) q[1];
rz(-0.64916237) q[2];
sx q[2];
rz(-1.1466951) q[2];
sx q[2];
rz(-2.8555388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.32375755) q[1];
sx q[1];
rz(-0.97684469) q[1];
sx q[1];
rz(-0.51333921) q[1];
rz(-pi) q[2];
x q[2];
rz(1.796631) q[3];
sx q[3];
rz(-2.2299181) q[3];
sx q[3];
rz(-0.82443217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85252243) q[2];
sx q[2];
rz(-1.728629) q[2];
sx q[2];
rz(-0.22656013) q[2];
rz(1.4551) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15014547) q[0];
sx q[0];
rz(-1.8658072) q[0];
sx q[0];
rz(0.62514296) q[0];
rz(1.217968) q[1];
sx q[1];
rz(-2.4324799) q[1];
sx q[1];
rz(1.2120754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0934132) q[0];
sx q[0];
rz(-2.163889) q[0];
sx q[0];
rz(-2.1593447) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3122968) q[2];
sx q[2];
rz(-2.2257559) q[2];
sx q[2];
rz(-2.0531246) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21704504) q[1];
sx q[1];
rz(-0.44483652) q[1];
sx q[1];
rz(-1.7052827) q[1];
rz(-pi) q[2];
rz(1.2792222) q[3];
sx q[3];
rz(-1.8599038) q[3];
sx q[3];
rz(-1.5856727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43442279) q[2];
sx q[2];
rz(-1.3753128) q[2];
sx q[2];
rz(-2.0909069) q[2];
rz(2.5896942) q[3];
sx q[3];
rz(-0.69010186) q[3];
sx q[3];
rz(1.8570073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5156749) q[0];
sx q[0];
rz(-2.3139625) q[0];
sx q[0];
rz(2.3232842) q[0];
rz(-2.3552409) q[1];
sx q[1];
rz(-2.4813589) q[1];
sx q[1];
rz(-2.9207041) q[1];
rz(3.1037504) q[2];
sx q[2];
rz(-1.2076245) q[2];
sx q[2];
rz(-2.5943499) q[2];
rz(3.0339387) q[3];
sx q[3];
rz(-2.4526377) q[3];
sx q[3];
rz(-2.4124877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
