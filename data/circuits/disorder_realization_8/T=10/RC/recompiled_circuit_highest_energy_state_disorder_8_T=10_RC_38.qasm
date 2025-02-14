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
rz(2.3249792) q[0];
sx q[0];
rz(-2.6007574) q[0];
sx q[0];
rz(1.9833366) q[0];
rz(0.090016063) q[1];
sx q[1];
rz(-2.6675192) q[1];
sx q[1];
rz(2.3064244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0220563) q[0];
sx q[0];
rz(-1.8827264) q[0];
sx q[0];
rz(0.74260143) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41423256) q[2];
sx q[2];
rz(-2.4576839) q[2];
sx q[2];
rz(1.4642844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1013284) q[1];
sx q[1];
rz(-2.6928386) q[1];
sx q[1];
rz(-0.59228102) q[1];
x q[2];
rz(-2.4477169) q[3];
sx q[3];
rz(-0.78842794) q[3];
sx q[3];
rz(0.83886787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0498672) q[2];
sx q[2];
rz(-2.2455402) q[2];
sx q[2];
rz(0.33207616) q[2];
rz(-0.24886985) q[3];
sx q[3];
rz(-1.1955465) q[3];
sx q[3];
rz(2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2700972) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(-0.53263295) q[0];
rz(3.0337785) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(2.8210988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1030514) q[0];
sx q[0];
rz(-1.6344995) q[0];
sx q[0];
rz(-2.2233783) q[0];
rz(1.8024496) q[2];
sx q[2];
rz(-1.4024874) q[2];
sx q[2];
rz(1.9274118) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7254843) q[1];
sx q[1];
rz(-1.7303161) q[1];
sx q[1];
rz(1.7249589) q[1];
rz(0.32933195) q[3];
sx q[3];
rz(-1.7747702) q[3];
sx q[3];
rz(-0.70131174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8104441) q[2];
sx q[2];
rz(-0.30511567) q[2];
sx q[2];
rz(-1.5020874) q[2];
rz(-2.8034927) q[3];
sx q[3];
rz(-2.2564042) q[3];
sx q[3];
rz(0.52687183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51760393) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(-2.7841618) q[0];
rz(0.39237157) q[1];
sx q[1];
rz(-0.79634276) q[1];
sx q[1];
rz(1.2145112) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8100963) q[0];
sx q[0];
rz(-1.4292681) q[0];
sx q[0];
rz(3.1241199) q[0];
rz(-2.8896595) q[2];
sx q[2];
rz(-2.5802543) q[2];
sx q[2];
rz(0.47719819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.743108) q[1];
sx q[1];
rz(-1.9581984) q[1];
sx q[1];
rz(-0.0031664567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5143993) q[3];
sx q[3];
rz(-0.58844756) q[3];
sx q[3];
rz(1.1882888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3452722) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(-2.7447682) q[2];
rz(-1.8848298) q[3];
sx q[3];
rz(-1.4182914) q[3];
sx q[3];
rz(0.61029148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20763718) q[0];
sx q[0];
rz(-1.7455245) q[0];
sx q[0];
rz(1.9130094) q[0];
rz(1.2127016) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(0.62087762) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1771496) q[0];
sx q[0];
rz(-1.7375792) q[0];
sx q[0];
rz(-2.0509999) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9758756) q[2];
sx q[2];
rz(-1.4172557) q[2];
sx q[2];
rz(0.64879791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5073237) q[1];
sx q[1];
rz(-0.36959999) q[1];
sx q[1];
rz(-2.3564767) q[1];
x q[2];
rz(0.38264783) q[3];
sx q[3];
rz(-1.6810025) q[3];
sx q[3];
rz(-2.4085542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2584201) q[2];
sx q[2];
rz(-1.1383388) q[2];
sx q[2];
rz(2.9849198) q[2];
rz(-0.91642085) q[3];
sx q[3];
rz(-1.0333002) q[3];
sx q[3];
rz(-2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85585344) q[0];
sx q[0];
rz(-1.9802977) q[0];
sx q[0];
rz(0.39749417) q[0];
rz(-0.022857895) q[1];
sx q[1];
rz(-0.50067478) q[1];
sx q[1];
rz(-2.4668677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200505) q[0];
sx q[0];
rz(-1.3617317) q[0];
sx q[0];
rz(-1.3305713) q[0];
rz(-pi) q[1];
rz(2.4910035) q[2];
sx q[2];
rz(-0.15531596) q[2];
sx q[2];
rz(2.3273205) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6276363) q[1];
sx q[1];
rz(-1.7412724) q[1];
sx q[1];
rz(-1.7254616) q[1];
x q[2];
rz(2.2915693) q[3];
sx q[3];
rz(-2.7394419) q[3];
sx q[3];
rz(-1.4141724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61458331) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(-0.69925365) q[2];
rz(-1.7953385) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(2.4301372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72494495) q[0];
sx q[0];
rz(-2.7236433) q[0];
sx q[0];
rz(2.7432192) q[0];
rz(-2.1669972) q[1];
sx q[1];
rz(-1.83788) q[1];
sx q[1];
rz(1.4385361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2319039) q[0];
sx q[0];
rz(-1.4657019) q[0];
sx q[0];
rz(-0.30334453) q[0];
rz(-1.6104524) q[2];
sx q[2];
rz(-2.8855814) q[2];
sx q[2];
rz(-0.91259749) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4230835) q[1];
sx q[1];
rz(-2.4197277) q[1];
sx q[1];
rz(1.9265429) q[1];
x q[2];
rz(1.8277728) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(-2.0723267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4633816) q[2];
sx q[2];
rz(-1.3621829) q[2];
sx q[2];
rz(-1.8360651) q[2];
rz(1.3156923) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(-0.8684043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.0618133) q[0];
sx q[0];
rz(-0.4158622) q[0];
sx q[0];
rz(1.6449991) q[0];
rz(2.1553701) q[1];
sx q[1];
rz(-1.58135) q[1];
sx q[1];
rz(2.644002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3357996) q[0];
sx q[0];
rz(-1.8247274) q[0];
sx q[0];
rz(3.0428314) q[0];
rz(0.87734434) q[2];
sx q[2];
rz(-2.3893271) q[2];
sx q[2];
rz(-1.8107506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0329602) q[1];
sx q[1];
rz(-1.7268306) q[1];
sx q[1];
rz(-1.969127) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4769745) q[3];
sx q[3];
rz(-1.854343) q[3];
sx q[3];
rz(2.97992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30279532) q[2];
sx q[2];
rz(-1.5668198) q[2];
sx q[2];
rz(2.8509169) q[2];
rz(0.21026462) q[3];
sx q[3];
rz(-2.1472774) q[3];
sx q[3];
rz(-2.1073585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57598376) q[0];
sx q[0];
rz(-1.9714332) q[0];
sx q[0];
rz(1.1822816) q[0];
rz(-2.9871509) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(-2.2363037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31667865) q[0];
sx q[0];
rz(-1.8435045) q[0];
sx q[0];
rz(0.012556745) q[0];
rz(-0.45640517) q[2];
sx q[2];
rz(-1.785136) q[2];
sx q[2];
rz(2.6021007) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0539848) q[1];
sx q[1];
rz(-1.7167257) q[1];
sx q[1];
rz(-0.75775679) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9117113) q[3];
sx q[3];
rz(-1.3671759) q[3];
sx q[3];
rz(-0.019817185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0514544) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(-0.95019379) q[2];
rz(3.0692302) q[3];
sx q[3];
rz(-1.7642998) q[3];
sx q[3];
rz(2.4322521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224834) q[0];
sx q[0];
rz(-2.1868732) q[0];
sx q[0];
rz(-2.0461653) q[0];
rz(1.0911881) q[1];
sx q[1];
rz(-2.2406816) q[1];
sx q[1];
rz(0.99517623) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93928775) q[0];
sx q[0];
rz(-1.5065333) q[0];
sx q[0];
rz(0.66575428) q[0];
rz(-pi) q[1];
rz(2.8666696) q[2];
sx q[2];
rz(-1.8951534) q[2];
sx q[2];
rz(-0.37516007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1171451) q[1];
sx q[1];
rz(-0.096344171) q[1];
sx q[1];
rz(0.49343719) q[1];
x q[2];
rz(-2.9561958) q[3];
sx q[3];
rz(-1.3747921) q[3];
sx q[3];
rz(-0.40962266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6301849) q[2];
sx q[2];
rz(-2.6738622) q[2];
sx q[2];
rz(-2.4616145) q[2];
rz(1.6857111) q[3];
sx q[3];
rz(-0.89434353) q[3];
sx q[3];
rz(1.9623914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7162914) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(-2.5262078) q[0];
rz(1.3353434) q[1];
sx q[1];
rz(-1.5348624) q[1];
sx q[1];
rz(-1.3478442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5292717) q[0];
sx q[0];
rz(-0.35824305) q[0];
sx q[0];
rz(-0.1310346) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.32424) q[2];
sx q[2];
rz(-2.0869227) q[2];
sx q[2];
rz(-1.1990666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19619689) q[1];
sx q[1];
rz(-0.41204231) q[1];
sx q[1];
rz(-0.41007385) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6928635) q[3];
sx q[3];
rz(-0.81840912) q[3];
sx q[3];
rz(0.19644745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8936257) q[2];
sx q[2];
rz(-1.3004356) q[2];
sx q[2];
rz(0.51353961) q[2];
rz(-2.9945471) q[3];
sx q[3];
rz(-2.5976318) q[3];
sx q[3];
rz(1.8769544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2697987) q[0];
sx q[0];
rz(-2.2525621) q[0];
sx q[0];
rz(2.9004108) q[0];
rz(-1.3006032) q[1];
sx q[1];
rz(-1.069297) q[1];
sx q[1];
rz(-0.020513608) q[1];
rz(2.9993771) q[2];
sx q[2];
rz(-1.1995865) q[2];
sx q[2];
rz(-2.8778278) q[2];
rz(2.0165689) q[3];
sx q[3];
rz(-2.1504938) q[3];
sx q[3];
rz(0.20636054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
