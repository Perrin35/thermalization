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
rz(-1.1582561) q[0];
rz(0.090016063) q[1];
sx q[1];
rz(-2.6675192) q[1];
sx q[1];
rz(-0.83516821) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0220563) q[0];
sx q[0];
rz(-1.2588663) q[0];
sx q[0];
rz(-0.74260143) q[0];
rz(2.5005241) q[2];
sx q[2];
rz(-1.3136697) q[2];
sx q[2];
rz(-2.7066305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1013284) q[1];
sx q[1];
rz(-0.44875408) q[1];
sx q[1];
rz(-2.5493116) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4477169) q[3];
sx q[3];
rz(-0.78842794) q[3];
sx q[3];
rz(-2.3027248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0917255) q[2];
sx q[2];
rz(-0.89605248) q[2];
sx q[2];
rz(-2.8095165) q[2];
rz(-0.24886985) q[3];
sx q[3];
rz(-1.1955465) q[3];
sx q[3];
rz(2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2700972) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(0.53263295) q[0];
rz(-3.0337785) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(-2.8210988) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6924393) q[0];
sx q[0];
rz(-2.4863613) q[0];
sx q[0];
rz(1.4661319) q[0];
rz(-2.2079289) q[2];
sx q[2];
rz(-2.8561428) q[2];
sx q[2];
rz(-2.1674726) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7832766) q[1];
sx q[1];
rz(-2.9202097) q[1];
sx q[1];
rz(-2.3795147) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56902253) q[3];
sx q[3];
rz(-0.38541615) q[3];
sx q[3];
rz(0.33447124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8104441) q[2];
sx q[2];
rz(-2.836477) q[2];
sx q[2];
rz(1.6395052) q[2];
rz(0.33809996) q[3];
sx q[3];
rz(-0.88518849) q[3];
sx q[3];
rz(2.6147208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6239887) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(2.7841618) q[0];
rz(2.7492211) q[1];
sx q[1];
rz(-0.79634276) q[1];
sx q[1];
rz(-1.2145112) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3833575) q[0];
sx q[0];
rz(-1.5880944) q[0];
sx q[0];
rz(1.4292468) q[0];
x q[1];
rz(2.8896595) q[2];
sx q[2];
rz(-2.5802543) q[2];
sx q[2];
rz(2.6643945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.17111546) q[1];
sx q[1];
rz(-1.5678645) q[1];
sx q[1];
rz(-1.1833925) q[1];
rz(-pi) q[2];
rz(1.5143993) q[3];
sx q[3];
rz(-0.58844756) q[3];
sx q[3];
rz(-1.9533039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3452722) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(-2.7447682) q[2];
rz(-1.8848298) q[3];
sx q[3];
rz(-1.7233012) q[3];
sx q[3];
rz(2.5313012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20763718) q[0];
sx q[0];
rz(-1.7455245) q[0];
sx q[0];
rz(1.2285832) q[0];
rz(-1.928891) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(0.62087762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1771496) q[0];
sx q[0];
rz(-1.4040134) q[0];
sx q[0];
rz(-2.0509999) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9758756) q[2];
sx q[2];
rz(-1.4172557) q[2];
sx q[2];
rz(0.64879791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5073237) q[1];
sx q[1];
rz(-0.36959999) q[1];
sx q[1];
rz(0.78511591) q[1];
rz(0.28811437) q[3];
sx q[3];
rz(-2.7441437) q[3];
sx q[3];
rz(2.0370874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88317251) q[2];
sx q[2];
rz(-1.1383388) q[2];
sx q[2];
rz(-2.9849198) q[2];
rz(2.2251718) q[3];
sx q[3];
rz(-2.1082924) q[3];
sx q[3];
rz(2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857392) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(0.39749417) q[0];
rz(0.022857895) q[1];
sx q[1];
rz(-0.50067478) q[1];
sx q[1];
rz(-0.674725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215421) q[0];
sx q[0];
rz(-1.7798609) q[0];
sx q[0];
rz(1.8110214) q[0];
rz(-pi) q[1];
rz(0.12395383) q[2];
sx q[2];
rz(-1.4769685) q[2];
sx q[2];
rz(-1.4012865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0832842) q[1];
sx q[1];
rz(-1.7232019) q[1];
sx q[1];
rz(0.17249523) q[1];
x q[2];
rz(0.85002331) q[3];
sx q[3];
rz(-2.7394419) q[3];
sx q[3];
rz(-1.7274203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61458331) q[2];
sx q[2];
rz(-2.535227) q[2];
sx q[2];
rz(-2.442339) q[2];
rz(1.7953385) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(0.7114555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72494495) q[0];
sx q[0];
rz(-2.7236433) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(2.1669972) q[1];
sx q[1];
rz(-1.3037126) q[1];
sx q[1];
rz(-1.7030565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2319039) q[0];
sx q[0];
rz(-1.6758907) q[0];
sx q[0];
rz(2.8382481) q[0];
rz(-3.1312156) q[2];
sx q[2];
rz(-1.8266018) q[2];
sx q[2];
rz(2.1880045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4230835) q[1];
sx q[1];
rz(-2.4197277) q[1];
sx q[1];
rz(-1.2150498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8277728) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(1.069266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4633816) q[2];
sx q[2];
rz(-1.7794098) q[2];
sx q[2];
rz(1.8360651) q[2];
rz(1.3156923) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(2.2731884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797794) q[0];
sx q[0];
rz(-2.7257305) q[0];
sx q[0];
rz(1.4965936) q[0];
rz(-0.98622259) q[1];
sx q[1];
rz(-1.58135) q[1];
sx q[1];
rz(2.644002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7114899) q[0];
sx q[0];
rz(-0.27207366) q[0];
sx q[0];
rz(-1.2077232) q[0];
x q[1];
rz(-0.94697081) q[2];
sx q[2];
rz(-1.1188036) q[2];
sx q[2];
rz(-0.78540451) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5274682) q[1];
sx q[1];
rz(-1.1775727) q[1];
sx q[1];
rz(0.16904633) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4769745) q[3];
sx q[3];
rz(-1.2872496) q[3];
sx q[3];
rz(-2.97992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8387973) q[2];
sx q[2];
rz(-1.5747728) q[2];
sx q[2];
rz(2.8509169) q[2];
rz(-0.21026462) q[3];
sx q[3];
rz(-2.1472774) q[3];
sx q[3];
rz(-1.0342342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57598376) q[0];
sx q[0];
rz(-1.1701595) q[0];
sx q[0];
rz(1.1822816) q[0];
rz(2.9871509) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(2.2363037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840928) q[0];
sx q[0];
rz(-1.5587036) q[0];
sx q[0];
rz(1.2980677) q[0];
x q[1];
rz(2.6851875) q[2];
sx q[2];
rz(-1.785136) q[2];
sx q[2];
rz(2.6021007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8108927) q[1];
sx q[1];
rz(-0.76892439) q[1];
sx q[1];
rz(-2.9309209) q[1];
rz(1.9117113) q[3];
sx q[3];
rz(-1.7744167) q[3];
sx q[3];
rz(3.1217755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0901383) q[2];
sx q[2];
rz(-1.5807512) q[2];
sx q[2];
rz(-0.95019379) q[2];
rz(-0.072362445) q[3];
sx q[3];
rz(-1.7642998) q[3];
sx q[3];
rz(2.4322521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41910928) q[0];
sx q[0];
rz(-0.95471946) q[0];
sx q[0];
rz(1.0954274) q[0];
rz(-1.0911881) q[1];
sx q[1];
rz(-0.90091101) q[1];
sx q[1];
rz(-2.1464164) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54992095) q[0];
sx q[0];
rz(-2.4732145) q[0];
sx q[0];
rz(3.0377798) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9068933) q[2];
sx q[2];
rz(-1.3105515) q[2];
sx q[2];
rz(2.0356095) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51982626) q[1];
sx q[1];
rz(-1.4859746) q[1];
sx q[1];
rz(-1.6165401) q[1];
rz(0.18539683) q[3];
sx q[3];
rz(-1.7668006) q[3];
sx q[3];
rz(-2.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5114078) q[2];
sx q[2];
rz(-2.6738622) q[2];
sx q[2];
rz(2.4616145) q[2];
rz(1.4558815) q[3];
sx q[3];
rz(-2.2472491) q[3];
sx q[3];
rz(-1.1792012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42530123) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(-2.5262078) q[0];
rz(1.8062493) q[1];
sx q[1];
rz(-1.6067303) q[1];
sx q[1];
rz(-1.3478442) q[1];
rz(-pi) q[2];
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
rz(-1.0546699) q[2];
sx q[2];
rz(-1.942526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19619689) q[1];
sx q[1];
rz(-0.41204231) q[1];
sx q[1];
rz(0.41007385) q[1];
rz(-pi) q[2];
rz(1.4487292) q[3];
sx q[3];
rz(-0.81840912) q[3];
sx q[3];
rz(-0.19644745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8936257) q[2];
sx q[2];
rz(-1.8411571) q[2];
sx q[2];
rz(-2.628053) q[2];
rz(-0.14704554) q[3];
sx q[3];
rz(-0.5439609) q[3];
sx q[3];
rz(1.8769544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87179398) q[0];
sx q[0];
rz(-2.2525621) q[0];
sx q[0];
rz(2.9004108) q[0];
rz(-1.8409894) q[1];
sx q[1];
rz(-2.0722957) q[1];
sx q[1];
rz(3.121079) q[1];
rz(1.1961436) q[2];
sx q[2];
rz(-1.7032663) q[2];
sx q[2];
rz(1.8864529) q[2];
rz(1.1250238) q[3];
sx q[3];
rz(-0.99109886) q[3];
sx q[3];
rz(-2.9352321) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
