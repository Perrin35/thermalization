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
rz(-2.3174602) q[0];
sx q[0];
rz(-1.1622518) q[0];
sx q[0];
rz(2.9293769) q[0];
rz(-2.4556887) q[1];
sx q[1];
rz(-0.75690126) q[1];
sx q[1];
rz(-0.25605717) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6997509) q[0];
sx q[0];
rz(-2.3810593) q[0];
sx q[0];
rz(-0.89450128) q[0];
x q[1];
rz(2.1317057) q[2];
sx q[2];
rz(-1.6445064) q[2];
sx q[2];
rz(0.91154237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4997811) q[1];
sx q[1];
rz(-2.8781673) q[1];
sx q[1];
rz(0.66536565) q[1];
rz(-pi) q[2];
rz(-0.30976661) q[3];
sx q[3];
rz(-2.081678) q[3];
sx q[3];
rz(2.8155308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2519007) q[2];
sx q[2];
rz(-2.8163741) q[2];
sx q[2];
rz(-1.088781) q[2];
rz(-0.47993663) q[3];
sx q[3];
rz(-1.9710541) q[3];
sx q[3];
rz(2.0282733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29555175) q[0];
sx q[0];
rz(-2.9132081) q[0];
sx q[0];
rz(-0.26902714) q[0];
rz(-0.61664063) q[1];
sx q[1];
rz(-0.65526217) q[1];
sx q[1];
rz(-2.2356967) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8347816) q[0];
sx q[0];
rz(-0.7894581) q[0];
sx q[0];
rz(-3.03699) q[0];
rz(2.0935161) q[2];
sx q[2];
rz(-2.0459896) q[2];
sx q[2];
rz(-1.1876729) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0695659) q[1];
sx q[1];
rz(-2.8588738) q[1];
sx q[1];
rz(2.7064884) q[1];
rz(2.5814186) q[3];
sx q[3];
rz(-1.2382714) q[3];
sx q[3];
rz(1.9412845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9783832) q[2];
sx q[2];
rz(-0.78028148) q[2];
sx q[2];
rz(-0.015446375) q[2];
rz(2.3695626) q[3];
sx q[3];
rz(-2.6431712) q[3];
sx q[3];
rz(-1.2482524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5323199) q[0];
sx q[0];
rz(-2.5072704) q[0];
sx q[0];
rz(-0.79737216) q[0];
rz(-2.4728921) q[1];
sx q[1];
rz(-2.1233605) q[1];
sx q[1];
rz(0.55577898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014975958) q[0];
sx q[0];
rz(-2.3117723) q[0];
sx q[0];
rz(1.2955342) q[0];
rz(-pi) q[1];
rz(3.0012461) q[2];
sx q[2];
rz(-2.2303584) q[2];
sx q[2];
rz(1.5857343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4517256) q[1];
sx q[1];
rz(-0.84459442) q[1];
sx q[1];
rz(1.5548585) q[1];
x q[2];
rz(-1.9706412) q[3];
sx q[3];
rz(-2.2388864) q[3];
sx q[3];
rz(-3.060241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7734311) q[2];
sx q[2];
rz(-1.1661466) q[2];
sx q[2];
rz(-1.9318761) q[2];
rz(3.0220253) q[3];
sx q[3];
rz(-2.5366294) q[3];
sx q[3];
rz(-0.97924489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46215737) q[0];
sx q[0];
rz(-1.698751) q[0];
sx q[0];
rz(2.6642642) q[0];
rz(0.10074549) q[1];
sx q[1];
rz(-2.0858177) q[1];
sx q[1];
rz(0.13564067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3326615) q[0];
sx q[0];
rz(-0.96600129) q[0];
sx q[0];
rz(0.51938341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2775947) q[2];
sx q[2];
rz(-1.141077) q[2];
sx q[2];
rz(-1.4336842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6203719) q[1];
sx q[1];
rz(-1.4209827) q[1];
sx q[1];
rz(0.42370875) q[1];
x q[2];
rz(-2.9075525) q[3];
sx q[3];
rz(-0.79667366) q[3];
sx q[3];
rz(0.19114574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9292844) q[2];
sx q[2];
rz(-1.6529275) q[2];
sx q[2];
rz(-2.6279214) q[2];
rz(-0.16289991) q[3];
sx q[3];
rz(-0.25560156) q[3];
sx q[3];
rz(1.1027078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61360079) q[0];
sx q[0];
rz(-3.1311212) q[0];
sx q[0];
rz(-2.7136059) q[0];
rz(-1.9635268) q[1];
sx q[1];
rz(-1.4617498) q[1];
sx q[1];
rz(-1.7286495) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21715054) q[0];
sx q[0];
rz(-0.61762327) q[0];
sx q[0];
rz(-0.87583603) q[0];
rz(-pi) q[1];
rz(-1.5287077) q[2];
sx q[2];
rz(-1.6070935) q[2];
sx q[2];
rz(-0.55755471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.481491) q[1];
sx q[1];
rz(-0.30687422) q[1];
sx q[1];
rz(2.9576388) q[1];
x q[2];
rz(2.4705162) q[3];
sx q[3];
rz(-0.6830712) q[3];
sx q[3];
rz(-0.67923123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54521769) q[2];
sx q[2];
rz(-1.7866106) q[2];
sx q[2];
rz(0.4893111) q[2];
rz(2.4423068) q[3];
sx q[3];
rz(-2.0977061) q[3];
sx q[3];
rz(1.4780686) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74772239) q[0];
sx q[0];
rz(-2.1339397) q[0];
sx q[0];
rz(-1.3442159) q[0];
rz(2.4463553) q[1];
sx q[1];
rz(-1.2691011) q[1];
sx q[1];
rz(0.42323798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7376331) q[0];
sx q[0];
rz(-0.18059093) q[0];
sx q[0];
rz(2.7904912) q[0];
rz(0.10144039) q[2];
sx q[2];
rz(-0.82221088) q[2];
sx q[2];
rz(1.6486096) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19638229) q[1];
sx q[1];
rz(-1.4221514) q[1];
sx q[1];
rz(-2.3988612) q[1];
x q[2];
rz(-0.18169348) q[3];
sx q[3];
rz(-1.5156997) q[3];
sx q[3];
rz(0.7167992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77373475) q[2];
sx q[2];
rz(-2.9025142) q[2];
sx q[2];
rz(1.9284922) q[2];
rz(-0.086094543) q[3];
sx q[3];
rz(-1.8910858) q[3];
sx q[3];
rz(1.3849753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64751476) q[0];
sx q[0];
rz(-2.8120742) q[0];
sx q[0];
rz(-0.20429097) q[0];
rz(0.34945166) q[1];
sx q[1];
rz(-1.5710257) q[1];
sx q[1];
rz(1.6250767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3649611) q[0];
sx q[0];
rz(-1.6775515) q[0];
sx q[0];
rz(-2.0364709) q[0];
x q[1];
rz(0.3150926) q[2];
sx q[2];
rz(-1.1865864) q[2];
sx q[2];
rz(-2.0496673) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.094665729) q[1];
sx q[1];
rz(-3.0311208) q[1];
sx q[1];
rz(-1.1383971) q[1];
x q[2];
rz(-1.5920226) q[3];
sx q[3];
rz(-2.5633865) q[3];
sx q[3];
rz(-2.8631722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1427631) q[2];
sx q[2];
rz(-2.5122354) q[2];
sx q[2];
rz(-2.2289185) q[2];
rz(0.015803745) q[3];
sx q[3];
rz(-1.088257) q[3];
sx q[3];
rz(-0.44939941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939746) q[0];
sx q[0];
rz(-2.0517218) q[0];
sx q[0];
rz(-2.3225978) q[0];
rz(0.35797572) q[1];
sx q[1];
rz(-2.0134108) q[1];
sx q[1];
rz(0.63041675) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6651745) q[0];
sx q[0];
rz(-3.0673084) q[0];
sx q[0];
rz(-1.3528385) q[0];
rz(-2.0840896) q[2];
sx q[2];
rz(-0.78864151) q[2];
sx q[2];
rz(3.0190133) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.158327) q[1];
sx q[1];
rz(-0.7559146) q[1];
sx q[1];
rz(-2.1793038) q[1];
rz(-2.181545) q[3];
sx q[3];
rz(-2.4195101) q[3];
sx q[3];
rz(-0.86632767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7716052) q[2];
sx q[2];
rz(-2.6321754) q[2];
sx q[2];
rz(0.39311692) q[2];
rz(0.76817948) q[3];
sx q[3];
rz(-0.79777515) q[3];
sx q[3];
rz(-1.2099077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2792252) q[0];
sx q[0];
rz(-0.51280642) q[0];
sx q[0];
rz(0.38238907) q[0];
rz(-2.6928316) q[1];
sx q[1];
rz(-0.02350137) q[1];
sx q[1];
rz(2.1429367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6210839) q[0];
sx q[0];
rz(-1.8077842) q[0];
sx q[0];
rz(-2.468796) q[0];
rz(-pi) q[1];
rz(-3.0566755) q[2];
sx q[2];
rz(-1.2341656) q[2];
sx q[2];
rz(1.6767329) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.081621165) q[1];
sx q[1];
rz(-1.27969) q[1];
sx q[1];
rz(-1.7779093) q[1];
rz(-pi) q[2];
rz(-1.6395929) q[3];
sx q[3];
rz(-2.2584887) q[3];
sx q[3];
rz(-0.53737924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0950332) q[2];
sx q[2];
rz(-0.79594374) q[2];
sx q[2];
rz(0.98540068) q[2];
rz(2.4424148) q[3];
sx q[3];
rz(-2.2907084) q[3];
sx q[3];
rz(-0.053475577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92626524) q[0];
sx q[0];
rz(-2.7113921) q[0];
sx q[0];
rz(-2.3660124) q[0];
rz(-2.8714478) q[1];
sx q[1];
rz(-1.4560478) q[1];
sx q[1];
rz(0.57033479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4264589) q[0];
sx q[0];
rz(-1.1021779) q[0];
sx q[0];
rz(-0.56632407) q[0];
rz(-pi) q[1];
rz(-3.0351972) q[2];
sx q[2];
rz(-1.0950076) q[2];
sx q[2];
rz(-1.758213) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.476984) q[1];
sx q[1];
rz(-0.23566569) q[1];
sx q[1];
rz(0.13500555) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4943284) q[3];
sx q[3];
rz(-0.88275331) q[3];
sx q[3];
rz(-0.23721671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9242085) q[2];
sx q[2];
rz(-1.9839828) q[2];
sx q[2];
rz(0.50979924) q[2];
rz(0.23586759) q[3];
sx q[3];
rz(-2.9678952) q[3];
sx q[3];
rz(0.97714669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8248642) q[0];
sx q[0];
rz(-1.7289799) q[0];
sx q[0];
rz(-1.254542) q[0];
rz(-0.53312373) q[1];
sx q[1];
rz(-1.4630547) q[1];
sx q[1];
rz(-2.0933082) q[1];
rz(-3.127373) q[2];
sx q[2];
rz(-1.1455677) q[2];
sx q[2];
rz(-1.3451411) q[2];
rz(2.6410136) q[3];
sx q[3];
rz(-1.4698613) q[3];
sx q[3];
rz(-2.1657497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
