OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7822617) q[0];
sx q[0];
rz(-1.2027272) q[0];
sx q[0];
rz(2.9666535) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0543047) q[2];
sx q[2];
rz(-0.44863551) q[2];
sx q[2];
rz(-1.0686312) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5056155) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(-0.14640267) q[1];
x q[2];
rz(-0.9790768) q[3];
sx q[3];
rz(-2.7032529) q[3];
sx q[3];
rz(2.0548267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(-2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(0.97066561) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(-0.81545365) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6204651) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(-0.64882664) q[0];
rz(1.843812) q[2];
sx q[2];
rz(-0.85999876) q[2];
sx q[2];
rz(-3.0836011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7533469) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(-1.5335598) q[1];
rz(-pi) q[2];
rz(-1.7933513) q[3];
sx q[3];
rz(-2.2046304) q[3];
sx q[3];
rz(-2.8702877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(-1.3300928) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(-2.1420746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.2721491) q[0];
sx q[0];
rz(-2.9850328) q[0];
rz(-0.91018422) q[2];
sx q[2];
rz(-2.0816457) q[2];
sx q[2];
rz(2.3597033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7175908) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(1.8600149) q[1];
rz(1.7129094) q[3];
sx q[3];
rz(-2.6283773) q[3];
sx q[3];
rz(-1.7538479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(2.8667563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4764458) q[0];
sx q[0];
rz(-0.11867141) q[0];
sx q[0];
rz(0.84400405) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5675797) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(-2.761063) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6064925) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(2.5938354) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12675385) q[3];
sx q[3];
rz(-2.2634014) q[3];
sx q[3];
rz(-0.17514378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.6246187) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8109587) q[0];
sx q[0];
rz(-1.9348382) q[0];
sx q[0];
rz(2.2383658) q[0];
rz(-pi) q[1];
rz(2.5324608) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(-2.6513211) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72488898) q[1];
sx q[1];
rz(-0.54425889) q[1];
sx q[1];
rz(-2.4328028) q[1];
rz(2.153271) q[3];
sx q[3];
rz(-2.0867996) q[3];
sx q[3];
rz(2.2367246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8020442) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(-0.91066796) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(0.34067672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7621988) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(1.5048774) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3855626) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(2.0331969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1813982) q[1];
sx q[1];
rz(-1.3409233) q[1];
sx q[1];
rz(-0.48744907) q[1];
x q[2];
rz(-0.86990279) q[3];
sx q[3];
rz(-0.60855908) q[3];
sx q[3];
rz(0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(-2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1691549) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(-3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-0.46494928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8220673) q[0];
sx q[0];
rz(-1.3789346) q[0];
sx q[0];
rz(-1.096154) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99636997) q[2];
sx q[2];
rz(-1.9949706) q[2];
sx q[2];
rz(3.0976354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.834224) q[1];
sx q[1];
rz(-2.8982179) q[1];
sx q[1];
rz(0.4537531) q[1];
rz(-pi) q[2];
rz(2.081359) q[3];
sx q[3];
rz(-0.052882346) q[3];
sx q[3];
rz(-2.3898861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(3.0632339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5057482) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23711726) q[2];
sx q[2];
rz(-1.8097005) q[2];
sx q[2];
rz(0.6616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2625418) q[1];
sx q[1];
rz(-1.3971551) q[1];
sx q[1];
rz(2.2294728) q[1];
rz(-pi) q[2];
rz(3.1073242) q[3];
sx q[3];
rz(-1.7274324) q[3];
sx q[3];
rz(0.53965118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3437929) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-2.267568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0747781) q[0];
sx q[0];
rz(-1.5794808) q[0];
sx q[0];
rz(-0.76807036) q[0];
rz(2.2888695) q[2];
sx q[2];
rz(-0.52041473) q[2];
sx q[2];
rz(2.4783217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6898432) q[1];
sx q[1];
rz(-2.6362231) q[1];
sx q[1];
rz(-1.7187353) q[1];
x q[2];
rz(1.9130575) q[3];
sx q[3];
rz(-2.1517793) q[3];
sx q[3];
rz(-0.57623219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(-1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2154685) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(1.8490851) q[0];
x q[1];
rz(2.5095021) q[2];
sx q[2];
rz(-2.9846016) q[2];
sx q[2];
rz(-1.1319515) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19503838) q[1];
sx q[1];
rz(-0.95609162) q[1];
sx q[1];
rz(2.9123995) q[1];
rz(-pi) q[2];
rz(-0.52272777) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(-0.47375351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5205004) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(-1.1307217) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(-0.48537985) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
