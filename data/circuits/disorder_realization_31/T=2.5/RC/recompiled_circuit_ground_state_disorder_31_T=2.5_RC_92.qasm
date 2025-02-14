OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.06935057) q[0];
sx q[0];
rz(-2.0877593) q[0];
sx q[0];
rz(-1.2487489) q[0];
rz(-1.2381923) q[1];
sx q[1];
rz(3.5817322) q[1];
sx q[1];
rz(11.323827) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71860367) q[0];
sx q[0];
rz(-1.1866633) q[0];
sx q[0];
rz(-2.5297013) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.473601) q[2];
sx q[2];
rz(-2.3665135) q[2];
sx q[2];
rz(0.72778406) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93672458) q[1];
sx q[1];
rz(-2.8191787) q[1];
sx q[1];
rz(-0.0017044981) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5322024) q[3];
sx q[3];
rz(-1.2838138) q[3];
sx q[3];
rz(3.0521986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.093546346) q[2];
sx q[2];
rz(-0.28154937) q[2];
sx q[2];
rz(1.2151037) q[2];
rz(1.18527) q[3];
sx q[3];
rz(-2.2351041) q[3];
sx q[3];
rz(-1.5989446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1854061) q[0];
sx q[0];
rz(-0.39919272) q[0];
sx q[0];
rz(-1.6831552) q[0];
rz(2.9996808) q[1];
sx q[1];
rz(-1.9344354) q[1];
sx q[1];
rz(-1.2123607) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1233009) q[0];
sx q[0];
rz(-0.93138501) q[0];
sx q[0];
rz(2.713504) q[0];
rz(-0.59213068) q[2];
sx q[2];
rz(-1.5730324) q[2];
sx q[2];
rz(-0.9809025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3691359) q[1];
sx q[1];
rz(-1.4765655) q[1];
sx q[1];
rz(-0.89558954) q[1];
rz(-pi) q[2];
rz(1.2207915) q[3];
sx q[3];
rz(-0.98549609) q[3];
sx q[3];
rz(3.0452951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.096574664) q[2];
sx q[2];
rz(-0.18288945) q[2];
sx q[2];
rz(2.1796687) q[2];
rz(-2.9436881) q[3];
sx q[3];
rz(-1.8353381) q[3];
sx q[3];
rz(-1.643868) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69979954) q[0];
sx q[0];
rz(-2.4213591) q[0];
sx q[0];
rz(2.0752456) q[0];
rz(-2.9329246) q[1];
sx q[1];
rz(-2.6228948) q[1];
sx q[1];
rz(-0.64812237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75877612) q[0];
sx q[0];
rz(-2.1471239) q[0];
sx q[0];
rz(-0.062950171) q[0];
x q[1];
rz(-0.49049536) q[2];
sx q[2];
rz(-1.6332262) q[2];
sx q[2];
rz(-1.708781) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0433139) q[1];
sx q[1];
rz(-2.2873291) q[1];
sx q[1];
rz(-2.067042) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6567635) q[3];
sx q[3];
rz(-0.86557612) q[3];
sx q[3];
rz(1.7253226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85731792) q[2];
sx q[2];
rz(-1.4633353) q[2];
sx q[2];
rz(0.54764444) q[2];
rz(1.0037496) q[3];
sx q[3];
rz(-2.818483) q[3];
sx q[3];
rz(-2.9799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169287) q[0];
sx q[0];
rz(-2.014761) q[0];
sx q[0];
rz(-0.3666077) q[0];
rz(-1.2855351) q[1];
sx q[1];
rz(-1.6211685) q[1];
sx q[1];
rz(2.3191648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3503338) q[0];
sx q[0];
rz(-2.9295577) q[0];
sx q[0];
rz(-2.7705482) q[0];
rz(1.8404615) q[2];
sx q[2];
rz(-1.7326151) q[2];
sx q[2];
rz(-1.2310892) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82717035) q[1];
sx q[1];
rz(-2.2232375) q[1];
sx q[1];
rz(0.62364044) q[1];
rz(2.8662258) q[3];
sx q[3];
rz(-1.5689225) q[3];
sx q[3];
rz(3.0089859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.503868) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(-1.5857504) q[2];
rz(1.2268892) q[3];
sx q[3];
rz(-0.63819686) q[3];
sx q[3];
rz(1.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1894839) q[0];
sx q[0];
rz(-1.1186849) q[0];
sx q[0];
rz(0.9504016) q[0];
rz(-2.9474126) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(2.8607184) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9683198) q[0];
sx q[0];
rz(-1.8283947) q[0];
sx q[0];
rz(0.93846847) q[0];
rz(-pi) q[1];
rz(2.9921586) q[2];
sx q[2];
rz(-1.583287) q[2];
sx q[2];
rz(-2.6246659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19843392) q[1];
sx q[1];
rz(-1.4869964) q[1];
sx q[1];
rz(1.4936624) q[1];
rz(-2.0858795) q[3];
sx q[3];
rz(-2.5722945) q[3];
sx q[3];
rz(-1.4610987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48996553) q[2];
sx q[2];
rz(-1.6941035) q[2];
sx q[2];
rz(-0.70072407) q[2];
rz(1.0939595) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(-2.3195364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1489498) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(-0.5067504) q[0];
rz(2.0172334) q[1];
sx q[1];
rz(-1.1686814) q[1];
sx q[1];
rz(2.7728424) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848727) q[0];
sx q[0];
rz(-0.73000008) q[0];
sx q[0];
rz(-1.2873136) q[0];
rz(-pi) q[1];
rz(-0.36823057) q[2];
sx q[2];
rz(-0.86129649) q[2];
sx q[2];
rz(1.8031098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1146093) q[1];
sx q[1];
rz(-2.8294246) q[1];
sx q[1];
rz(-1.8725558) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7140824) q[3];
sx q[3];
rz(-1.1509794) q[3];
sx q[3];
rz(1.1838208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2345978) q[2];
sx q[2];
rz(-2.5556421) q[2];
sx q[2];
rz(2.9242945) q[2];
rz(-1.6654738) q[3];
sx q[3];
rz(-2.3264591) q[3];
sx q[3];
rz(-2.7998717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9076964) q[0];
sx q[0];
rz(-0.32231575) q[0];
sx q[0];
rz(2.3436558) q[0];
rz(1.7012874) q[1];
sx q[1];
rz(-1.0402352) q[1];
sx q[1];
rz(0.61680102) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42617048) q[0];
sx q[0];
rz(-0.93807332) q[0];
sx q[0];
rz(-1.4769082) q[0];
rz(-pi) q[1];
rz(-1.2347414) q[2];
sx q[2];
rz(-1.0288887) q[2];
sx q[2];
rz(2.0481179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4492256) q[1];
sx q[1];
rz(-1.7150075) q[1];
sx q[1];
rz(1.3647791) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2374898) q[3];
sx q[3];
rz(-2.1573503) q[3];
sx q[3];
rz(1.3400638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1282318) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(-0.70179233) q[2];
rz(-0.93891406) q[3];
sx q[3];
rz(-2.2434668) q[3];
sx q[3];
rz(-1.8963337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.668648) q[0];
sx q[0];
rz(-0.15707792) q[0];
sx q[0];
rz(-3.0182086) q[0];
rz(2.3911047) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(-2.1231245) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7582507) q[0];
sx q[0];
rz(-1.4660379) q[0];
sx q[0];
rz(-1.3488171) q[0];
rz(-pi) q[1];
rz(2.6776601) q[2];
sx q[2];
rz(-1.1008769) q[2];
sx q[2];
rz(-1.7352833) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2606517) q[1];
sx q[1];
rz(-2.6545432) q[1];
sx q[1];
rz(-3.0609291) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6945778) q[3];
sx q[3];
rz(-0.63237337) q[3];
sx q[3];
rz(-3.0770709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6259152) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(0.53543004) q[2];
rz(-1.2299906) q[3];
sx q[3];
rz(-1.001469) q[3];
sx q[3];
rz(3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5984421) q[0];
sx q[0];
rz(-0.66023985) q[0];
sx q[0];
rz(1.862233) q[0];
rz(1.2641501) q[1];
sx q[1];
rz(-1.8708355) q[1];
sx q[1];
rz(-2.989891) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710508) q[0];
sx q[0];
rz(-0.40239516) q[0];
sx q[0];
rz(2.9950525) q[0];
rz(0.13387605) q[2];
sx q[2];
rz(-1.5508442) q[2];
sx q[2];
rz(0.48568113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98916228) q[1];
sx q[1];
rz(-2.1235211) q[1];
sx q[1];
rz(-0.3684045) q[1];
rz(-pi) q[2];
rz(0.72355481) q[3];
sx q[3];
rz(-1.9949942) q[3];
sx q[3];
rz(-1.3217317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5558527) q[2];
sx q[2];
rz(-0.93583217) q[2];
sx q[2];
rz(1.2305416) q[2];
rz(1.3452283) q[3];
sx q[3];
rz(-1.7788922) q[3];
sx q[3];
rz(-1.2983373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718488) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(0.44542435) q[0];
rz(-2.716966) q[1];
sx q[1];
rz(-1.5698965) q[1];
sx q[1];
rz(2.5158023) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9827756) q[0];
sx q[0];
rz(-2.3390769) q[0];
sx q[0];
rz(-2.3071852) q[0];
x q[1];
rz(-0.16846229) q[2];
sx q[2];
rz(-2.0474153) q[2];
sx q[2];
rz(2.1229975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6546331) q[1];
sx q[1];
rz(-0.59096293) q[1];
sx q[1];
rz(-2.0396359) q[1];
rz(-pi) q[2];
rz(-0.15066093) q[3];
sx q[3];
rz(-1.6087349) q[3];
sx q[3];
rz(0.67478363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9618535) q[2];
sx q[2];
rz(-1.0280131) q[2];
sx q[2];
rz(-1.6801768) q[2];
rz(3.0840868) q[3];
sx q[3];
rz(-1.688136) q[3];
sx q[3];
rz(-0.97829515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.7326603) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.59857359) q[0];
rz(-0.21845017) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(-0.54616164) q[2];
sx q[2];
rz(-1.6424137) q[2];
sx q[2];
rz(0.066802468) q[2];
rz(0.784008) q[3];
sx q[3];
rz(-1.5386469) q[3];
sx q[3];
rz(0.28924573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
