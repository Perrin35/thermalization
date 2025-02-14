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
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(1.3527704) q[0];
rz(-1.9569205) q[1];
sx q[1];
rz(-2.4689622) q[1];
sx q[1];
rz(-1.5157359) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0294581) q[0];
sx q[0];
rz(-1.7687651) q[0];
sx q[0];
rz(-0.1898808) q[0];
rz(2.3358222) q[2];
sx q[2];
rz(-1.9293074) q[2];
sx q[2];
rz(1.9961949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3876311) q[1];
sx q[1];
rz(-1.921614) q[1];
sx q[1];
rz(-2.6430348) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58324849) q[3];
sx q[3];
rz(-2.6494559) q[3];
sx q[3];
rz(1.0517091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0036014) q[2];
sx q[2];
rz(-3.0333952) q[2];
sx q[2];
rz(0.18599621) q[2];
rz(-3.1229535) q[3];
sx q[3];
rz(-1.847307) q[3];
sx q[3];
rz(-1.0703465) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828534) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(-2.8232316) q[0];
rz(0.63703498) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(2.6723518) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97142727) q[0];
sx q[0];
rz(-1.3999456) q[0];
sx q[0];
rz(0.20194443) q[0];
rz(-pi) q[1];
rz(3.0576433) q[2];
sx q[2];
rz(-1.9453887) q[2];
sx q[2];
rz(-0.89886802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.138006) q[1];
sx q[1];
rz(-0.39630656) q[1];
sx q[1];
rz(2.519964) q[1];
rz(3.1358767) q[3];
sx q[3];
rz(-0.88086956) q[3];
sx q[3];
rz(1.6315414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3176754) q[2];
sx q[2];
rz(-2.9313512) q[2];
sx q[2];
rz(2.4483185) q[2];
rz(-1.3200101) q[3];
sx q[3];
rz(-1.8157248) q[3];
sx q[3];
rz(-2.0461953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005118) q[0];
sx q[0];
rz(-0.63758272) q[0];
sx q[0];
rz(2.6620423) q[0];
rz(2.6374822) q[1];
sx q[1];
rz(-1.9934374) q[1];
sx q[1];
rz(3.0832916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2445912) q[0];
sx q[0];
rz(-2.069733) q[0];
sx q[0];
rz(-1.8656436) q[0];
rz(-pi) q[1];
rz(-0.68966663) q[2];
sx q[2];
rz(-1.153077) q[2];
sx q[2];
rz(2.6725685) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41861288) q[1];
sx q[1];
rz(-1.449391) q[1];
sx q[1];
rz(2.1239807) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.635635) q[3];
sx q[3];
rz(-1.3383597) q[3];
sx q[3];
rz(2.1955804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7551859) q[2];
sx q[2];
rz(-1.9599954) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(1.6218328) q[3];
sx q[3];
rz(-0.45791364) q[3];
sx q[3];
rz(-2.7601385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.195049) q[0];
sx q[0];
rz(-0.46849546) q[0];
sx q[0];
rz(3.0750437) q[0];
rz(-1.6171803) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(-2.4066511) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682869) q[0];
sx q[0];
rz(-2.5096008) q[0];
sx q[0];
rz(0.39936622) q[0];
x q[1];
rz(2.8822172) q[2];
sx q[2];
rz(-0.7470567) q[2];
sx q[2];
rz(-2.7310028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1030324) q[1];
sx q[1];
rz(-1.2913229) q[1];
sx q[1];
rz(-2.7725262) q[1];
rz(-2.9959277) q[3];
sx q[3];
rz(-2.2861135) q[3];
sx q[3];
rz(0.13038929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28896004) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(-0.10870474) q[2];
rz(-0.92681256) q[3];
sx q[3];
rz(-2.1709397) q[3];
sx q[3];
rz(-1.5638117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0842593) q[0];
sx q[0];
rz(-2.5045392) q[0];
sx q[0];
rz(2.0507226) q[0];
rz(-0.57251656) q[1];
sx q[1];
rz(-1.1895836) q[1];
sx q[1];
rz(-0.82430878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6466581) q[0];
sx q[0];
rz(-1.3368155) q[0];
sx q[0];
rz(-1.184638) q[0];
rz(0.53218158) q[2];
sx q[2];
rz(-2.144118) q[2];
sx q[2];
rz(-0.5328726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5368176) q[1];
sx q[1];
rz(-2.1120434) q[1];
sx q[1];
rz(2.6763335) q[1];
x q[2];
rz(-2.0044672) q[3];
sx q[3];
rz(-1.9423565) q[3];
sx q[3];
rz(-2.5783758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6606286) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(-2.8503897) q[2];
rz(-0.63699841) q[3];
sx q[3];
rz(-1.6773418) q[3];
sx q[3];
rz(1.252482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7630735) q[0];
sx q[0];
rz(-0.28994361) q[0];
sx q[0];
rz(-0.13105233) q[0];
rz(-1.1669) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(3.1189721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2289091) q[0];
sx q[0];
rz(-2.1118409) q[0];
sx q[0];
rz(1.9846657) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51949595) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(2.1313388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9359259) q[1];
sx q[1];
rz(-2.4956411) q[1];
sx q[1];
rz(-1.0420858) q[1];
x q[2];
rz(1.3778481) q[3];
sx q[3];
rz(-1.7304599) q[3];
sx q[3];
rz(2.1482688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8733946) q[2];
sx q[2];
rz(-0.70316535) q[2];
sx q[2];
rz(-1.084729) q[2];
rz(1.9966513) q[3];
sx q[3];
rz(-1.2866674) q[3];
sx q[3];
rz(1.0219319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5893843) q[0];
sx q[0];
rz(-1.6227868) q[0];
sx q[0];
rz(-2.6680706) q[0];
rz(-1.9206958) q[1];
sx q[1];
rz(-0.41243204) q[1];
sx q[1];
rz(-1.8730877) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3800156) q[0];
sx q[0];
rz(-0.76671874) q[0];
sx q[0];
rz(-0.79464998) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2008197) q[2];
sx q[2];
rz(-0.93859276) q[2];
sx q[2];
rz(-0.8559627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2655609) q[1];
sx q[1];
rz(-0.68194333) q[1];
sx q[1];
rz(-0.69721402) q[1];
rz(1.394519) q[3];
sx q[3];
rz(-0.44197861) q[3];
sx q[3];
rz(0.37261841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1641757) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(0.24641985) q[2];
rz(-0.94270802) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(-2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3281658) q[0];
sx q[0];
rz(-1.270371) q[0];
sx q[0];
rz(-3.0810007) q[0];
rz(1.5024441) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(-1.5938119) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0133724) q[0];
sx q[0];
rz(-0.0039847535) q[0];
sx q[0];
rz(0.26040034) q[0];
rz(2.8558989) q[2];
sx q[2];
rz(-2.5473352) q[2];
sx q[2];
rz(-2.2069634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8503849) q[1];
sx q[1];
rz(-0.7813079) q[1];
sx q[1];
rz(2.2154097) q[1];
x q[2];
rz(1.1154864) q[3];
sx q[3];
rz(-0.90789617) q[3];
sx q[3];
rz(0.864322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8760406) q[2];
sx q[2];
rz(-2.2206842) q[2];
sx q[2];
rz(1.457224) q[2];
rz(-0.11988104) q[3];
sx q[3];
rz(-2.9370152) q[3];
sx q[3];
rz(-1.0286819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(2.0062398) q[0];
rz(-0.10336939) q[1];
sx q[1];
rz(-1.7366948) q[1];
sx q[1];
rz(2.1939383) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4554868) q[0];
sx q[0];
rz(-2.4999833) q[0];
sx q[0];
rz(-3.0296586) q[0];
rz(-pi) q[1];
rz(1.0815166) q[2];
sx q[2];
rz(-0.7502509) q[2];
sx q[2];
rz(-2.3231151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3909437) q[1];
sx q[1];
rz(-0.90439561) q[1];
sx q[1];
rz(2.8007068) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71368788) q[3];
sx q[3];
rz(-2.2550224) q[3];
sx q[3];
rz(2.8005176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7171628) q[2];
sx q[2];
rz(-2.5134176) q[2];
sx q[2];
rz(2.2372645) q[2];
rz(1.8782015) q[3];
sx q[3];
rz(-1.2369913) q[3];
sx q[3];
rz(-0.5164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(-0.98536056) q[0];
rz(-2.789978) q[1];
sx q[1];
rz(-2.1814587) q[1];
sx q[1];
rz(-2.7947289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540005) q[0];
sx q[0];
rz(-1.3525602) q[0];
sx q[0];
rz(-2.2130593) q[0];
x q[1];
rz(2.0167256) q[2];
sx q[2];
rz(-1.4691778) q[2];
sx q[2];
rz(-1.1559238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0488009) q[1];
sx q[1];
rz(-1.4818703) q[1];
sx q[1];
rz(-0.6912749) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7010491) q[3];
sx q[3];
rz(-1.2367354) q[3];
sx q[3];
rz(1.3584709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88006192) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(-3.0066709) q[2];
rz(3.0123582) q[3];
sx q[3];
rz(-0.40004572) q[3];
sx q[3];
rz(1.038704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018910949) q[0];
sx q[0];
rz(-1.8862579) q[0];
sx q[0];
rz(0.12611623) q[0];
rz(0.46161721) q[1];
sx q[1];
rz(-1.4001662) q[1];
sx q[1];
rz(-2.9495159) q[1];
rz(0.89712894) q[2];
sx q[2];
rz(-1.9403602) q[2];
sx q[2];
rz(-0.42862949) q[2];
rz(-1.73883) q[3];
sx q[3];
rz(-1.8258655) q[3];
sx q[3];
rz(-2.0139555) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
