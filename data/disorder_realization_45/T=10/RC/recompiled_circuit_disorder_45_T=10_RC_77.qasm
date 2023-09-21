OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(-0.64136139) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(1.3134365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9986721) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(-1.8510173) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2506867) q[1];
sx q[1];
rz(-2.2911934) q[1];
sx q[1];
rz(2.5198063) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9786505) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(-0.70139635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.550094) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-0.47810289) q[2];
rz(1.452662) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(2.5959192) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1801382) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(2.1226728) q[0];
rz(-1.4787176) q[1];
sx q[1];
rz(-0.61518413) q[1];
sx q[1];
rz(-0.63308024) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4396673) q[0];
sx q[0];
rz(-0.81322008) q[0];
sx q[0];
rz(-1.4600091) q[0];
rz(-pi) q[1];
rz(2.1585629) q[2];
sx q[2];
rz(-2.2289742) q[2];
sx q[2];
rz(-0.36511974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2411023) q[1];
sx q[1];
rz(-1.2304658) q[1];
sx q[1];
rz(1.9963005) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8301864) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(-0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(3.0241372) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85787073) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-2.03405) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(-0.69082469) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72969499) q[0];
sx q[0];
rz(-1.4955965) q[0];
sx q[0];
rz(0.15326432) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3356528) q[2];
sx q[2];
rz(-2.3217839) q[2];
sx q[2];
rz(-2.2355134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.01268) q[1];
sx q[1];
rz(-1.7528273) q[1];
sx q[1];
rz(0.026489594) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3986474) q[3];
sx q[3];
rz(-2.6446614) q[3];
sx q[3];
rz(0.039615354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(-1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(-0.52247125) q[0];
rz(2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(2.5879588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.722256) q[0];
sx q[0];
rz(-2.0291078) q[0];
sx q[0];
rz(2.3769747) q[0];
rz(0.81981084) q[2];
sx q[2];
rz(-0.95633436) q[2];
sx q[2];
rz(-2.7053506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0370889) q[1];
sx q[1];
rz(-1.4375086) q[1];
sx q[1];
rz(0.78661211) q[1];
rz(-pi) q[2];
rz(-3.1133075) q[3];
sx q[3];
rz(-0.79704282) q[3];
sx q[3];
rz(-0.68238168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5840977) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(-0.90302145) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(-1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49098) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(0.18403149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2011916) q[0];
sx q[0];
rz(-1.860025) q[0];
sx q[0];
rz(1.679923) q[0];
x q[1];
rz(-0.36501546) q[2];
sx q[2];
rz(-1.8301617) q[2];
sx q[2];
rz(2.5600381) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6726471) q[1];
sx q[1];
rz(-0.99860672) q[1];
sx q[1];
rz(-1.0020301) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22589485) q[3];
sx q[3];
rz(-0.73879209) q[3];
sx q[3];
rz(-1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(-0.42282894) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(2.0146577) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891156) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(0.24041644) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(2.9930847) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0200955) q[0];
sx q[0];
rz(-1.7317061) q[0];
sx q[0];
rz(-2.4077329) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4609023) q[2];
sx q[2];
rz(-1.7981148) q[2];
sx q[2];
rz(-2.0109039) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0420694) q[1];
sx q[1];
rz(-1.6798786) q[1];
sx q[1];
rz(-1.8314349) q[1];
rz(-pi) q[2];
rz(-2.6753747) q[3];
sx q[3];
rz(-1.2122452) q[3];
sx q[3];
rz(0.31043226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85577661) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(2.6532069) q[2];
rz(-0.48940247) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69328904) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(-1.12524) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8669406) q[0];
sx q[0];
rz(-0.38892239) q[0];
sx q[0];
rz(-1.8453983) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2721887) q[2];
sx q[2];
rz(-1.4223137) q[2];
sx q[2];
rz(2.6872925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69562558) q[1];
sx q[1];
rz(-1.274316) q[1];
sx q[1];
rz(0.28709025) q[1];
x q[2];
rz(-1.0057698) q[3];
sx q[3];
rz(-2.1566448) q[3];
sx q[3];
rz(-2.1573531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0685048) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(-1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97061625) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(-1.9203141) q[1];
sx q[1];
rz(-1.5131283) q[1];
sx q[1];
rz(2.1309526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686417) q[0];
sx q[0];
rz(-0.7612345) q[0];
sx q[0];
rz(1.8038473) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8433617) q[2];
sx q[2];
rz(-2.0474307) q[2];
sx q[2];
rz(0.44520608) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9606564) q[1];
sx q[1];
rz(-1.8005383) q[1];
sx q[1];
rz(-0.22682637) q[1];
x q[2];
rz(2.9270494) q[3];
sx q[3];
rz(-2.7326267) q[3];
sx q[3];
rz(1.4734801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(2.5726035) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(0.49474299) q[0];
rz(0.61839473) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(3.0659952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8537124) q[0];
sx q[0];
rz(-1.8822077) q[0];
sx q[0];
rz(2.1561554) q[0];
rz(-1.7463023) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(-2.4621071) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1765715) q[1];
sx q[1];
rz(-2.2551943) q[1];
sx q[1];
rz(1.4375163) q[1];
rz(-pi) q[2];
rz(-1.3887614) q[3];
sx q[3];
rz(-0.958003) q[3];
sx q[3];
rz(-2.777367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8081234) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(-1.8010275) q[2];
rz(-2.8373485) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-2.9647968) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84687418) q[0];
sx q[0];
rz(-0.28781578) q[0];
sx q[0];
rz(-0.78886445) q[0];
rz(-pi) q[1];
rz(-2.6920325) q[2];
sx q[2];
rz(-0.170378) q[2];
sx q[2];
rz(0.61194387) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2446574) q[1];
sx q[1];
rz(-1.999563) q[1];
sx q[1];
rz(2.7662617) q[1];
rz(-pi) q[2];
rz(-1.8701843) q[3];
sx q[3];
rz(-0.96794879) q[3];
sx q[3];
rz(0.30764461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(2.2484696) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72398913) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(3.1148615) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-2.2970207) q[2];
sx q[2];
rz(-1.9324586) q[2];
sx q[2];
rz(2.4697138) q[2];
rz(2.0704913) q[3];
sx q[3];
rz(-1.0715967) q[3];
sx q[3];
rz(-1.788492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];