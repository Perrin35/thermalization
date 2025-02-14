OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.42263) q[0];
sx q[0];
rz(-2.842272) q[0];
sx q[0];
rz(2.646995) q[0];
rz(-1.9994796) q[1];
sx q[1];
rz(-2.1358868) q[1];
sx q[1];
rz(-1.1297273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0690735) q[0];
sx q[0];
rz(-1.4125709) q[0];
sx q[0];
rz(2.4618966) q[0];
x q[1];
rz(-0.45290516) q[2];
sx q[2];
rz(-1.1158841) q[2];
sx q[2];
rz(-2.186113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8985624) q[1];
sx q[1];
rz(-0.77706017) q[1];
sx q[1];
rz(-1.206336) q[1];
rz(1.0315597) q[3];
sx q[3];
rz(-2.4162216) q[3];
sx q[3];
rz(1.6417208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31892458) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(-0.85282105) q[2];
rz(-1.8114113) q[3];
sx q[3];
rz(-2.4460402) q[3];
sx q[3];
rz(-3.0947963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80832076) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(1.464123) q[0];
rz(-1.4785712) q[1];
sx q[1];
rz(-1.1590978) q[1];
sx q[1];
rz(-1.9333855) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0639275) q[0];
sx q[0];
rz(-2.0443235) q[0];
sx q[0];
rz(2.5350476) q[0];
rz(-pi) q[1];
rz(-1.0958395) q[2];
sx q[2];
rz(-1.37687) q[2];
sx q[2];
rz(-0.79481193) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18259163) q[1];
sx q[1];
rz(-1.6920857) q[1];
sx q[1];
rz(1.8392929) q[1];
x q[2];
rz(-0.34388108) q[3];
sx q[3];
rz(-0.47894997) q[3];
sx q[3];
rz(-0.44287455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6055484) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(-1.4844683) q[2];
rz(3.1260955) q[3];
sx q[3];
rz(-1.9280547) q[3];
sx q[3];
rz(-1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.991796) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(0.27161828) q[0];
rz(-2.244921) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(-1.8439878) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54416849) q[0];
sx q[0];
rz(-2.6784458) q[0];
sx q[0];
rz(-1.5933655) q[0];
x q[1];
rz(2.6372725) q[2];
sx q[2];
rz(-1.2337304) q[2];
sx q[2];
rz(-0.45318174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43186346) q[1];
sx q[1];
rz(-1.6846906) q[1];
sx q[1];
rz(1.8191871) q[1];
x q[2];
rz(3.1362757) q[3];
sx q[3];
rz(-0.76197366) q[3];
sx q[3];
rz(3.0865106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74035949) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(-0.05376251) q[2];
rz(1.2939804) q[3];
sx q[3];
rz(-2.3432178) q[3];
sx q[3];
rz(0.23538858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2006705) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(1.4642375) q[0];
rz(-1.1306521) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(3.0272223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76246156) q[0];
sx q[0];
rz(-1.987559) q[0];
sx q[0];
rz(-2.5083191) q[0];
rz(-pi) q[1];
rz(-2.4855544) q[2];
sx q[2];
rz(-1.1636574) q[2];
sx q[2];
rz(0.40358692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3614157) q[1];
sx q[1];
rz(-0.89539278) q[1];
sx q[1];
rz(-0.03390582) q[1];
x q[2];
rz(2.9435677) q[3];
sx q[3];
rz(-1.5772444) q[3];
sx q[3];
rz(-0.40606582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6681246) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(3.065897) q[2];
rz(-2.7068052) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(-3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39744034) q[0];
sx q[0];
rz(-1.8295153) q[0];
sx q[0];
rz(0.42006668) q[0];
rz(0.21332598) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(1.2423645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7540589) q[0];
sx q[0];
rz(-2.2926169) q[0];
sx q[0];
rz(-0.57768627) q[0];
x q[1];
rz(-1.2790643) q[2];
sx q[2];
rz(-1.9856678) q[2];
sx q[2];
rz(-1.8366829) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8458357) q[1];
sx q[1];
rz(-1.6713665) q[1];
sx q[1];
rz(0.015458903) q[1];
rz(2.9202634) q[3];
sx q[3];
rz(-1.6175272) q[3];
sx q[3];
rz(-0.58906259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3199978) q[2];
sx q[2];
rz(-0.91511202) q[2];
sx q[2];
rz(2.7670009) q[2];
rz(1.0125259) q[3];
sx q[3];
rz(-0.39207021) q[3];
sx q[3];
rz(-2.4969126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11151611) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(-0.43055713) q[0];
rz(2.997609) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(2.5103501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66723204) q[0];
sx q[0];
rz(-0.96275126) q[0];
sx q[0];
rz(1.505018) q[0];
rz(-pi) q[1];
rz(1.9912791) q[2];
sx q[2];
rz(-1.704272) q[2];
sx q[2];
rz(-1.6949289) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5908554) q[1];
sx q[1];
rz(-0.8181347) q[1];
sx q[1];
rz(-0.70597752) q[1];
rz(0.9489551) q[3];
sx q[3];
rz(-2.6027711) q[3];
sx q[3];
rz(-2.1486189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1708019) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(-1.7010752) q[2];
rz(-1.7801646) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.68706566) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(-2.2173296) q[0];
rz(-0.29280064) q[1];
sx q[1];
rz(-1.9127138) q[1];
sx q[1];
rz(2.3407095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070904254) q[0];
sx q[0];
rz(-1.3137806) q[0];
sx q[0];
rz(-1.9563673) q[0];
x q[1];
rz(2.7960294) q[2];
sx q[2];
rz(-2.2032732) q[2];
sx q[2];
rz(-1.3407624) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.427721) q[1];
sx q[1];
rz(-1.6000976) q[1];
sx q[1];
rz(-3.0364365) q[1];
rz(-pi) q[2];
rz(3.1195779) q[3];
sx q[3];
rz(-1.5436633) q[3];
sx q[3];
rz(-0.61329816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.822829) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(1.2804383) q[2];
rz(2.473623) q[3];
sx q[3];
rz(-2.6536055) q[3];
sx q[3];
rz(-0.41332301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845344) q[0];
sx q[0];
rz(-1.2385383) q[0];
sx q[0];
rz(-1.9588233) q[0];
rz(-0.23712748) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(3.0822486) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6366258) q[0];
sx q[0];
rz(-1.6904181) q[0];
sx q[0];
rz(1.8455532) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4164532) q[2];
sx q[2];
rz(-1.1584917) q[2];
sx q[2];
rz(-2.7453842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13798156) q[1];
sx q[1];
rz(-1.4661769) q[1];
sx q[1];
rz(-1.305718) q[1];
rz(-pi) q[2];
rz(1.8944727) q[3];
sx q[3];
rz(-1.4973065) q[3];
sx q[3];
rz(2.2761114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5817029) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(-1.7891368) q[2];
rz(-1.7743568) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(2.2860315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(2.6819041) q[0];
rz(-0.028060878) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(-1.185816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677886) q[0];
sx q[0];
rz(-1.8361183) q[0];
sx q[0];
rz(-1.5242884) q[0];
rz(0.90699754) q[2];
sx q[2];
rz(-3.0131648) q[2];
sx q[2];
rz(-2.3333486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.016077135) q[1];
sx q[1];
rz(-2.0488157) q[1];
sx q[1];
rz(2.1090871) q[1];
x q[2];
rz(0.61798685) q[3];
sx q[3];
rz(-1.8027961) q[3];
sx q[3];
rz(-0.31832507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49843732) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(0.15677491) q[2];
rz(0.59569851) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(-3.116385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6237685) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(-0.18381707) q[0];
rz(-0.078016438) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-2.877291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.030609) q[0];
sx q[0];
rz(-1.8389337) q[0];
sx q[0];
rz(0.58695745) q[0];
x q[1];
rz(-1.043522) q[2];
sx q[2];
rz(-0.95607483) q[2];
sx q[2];
rz(-1.5103024) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1127013) q[1];
sx q[1];
rz(-1.3910146) q[1];
sx q[1];
rz(1.4110908) q[1];
rz(-pi) q[2];
rz(0.49510689) q[3];
sx q[3];
rz(-1.4864085) q[3];
sx q[3];
rz(1.6364678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(-1.7372355) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(3.0586045) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(0.71113853) q[1];
sx q[1];
rz(-1.8322721) q[1];
sx q[1];
rz(-2.4087404) q[1];
rz(-0.044471519) q[2];
sx q[2];
rz(-1.0723249) q[2];
sx q[2];
rz(1.9981801) q[2];
rz(2.5855999) q[3];
sx q[3];
rz(-1.0992194) q[3];
sx q[3];
rz(1.9627375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
