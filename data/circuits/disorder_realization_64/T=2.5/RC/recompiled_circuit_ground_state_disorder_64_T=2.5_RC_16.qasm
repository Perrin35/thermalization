OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.46848133) q[0];
sx q[0];
rz(-0.10804478) q[0];
sx q[0];
rz(-2.151902) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(1.4563814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0385168) q[0];
sx q[0];
rz(-1.4355735) q[0];
sx q[0];
rz(-1.9413906) q[0];
x q[1];
rz(1.5338495) q[2];
sx q[2];
rz(-1.4375028) q[2];
sx q[2];
rz(0.61090602) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65653458) q[1];
sx q[1];
rz(-2.5042201) q[1];
sx q[1];
rz(-2.0591048) q[1];
x q[2];
rz(1.637804) q[3];
sx q[3];
rz(-0.57811224) q[3];
sx q[3];
rz(1.029795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1504537) q[2];
sx q[2];
rz(-0.011971124) q[2];
sx q[2];
rz(-2.0963304) q[2];
rz(0.996905) q[3];
sx q[3];
rz(-3.1361339) q[3];
sx q[3];
rz(-1.8256942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.588722) q[0];
sx q[0];
rz(-1.9010239) q[0];
sx q[0];
rz(-1.7887822) q[0];
rz(-0.040933985) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(1.5418242) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1800265) q[0];
sx q[0];
rz(-1.7124023) q[0];
sx q[0];
rz(1.4292634) q[0];
rz(-pi) q[1];
rz(-1.5495318) q[2];
sx q[2];
rz(-1.5545115) q[2];
sx q[2];
rz(-0.88693888) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5168851) q[1];
sx q[1];
rz(-1.5398281) q[1];
sx q[1];
rz(0.9743486) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1349234) q[3];
sx q[3];
rz(-0.65186497) q[3];
sx q[3];
rz(1.3019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6359977) q[2];
sx q[2];
rz(-3.1090241) q[2];
sx q[2];
rz(-0.51844281) q[2];
rz(-2.017766) q[3];
sx q[3];
rz(-0.80945194) q[3];
sx q[3];
rz(-2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951303) q[0];
sx q[0];
rz(-3.0914682) q[0];
sx q[0];
rz(-1.5396402) q[0];
rz(0.72499544) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(-0.67951387) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5711229) q[0];
sx q[0];
rz(-0.75665604) q[0];
sx q[0];
rz(0.65780145) q[0];
rz(-3.1199961) q[2];
sx q[2];
rz(-2.9271002) q[2];
sx q[2];
rz(-0.15211764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5522436) q[1];
sx q[1];
rz(-1.2584983) q[1];
sx q[1];
rz(0.27321385) q[1];
rz(-2.1587426) q[3];
sx q[3];
rz(-1.3508537) q[3];
sx q[3];
rz(-1.4877121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9069549) q[2];
sx q[2];
rz(-0.24181557) q[2];
sx q[2];
rz(2.6406636) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-0.026345043) q[3];
sx q[3];
rz(-2.946089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86922115) q[0];
sx q[0];
rz(-3.0932194) q[0];
sx q[0];
rz(2.3424171) q[0];
rz(-0.29798206) q[1];
sx q[1];
rz(-2.8831392) q[1];
sx q[1];
rz(0.90202773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.018343) q[0];
sx q[0];
rz(-0.76597528) q[0];
sx q[0];
rz(0.19874707) q[0];
x q[1];
rz(-0.24803411) q[2];
sx q[2];
rz(-2.5861861) q[2];
sx q[2];
rz(-2.390304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2711401) q[1];
sx q[1];
rz(-1.5612523) q[1];
sx q[1];
rz(-1.4231667) q[1];
rz(-1.4933735) q[3];
sx q[3];
rz(-1.5158049) q[3];
sx q[3];
rz(-2.1907326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67348376) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(-1.4704977) q[2];
rz(-2.9521613) q[3];
sx q[3];
rz(-0.048308689) q[3];
sx q[3];
rz(0.76907492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9162132) q[0];
sx q[0];
rz(-3.0140641) q[0];
sx q[0];
rz(-2.7498229) q[0];
rz(-2.0562992) q[1];
sx q[1];
rz(-0.009805209) q[1];
sx q[1];
rz(-2.772803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2854239) q[0];
sx q[0];
rz(-2.4510181) q[0];
sx q[0];
rz(0.57764932) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1541751) q[2];
sx q[2];
rz(-1.3236127) q[2];
sx q[2];
rz(-2.8122623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49807355) q[1];
sx q[1];
rz(-0.0066702492) q[1];
sx q[1];
rz(-1.6825063) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4318352) q[3];
sx q[3];
rz(-1.7492948) q[3];
sx q[3];
rz(2.3861726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58086777) q[2];
sx q[2];
rz(-0.11595011) q[2];
sx q[2];
rz(-1.3690534) q[2];
rz(-0.12618682) q[3];
sx q[3];
rz(-2.7044665) q[3];
sx q[3];
rz(-1.060846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5419902) q[0];
sx q[0];
rz(-2.3805711) q[0];
sx q[0];
rz(-1.5916995) q[0];
rz(0.97269336) q[1];
sx q[1];
rz(-0.21316554) q[1];
sx q[1];
rz(1.0290283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21287928) q[0];
sx q[0];
rz(-2.166011) q[0];
sx q[0];
rz(-1.2052324) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9036361) q[2];
sx q[2];
rz(-1.4334442) q[2];
sx q[2];
rz(-0.99983012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8464079) q[1];
sx q[1];
rz(-1.666953) q[1];
sx q[1];
rz(3.1411562) q[1];
x q[2];
rz(0.036110445) q[3];
sx q[3];
rz(-1.7381769) q[3];
sx q[3];
rz(3.0073364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.782393) q[2];
sx q[2];
rz(-1.426037) q[2];
sx q[2];
rz(-0.4314118) q[2];
rz(-2.1443478) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(-0.55716151) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3417086) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(2.0836015) q[0];
rz(-2.3175088) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(0.81738671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0426328) q[0];
sx q[0];
rz(-1.3606679) q[0];
sx q[0];
rz(-1.9779343) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5638875) q[2];
sx q[2];
rz(-1.5525609) q[2];
sx q[2];
rz(1.4438786) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89725607) q[1];
sx q[1];
rz(-0.24459363) q[1];
sx q[1];
rz(1.4084498) q[1];
rz(-2.1002186) q[3];
sx q[3];
rz(-2.2499871) q[3];
sx q[3];
rz(-0.27959945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.151513) q[2];
sx q[2];
rz(-1.9033056) q[2];
sx q[2];
rz(-1.6286758) q[2];
rz(0.40142909) q[3];
sx q[3];
rz(-3.1135058) q[3];
sx q[3];
rz(2.1872971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7734739) q[0];
sx q[0];
rz(-0.064726949) q[0];
sx q[0];
rz(-1.3637967) q[0];
rz(-3.0996481) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(1.0036453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.100228) q[0];
sx q[0];
rz(-1.6637088) q[0];
sx q[0];
rz(0.63444699) q[0];
rz(-1.9081536) q[2];
sx q[2];
rz(-1.5784831) q[2];
sx q[2];
rz(-0.53042646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3253206) q[1];
sx q[1];
rz(-0.31041103) q[1];
sx q[1];
rz(-2.2903633) q[1];
x q[2];
rz(-0.65559994) q[3];
sx q[3];
rz(-1.3274945) q[3];
sx q[3];
rz(2.6213364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7046788) q[2];
sx q[2];
rz(-0.059160058) q[2];
sx q[2];
rz(-1.7436279) q[2];
rz(-2.6513903) q[3];
sx q[3];
rz(-0.043488113) q[3];
sx q[3];
rz(1.9628261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040084664) q[0];
sx q[0];
rz(-3.0136216) q[0];
sx q[0];
rz(2.9301933) q[0];
rz(2.0231817) q[1];
sx q[1];
rz(-0.011592955) q[1];
sx q[1];
rz(1.5107907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.495962) q[0];
sx q[0];
rz(-1.8798774) q[0];
sx q[0];
rz(-1.2675955) q[0];
rz(-pi) q[1];
rz(-0.62497834) q[2];
sx q[2];
rz(-1.9328261) q[2];
sx q[2];
rz(-3.0292567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80434752) q[1];
sx q[1];
rz(-1.7253863) q[1];
sx q[1];
rz(1.5534326) q[1];
rz(0.49064891) q[3];
sx q[3];
rz(-0.18814859) q[3];
sx q[3];
rz(0.33056459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7819034) q[2];
sx q[2];
rz(-3.1236881) q[2];
sx q[2];
rz(0.27931279) q[2];
rz(2.277788) q[3];
sx q[3];
rz(-0.0044718663) q[3];
sx q[3];
rz(0.68252501) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.913468) q[0];
sx q[0];
rz(-2.0371912) q[0];
sx q[0];
rz(-1.3155235) q[0];
rz(-2.3663991) q[1];
sx q[1];
rz(-0.23861353) q[1];
sx q[1];
rz(1.388789) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42487803) q[0];
sx q[0];
rz(-1.2880863) q[0];
sx q[0];
rz(1.7260051) q[0];
rz(-pi) q[1];
rz(-2.7371762) q[2];
sx q[2];
rz(-1.7787654) q[2];
sx q[2];
rz(0.31851092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57143104) q[1];
sx q[1];
rz(-1.5699016) q[1];
sx q[1];
rz(-1.5714297) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78755112) q[3];
sx q[3];
rz(-0.56122223) q[3];
sx q[3];
rz(2.7360327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76643252) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(1.6895705) q[2];
rz(-0.25894138) q[3];
sx q[3];
rz(-0.18766923) q[3];
sx q[3];
rz(1.1567206) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6257085) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(-1.6547849) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(3.1080053) q[2];
sx q[2];
rz(-2.941249) q[2];
sx q[2];
rz(-2.8810101) q[2];
rz(1.9176626) q[3];
sx q[3];
rz(-1.5462688) q[3];
sx q[3];
rz(1.605223) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
