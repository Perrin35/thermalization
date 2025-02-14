OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6731113) q[0];
sx q[0];
rz(-3.0335479) q[0];
sx q[0];
rz(-0.98969069) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(1.4563814) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8016865) q[0];
sx q[0];
rz(-2.7481724) q[0];
sx q[0];
rz(1.2114458) q[0];
x q[1];
rz(-3.0082092) q[2];
sx q[2];
rz(-1.6074153) q[2];
sx q[2];
rz(-2.1866148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0691245) q[1];
sx q[1];
rz(-2.1241423) q[1];
sx q[1];
rz(2.8072559) q[1];
x q[2];
rz(-0.043660284) q[3];
sx q[3];
rz(-0.99414765) q[3];
sx q[3];
rz(0.94983627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1504537) q[2];
sx q[2];
rz(-0.011971124) q[2];
sx q[2];
rz(2.0963304) q[2];
rz(0.996905) q[3];
sx q[3];
rz(-3.1361339) q[3];
sx q[3];
rz(-1.8256942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.588722) q[0];
sx q[0];
rz(-1.2405688) q[0];
sx q[0];
rz(-1.3528104) q[0];
rz(-0.040933985) q[1];
sx q[1];
rz(-1.2177688) q[1];
sx q[1];
rz(1.5997684) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41087655) q[0];
sx q[0];
rz(-1.4306895) q[0];
sx q[0];
rz(2.9985757) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91716839) q[2];
sx q[2];
rz(-0.02678314) q[2];
sx q[2];
rz(1.8042804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92506986) q[1];
sx q[1];
rz(-2.1669186) q[1];
sx q[1];
rz(0.03742569) q[1];
rz(-pi) q[2];
rz(-1.0066693) q[3];
sx q[3];
rz(-2.4897277) q[3];
sx q[3];
rz(1.3019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6359977) q[2];
sx q[2];
rz(-3.1090241) q[2];
sx q[2];
rz(-0.51844281) q[2];
rz(-2.017766) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951303) q[0];
sx q[0];
rz(-0.050124425) q[0];
sx q[0];
rz(1.5396402) q[0];
rz(0.72499544) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(2.4620788) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8963663) q[0];
sx q[0];
rz(-0.99649444) q[0];
sx q[0];
rz(-1.0473053) q[0];
x q[1];
rz(1.5660921) q[2];
sx q[2];
rz(-1.3563547) q[2];
sx q[2];
rz(0.1300148) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.589349) q[1];
sx q[1];
rz(-1.2584983) q[1];
sx q[1];
rz(-2.8683788) q[1];
rz(-1.1876538) q[3];
sx q[3];
rz(-2.5184298) q[3];
sx q[3];
rz(2.9083657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2346377) q[2];
sx q[2];
rz(-0.24181557) q[2];
sx q[2];
rz(0.50092906) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-3.1152476) q[3];
sx q[3];
rz(2.946089) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86922115) q[0];
sx q[0];
rz(-0.048373241) q[0];
sx q[0];
rz(2.3424171) q[0];
rz(-0.29798206) q[1];
sx q[1];
rz(-2.8831392) q[1];
sx q[1];
rz(-2.2395649) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74583861) q[0];
sx q[0];
rz(-0.82358783) q[0];
sx q[0];
rz(-1.7584778) q[0];
rz(-pi) q[1];
rz(0.54157864) q[2];
sx q[2];
rz(-1.4409833) q[2];
sx q[2];
rz(2.1101348) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2711401) q[1];
sx q[1];
rz(-1.5612523) q[1];
sx q[1];
rz(1.4231667) q[1];
rz(-pi) q[2];
rz(1.6482192) q[3];
sx q[3];
rz(-1.6257878) q[3];
sx q[3];
rz(-0.95086004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67348376) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(1.4704977) q[2];
rz(-0.18943131) q[3];
sx q[3];
rz(-0.048308689) q[3];
sx q[3];
rz(2.3725177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9162132) q[0];
sx q[0];
rz(-0.12752859) q[0];
sx q[0];
rz(-0.39176971) q[0];
rz(-1.0852934) q[1];
sx q[1];
rz(-0.009805209) q[1];
sx q[1];
rz(2.772803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2854239) q[0];
sx q[0];
rz(-0.69057451) q[0];
sx q[0];
rz(2.5639433) q[0];
rz(-pi) q[1];
x q[1];
rz(2.872345) q[2];
sx q[2];
rz(-1.1675861) q[2];
sx q[2];
rz(1.1336114) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5318067) q[1];
sx q[1];
rz(-1.5641677) q[1];
sx q[1];
rz(-3.1408491) q[1];
rz(-2.4318352) q[3];
sx q[3];
rz(-1.7492948) q[3];
sx q[3];
rz(-0.75542007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58086777) q[2];
sx q[2];
rz(-0.11595011) q[2];
sx q[2];
rz(-1.3690534) q[2];
rz(-0.12618682) q[3];
sx q[3];
rz(-0.43712619) q[3];
sx q[3];
rz(1.060846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5419902) q[0];
sx q[0];
rz(-2.3805711) q[0];
sx q[0];
rz(-1.5498932) q[0];
rz(0.97269336) q[1];
sx q[1];
rz(-2.9284271) q[1];
sx q[1];
rz(-1.0290283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3297076) q[0];
sx q[0];
rz(-2.4548479) q[0];
sx q[0];
rz(-0.48573692) q[0];
x q[1];
rz(-2.9036361) q[2];
sx q[2];
rz(-1.4334442) q[2];
sx q[2];
rz(-0.99983012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8659391) q[1];
sx q[1];
rz(-1.5712308) q[1];
sx q[1];
rz(-1.666953) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3602882) q[3];
sx q[3];
rz(-2.9703968) q[3];
sx q[3];
rz(3.0623113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35919967) q[2];
sx q[2];
rz(-1.426037) q[2];
sx q[2];
rz(0.4314118) q[2];
rz(-0.99724489) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3417086) q[0];
sx q[0];
rz(-2.6567617) q[0];
sx q[0];
rz(1.0579911) q[0];
rz(2.3175088) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(2.3242059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0426328) q[0];
sx q[0];
rz(-1.7809247) q[0];
sx q[0];
rz(1.9779343) q[0];
rz(-pi) q[1];
rz(0.018235834) q[2];
sx q[2];
rz(-1.577704) q[2];
sx q[2];
rz(-3.0148009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73002023) q[1];
sx q[1];
rz(-1.8121094) q[1];
sx q[1];
rz(0.040320591) q[1];
rz(-pi) q[2];
rz(2.1002186) q[3];
sx q[3];
rz(-2.2499871) q[3];
sx q[3];
rz(0.27959945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.151513) q[2];
sx q[2];
rz(-1.9033056) q[2];
sx q[2];
rz(-1.6286758) q[2];
rz(2.7401636) q[3];
sx q[3];
rz(-0.028086834) q[3];
sx q[3];
rz(-0.95429558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3681188) q[0];
sx q[0];
rz(-0.064726949) q[0];
sx q[0];
rz(-1.3637967) q[0];
rz(3.0996481) q[1];
sx q[1];
rz(-3.0105803) q[1];
sx q[1];
rz(-2.1379474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59597385) q[0];
sx q[0];
rz(-2.501308) q[0];
sx q[0];
rz(0.15592928) q[0];
rz(-pi) q[1];
rz(-1.233439) q[2];
sx q[2];
rz(-1.5784831) q[2];
sx q[2];
rz(0.53042646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4499047) q[1];
sx q[1];
rz(-1.7734911) q[1];
sx q[1];
rz(-1.8075289) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2673301) q[3];
sx q[3];
rz(-0.93765536) q[3];
sx q[3];
rz(-1.9078524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7046788) q[2];
sx q[2];
rz(-0.059160058) q[2];
sx q[2];
rz(-1.3979647) q[2];
rz(-2.6513903) q[3];
sx q[3];
rz(-3.0981045) q[3];
sx q[3];
rz(1.1787666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.101508) q[0];
sx q[0];
rz(-0.12797102) q[0];
sx q[0];
rz(-0.21139938) q[0];
rz(-1.118411) q[1];
sx q[1];
rz(-3.1299997) q[1];
sx q[1];
rz(-1.5107907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4450685) q[0];
sx q[0];
rz(-0.42953209) q[0];
sx q[0];
rz(2.3897445) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5166143) q[2];
sx q[2];
rz(-1.9328261) q[2];
sx q[2];
rz(0.11233594) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69204077) q[1];
sx q[1];
rz(-0.15555432) q[1];
sx q[1];
rz(-0.1109619) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9752067) q[3];
sx q[3];
rz(-1.659044) q[3];
sx q[3];
rz(1.4180753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3596892) q[2];
sx q[2];
rz(-3.1236881) q[2];
sx q[2];
rz(2.8622799) q[2];
rz(2.277788) q[3];
sx q[3];
rz(-0.0044718663) q[3];
sx q[3];
rz(-2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.913468) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(1.3155235) q[0];
rz(0.77519351) q[1];
sx q[1];
rz(-0.23861353) q[1];
sx q[1];
rz(-1.7528037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9520541) q[0];
sx q[0];
rz(-1.4217958) q[0];
sx q[0];
rz(0.28596626) q[0];
rz(-pi) q[1];
rz(-2.7371762) q[2];
sx q[2];
rz(-1.7787654) q[2];
sx q[2];
rz(0.31851092) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57143104) q[1];
sx q[1];
rz(-1.5716911) q[1];
sx q[1];
rz(1.5714297) q[1];
rz(-pi) q[2];
rz(-0.41749145) q[3];
sx q[3];
rz(-1.9575099) q[3];
sx q[3];
rz(-1.8698805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76643252) q[2];
sx q[2];
rz(-0.015724026) q[2];
sx q[2];
rz(1.4520221) q[2];
rz(0.25894138) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(-1.9848721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6257085) q[0];
sx q[0];
rz(-0.71737552) q[0];
sx q[0];
rz(1.4364568) q[0];
rz(1.4868078) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(1.5639772) q[2];
sx q[2];
rz(-1.3705672) q[2];
sx q[2];
rz(-2.9152824) q[2];
rz(-0.026080118) q[3];
sx q[3];
rz(-1.2240388) q[3];
sx q[3];
rz(0.025561533) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
