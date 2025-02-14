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
rz(0.98969069) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(-1.6852112) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3399062) q[0];
sx q[0];
rz(-2.7481724) q[0];
sx q[0];
rz(-1.2114458) q[0];
rz(-pi) q[1];
rz(-0.13338344) q[2];
sx q[2];
rz(-1.6074153) q[2];
sx q[2];
rz(2.1866148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4850581) q[1];
sx q[1];
rz(-0.63737255) q[1];
sx q[1];
rz(-2.0591048) q[1];
rz(0.043660284) q[3];
sx q[3];
rz(-2.147445) q[3];
sx q[3];
rz(-2.1917564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9911389) q[2];
sx q[2];
rz(-0.011971124) q[2];
sx q[2];
rz(2.0963304) q[2];
rz(0.996905) q[3];
sx q[3];
rz(-0.0054587047) q[3];
sx q[3];
rz(-1.3158984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.588722) q[0];
sx q[0];
rz(-1.9010239) q[0];
sx q[0];
rz(1.7887822) q[0];
rz(-3.1006587) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(-1.5418242) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41087655) q[0];
sx q[0];
rz(-1.4306895) q[0];
sx q[0];
rz(-0.1430169) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2244243) q[2];
sx q[2];
rz(-3.1148095) q[2];
sx q[2];
rz(1.3373122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5168851) q[1];
sx q[1];
rz(-1.5398281) q[1];
sx q[1];
rz(-2.1672441) q[1];
rz(-pi) q[2];
rz(-0.38741855) q[3];
sx q[3];
rz(-1.0325047) q[3];
sx q[3];
rz(1.1673934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.505595) q[2];
sx q[2];
rz(-3.1090241) q[2];
sx q[2];
rz(2.6231498) q[2];
rz(2.017766) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(0.43784416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464624) q[0];
sx q[0];
rz(-3.0914682) q[0];
sx q[0];
rz(-1.5396402) q[0];
rz(0.72499544) q[1];
sx q[1];
rz(-3.1097737) q[1];
sx q[1];
rz(-2.4620788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2452263) q[0];
sx q[0];
rz(-2.1450982) q[0];
sx q[0];
rz(-2.0942874) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9271487) q[2];
sx q[2];
rz(-1.5661998) q[2];
sx q[2];
rz(1.4397804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0371548) q[1];
sx q[1];
rz(-1.8304811) q[1];
sx q[1];
rz(-1.2472769) q[1];
x q[2];
rz(2.1587426) q[3];
sx q[3];
rz(-1.3508537) q[3];
sx q[3];
rz(1.4877121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9069549) q[2];
sx q[2];
rz(-2.8997771) q[2];
sx q[2];
rz(-2.6406636) q[2];
rz(2.6856954) q[3];
sx q[3];
rz(-0.026345043) q[3];
sx q[3];
rz(2.946089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2723715) q[0];
sx q[0];
rz(-3.0932194) q[0];
sx q[0];
rz(0.79917556) q[0];
rz(-0.29798206) q[1];
sx q[1];
rz(-0.25845343) q[1];
sx q[1];
rz(-0.90202773) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1232496) q[0];
sx q[0];
rz(-2.3756174) q[0];
sx q[0];
rz(2.9428456) q[0];
rz(-pi) q[1];
rz(-1.7219826) q[2];
sx q[2];
rz(-2.107321) q[2];
sx q[2];
rz(-0.46162185) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9060256) q[1];
sx q[1];
rz(-0.14793554) q[1];
sx q[1];
rz(-1.6355913) q[1];
rz(-pi) q[2];
rz(-0.055156339) q[3];
sx q[3];
rz(-1.4934908) q[3];
sx q[3];
rz(-0.62420024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67348376) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(1.4704977) q[2];
rz(0.18943131) q[3];
sx q[3];
rz(-3.093284) q[3];
sx q[3];
rz(-0.76907492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.9162132) q[0];
sx q[0];
rz(-3.0140641) q[0];
sx q[0];
rz(-2.7498229) q[0];
rz(-2.0562992) q[1];
sx q[1];
rz(-0.009805209) q[1];
sx q[1];
rz(0.36878961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2854239) q[0];
sx q[0];
rz(-2.4510181) q[0];
sx q[0];
rz(-2.5639433) q[0];
rz(-pi) q[1];
rz(-2.872345) q[2];
sx q[2];
rz(-1.1675861) q[2];
sx q[2];
rz(2.0079812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5318067) q[1];
sx q[1];
rz(-1.5641677) q[1];
sx q[1];
rz(3.1408491) q[1];
rz(-1.3372793) q[3];
sx q[3];
rz(-0.87461014) q[3];
sx q[3];
rz(0.96674572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5607249) q[2];
sx q[2];
rz(-3.0256425) q[2];
sx q[2];
rz(1.7725393) q[2];
rz(-0.12618682) q[3];
sx q[3];
rz(-0.43712619) q[3];
sx q[3];
rz(1.060846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5419902) q[0];
sx q[0];
rz(-0.76102155) q[0];
sx q[0];
rz(-1.5498932) q[0];
rz(0.97269336) q[1];
sx q[1];
rz(-0.21316554) q[1];
sx q[1];
rz(-2.1125643) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9950819) q[0];
sx q[0];
rz(-1.8712988) q[0];
sx q[0];
rz(-2.5142558) q[0];
x q[1];
rz(1.712079) q[2];
sx q[2];
rz(-1.3351235) q[2];
sx q[2];
rz(-0.53776803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29063836) q[1];
sx q[1];
rz(-3.045435) q[1];
sx q[1];
rz(1.566271) q[1];
rz(-pi) q[2];
rz(3.1054822) q[3];
sx q[3];
rz(-1.7381769) q[3];
sx q[3];
rz(-3.0073364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.782393) q[2];
sx q[2];
rz(-1.7155557) q[2];
sx q[2];
rz(2.7101809) q[2];
rz(2.1443478) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3417086) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(-2.0836015) q[0];
rz(2.3175088) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(-0.81738671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0213172) q[0];
sx q[0];
rz(-2.6861354) q[0];
sx q[0];
rz(2.0648455) q[0];
rz(-0.36211966) q[2];
sx q[2];
rz(-3.1220925) q[2];
sx q[2];
rz(1.8060613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83113545) q[1];
sx q[1];
rz(-1.609948) q[1];
sx q[1];
rz(-1.3292946) q[1];
x q[2];
rz(-2.582586) q[3];
sx q[3];
rz(-0.83448258) q[3];
sx q[3];
rz(1.0295537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99007964) q[2];
sx q[2];
rz(-1.9033056) q[2];
sx q[2];
rz(-1.6286758) q[2];
rz(-0.40142909) q[3];
sx q[3];
rz(-0.028086834) q[3];
sx q[3];
rz(2.1872971) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7734739) q[0];
sx q[0];
rz(-0.064726949) q[0];
sx q[0];
rz(-1.7777959) q[0];
rz(-3.0996481) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(-2.1379474) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7391989) q[0];
sx q[0];
rz(-0.93951997) q[0];
sx q[0];
rz(1.4556134) q[0];
rz(-pi) q[1];
rz(-1.9081536) q[2];
sx q[2];
rz(-1.5784831) q[2];
sx q[2];
rz(-0.53042646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81627203) q[1];
sx q[1];
rz(-2.8311816) q[1];
sx q[1];
rz(2.2903633) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3866638) q[3];
sx q[3];
rz(-0.69299504) q[3];
sx q[3];
rz(2.3946144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4369138) q[2];
sx q[2];
rz(-3.0824326) q[2];
sx q[2];
rz(1.3979647) q[2];
rz(2.6513903) q[3];
sx q[3];
rz(-0.043488113) q[3];
sx q[3];
rz(-1.9628261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040084664) q[0];
sx q[0];
rz(-3.0136216) q[0];
sx q[0];
rz(0.21139938) q[0];
rz(-1.118411) q[1];
sx q[1];
rz(-0.011592955) q[1];
sx q[1];
rz(1.5107907) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6965242) q[0];
sx q[0];
rz(-2.7120606) q[0];
sx q[0];
rz(-2.3897445) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62497834) q[2];
sx q[2];
rz(-1.9328261) q[2];
sx q[2];
rz(-3.0292567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76912266) q[1];
sx q[1];
rz(-1.5536397) q[1];
sx q[1];
rz(0.15461289) q[1];
rz(-pi) q[2];
rz(-1.4813194) q[3];
sx q[3];
rz(-1.4050639) q[3];
sx q[3];
rz(0.16752088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3596892) q[2];
sx q[2];
rz(-0.017904559) q[2];
sx q[2];
rz(-2.8622799) q[2];
rz(0.8638047) q[3];
sx q[3];
rz(-3.1371208) q[3];
sx q[3];
rz(0.68252501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.913468) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(1.3155235) q[0];
rz(-0.77519351) q[1];
sx q[1];
rz(-2.9029791) q[1];
sx q[1];
rz(-1.7528037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42487803) q[0];
sx q[0];
rz(-1.8535063) q[0];
sx q[0];
rz(1.7260051) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3451683) q[2];
sx q[2];
rz(-1.9660082) q[2];
sx q[2];
rz(1.164142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9541746) q[1];
sx q[1];
rz(-0.0010962195) q[1];
sx q[1];
rz(-2.5256059) q[1];
x q[2];
rz(2.7241012) q[3];
sx q[3];
rz(-1.1840828) q[3];
sx q[3];
rz(-1.2717122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3751601) q[2];
sx q[2];
rz(-0.015724026) q[2];
sx q[2];
rz(1.6895705) q[2];
rz(-0.25894138) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(1.9848721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5158841) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(1.4868078) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(-0.20023365) q[2];
sx q[2];
rz(-1.5774792) q[2];
sx q[2];
rz(1.7984628) q[2];
rz(1.4987569) q[3];
sx q[3];
rz(-2.7938953) q[3];
sx q[3];
rz(0.10216879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
