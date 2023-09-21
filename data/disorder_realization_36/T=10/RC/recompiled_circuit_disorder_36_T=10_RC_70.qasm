OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5776423) q[0];
sx q[0];
rz(-2.2221774) q[0];
sx q[0];
rz(1.7786068) q[0];
x q[1];
rz(2.0787813) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(0.21131549) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2540993) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(-1.1521794) q[1];
rz(-pi) q[2];
rz(-2.4926315) q[3];
sx q[3];
rz(-1.0421841) q[3];
sx q[3];
rz(1.987207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802404) q[0];
sx q[0];
rz(-0.70525673) q[0];
sx q[0];
rz(-0.7028701) q[0];
rz(0.77180736) q[2];
sx q[2];
rz(-0.1856005) q[2];
sx q[2];
rz(1.1112569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.37296346) q[1];
sx q[1];
rz(-2.0634723) q[1];
sx q[1];
rz(-1.3586587) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2250697) q[3];
sx q[3];
rz(-2.407981) q[3];
sx q[3];
rz(2.3678126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-0.4804002) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95124328) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074791) q[0];
sx q[0];
rz(-2.0322324) q[0];
sx q[0];
rz(-3.1397318) q[0];
rz(-pi) q[1];
rz(3.0263607) q[2];
sx q[2];
rz(-0.92667246) q[2];
sx q[2];
rz(2.6098721) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.31049) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(2.8527841) q[1];
rz(-pi) q[2];
rz(-2.9680786) q[3];
sx q[3];
rz(-1.1713235) q[3];
sx q[3];
rz(0.80843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(-0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(-0.048359811) q[0];
rz(-0.16391779) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(-1.4455459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3460798) q[0];
sx q[0];
rz(-0.68329408) q[0];
sx q[0];
rz(-2.8542095) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7518026) q[2];
sx q[2];
rz(-2.9082657) q[2];
sx q[2];
rz(-2.0248272) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1452892) q[1];
sx q[1];
rz(-0.22876303) q[1];
sx q[1];
rz(0.28782515) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0914145) q[3];
sx q[3];
rz(-1.319066) q[3];
sx q[3];
rz(-0.43581918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.6131489) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4399453) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(-1.8540927) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094164205) q[0];
sx q[0];
rz(-1.578997) q[0];
sx q[0];
rz(-0.8568944) q[0];
x q[1];
rz(2.6536077) q[2];
sx q[2];
rz(-2.2099566) q[2];
sx q[2];
rz(1.1346863) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0336049) q[1];
sx q[1];
rz(-1.0228979) q[1];
sx q[1];
rz(-0.69544905) q[1];
rz(-pi) q[2];
rz(0.55027996) q[3];
sx q[3];
rz(-1.0770505) q[3];
sx q[3];
rz(-0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.95191082) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(0.50950766) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(0.12621005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.135658) q[0];
sx q[0];
rz(-1.6023484) q[0];
sx q[0];
rz(0.21261442) q[0];
rz(-pi) q[1];
rz(-0.77504471) q[2];
sx q[2];
rz(-1.2942874) q[2];
sx q[2];
rz(1.0724049) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0969442) q[1];
sx q[1];
rz(-1.7587887) q[1];
sx q[1];
rz(-2.2915927) q[1];
rz(-1.0522271) q[3];
sx q[3];
rz(-2.7675981) q[3];
sx q[3];
rz(2.3911589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90298992) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(-2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0347621) q[0];
sx q[0];
rz(-0.97495279) q[0];
sx q[0];
rz(-1.725127) q[0];
rz(2.9497629) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(1.1508458) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0911078) q[1];
sx q[1];
rz(-2.9851966) q[1];
sx q[1];
rz(2.1662103) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5114602) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(0.4160479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(0.5665468) q[3];
sx q[3];
rz(-2.6930801) q[3];
sx q[3];
rz(-2.6678273) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375181) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-2.2699845) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5985142) q[0];
sx q[0];
rz(-0.075965479) q[0];
sx q[0];
rz(0.11673467) q[0];
rz(-pi) q[1];
rz(-1.7224738) q[2];
sx q[2];
rz(-1.9722087) q[2];
sx q[2];
rz(2.7733208) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.826088) q[1];
sx q[1];
rz(-1.2102038) q[1];
sx q[1];
rz(0.70198595) q[1];
rz(-1.7730764) q[3];
sx q[3];
rz(-2.5296122) q[3];
sx q[3];
rz(0.66093854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-1.139572) q[2];
rz(1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(-0.16601673) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4964543) q[0];
sx q[0];
rz(-1.5418515) q[0];
sx q[0];
rz(1.4883947) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7622856) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(-2.7188403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2251687) q[1];
sx q[1];
rz(-1.4125707) q[1];
sx q[1];
rz(-1.2374864) q[1];
rz(-2.7306261) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(0.70511234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-2.7900556) q[2];
rz(-1.0567788) q[3];
sx q[3];
rz(-2.6119699) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(-1.836401) q[0];
rz(0.38048831) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(2.8881853) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6127374) q[0];
sx q[0];
rz(-0.40898541) q[0];
sx q[0];
rz(-1.3091875) q[0];
rz(-pi) q[1];
rz(2.017574) q[2];
sx q[2];
rz(-2.7687216) q[2];
sx q[2];
rz(0.98137059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0512143) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(1.8263837) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5305383) q[3];
sx q[3];
rz(-1.2130034) q[3];
sx q[3];
rz(-0.65294453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9713365) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(2.4304216) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-1.3118369) q[2];
sx q[2];
rz(-1.2987518) q[2];
sx q[2];
rz(-2.1546298) q[2];
rz(2.3861804) q[3];
sx q[3];
rz(-2.073954) q[3];
sx q[3];
rz(3.0975773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
