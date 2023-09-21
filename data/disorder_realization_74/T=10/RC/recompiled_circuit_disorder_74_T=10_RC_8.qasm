OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(-0.98288012) q[0];
sx q[0];
rz(1.0097526) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(-0.90254012) q[1];
sx q[1];
rz(-1.8523822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6362808) q[0];
sx q[0];
rz(-1.8587776) q[0];
sx q[0];
rz(-2.050839) q[0];
rz(0.024359811) q[2];
sx q[2];
rz(-1.5775541) q[2];
sx q[2];
rz(-2.9565405) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9413486) q[1];
sx q[1];
rz(-1.3625047) q[1];
sx q[1];
rz(1.8336481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3025769) q[3];
sx q[3];
rz(-1.5032288) q[3];
sx q[3];
rz(0.64642954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5499128) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(0.13554779) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-2.2048075) q[0];
sx q[0];
rz(-2.360789) q[0];
rz(-2.8813598) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.388893) q[0];
sx q[0];
rz(-1.7912731) q[0];
sx q[0];
rz(1.7874784) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8333089) q[2];
sx q[2];
rz(-1.1051902) q[2];
sx q[2];
rz(-0.75454933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0158851) q[1];
sx q[1];
rz(-1.410555) q[1];
sx q[1];
rz(-1.4829774) q[1];
x q[2];
rz(2.1900858) q[3];
sx q[3];
rz(-0.12976876) q[3];
sx q[3];
rz(1.4459923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.6080674) q[2];
sx q[2];
rz(2.1874645) q[2];
rz(-1.4387087) q[3];
sx q[3];
rz(-1.3043159) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057782877) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(-2.4457248) q[0];
rz(-0.35481915) q[1];
sx q[1];
rz(-1.6638919) q[1];
sx q[1];
rz(-2.9755039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3375219) q[0];
sx q[0];
rz(-1.7162168) q[0];
sx q[0];
rz(-1.8619596) q[0];
rz(-pi) q[1];
rz(2.4737349) q[2];
sx q[2];
rz(-0.72615004) q[2];
sx q[2];
rz(-2.6454676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6925466) q[1];
sx q[1];
rz(-2.0448301) q[1];
sx q[1];
rz(-1.0707335) q[1];
rz(2.8852799) q[3];
sx q[3];
rz(-0.87402065) q[3];
sx q[3];
rz(3.0031406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7819536) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-0.98177838) q[2];
rz(2.0180457) q[3];
sx q[3];
rz(-0.26505622) q[3];
sx q[3];
rz(1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05900255) q[0];
sx q[0];
rz(-0.95472097) q[0];
sx q[0];
rz(-2.847805) q[0];
rz(2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(0.77484432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9967277) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(-3.0680455) q[0];
rz(-1.5261126) q[2];
sx q[2];
rz(-1.9384906) q[2];
sx q[2];
rz(-1.9862663) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2312154) q[1];
sx q[1];
rz(-0.41035715) q[1];
sx q[1];
rz(1.3740205) q[1];
x q[2];
rz(-0.1713486) q[3];
sx q[3];
rz(-1.8257986) q[3];
sx q[3];
rz(1.9920497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1321156) q[2];
sx q[2];
rz(-1.0813794) q[2];
sx q[2];
rz(-0.57007989) q[2];
rz(-1.8360957) q[3];
sx q[3];
rz(-2.2142742) q[3];
sx q[3];
rz(-1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(2.9602125) q[0];
rz(-0.60896215) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(-3.0338874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2337423) q[0];
sx q[0];
rz(-1.6630917) q[0];
sx q[0];
rz(-1.2170296) q[0];
x q[1];
rz(-0.92685076) q[2];
sx q[2];
rz(-1.1139718) q[2];
sx q[2];
rz(-0.18075519) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4957461) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(-0.87162019) q[1];
rz(-2.0884573) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(-0.79013463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.084215) q[2];
sx q[2];
rz(-1.4930909) q[2];
sx q[2];
rz(-2.8520612) q[2];
rz(-0.036751898) q[3];
sx q[3];
rz(-0.74220243) q[3];
sx q[3];
rz(-3.1307722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.058379563) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(-1.6311197) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(1.0169792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6769584) q[0];
sx q[0];
rz(-0.58566739) q[0];
sx q[0];
rz(2.2023489) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9404042) q[2];
sx q[2];
rz(-1.5941465) q[2];
sx q[2];
rz(-2.3692998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.444126) q[1];
sx q[1];
rz(-1.0870458) q[1];
sx q[1];
rz(-2.6078348) q[1];
x q[2];
rz(2.6050657) q[3];
sx q[3];
rz(-1.121472) q[3];
sx q[3];
rz(-2.7262053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(0.63465676) q[2];
rz(2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(-2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909661) q[0];
sx q[0];
rz(-2.4996596) q[0];
sx q[0];
rz(0.84386688) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(1.2197781) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0735002) q[0];
sx q[0];
rz(-1.6627899) q[0];
sx q[0];
rz(-1.5324701) q[0];
rz(-0.90791038) q[2];
sx q[2];
rz(-2.6320576) q[2];
sx q[2];
rz(0.31928911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0724832) q[1];
sx q[1];
rz(-1.3246037) q[1];
sx q[1];
rz(-2.2407131) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7759833) q[3];
sx q[3];
rz(-2.725051) q[3];
sx q[3];
rz(-2.6147571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5119778) q[2];
sx q[2];
rz(-2.2045366) q[2];
sx q[2];
rz(0.31663319) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.7946295) q[3];
sx q[3];
rz(-0.75604701) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0379631) q[0];
sx q[0];
rz(-1.0873955) q[0];
sx q[0];
rz(-0.33139247) q[0];
rz(1.9695075) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(-2.7323639) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9650139) q[0];
sx q[0];
rz(-1.2897549) q[0];
sx q[0];
rz(2.2289508) q[0];
x q[1];
rz(-0.51034575) q[2];
sx q[2];
rz(-0.58594698) q[2];
sx q[2];
rz(0.79615359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48620263) q[1];
sx q[1];
rz(-1.7230002) q[1];
sx q[1];
rz(-2.5712719) q[1];
rz(-pi) q[2];
rz(0.51729047) q[3];
sx q[3];
rz(-1.7833157) q[3];
sx q[3];
rz(-1.1503435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1428712) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(-0.55465737) q[2];
rz(-2.5000642) q[3];
sx q[3];
rz(-2.9252164) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759578) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(-1.0558646) q[1];
sx q[1];
rz(-1.5664342) q[1];
sx q[1];
rz(1.6329637) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976534) q[0];
sx q[0];
rz(-0.88730747) q[0];
sx q[0];
rz(2.9336714) q[0];
x q[1];
rz(-0.091392322) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(2.6932655) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29533169) q[1];
sx q[1];
rz(-1.2519072) q[1];
sx q[1];
rz(-1.0237414) q[1];
x q[2];
rz(-1.5040594) q[3];
sx q[3];
rz(-2.4373694) q[3];
sx q[3];
rz(-0.080554068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14279723) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(-2.7049086) q[2];
rz(1.3302594) q[3];
sx q[3];
rz(-0.67970777) q[3];
sx q[3];
rz(-2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10903877) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(0.56525266) q[0];
rz(0.25282192) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.3814829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917282) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(1.6185332) q[0];
rz(-pi) q[1];
rz(1.6818468) q[2];
sx q[2];
rz(-2.0425218) q[2];
sx q[2];
rz(1.9327088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6773721) q[1];
sx q[1];
rz(-0.83253011) q[1];
sx q[1];
rz(-0.99454753) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28678203) q[3];
sx q[3];
rz(-0.87773318) q[3];
sx q[3];
rz(-1.1657438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6910203) q[2];
sx q[2];
rz(-0.82442966) q[2];
sx q[2];
rz(2.5449469) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-0.89533007) q[3];
sx q[3];
rz(0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(0.045885432) q[2];
sx q[2];
rz(-1.3165717) q[2];
sx q[2];
rz(-1.7236621) q[2];
rz(1.4894555) q[3];
sx q[3];
rz(-1.0621536) q[3];
sx q[3];
rz(-2.5928706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
