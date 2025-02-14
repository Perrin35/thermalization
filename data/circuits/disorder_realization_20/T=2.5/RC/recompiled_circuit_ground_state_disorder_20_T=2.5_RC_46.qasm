OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29937509) q[0];
sx q[0];
rz(-2.8111281) q[0];
sx q[0];
rz(2.0781031) q[0];
rz(-0.039634135) q[1];
sx q[1];
rz(-0.57365817) q[1];
sx q[1];
rz(2.0613476) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325057) q[0];
sx q[0];
rz(-1.5686146) q[0];
sx q[0];
rz(-2.0337142) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6559907) q[2];
sx q[2];
rz(-2.9692215) q[2];
sx q[2];
rz(-1.048798) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0277522) q[1];
sx q[1];
rz(-0.58958399) q[1];
sx q[1];
rz(-2.6623308) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64452602) q[3];
sx q[3];
rz(-1.9686507) q[3];
sx q[3];
rz(0.81105876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0237191) q[2];
sx q[2];
rz(-1.7261852) q[2];
sx q[2];
rz(2.7685557) q[2];
rz(-0.55752623) q[3];
sx q[3];
rz(-1.6256465) q[3];
sx q[3];
rz(-2.663234) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28906223) q[0];
sx q[0];
rz(-1.7076778) q[0];
sx q[0];
rz(2.1321994) q[0];
rz(2.8139662) q[1];
sx q[1];
rz(-2.8150924) q[1];
sx q[1];
rz(-2.0351298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94325667) q[0];
sx q[0];
rz(-3.0064764) q[0];
sx q[0];
rz(0.19395239) q[0];
rz(1.5861139) q[2];
sx q[2];
rz(-2.2557182) q[2];
sx q[2];
rz(-0.1372125) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4047186) q[1];
sx q[1];
rz(-2.0366743) q[1];
sx q[1];
rz(2.2437255) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7296612) q[3];
sx q[3];
rz(-1.3917482) q[3];
sx q[3];
rz(-1.826394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24836765) q[2];
sx q[2];
rz(-1.1838341) q[2];
sx q[2];
rz(0.66967213) q[2];
rz(1.9855965) q[3];
sx q[3];
rz(-2.7892734) q[3];
sx q[3];
rz(-2.7922351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356693) q[0];
sx q[0];
rz(-2.7140129) q[0];
sx q[0];
rz(-1.8100354) q[0];
rz(1.9412899) q[1];
sx q[1];
rz(-0.2526865) q[1];
sx q[1];
rz(2.4694064) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12608317) q[0];
sx q[0];
rz(-2.042217) q[0];
sx q[0];
rz(-1.7065364) q[0];
x q[1];
rz(1.9750017) q[2];
sx q[2];
rz(-1.4160755) q[2];
sx q[2];
rz(-0.2153309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7461639) q[1];
sx q[1];
rz(-0.47074598) q[1];
sx q[1];
rz(2.7676299) q[1];
rz(-pi) q[2];
x q[2];
rz(2.658031) q[3];
sx q[3];
rz(-3.0322595) q[3];
sx q[3];
rz(3.1121174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65472284) q[2];
sx q[2];
rz(-1.9630311) q[2];
sx q[2];
rz(0.85050026) q[2];
rz(-0.9507829) q[3];
sx q[3];
rz(-2.6362004) q[3];
sx q[3];
rz(-2.2397485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.013407) q[0];
sx q[0];
rz(-0.81939092) q[0];
sx q[0];
rz(2.4017781) q[0];
rz(-0.20135227) q[1];
sx q[1];
rz(-1.9722152) q[1];
sx q[1];
rz(-0.22612017) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0126919) q[0];
sx q[0];
rz(-1.6058815) q[0];
sx q[0];
rz(-3.1265774) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9490384) q[2];
sx q[2];
rz(-1.8762445) q[2];
sx q[2];
rz(-2.7482207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12783957) q[1];
sx q[1];
rz(-2.3495449) q[1];
sx q[1];
rz(2.386078) q[1];
rz(2.8299895) q[3];
sx q[3];
rz(-2.7392276) q[3];
sx q[3];
rz(-0.14801797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7390593) q[2];
sx q[2];
rz(-2.2465574) q[2];
sx q[2];
rz(2.1441937) q[2];
rz(-2.7282257) q[3];
sx q[3];
rz(-1.2099384) q[3];
sx q[3];
rz(1.0443784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.0053552) q[0];
sx q[0];
rz(-0.69664609) q[0];
sx q[0];
rz(0.053939017) q[0];
rz(-0.18221642) q[1];
sx q[1];
rz(-2.2257664) q[1];
sx q[1];
rz(-0.30311662) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5526472) q[0];
sx q[0];
rz(-0.56845462) q[0];
sx q[0];
rz(-2.2442152) q[0];
rz(-pi) q[1];
rz(-0.19333573) q[2];
sx q[2];
rz(-0.92776644) q[2];
sx q[2];
rz(-1.0372127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7531277) q[1];
sx q[1];
rz(-1.376994) q[1];
sx q[1];
rz(1.0419921) q[1];
rz(-3.1040499) q[3];
sx q[3];
rz(-0.75858603) q[3];
sx q[3];
rz(-0.19877082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4431241) q[2];
sx q[2];
rz(-1.2681401) q[2];
sx q[2];
rz(2.7847086) q[2];
rz(-2.3667864) q[3];
sx q[3];
rz(-1.509343) q[3];
sx q[3];
rz(-0.43865144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2543432) q[0];
sx q[0];
rz(-1.9312504) q[0];
sx q[0];
rz(-0.78897011) q[0];
rz(-1.3428768) q[1];
sx q[1];
rz(-1.7306381) q[1];
sx q[1];
rz(-0.79997396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3244902) q[0];
sx q[0];
rz(-0.68590012) q[0];
sx q[0];
rz(1.5490129) q[0];
rz(-0.87514295) q[2];
sx q[2];
rz(-1.4261386) q[2];
sx q[2];
rz(-0.62334594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8084849) q[1];
sx q[1];
rz(-1.4623545) q[1];
sx q[1];
rz(-0.37745775) q[1];
rz(-pi) q[2];
rz(2.7423032) q[3];
sx q[3];
rz(-1.1762816) q[3];
sx q[3];
rz(2.7189915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0903025) q[2];
sx q[2];
rz(-2.3919969) q[2];
sx q[2];
rz(-1.7898111) q[2];
rz(-0.46818647) q[3];
sx q[3];
rz(-1.0951833) q[3];
sx q[3];
rz(0.93402544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4139597) q[0];
sx q[0];
rz(-0.41828823) q[0];
sx q[0];
rz(-3.1391414) q[0];
rz(-0.0062423627) q[1];
sx q[1];
rz(-2.7793482) q[1];
sx q[1];
rz(-2.7856766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.064156) q[0];
sx q[0];
rz(-1.9165123) q[0];
sx q[0];
rz(-0.92597175) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8121901) q[2];
sx q[2];
rz(-0.74879941) q[2];
sx q[2];
rz(1.8210653) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.351928) q[1];
sx q[1];
rz(-1.8311336) q[1];
sx q[1];
rz(0.89485333) q[1];
x q[2];
rz(0.53419279) q[3];
sx q[3];
rz(-1.4823128) q[3];
sx q[3];
rz(-0.26966394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1459085) q[2];
sx q[2];
rz(-1.9084787) q[2];
sx q[2];
rz(-1.0710867) q[2];
rz(0.60761014) q[3];
sx q[3];
rz(-1.1039609) q[3];
sx q[3];
rz(1.5266017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0869658) q[0];
sx q[0];
rz(-2.2038951) q[0];
sx q[0];
rz(-1.1338393) q[0];
rz(-2.0896104) q[1];
sx q[1];
rz(-1.4583505) q[1];
sx q[1];
rz(2.4005311) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4008557) q[0];
sx q[0];
rz(-1.0154402) q[0];
sx q[0];
rz(-1.8779138) q[0];
rz(-0.5036139) q[2];
sx q[2];
rz(-2.1237649) q[2];
sx q[2];
rz(1.210975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3723124) q[1];
sx q[1];
rz(-2.0959294) q[1];
sx q[1];
rz(-1.5430081) q[1];
rz(-0.94101936) q[3];
sx q[3];
rz(-1.2411323) q[3];
sx q[3];
rz(2.8235112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97003254) q[2];
sx q[2];
rz(-0.67983183) q[2];
sx q[2];
rz(2.6762834) q[2];
rz(-2.4252637) q[3];
sx q[3];
rz(-1.6445487) q[3];
sx q[3];
rz(1.1819476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5301836) q[0];
sx q[0];
rz(-0.94785988) q[0];
sx q[0];
rz(1.1780257) q[0];
rz(1.0097965) q[1];
sx q[1];
rz(-1.806908) q[1];
sx q[1];
rz(0.10733265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7032181) q[0];
sx q[0];
rz(-2.1461663) q[0];
sx q[0];
rz(-2.1557425) q[0];
rz(-1.7221409) q[2];
sx q[2];
rz(-0.79337315) q[2];
sx q[2];
rz(1.8167439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33683981) q[1];
sx q[1];
rz(-0.64267413) q[1];
sx q[1];
rz(-0.39037946) q[1];
x q[2];
rz(2.168591) q[3];
sx q[3];
rz(-1.8174371) q[3];
sx q[3];
rz(-1.1688978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44059077) q[2];
sx q[2];
rz(-0.89459449) q[2];
sx q[2];
rz(0.44006285) q[2];
rz(-0.20982404) q[3];
sx q[3];
rz(-1.194229) q[3];
sx q[3];
rz(-2.8275209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1485727) q[0];
sx q[0];
rz(-0.58605376) q[0];
sx q[0];
rz(1.7589737) q[0];
rz(2.4644409) q[1];
sx q[1];
rz(-1.4644198) q[1];
sx q[1];
rz(2.009353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099927946) q[0];
sx q[0];
rz(-1.7986725) q[0];
sx q[0];
rz(-2.9128338) q[0];
rz(-1.3034362) q[2];
sx q[2];
rz(-2.3472381) q[2];
sx q[2];
rz(1.8092138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5508278) q[1];
sx q[1];
rz(-0.71749291) q[1];
sx q[1];
rz(-0.036542459) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4200045) q[3];
sx q[3];
rz(-0.81755985) q[3];
sx q[3];
rz(-2.8374654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18572346) q[2];
sx q[2];
rz(-2.393554) q[2];
sx q[2];
rz(1.424074) q[2];
rz(2.2636223) q[3];
sx q[3];
rz(-2.9214171) q[3];
sx q[3];
rz(-0.018208114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743185) q[0];
sx q[0];
rz(-1.6813288) q[0];
sx q[0];
rz(1.9718476) q[0];
rz(-2.7727903) q[1];
sx q[1];
rz(-1.7316876) q[1];
sx q[1];
rz(0.015451886) q[1];
rz(-2.5813492) q[2];
sx q[2];
rz(-0.50749736) q[2];
sx q[2];
rz(0.56847405) q[2];
rz(0.54184171) q[3];
sx q[3];
rz(-2.4892157) q[3];
sx q[3];
rz(2.6996725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
