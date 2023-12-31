OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(3.9305384) q[1];
sx q[1];
rz(9.5126704) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1577507) q[0];
sx q[0];
rz(-1.112126) q[0];
sx q[0];
rz(-2.0660603) q[0];
rz(-pi) q[1];
rz(-0.90859969) q[2];
sx q[2];
rz(-0.63630051) q[2];
sx q[2];
rz(-1.5719906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3808448) q[1];
sx q[1];
rz(-1.5729841) q[1];
sx q[1];
rz(2.6373177) q[1];
rz(-pi) q[2];
rz(-3.1357364) q[3];
sx q[3];
rz(-1.8294888) q[3];
sx q[3];
rz(0.70219016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0007881) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(2.5722356) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5434718) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.233261) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9995995) q[0];
sx q[0];
rz(-0.90786952) q[0];
sx q[0];
rz(-0.87447383) q[0];
rz(-pi) q[1];
rz(2.9821175) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(-0.28123873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.711449) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(1.7769017) q[1];
x q[2];
rz(-1.3319098) q[3];
sx q[3];
rz(-2.2242862) q[3];
sx q[3];
rz(1.7042004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(2.3114752) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8925791) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(-2.3705204) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(3.070389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6547346) q[0];
sx q[0];
rz(-1.3329957) q[0];
sx q[0];
rz(0.20240692) q[0];
rz(-1.0686915) q[2];
sx q[2];
rz(-0.50902589) q[2];
sx q[2];
rz(-0.44391649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0752807) q[1];
sx q[1];
rz(-1.4367141) q[1];
sx q[1];
rz(0.65384298) q[1];
rz(0.86665385) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(-1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(0.52345792) q[2];
rz(-2.8097025) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.7317384) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(-0.82558924) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.4473787) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7839514) q[0];
sx q[0];
rz(-2.8581627) q[0];
sx q[0];
rz(-2.1901603) q[0];
x q[1];
rz(-1.5201735) q[2];
sx q[2];
rz(-0.17855893) q[2];
sx q[2];
rz(-2.6176917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96781603) q[1];
sx q[1];
rz(-0.88138352) q[1];
sx q[1];
rz(-2.9348228) q[1];
x q[2];
rz(2.5656384) q[3];
sx q[3];
rz(-2.7154657) q[3];
sx q[3];
rz(-2.0800316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(-0.14492598) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99575627) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-2.4575535) q[0];
rz(-1.1902635) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(-1.9285944) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66515231) q[0];
sx q[0];
rz(-1.1993009) q[0];
sx q[0];
rz(1.2907589) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91765399) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(-1.7078924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0141543) q[1];
sx q[1];
rz(-1.6886097) q[1];
sx q[1];
rz(0.25728667) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9641987) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(0.28211668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-2.5693494) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(-0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1155788) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(-0.28636006) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(0.58475959) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50305784) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(1.1919341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5445968) q[2];
sx q[2];
rz(-2.1850586) q[2];
sx q[2];
rz(2.7363077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4510182) q[1];
sx q[1];
rz(-1.8585397) q[1];
sx q[1];
rz(0.87177999) q[1];
rz(-pi) q[2];
rz(0.16634511) q[3];
sx q[3];
rz(-1.5208941) q[3];
sx q[3];
rz(0.14453105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(2.1785054) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0528089) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(3.0896297) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(0.60639492) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.279778) q[0];
sx q[0];
rz(-1.8068411) q[0];
sx q[0];
rz(-2.1104913) q[0];
rz(-pi) q[1];
rz(2.0386001) q[2];
sx q[2];
rz(-0.71873795) q[2];
sx q[2];
rz(-1.4008092) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4013306) q[1];
sx q[1];
rz(-1.6832502) q[1];
sx q[1];
rz(-2.3349891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9803195) q[3];
sx q[3];
rz(-1.2425353) q[3];
sx q[3];
rz(1.8144153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(-2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(2.3876277) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(0.22769134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4906625) q[0];
sx q[0];
rz(-0.41914808) q[0];
sx q[0];
rz(1.3952257) q[0];
rz(-0.44616206) q[2];
sx q[2];
rz(-1.9036306) q[2];
sx q[2];
rz(-1.7762426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5880809) q[1];
sx q[1];
rz(-1.7953824) q[1];
sx q[1];
rz(-2.5976546) q[1];
rz(0.065683059) q[3];
sx q[3];
rz(-2.5816397) q[3];
sx q[3];
rz(-0.63510676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(2.3030346) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(-2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39695981) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.9931591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52590695) q[0];
sx q[0];
rz(-1.4988168) q[0];
sx q[0];
rz(0.023948897) q[0];
rz(-pi) q[1];
rz(-2.7788413) q[2];
sx q[2];
rz(-2.0667549) q[2];
sx q[2];
rz(2.9438058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0828404) q[1];
sx q[1];
rz(-2.0420923) q[1];
sx q[1];
rz(-2.6575762) q[1];
rz(2.7896342) q[3];
sx q[3];
rz(-1.6037914) q[3];
sx q[3];
rz(-1.9381423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0917197) q[2];
sx q[2];
rz(-1.0174948) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(-2.9459279) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739928) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(2.4018438) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(-0.69828066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65933278) q[0];
sx q[0];
rz(-1.8730622) q[0];
sx q[0];
rz(-1.4955826) q[0];
rz(1.8928705) q[2];
sx q[2];
rz(-3.0634355) q[2];
sx q[2];
rz(-2.8305778) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0776978) q[1];
sx q[1];
rz(-0.5118891) q[1];
sx q[1];
rz(1.1179395) q[1];
x q[2];
rz(-2.1383168) q[3];
sx q[3];
rz(-1.7756117) q[3];
sx q[3];
rz(1.8703465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(0.84021604) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-0.60224709) q[2];
sx q[2];
rz(-1.8613653) q[2];
sx q[2];
rz(-2.0801434) q[2];
rz(-2.9914231) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
