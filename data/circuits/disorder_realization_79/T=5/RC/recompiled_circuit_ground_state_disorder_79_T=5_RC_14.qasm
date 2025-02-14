OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.76685846) q[0];
sx q[0];
rz(2.2280966) q[0];
sx q[0];
rz(10.812617) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(-1.8228276) q[1];
sx q[1];
rz(0.28611046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494516) q[0];
sx q[0];
rz(-2.0311264) q[0];
sx q[0];
rz(1.597531) q[0];
x q[1];
rz(-0.95979664) q[2];
sx q[2];
rz(-1.6382484) q[2];
sx q[2];
rz(2.5416201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0148791) q[1];
sx q[1];
rz(-2.5225547) q[1];
sx q[1];
rz(-0.45780413) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6872365) q[3];
sx q[3];
rz(-0.44346657) q[3];
sx q[3];
rz(-0.93772353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.046595786) q[2];
sx q[2];
rz(-1.4904212) q[2];
sx q[2];
rz(2.4563834) q[2];
rz(-1.6738711) q[3];
sx q[3];
rz(-2.0592212) q[3];
sx q[3];
rz(0.26836029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14107038) q[0];
sx q[0];
rz(-1.594161) q[0];
sx q[0];
rz(-0.94456124) q[0];
rz(-3.0682849) q[1];
sx q[1];
rz(-2.3088375) q[1];
sx q[1];
rz(-1.8836969) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8735786) q[0];
sx q[0];
rz(-2.0248381) q[0];
sx q[0];
rz(-2.193903) q[0];
x q[1];
rz(1.3665359) q[2];
sx q[2];
rz(-2.0957207) q[2];
sx q[2];
rz(-1.7634942) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98698178) q[1];
sx q[1];
rz(-1.4263337) q[1];
sx q[1];
rz(-2.7491436) q[1];
x q[2];
rz(-1.5242833) q[3];
sx q[3];
rz(-1.6078954) q[3];
sx q[3];
rz(2.056289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32626095) q[2];
sx q[2];
rz(-1.7504642) q[2];
sx q[2];
rz(1.0247256) q[2];
rz(2.4552086) q[3];
sx q[3];
rz(-0.73203433) q[3];
sx q[3];
rz(2.2288403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82655418) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(2.3231373) q[0];
rz(-1.2089027) q[1];
sx q[1];
rz(-1.6345638) q[1];
sx q[1];
rz(-0.85982972) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6085109) q[0];
sx q[0];
rz(-1.2405335) q[0];
sx q[0];
rz(0.94655605) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0947862) q[2];
sx q[2];
rz(-0.90887672) q[2];
sx q[2];
rz(-2.2713585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26759717) q[1];
sx q[1];
rz(-1.1073879) q[1];
sx q[1];
rz(-0.74041669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1016261) q[3];
sx q[3];
rz(-0.2519603) q[3];
sx q[3];
rz(-2.2445238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7439421) q[2];
sx q[2];
rz(-0.33471477) q[2];
sx q[2];
rz(-3.0261377) q[2];
rz(-2.7502821) q[3];
sx q[3];
rz(-2.1152928) q[3];
sx q[3];
rz(1.9976043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35599577) q[0];
sx q[0];
rz(-1.1859256) q[0];
sx q[0];
rz(-0.24056973) q[0];
rz(1.3661512) q[1];
sx q[1];
rz(-1.4931449) q[1];
sx q[1];
rz(-1.2496525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24883305) q[0];
sx q[0];
rz(-2.2035778) q[0];
sx q[0];
rz(2.8248293) q[0];
rz(-pi) q[1];
rz(-2.0555105) q[2];
sx q[2];
rz(-2.4668985) q[2];
sx q[2];
rz(0.36568588) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74339691) q[1];
sx q[1];
rz(-0.92602387) q[1];
sx q[1];
rz(-1.2750687) q[1];
x q[2];
rz(-0.68393882) q[3];
sx q[3];
rz(-1.9762282) q[3];
sx q[3];
rz(2.8059174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8852641) q[2];
sx q[2];
rz(-2.5958131) q[2];
sx q[2];
rz(-1.7735927) q[2];
rz(2.8374425) q[3];
sx q[3];
rz(-1.5757898) q[3];
sx q[3];
rz(2.4450541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.91365415) q[0];
sx q[0];
rz(-0.21622394) q[0];
sx q[0];
rz(0.53713334) q[0];
rz(2.3917603) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(2.4044663) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43057399) q[0];
sx q[0];
rz(-1.6665742) q[0];
sx q[0];
rz(-1.6349313) q[0];
x q[1];
rz(-1.7456878) q[2];
sx q[2];
rz(-1.6564257) q[2];
sx q[2];
rz(-2.0599287) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3143718) q[1];
sx q[1];
rz(-1.6984358) q[1];
sx q[1];
rz(2.8807544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6629822) q[3];
sx q[3];
rz(-1.2488197) q[3];
sx q[3];
rz(-0.76065242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4445112) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(0.14661655) q[2];
rz(-1.4491436) q[3];
sx q[3];
rz(-1.2796947) q[3];
sx q[3];
rz(-0.52391887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61831063) q[0];
sx q[0];
rz(-1.016541) q[0];
sx q[0];
rz(-1.7200394) q[0];
rz(-1.8824185) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(-2.7377985) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84995364) q[0];
sx q[0];
rz(-1.2899512) q[0];
sx q[0];
rz(2.3333059) q[0];
rz(-pi) q[1];
rz(2.8844599) q[2];
sx q[2];
rz(-1.3916573) q[2];
sx q[2];
rz(0.15878294) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6359977) q[1];
sx q[1];
rz(-2.2790716) q[1];
sx q[1];
rz(-0.8188894) q[1];
rz(0.29098068) q[3];
sx q[3];
rz(-2.432219) q[3];
sx q[3];
rz(0.13252276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28973618) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(2.9750321) q[2];
rz(-1.5038495) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(-0.09854266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2378167) q[0];
sx q[0];
rz(-2.7523968) q[0];
sx q[0];
rz(0.87752262) q[0];
rz(2.1794043) q[1];
sx q[1];
rz(-2.018237) q[1];
sx q[1];
rz(-0.35433623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0838881) q[0];
sx q[0];
rz(-1.5530927) q[0];
sx q[0];
rz(2.903865) q[0];
rz(-1.7155859) q[2];
sx q[2];
rz(-1.5476523) q[2];
sx q[2];
rz(2.1320081) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52088941) q[1];
sx q[1];
rz(-2.1701309) q[1];
sx q[1];
rz(0.66461892) q[1];
rz(0.13451666) q[3];
sx q[3];
rz(-2.627774) q[3];
sx q[3];
rz(1.5863933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3366036) q[2];
sx q[2];
rz(-0.73770928) q[2];
sx q[2];
rz(-2.0965516) q[2];
rz(2.7392144) q[3];
sx q[3];
rz(-1.9741524) q[3];
sx q[3];
rz(-1.4761402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43110111) q[0];
sx q[0];
rz(-0.11678188) q[0];
sx q[0];
rz(-0.85103273) q[0];
rz(-2.7878413) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(-2.2869349) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8847722) q[0];
sx q[0];
rz(-1.2899634) q[0];
sx q[0];
rz(1.3815085) q[0];
x q[1];
rz(1.7206011) q[2];
sx q[2];
rz(-1.8179107) q[2];
sx q[2];
rz(0.30943662) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51616299) q[1];
sx q[1];
rz(-0.99402797) q[1];
sx q[1];
rz(1.7712084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1885957) q[3];
sx q[3];
rz(-2.7321987) q[3];
sx q[3];
rz(-2.2379217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9510368) q[2];
sx q[2];
rz(-0.64728105) q[2];
sx q[2];
rz(-0.62038842) q[2];
rz(-2.9151211) q[3];
sx q[3];
rz(-1.3969996) q[3];
sx q[3];
rz(-2.5605719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1802375) q[0];
sx q[0];
rz(-1.9397475) q[0];
sx q[0];
rz(-3.1264547) q[0];
rz(-2.119428) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(-1.9415564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43754345) q[0];
sx q[0];
rz(-1.3945197) q[0];
sx q[0];
rz(2.2793819) q[0];
rz(-pi) q[1];
rz(1.6300417) q[2];
sx q[2];
rz(-0.97840913) q[2];
sx q[2];
rz(-2.5424344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2150458) q[1];
sx q[1];
rz(-2.6290942) q[1];
sx q[1];
rz(1.6319176) q[1];
rz(-2.3288378) q[3];
sx q[3];
rz(-0.65633196) q[3];
sx q[3];
rz(-1.1876653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6156561) q[2];
sx q[2];
rz(-0.91292149) q[2];
sx q[2];
rz(-0.23942648) q[2];
rz(-3.0827403) q[3];
sx q[3];
rz(-0.81164304) q[3];
sx q[3];
rz(-2.8656901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(-1.9489768) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(2.1666727) q[0];
rz(0.67543593) q[1];
sx q[1];
rz(-0.91921872) q[1];
sx q[1];
rz(-0.98178896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73327195) q[0];
sx q[0];
rz(-1.3106723) q[0];
sx q[0];
rz(0.97113804) q[0];
x q[1];
rz(1.6390332) q[2];
sx q[2];
rz(-1.1600798) q[2];
sx q[2];
rz(-1.0194743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60335873) q[1];
sx q[1];
rz(-2.3673944) q[1];
sx q[1];
rz(1.1959568) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4922139) q[3];
sx q[3];
rz(-2.693697) q[3];
sx q[3];
rz(-0.092002951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56741095) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(-1.9166454) q[2];
rz(2.775906) q[3];
sx q[3];
rz(-2.1920125) q[3];
sx q[3];
rz(1.1398116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0443307) q[0];
sx q[0];
rz(-1.8987569) q[0];
sx q[0];
rz(1.4415997) q[0];
rz(0.080009566) q[1];
sx q[1];
rz(-1.3792104) q[1];
sx q[1];
rz(3.1389799) q[1];
rz(0.44543191) q[2];
sx q[2];
rz(-1.6966698) q[2];
sx q[2];
rz(2.0296283) q[2];
rz(-1.5936218) q[3];
sx q[3];
rz(-1.1183839) q[3];
sx q[3];
rz(-0.090261685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
