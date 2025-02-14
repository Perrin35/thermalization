OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(1.3774435) q[0];
sx q[0];
rz(5.8537771) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(-1.5996999) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5326951) q[0];
sx q[0];
rz(-2.7809394) q[0];
sx q[0];
rz(2.2171573) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6060271) q[2];
sx q[2];
rz(-1.7144979) q[2];
sx q[2];
rz(1.8409539) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0860025) q[1];
sx q[1];
rz(-1.1430972) q[1];
sx q[1];
rz(1.5375288) q[1];
rz(-pi) q[2];
rz(1.8640932) q[3];
sx q[3];
rz(-1.9175005) q[3];
sx q[3];
rz(1.6270005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30449197) q[2];
sx q[2];
rz(-1.9383177) q[2];
sx q[2];
rz(-0.692918) q[2];
rz(-0.20017008) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(1.1339124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3226586) q[0];
sx q[0];
rz(-0.19911961) q[0];
sx q[0];
rz(-2.0767427) q[0];
rz(0.62243593) q[1];
sx q[1];
rz(-1.348) q[1];
sx q[1];
rz(-0.49076864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9483295) q[0];
sx q[0];
rz(-1.5800467) q[0];
sx q[0];
rz(3.1178283) q[0];
x q[1];
rz(-2.6726637) q[2];
sx q[2];
rz(-0.46655077) q[2];
sx q[2];
rz(-2.5836437) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.307617) q[1];
sx q[1];
rz(-2.3217891) q[1];
sx q[1];
rz(-1.7012672) q[1];
rz(-pi) q[2];
rz(1.2964046) q[3];
sx q[3];
rz(-0.17582045) q[3];
sx q[3];
rz(0.71064132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48851594) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(-2.901279) q[2];
rz(-2.6416685) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(-1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2568473) q[0];
sx q[0];
rz(-2.186543) q[0];
sx q[0];
rz(-2.1599059) q[0];
rz(2.6495972) q[1];
sx q[1];
rz(-1.0849846) q[1];
sx q[1];
rz(-1.5031776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5444191) q[0];
sx q[0];
rz(-1.0962631) q[0];
sx q[0];
rz(-2.3719792) q[0];
rz(-pi) q[1];
rz(-3.0122224) q[2];
sx q[2];
rz(-2.7384842) q[2];
sx q[2];
rz(1.0327686) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7042289) q[1];
sx q[1];
rz(-1.2384602) q[1];
sx q[1];
rz(1.9802914) q[1];
rz(-1.668456) q[3];
sx q[3];
rz(-1.5499695) q[3];
sx q[3];
rz(1.5690924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65695277) q[2];
sx q[2];
rz(-2.6159365) q[2];
sx q[2];
rz(3.0204115) q[2];
rz(-1.0080522) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(1.0813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0722395) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(0.51373154) q[0];
rz(0.95798245) q[1];
sx q[1];
rz(-1.1424516) q[1];
sx q[1];
rz(-1.6987919) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49794562) q[0];
sx q[0];
rz(-1.7352603) q[0];
sx q[0];
rz(2.7269159) q[0];
x q[1];
rz(1.3495693) q[2];
sx q[2];
rz(-1.0509247) q[2];
sx q[2];
rz(-1.9398361) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5123547) q[1];
sx q[1];
rz(-1.4886453) q[1];
sx q[1];
rz(-2.1869529) q[1];
rz(-0.6616627) q[3];
sx q[3];
rz(-1.6473624) q[3];
sx q[3];
rz(-0.14684248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5865667) q[2];
sx q[2];
rz(-0.96379605) q[2];
sx q[2];
rz(-1.8761934) q[2];
rz(-1.1749367) q[3];
sx q[3];
rz(-0.73052162) q[3];
sx q[3];
rz(1.9780212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2555399) q[0];
sx q[0];
rz(-0.20296725) q[0];
sx q[0];
rz(-1.8170005) q[0];
rz(2.1728204) q[1];
sx q[1];
rz(-1.6561534) q[1];
sx q[1];
rz(-1.0383777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7234112) q[0];
sx q[0];
rz(-1.9514284) q[0];
sx q[0];
rz(1.3971973) q[0];
x q[1];
rz(-0.81401396) q[2];
sx q[2];
rz(-0.90384877) q[2];
sx q[2];
rz(0.89165724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5165625) q[1];
sx q[1];
rz(-1.5296401) q[1];
sx q[1];
rz(-2.7709318) q[1];
x q[2];
rz(2.5143751) q[3];
sx q[3];
rz(-0.76290799) q[3];
sx q[3];
rz(1.7064394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5021299) q[2];
sx q[2];
rz(-0.30439964) q[2];
sx q[2];
rz(2.2501865) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.6311878) q[3];
sx q[3];
rz(0.6663028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0797743) q[0];
sx q[0];
rz(-0.49538716) q[0];
sx q[0];
rz(2.1451982) q[0];
rz(2.2415316) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(-0.17734227) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3729438) q[0];
sx q[0];
rz(-1.5705612) q[0];
sx q[0];
rz(1.7223486) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93219535) q[2];
sx q[2];
rz(-1.5538486) q[2];
sx q[2];
rz(-2.058896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4882433) q[1];
sx q[1];
rz(-2.0470139) q[1];
sx q[1];
rz(0.32250065) q[1];
rz(-pi) q[2];
rz(-3.0445339) q[3];
sx q[3];
rz(-2.0776111) q[3];
sx q[3];
rz(2.7926796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29048723) q[2];
sx q[2];
rz(-1.1815716) q[2];
sx q[2];
rz(-2.8372852) q[2];
rz(-2.815222) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(2.2296026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07311634) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(-0.9325183) q[0];
rz(-3.0724691) q[1];
sx q[1];
rz(-0.2176452) q[1];
sx q[1];
rz(1.0401915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067834082) q[0];
sx q[0];
rz(-2.375575) q[0];
sx q[0];
rz(3.104408) q[0];
rz(-pi) q[1];
rz(-0.97904737) q[2];
sx q[2];
rz(-0.30610105) q[2];
sx q[2];
rz(0.09397587) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8833163) q[1];
sx q[1];
rz(-2.3280716) q[1];
sx q[1];
rz(-1.798428) q[1];
rz(-pi) q[2];
rz(-1.0024985) q[3];
sx q[3];
rz(-2.1868984) q[3];
sx q[3];
rz(1.7006066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13504623) q[2];
sx q[2];
rz(-2.4032205) q[2];
sx q[2];
rz(0.73299232) q[2];
rz(-2.0586355) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(-2.1347031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5873544) q[0];
sx q[0];
rz(-3.0576958) q[0];
sx q[0];
rz(-0.49302897) q[0];
rz(-2.9250277) q[1];
sx q[1];
rz(-2.000587) q[1];
sx q[1];
rz(-2.1176178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42070779) q[0];
sx q[0];
rz(-2.5455089) q[0];
sx q[0];
rz(-3.1203624) q[0];
rz(-2.8220213) q[2];
sx q[2];
rz(-0.34863483) q[2];
sx q[2];
rz(-0.92952585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7143284) q[1];
sx q[1];
rz(-1.9275394) q[1];
sx q[1];
rz(-1.4209695) q[1];
x q[2];
rz(1.1795565) q[3];
sx q[3];
rz(-0.68022281) q[3];
sx q[3];
rz(-2.8381951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0563125) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(2.4658266) q[2];
rz(2.7285649) q[3];
sx q[3];
rz(-1.690381) q[3];
sx q[3];
rz(-2.5954424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6290879) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(3.107048) q[0];
rz(-0.054556219) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(-2.3202855) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9050582) q[0];
sx q[0];
rz(-0.58249677) q[0];
sx q[0];
rz(0.18180099) q[0];
rz(-0.088772687) q[2];
sx q[2];
rz(-1.5800416) q[2];
sx q[2];
rz(-0.68999664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37335884) q[1];
sx q[1];
rz(-1.4211417) q[1];
sx q[1];
rz(1.6878769) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48252941) q[3];
sx q[3];
rz(-0.97597968) q[3];
sx q[3];
rz(-1.1728668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7811232) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(-0.13295573) q[2];
rz(0.41108701) q[3];
sx q[3];
rz(-1.5186331) q[3];
sx q[3];
rz(-1.0857922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046831176) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(-1.8827615) q[0];
rz(1.5065441) q[1];
sx q[1];
rz(-0.656744) q[1];
sx q[1];
rz(-0.81319317) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55594873) q[0];
sx q[0];
rz(-0.085594479) q[0];
sx q[0];
rz(2.7187347) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4237464) q[2];
sx q[2];
rz(-0.67710256) q[2];
sx q[2];
rz(2.8092217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.606914) q[1];
sx q[1];
rz(-2.3707304) q[1];
sx q[1];
rz(0.67463629) q[1];
rz(0.2701668) q[3];
sx q[3];
rz(-2.674683) q[3];
sx q[3];
rz(-1.2845357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4930111) q[2];
sx q[2];
rz(-2.0040671) q[2];
sx q[2];
rz(-2.0900334) q[2];
rz(2.1227396) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(0.65783182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9857585) q[0];
sx q[0];
rz(-2.4875165) q[0];
sx q[0];
rz(-2.2610337) q[0];
rz(-1.3195994) q[1];
sx q[1];
rz(-1.9965912) q[1];
sx q[1];
rz(0.051864787) q[1];
rz(-0.32991275) q[2];
sx q[2];
rz(-2.920426) q[2];
sx q[2];
rz(-2.5830808) q[2];
rz(-1.6307835) q[3];
sx q[3];
rz(-1.7675478) q[3];
sx q[3];
rz(0.02789733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
