OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8837589) q[0];
sx q[0];
rz(-0.45816445) q[0];
sx q[0];
rz(-2.5786361) q[0];
rz(1.9986721) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(1.8510173) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3798843) q[1];
sx q[1];
rz(-1.1176846) q[1];
sx q[1];
rz(-0.74700991) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53545973) q[3];
sx q[3];
rz(-2.4389431) q[3];
sx q[3];
rz(-0.040576064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59149867) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(0.47810289) q[2];
rz(-1.6889307) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(-12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(2.1226728) q[0];
rz(1.6628751) q[1];
sx q[1];
rz(-0.61518413) q[1];
sx q[1];
rz(2.5085124) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86242005) q[0];
sx q[0];
rz(-2.3775568) q[0];
sx q[0];
rz(-3.0252302) q[0];
rz(-pi) q[1];
rz(2.1585629) q[2];
sx q[2];
rz(-0.9126185) q[2];
sx q[2];
rz(0.36511974) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96453297) q[1];
sx q[1];
rz(-2.6033083) q[1];
sx q[1];
rz(-0.86175491) q[1];
x q[2];
rz(0.88926104) q[3];
sx q[3];
rz(-1.325843) q[3];
sx q[3];
rz(-0.8196866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35976609) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(3.0241372) q[2];
rz(-0.30101267) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(-1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-2.03405) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(-0.69082469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3120964) q[0];
sx q[0];
rz(-1.4179686) q[0];
sx q[0];
rz(-1.6468847) q[0];
rz(-pi) q[1];
rz(0.76491852) q[2];
sx q[2];
rz(-1.3996482) q[2];
sx q[2];
rz(-2.314832) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9835811) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(1.4278825) q[1];
x q[2];
rz(0.74294528) q[3];
sx q[3];
rz(-0.49693123) q[3];
sx q[3];
rz(-0.039615354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(3.0351191) q[2];
rz(-1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0687662) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(-0.52247125) q[0];
rz(2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(-0.55363384) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5579917) q[0];
sx q[0];
rz(-0.86692536) q[0];
sx q[0];
rz(0.61917275) q[0];
rz(-pi) q[1];
rz(2.373898) q[2];
sx q[2];
rz(-0.97937095) q[2];
sx q[2];
rz(-2.5007345) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1045038) q[1];
sx q[1];
rz(-1.7040841) q[1];
sx q[1];
rz(-2.3549805) q[1];
rz(2.3447498) q[3];
sx q[3];
rz(-1.5505655) q[3];
sx q[3];
rz(-0.86864558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.557495) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(2.6468357) q[2];
rz(0.90302145) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(1.8516699) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49098) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(1.0850798) q[0];
rz(2.1814573) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(0.18403149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2011916) q[0];
sx q[0];
rz(-1.2815676) q[0];
sx q[0];
rz(1.679923) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2940302) q[2];
sx q[2];
rz(-1.2185316) q[2];
sx q[2];
rz(0.89154746) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3369493) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(2.4451392) q[1];
rz(-pi) q[2];
rz(0.22589485) q[3];
sx q[3];
rz(-2.4028006) q[3];
sx q[3];
rz(-1.5238638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90157834) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(-2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(2.9011762) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(2.9930847) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69420544) q[0];
sx q[0];
rz(-2.2930657) q[0];
sx q[0];
rz(-1.7859875) q[0];
x q[1];
rz(-0.22865061) q[2];
sx q[2];
rz(-1.4637404) q[2];
sx q[2];
rz(2.6766237) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2827742) q[1];
sx q[1];
rz(-2.8595279) q[1];
sx q[1];
rz(-1.9726522) q[1];
rz(1.968019) q[3];
sx q[3];
rz(-2.0052611) q[3];
sx q[3];
rz(-1.43515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.285816) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(-2.6532069) q[2];
rz(2.6521902) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(-1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4483036) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(-0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(2.0163527) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624807) q[0];
sx q[0];
rz(-1.1971803) q[0];
sx q[0];
rz(3.0309249) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15525012) q[2];
sx q[2];
rz(-1.2755738) q[2];
sx q[2];
rz(1.0709907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69562558) q[1];
sx q[1];
rz(-1.8672767) q[1];
sx q[1];
rz(-0.28709025) q[1];
x q[2];
rz(-0.66594395) q[3];
sx q[3];
rz(-1.1083372) q[3];
sx q[3];
rz(0.24941128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(2.4970064) q[2];
rz(-1.6623496) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.4054327) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97061625) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(1.53565) q[0];
rz(-1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-1.01064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686417) q[0];
sx q[0];
rz(-2.3803582) q[0];
sx q[0];
rz(-1.3377454) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0881565) q[2];
sx q[2];
rz(-2.5854955) q[2];
sx q[2];
rz(0.14511395) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.752754) q[1];
sx q[1];
rz(-2.8201582) q[1];
sx q[1];
rz(-2.3366117) q[1];
rz(-pi) q[2];
rz(-1.478785) q[3];
sx q[3];
rz(-1.9698471) q[3];
sx q[3];
rz(-1.9012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(-2.4475205) q[2];
rz(-0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(0.49474299) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(0.075597413) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2878802) q[0];
sx q[0];
rz(-1.8822077) q[0];
sx q[0];
rz(-2.1561554) q[0];
rz(-1.8972626) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(-1.2125804) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52121431) q[1];
sx q[1];
rz(-1.4676536) q[1];
sx q[1];
rz(0.68876536) q[1];
rz(-0.62065403) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8081234) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(-1.8010275) q[2];
rz(2.8373485) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2952404) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-0.17679581) q[0];
rz(-1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(2.7005844) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037017578) q[0];
sx q[0];
rz(-1.772176) q[0];
sx q[0];
rz(-1.7778394) q[0];
rz(-2.987864) q[2];
sx q[2];
rz(-1.6445465) q[2];
sx q[2];
rz(1.4027301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8969352) q[1];
sx q[1];
rz(-1.999563) q[1];
sx q[1];
rz(-0.37533092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40481683) q[3];
sx q[3];
rz(-2.4768618) q[3];
sx q[3];
rz(-0.80617314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(2.8397172) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(-3.1148615) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(2.6731861) q[2];
sx q[2];
rz(-0.9006587) q[2];
sx q[2];
rz(-2.547154) q[2];
rz(2.4206352) q[3];
sx q[3];
rz(-0.69098916) q[3];
sx q[3];
rz(0.50222764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
