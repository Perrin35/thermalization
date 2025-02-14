OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12914817) q[0];
sx q[0];
rz(-1.669786) q[0];
sx q[0];
rz(0.9653402) q[0];
rz(0.020429285) q[1];
sx q[1];
rz(-1.483622) q[1];
sx q[1];
rz(-1.3774011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55104461) q[0];
sx q[0];
rz(-0.88860308) q[0];
sx q[0];
rz(2.5709573) q[0];
x q[1];
rz(0.032261511) q[2];
sx q[2];
rz(-0.92515495) q[2];
sx q[2];
rz(1.0021068) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5038576) q[1];
sx q[1];
rz(-1.4252932) q[1];
sx q[1];
rz(2.9521431) q[1];
x q[2];
rz(3.1226452) q[3];
sx q[3];
rz(-0.80120443) q[3];
sx q[3];
rz(-2.612118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2237902) q[2];
sx q[2];
rz(-1.6518355) q[2];
sx q[2];
rz(1.5000783) q[2];
rz(2.3598059) q[3];
sx q[3];
rz(-0.89038554) q[3];
sx q[3];
rz(2.3894943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016945275) q[0];
sx q[0];
rz(-2.3638159) q[0];
sx q[0];
rz(-2.1710904) q[0];
rz(1.8151201) q[1];
sx q[1];
rz(-1.3423723) q[1];
sx q[1];
rz(0.91189799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058033179) q[0];
sx q[0];
rz(-0.61249295) q[0];
sx q[0];
rz(-2.8508194) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17120338) q[2];
sx q[2];
rz(-1.7641626) q[2];
sx q[2];
rz(-0.69583508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9262097) q[1];
sx q[1];
rz(-1.5299986) q[1];
sx q[1];
rz(0.87017362) q[1];
x q[2];
rz(-1.5390057) q[3];
sx q[3];
rz(-2.6864751) q[3];
sx q[3];
rz(-1.8546752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9087002) q[2];
sx q[2];
rz(-2.0286655) q[2];
sx q[2];
rz(-2.1523037) q[2];
rz(-0.42012897) q[3];
sx q[3];
rz(-1.0003072) q[3];
sx q[3];
rz(-2.4205128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99783889) q[0];
sx q[0];
rz(-1.9793352) q[0];
sx q[0];
rz(1.0379399) q[0];
rz(-2.2833917) q[1];
sx q[1];
rz(-0.52640262) q[1];
sx q[1];
rz(0.16955489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114544) q[0];
sx q[0];
rz(-1.4272604) q[0];
sx q[0];
rz(3.1375196) q[0];
x q[1];
rz(2.2409147) q[2];
sx q[2];
rz(-0.49374106) q[2];
sx q[2];
rz(2.7639219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.922561) q[1];
sx q[1];
rz(-1.2777849) q[1];
sx q[1];
rz(-2.0907563) q[1];
rz(-1.9603332) q[3];
sx q[3];
rz(-2.7523605) q[3];
sx q[3];
rz(1.8336982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9618591) q[2];
sx q[2];
rz(-0.51859513) q[2];
sx q[2];
rz(-0.40599424) q[2];
rz(2.3640442) q[3];
sx q[3];
rz(-0.33027875) q[3];
sx q[3];
rz(-2.3269261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4985713) q[0];
sx q[0];
rz(-1.850147) q[0];
sx q[0];
rz(-2.2136069) q[0];
rz(3.0827177) q[1];
sx q[1];
rz(-1.4480271) q[1];
sx q[1];
rz(-1.4307129) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3979523) q[0];
sx q[0];
rz(-0.40983202) q[0];
sx q[0];
rz(-1.3751956) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7785382) q[2];
sx q[2];
rz(-1.7573234) q[2];
sx q[2];
rz(-0.26939738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3763172) q[1];
sx q[1];
rz(-2.5936454) q[1];
sx q[1];
rz(-1.1330695) q[1];
x q[2];
rz(2.1998422) q[3];
sx q[3];
rz(-1.3247196) q[3];
sx q[3];
rz(-1.4914644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6835988) q[2];
sx q[2];
rz(-1.6432089) q[2];
sx q[2];
rz(-1.9963473) q[2];
rz(-2.2942719) q[3];
sx q[3];
rz(-1.3342131) q[3];
sx q[3];
rz(1.639521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97750807) q[0];
sx q[0];
rz(-1.9590398) q[0];
sx q[0];
rz(-2.7567647) q[0];
rz(-2.341914) q[1];
sx q[1];
rz(-2.622602) q[1];
sx q[1];
rz(1.5934561) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6131957) q[0];
sx q[0];
rz(-1.3269182) q[0];
sx q[0];
rz(2.8923678) q[0];
x q[1];
rz(-1.1742915) q[2];
sx q[2];
rz(-1.7616211) q[2];
sx q[2];
rz(-0.79939524) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56866327) q[1];
sx q[1];
rz(-1.8844114) q[1];
sx q[1];
rz(0.36797087) q[1];
rz(2.1372651) q[3];
sx q[3];
rz(-0.62823717) q[3];
sx q[3];
rz(2.7816176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3758214) q[2];
sx q[2];
rz(-2.4479726) q[2];
sx q[2];
rz(-0.080246298) q[2];
rz(-0.56096983) q[3];
sx q[3];
rz(-1.4813981) q[3];
sx q[3];
rz(-2.5188353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65492594) q[0];
sx q[0];
rz(-0.66513649) q[0];
sx q[0];
rz(1.2605865) q[0];
rz(-0.27443019) q[1];
sx q[1];
rz(-1.471289) q[1];
sx q[1];
rz(2.4226277) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711766) q[0];
sx q[0];
rz(-2.377802) q[0];
sx q[0];
rz(2.0877354) q[0];
x q[1];
rz(-0.8242714) q[2];
sx q[2];
rz(-1.5780515) q[2];
sx q[2];
rz(1.5430918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8985892) q[1];
sx q[1];
rz(-0.92783992) q[1];
sx q[1];
rz(-1.7220294) q[1];
x q[2];
rz(0.218504) q[3];
sx q[3];
rz(-1.9338622) q[3];
sx q[3];
rz(2.4485112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4869953) q[2];
sx q[2];
rz(-1.2664814) q[2];
sx q[2];
rz(0.60301644) q[2];
rz(1.6675789) q[3];
sx q[3];
rz(-2.1173729) q[3];
sx q[3];
rz(-2.730864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5176158) q[0];
sx q[0];
rz(-1.5331974) q[0];
sx q[0];
rz(-0.18727592) q[0];
rz(-2.1693443) q[1];
sx q[1];
rz(-0.15521237) q[1];
sx q[1];
rz(2.7211199) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3492476) q[0];
sx q[0];
rz(-0.92206832) q[0];
sx q[0];
rz(-2.1029841) q[0];
rz(-pi) q[1];
rz(-0.98780178) q[2];
sx q[2];
rz(-0.62623411) q[2];
sx q[2];
rz(0.69863182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9686465) q[1];
sx q[1];
rz(-1.1889646) q[1];
sx q[1];
rz(0.29995014) q[1];
rz(-2.1450348) q[3];
sx q[3];
rz(-1.7645932) q[3];
sx q[3];
rz(-0.36239788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.04574) q[2];
sx q[2];
rz(-1.1008215) q[2];
sx q[2];
rz(-2.8908758) q[2];
rz(-1.8935253) q[3];
sx q[3];
rz(-1.7496795) q[3];
sx q[3];
rz(-2.5417476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8641758) q[0];
sx q[0];
rz(-0.59159788) q[0];
sx q[0];
rz(-2.199882) q[0];
rz(-0.038453728) q[1];
sx q[1];
rz(-1.4310623) q[1];
sx q[1];
rz(-0.11631913) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424902) q[0];
sx q[0];
rz(-2.0138055) q[0];
sx q[0];
rz(-0.82216121) q[0];
rz(-pi) q[1];
rz(-0.88960008) q[2];
sx q[2];
rz(-1.7187013) q[2];
sx q[2];
rz(-3.0513632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.460187) q[1];
sx q[1];
rz(-1.3906154) q[1];
sx q[1];
rz(-0.26162511) q[1];
x q[2];
rz(-1.4606455) q[3];
sx q[3];
rz(-1.9578736) q[3];
sx q[3];
rz(0.36110669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2945127) q[2];
sx q[2];
rz(-2.3307762) q[2];
sx q[2];
rz(-2.1152367) q[2];
rz(1.2096842) q[3];
sx q[3];
rz(-2.1116202) q[3];
sx q[3];
rz(2.1799555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66861361) q[0];
sx q[0];
rz(-1.0740148) q[0];
sx q[0];
rz(-2.7630254) q[0];
rz(1.1096795) q[1];
sx q[1];
rz(-0.30866426) q[1];
sx q[1];
rz(-3.1139156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8016114) q[0];
sx q[0];
rz(-1.9343573) q[0];
sx q[0];
rz(-2.5021863) q[0];
x q[1];
rz(-0.1849298) q[2];
sx q[2];
rz(-1.3580631) q[2];
sx q[2];
rz(0.30480584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8308786) q[1];
sx q[1];
rz(-0.40604499) q[1];
sx q[1];
rz(2.0534404) q[1];
rz(-2.1975193) q[3];
sx q[3];
rz(-0.51261307) q[3];
sx q[3];
rz(-0.98820283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8370342) q[2];
sx q[2];
rz(-0.95418945) q[2];
sx q[2];
rz(0.41998106) q[2];
rz(0.84152451) q[3];
sx q[3];
rz(-2.1244815) q[3];
sx q[3];
rz(-0.37874547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3751752) q[0];
sx q[0];
rz(-0.11688047) q[0];
sx q[0];
rz(-2.8564659) q[0];
rz(2.9410703) q[1];
sx q[1];
rz(-1.9384117) q[1];
sx q[1];
rz(-2.3060422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1932521) q[0];
sx q[0];
rz(-0.86328816) q[0];
sx q[0];
rz(2.18594) q[0];
rz(-3.1329536) q[2];
sx q[2];
rz(-0.64707478) q[2];
sx q[2];
rz(-1.6948014) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0802059) q[1];
sx q[1];
rz(-2.5701984) q[1];
sx q[1];
rz(0.6462884) q[1];
x q[2];
rz(0.16160139) q[3];
sx q[3];
rz(-0.23862632) q[3];
sx q[3];
rz(-2.0493226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42651132) q[2];
sx q[2];
rz(-2.0501523) q[2];
sx q[2];
rz(-0.69768989) q[2];
rz(0.52608025) q[3];
sx q[3];
rz(-0.79280058) q[3];
sx q[3];
rz(1.5066159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1220916) q[0];
sx q[0];
rz(-2.2631336) q[0];
sx q[0];
rz(2.7418131) q[0];
rz(-1.1493692) q[1];
sx q[1];
rz(-0.72723564) q[1];
sx q[1];
rz(-0.74692187) q[1];
rz(-2.2483038) q[2];
sx q[2];
rz(-1.2947686) q[2];
sx q[2];
rz(-0.99560621) q[2];
rz(-2.5935843) q[3];
sx q[3];
rz(-1.8262134) q[3];
sx q[3];
rz(-1.2977737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
