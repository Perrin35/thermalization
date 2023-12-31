OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(0.8843511) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8774672) q[0];
sx q[0];
rz(-1.929951) q[0];
sx q[0];
rz(-2.9496664) q[0];
x q[1];
rz(-2.5764478) q[2];
sx q[2];
rz(-1.5663212) q[2];
sx q[2];
rz(2.7067513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.002418) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-2.8789218) q[1];
x q[2];
rz(2.7847071) q[3];
sx q[3];
rz(-2.2483453) q[3];
sx q[3];
rz(0.68912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-2.4856429) q[2];
rz(1.9338699) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(-2.1470127) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(-2.9128089) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(-0.0072335009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418966) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(-1.9641563) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8480808) q[2];
sx q[2];
rz(-1.8732757) q[2];
sx q[2];
rz(2.7345865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0821973) q[1];
sx q[1];
rz(-2.7121183) q[1];
sx q[1];
rz(-0.55179623) q[1];
rz(-pi) q[2];
rz(-2.238027) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(0.29153338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5269512) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(3.0200322) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(0.30383032) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6737297) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(-1.7720222) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(-2.8799768) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8825892) q[0];
sx q[0];
rz(-0.020375229) q[0];
sx q[0];
rz(2.0632319) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4357655) q[2];
sx q[2];
rz(-2.3927852) q[2];
sx q[2];
rz(2.7111862) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6690327) q[1];
sx q[1];
rz(-1.4497888) q[1];
sx q[1];
rz(-2.3585412) q[1];
rz(-pi) q[2];
rz(0.055734169) q[3];
sx q[3];
rz(-0.51364964) q[3];
sx q[3];
rz(1.0518215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4613142) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(-2.938081) q[2];
rz(-2.2198548) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.8932231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0715863) q[0];
sx q[0];
rz(-1.4741815) q[0];
sx q[0];
rz(0.065652547) q[0];
rz(-pi) q[1];
rz(-0.7123956) q[2];
sx q[2];
rz(-2.8678896) q[2];
sx q[2];
rz(-1.7143539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8716988) q[1];
sx q[1];
rz(-1.4596241) q[1];
sx q[1];
rz(-2.486869) q[1];
rz(-pi) q[2];
rz(2.4183211) q[3];
sx q[3];
rz(-1.7687106) q[3];
sx q[3];
rz(2.7151782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1127597) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(-0.83703414) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(2.2391589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31061253) q[0];
sx q[0];
rz(-1.6164301) q[0];
sx q[0];
rz(-2.5020585) q[0];
rz(-0.4544223) q[2];
sx q[2];
rz(-1.7917969) q[2];
sx q[2];
rz(-2.1511252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9296226) q[1];
sx q[1];
rz(-0.8319444) q[1];
sx q[1];
rz(-2.3421939) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5794472) q[3];
sx q[3];
rz(-0.31673613) q[3];
sx q[3];
rz(2.1474311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19501413) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(2.6384171) q[0];
rz(-1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(-0.20176372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70226442) q[0];
sx q[0];
rz(-1.4685255) q[0];
sx q[0];
rz(3.0754473) q[0];
rz(-pi) q[1];
rz(2.6642338) q[2];
sx q[2];
rz(-1.9743894) q[2];
sx q[2];
rz(1.8237643) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7355431) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(-0.72592782) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4089036) q[3];
sx q[3];
rz(-1.9273888) q[3];
sx q[3];
rz(2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21489828) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(-2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8478407) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-0.057549495) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38521117) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(-2.5368607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96037453) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(-2.9396217) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.645694) q[1];
sx q[1];
rz(-1.8897448) q[1];
sx q[1];
rz(-2.0046528) q[1];
x q[2];
rz(1.058217) q[3];
sx q[3];
rz(-1.3430542) q[3];
sx q[3];
rz(-0.77373576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-0.28392917) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2487508) q[0];
sx q[0];
rz(-2.1584956) q[0];
sx q[0];
rz(-2.6177004) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5416652) q[2];
sx q[2];
rz(-1.6919961) q[2];
sx q[2];
rz(-1.3231414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.200951) q[1];
sx q[1];
rz(-1.1108228) q[1];
sx q[1];
rz(0.56401395) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1030032) q[3];
sx q[3];
rz(-0.66017294) q[3];
sx q[3];
rz(0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.930687) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(2.912345) q[0];
rz(-0.30300888) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.680826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0871353) q[0];
sx q[0];
rz(-1.5883049) q[0];
sx q[0];
rz(2.8580335) q[0];
rz(-pi) q[1];
rz(-2.7621208) q[2];
sx q[2];
rz(-1.4037637) q[2];
sx q[2];
rz(0.52969474) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2158828) q[1];
sx q[1];
rz(-0.52535086) q[1];
sx q[1];
rz(-0.044563091) q[1];
rz(-pi) q[2];
rz(2.7542354) q[3];
sx q[3];
rz(-0.82327561) q[3];
sx q[3];
rz(0.92126095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.573695) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(0.38063231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.80523) q[0];
sx q[0];
rz(-0.21348937) q[0];
sx q[0];
rz(-1.2345033) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74110465) q[2];
sx q[2];
rz(-2.0217102) q[2];
sx q[2];
rz(1.2620743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3818647) q[1];
sx q[1];
rz(-1.3761531) q[1];
sx q[1];
rz(2.0774283) q[1];
rz(-2.5638644) q[3];
sx q[3];
rz(-1.2168222) q[3];
sx q[3];
rz(-2.3059394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(2.005119) q[2];
rz(-3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3363591) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-1.0271172) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(-1.3409875) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(2.9410578) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
