OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(-0.12726769) q[0];
sx q[0];
rz(0.98841086) q[0];
rz(-0.53108162) q[1];
sx q[1];
rz(3.4903033) q[1];
sx q[1];
rz(11.217584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006346) q[0];
sx q[0];
rz(-1.4057584) q[0];
sx q[0];
rz(-1.4532695) q[0];
rz(-pi) q[1];
rz(0.26308194) q[2];
sx q[2];
rz(-1.4882659) q[2];
sx q[2];
rz(-2.6334327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8475436) q[1];
sx q[1];
rz(-2.2806892) q[1];
sx q[1];
rz(-2.3178046) q[1];
rz(2.3380321) q[3];
sx q[3];
rz(-2.8674012) q[3];
sx q[3];
rz(0.95925946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-0.72221243) q[2];
sx q[2];
rz(0.34525004) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-2.6510986) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17125601) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(-0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(-2.125724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.779218) q[0];
sx q[0];
rz(-1.9694766) q[0];
sx q[0];
rz(-2.6108512) q[0];
rz(-0.83769669) q[2];
sx q[2];
rz(-2.838755) q[2];
sx q[2];
rz(2.8267415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.42900436) q[1];
sx q[1];
rz(-1.9468369) q[1];
sx q[1];
rz(1.7141378) q[1];
x q[2];
rz(-0.85605551) q[3];
sx q[3];
rz(-2.5119009) q[3];
sx q[3];
rz(0.39470181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(1.1091728) q[2];
rz(0.43168133) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(0.88009673) q[0];
rz(-1.8799211) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(2.9601011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3767552) q[0];
sx q[0];
rz(-1.1163045) q[0];
sx q[0];
rz(2.6548813) q[0];
x q[1];
rz(2.491465) q[2];
sx q[2];
rz(-2.1047154) q[2];
sx q[2];
rz(-0.18661737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.03141244) q[1];
sx q[1];
rz(-1.8437779) q[1];
sx q[1];
rz(-0.024755342) q[1];
rz(-pi) q[2];
rz(-2.6604026) q[3];
sx q[3];
rz(-1.378873) q[3];
sx q[3];
rz(-0.47062518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(-2.4530607) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(-2.9825488) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6200136) q[0];
sx q[0];
rz(-2.0984681) q[0];
sx q[0];
rz(2.2844727) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7646388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(2.8958547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(3.1161518) q[1];
rz(-pi) q[2];
rz(2.6166797) q[3];
sx q[3];
rz(-1.3789) q[3];
sx q[3];
rz(-1.6247941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8461385) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(-2.7159178) q[2];
rz(3.1221636) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46538019) q[0];
sx q[0];
rz(-3.1349482) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(-0.68583268) q[1];
sx q[1];
rz(-2.1452955) q[1];
sx q[1];
rz(-0.54840666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0869285) q[0];
sx q[0];
rz(-1.6097277) q[0];
sx q[0];
rz(0.995308) q[0];
rz(-0.19210179) q[2];
sx q[2];
rz(-2.3641799) q[2];
sx q[2];
rz(-2.8156413) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.027443176) q[1];
sx q[1];
rz(-0.47361923) q[1];
sx q[1];
rz(-1.7042392) q[1];
rz(1.3324276) q[3];
sx q[3];
rz(-2.7276464) q[3];
sx q[3];
rz(0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.4667286) q[2];
rz(-2.0560125) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(-0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.1941432) q[0];
sx q[0];
rz(2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(0.99463314) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1625527) q[0];
sx q[0];
rz(-1.0096706) q[0];
sx q[0];
rz(-2.8812376) q[0];
rz(-pi) q[1];
rz(0.0066130916) q[2];
sx q[2];
rz(-2.6380739) q[2];
sx q[2];
rz(0.40724444) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99018807) q[1];
sx q[1];
rz(-1.2328887) q[1];
sx q[1];
rz(-0.0068454725) q[1];
rz(-1.1676222) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(2.6350104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2991335) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(-2.3343202) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(-2.2454967) q[0];
rz(0.90944666) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(-0.075008579) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1589926) q[0];
sx q[0];
rz(-1.2121965) q[0];
sx q[0];
rz(0.067111777) q[0];
x q[1];
rz(-0.54589097) q[2];
sx q[2];
rz(-2.2635169) q[2];
sx q[2];
rz(0.32165124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1753462) q[1];
sx q[1];
rz(-0.65755492) q[1];
sx q[1];
rz(1.6772126) q[1];
x q[2];
rz(-2.6610664) q[3];
sx q[3];
rz(-0.98283813) q[3];
sx q[3];
rz(0.5476391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16671495) q[2];
sx q[2];
rz(-1.1324984) q[2];
sx q[2];
rz(-1.2274851) q[2];
rz(3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1535004) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(-0.41326997) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-2.4024898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9773213) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(3.0493899) q[0];
rz(-pi) q[1];
rz(1.1147538) q[2];
sx q[2];
rz(-1.5653492) q[2];
sx q[2];
rz(-1.3378439) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.821927) q[1];
sx q[1];
rz(-0.85532665) q[1];
sx q[1];
rz(-2.9294076) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3063752) q[3];
sx q[3];
rz(-1.5167055) q[3];
sx q[3];
rz(-2.5097178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83453137) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(-1.1018264) q[2];
rz(-0.39673355) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(-0.81469369) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6535764) q[0];
sx q[0];
rz(-2.5119913) q[0];
sx q[0];
rz(0.6189515) q[0];
rz(1.0194107) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(2.6143262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4370678) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(1.8591451) q[0];
rz(-2.1617266) q[2];
sx q[2];
rz(-0.96989378) q[2];
sx q[2];
rz(0.88496937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4317324) q[1];
sx q[1];
rz(-2.0159971) q[1];
sx q[1];
rz(-2.6688337) q[1];
rz(-pi) q[2];
rz(0.91068565) q[3];
sx q[3];
rz(-1.761697) q[3];
sx q[3];
rz(1.1235352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28425372) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(2.2100892) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(1.5739937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60944027) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(2.2562064) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(0.17818174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70885926) q[0];
sx q[0];
rz(-1.0646001) q[0];
sx q[0];
rz(-0.83336713) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4536132) q[2];
sx q[2];
rz(-1.2974206) q[2];
sx q[2];
rz(0.84699398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1106995) q[1];
sx q[1];
rz(-2.5712214) q[1];
sx q[1];
rz(-1.13899) q[1];
rz(-pi) q[2];
x q[2];
rz(1.448477) q[3];
sx q[3];
rz(-1.9983851) q[3];
sx q[3];
rz(1.7784255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(-0.12100425) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9085893) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(-1.9227149) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(-1.490996) q[2];
sx q[2];
rz(-0.36016338) q[2];
sx q[2];
rz(2.2463837) q[2];
rz(2.8015295) q[3];
sx q[3];
rz(-2.1652514) q[3];
sx q[3];
rz(1.3596331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
