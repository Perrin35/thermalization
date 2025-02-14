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
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(1.8486899) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(-1.9424865) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0082675) q[0];
sx q[0];
rz(-3.016351) q[0];
sx q[0];
rz(-2.6175346) q[0];
x q[1];
rz(-0.22817518) q[2];
sx q[2];
rz(-0.68636218) q[2];
sx q[2];
rz(0.31191269) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4292214) q[1];
sx q[1];
rz(-1.0069798) q[1];
sx q[1];
rz(0.22937201) q[1];
rz(-pi) q[2];
rz(0.88319234) q[3];
sx q[3];
rz(-0.84005648) q[3];
sx q[3];
rz(2.9321364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6370411) q[2];
sx q[2];
rz(-1.4013638) q[2];
sx q[2];
rz(1.6999647) q[2];
rz(-2.1291034) q[3];
sx q[3];
rz(-1.1463405) q[3];
sx q[3];
rz(0.47154021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3016894) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(2.1121693) q[0];
rz(-1.3599716) q[1];
sx q[1];
rz(-1.9554892) q[1];
sx q[1];
rz(-0.50672466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27924867) q[0];
sx q[0];
rz(-1.5900668) q[0];
sx q[0];
rz(-1.6750245) q[0];
rz(-2.7756491) q[2];
sx q[2];
rz(-2.2323713) q[2];
sx q[2];
rz(-0.55398527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9190127) q[1];
sx q[1];
rz(-1.5052541) q[1];
sx q[1];
rz(-2.1497141) q[1];
x q[2];
rz(-1.2348753) q[3];
sx q[3];
rz(-0.75390076) q[3];
sx q[3];
rz(-0.17673211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3158675) q[2];
sx q[2];
rz(-2.5434912) q[2];
sx q[2];
rz(-0.21710795) q[2];
rz(-1.0202967) q[3];
sx q[3];
rz(-1.812457) q[3];
sx q[3];
rz(0.98062688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8751136) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-0.76817051) q[0];
rz(0.29113302) q[1];
sx q[1];
rz(-1.7765287) q[1];
sx q[1];
rz(-1.9545782) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087045036) q[0];
sx q[0];
rz(-1.7816418) q[0];
sx q[0];
rz(1.6561942) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6441543) q[2];
sx q[2];
rz(-1.6147476) q[2];
sx q[2];
rz(-2.8188561) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87499133) q[1];
sx q[1];
rz(-0.35540798) q[1];
sx q[1];
rz(-0.0796109) q[1];
rz(-pi) q[2];
rz(0.95245016) q[3];
sx q[3];
rz(-2.420331) q[3];
sx q[3];
rz(-2.9968468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8698296) q[2];
sx q[2];
rz(-0.84438476) q[2];
sx q[2];
rz(-0.96770206) q[2];
rz(-0.23396954) q[3];
sx q[3];
rz(-0.7015737) q[3];
sx q[3];
rz(-0.35681891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02803239) q[0];
sx q[0];
rz(-1.5505294) q[0];
sx q[0];
rz(2.0618942) q[0];
rz(0.30403852) q[1];
sx q[1];
rz(-1.9875151) q[1];
sx q[1];
rz(-1.8255723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9213778) q[0];
sx q[0];
rz(-2.1002033) q[0];
sx q[0];
rz(1.7142943) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.082985254) q[2];
sx q[2];
rz(-0.97681649) q[2];
sx q[2];
rz(1.9165438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5583937) q[1];
sx q[1];
rz(-2.2474562) q[1];
sx q[1];
rz(0.79460245) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.06240978) q[3];
sx q[3];
rz(-1.1452598) q[3];
sx q[3];
rz(-1.1264914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6198081) q[2];
sx q[2];
rz(-1.2336171) q[2];
sx q[2];
rz(0.141315) q[2];
rz(-2.5610793) q[3];
sx q[3];
rz(-1.4760735) q[3];
sx q[3];
rz(2.9132402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88437951) q[0];
sx q[0];
rz(-0.79638201) q[0];
sx q[0];
rz(1.4759395) q[0];
rz(0.93302226) q[1];
sx q[1];
rz(-2.4304183) q[1];
sx q[1];
rz(-1.3915541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0071268) q[0];
sx q[0];
rz(-2.3023805) q[0];
sx q[0];
rz(-2.9842019) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4165164) q[2];
sx q[2];
rz(-1.1517467) q[2];
sx q[2];
rz(-1.2487457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4310225) q[1];
sx q[1];
rz(-0.52174134) q[1];
sx q[1];
rz(1.6304593) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0380596) q[3];
sx q[3];
rz(-1.2838072) q[3];
sx q[3];
rz(0.56515938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89620227) q[2];
sx q[2];
rz(-2.7707477) q[2];
sx q[2];
rz(0.24946269) q[2];
rz(1.6438515) q[3];
sx q[3];
rz(-1.7353053) q[3];
sx q[3];
rz(0.41516414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.45667085) q[0];
sx q[0];
rz(-2.6455854) q[0];
sx q[0];
rz(1.3859092) q[0];
rz(-1.8798401) q[1];
sx q[1];
rz(-1.933814) q[1];
sx q[1];
rz(2.6127846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.945707) q[0];
sx q[0];
rz(-1.5778827) q[0];
sx q[0];
rz(-2.0160227) q[0];
rz(-2.5031935) q[2];
sx q[2];
rz(-2.0631472) q[2];
sx q[2];
rz(-1.3251208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4836135) q[1];
sx q[1];
rz(-0.93320642) q[1];
sx q[1];
rz(-0.083483551) q[1];
x q[2];
rz(2.9824419) q[3];
sx q[3];
rz(-2.4503631) q[3];
sx q[3];
rz(-2.7214512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1888107) q[2];
sx q[2];
rz(-0.47241259) q[2];
sx q[2];
rz(1.3713651) q[2];
rz(-2.93907) q[3];
sx q[3];
rz(-1.4684418) q[3];
sx q[3];
rz(-1.9523581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33223575) q[0];
sx q[0];
rz(-1.7928596) q[0];
sx q[0];
rz(0.15383823) q[0];
rz(-0.73506749) q[1];
sx q[1];
rz(-2.0497597) q[1];
sx q[1];
rz(0.14911266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47220889) q[0];
sx q[0];
rz(-2.4824004) q[0];
sx q[0];
rz(-1.7071595) q[0];
x q[1];
rz(-2.0015099) q[2];
sx q[2];
rz(-0.39719279) q[2];
sx q[2];
rz(0.45969648) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2231517) q[1];
sx q[1];
rz(-1.0936803) q[1];
sx q[1];
rz(2.3972307) q[1];
x q[2];
rz(-0.090773067) q[3];
sx q[3];
rz(-2.0730264) q[3];
sx q[3];
rz(1.2365562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6814802) q[2];
sx q[2];
rz(-1.511938) q[2];
sx q[2];
rz(0.45664772) q[2];
rz(-1.7234507) q[3];
sx q[3];
rz(-0.8816312) q[3];
sx q[3];
rz(-0.68896967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6316471) q[0];
sx q[0];
rz(-0.25196415) q[0];
sx q[0];
rz(0.052635996) q[0];
rz(1.8317728) q[1];
sx q[1];
rz(-0.24886623) q[1];
sx q[1];
rz(1.9356669) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.238097) q[0];
sx q[0];
rz(-0.66199979) q[0];
sx q[0];
rz(0.9107301) q[0];
rz(-pi) q[1];
rz(-0.3605901) q[2];
sx q[2];
rz(-2.5506488) q[2];
sx q[2];
rz(-0.49345582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5284023) q[1];
sx q[1];
rz(-0.90189108) q[1];
sx q[1];
rz(1.7714798) q[1];
x q[2];
rz(-2.9067791) q[3];
sx q[3];
rz(-2.1800332) q[3];
sx q[3];
rz(-1.2138194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.84411088) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(0.55997509) q[2];
rz(2.0848134) q[3];
sx q[3];
rz(-2.4323075) q[3];
sx q[3];
rz(-0.81282508) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63447222) q[0];
sx q[0];
rz(-1.3993323) q[0];
sx q[0];
rz(1.5647474) q[0];
rz(1.5693846) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(-0.80894583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8485238) q[0];
sx q[0];
rz(-1.5331755) q[0];
sx q[0];
rz(2.6462808) q[0];
rz(-pi) q[1];
rz(-2.7465863) q[2];
sx q[2];
rz(-2.2904636) q[2];
sx q[2];
rz(1.0915983) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8802277) q[1];
sx q[1];
rz(-2.5305809) q[1];
sx q[1];
rz(-1.548012) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.06639938) q[3];
sx q[3];
rz(-1.3940991) q[3];
sx q[3];
rz(1.599358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6605777) q[2];
sx q[2];
rz(-2.1889841) q[2];
sx q[2];
rz(-3.0926404) q[2];
rz(1.8598716) q[3];
sx q[3];
rz(-2.6661524) q[3];
sx q[3];
rz(-1.6767282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(0.82185164) q[0];
sx q[0];
rz(-2.8806683) q[0];
sx q[0];
rz(-1.9191746) q[0];
rz(-1.6659196) q[1];
sx q[1];
rz(-1.2011352) q[1];
sx q[1];
rz(-0.45752057) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0973061) q[0];
sx q[0];
rz(-0.73839085) q[0];
sx q[0];
rz(3.0074658) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2234688) q[2];
sx q[2];
rz(-0.61033536) q[2];
sx q[2];
rz(0.68896919) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.690722) q[1];
sx q[1];
rz(-1.2216946) q[1];
sx q[1];
rz(1.786475) q[1];
x q[2];
rz(-0.68660929) q[3];
sx q[3];
rz(-1.2401738) q[3];
sx q[3];
rz(-2.7329684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3931302) q[2];
sx q[2];
rz(-2.8905383) q[2];
sx q[2];
rz(-1.040323) q[2];
rz(2.6535502) q[3];
sx q[3];
rz(-1.4405684) q[3];
sx q[3];
rz(-0.53781167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77448612) q[0];
sx q[0];
rz(-1.0716866) q[0];
sx q[0];
rz(-2.5197784) q[0];
rz(-1.294301) q[1];
sx q[1];
rz(-0.58302561) q[1];
sx q[1];
rz(2.2517712) q[1];
rz(-1.1414083) q[2];
sx q[2];
rz(-0.63776897) q[2];
sx q[2];
rz(-2.706131) q[2];
rz(0.46066416) q[3];
sx q[3];
rz(-0.74429269) q[3];
sx q[3];
rz(-0.85891354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
