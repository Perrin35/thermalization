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
rz(1.0267216) q[0];
sx q[0];
rz(10.717681) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(8.6051056) q[1];
sx q[1];
rz(17.65045) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0831868) q[0];
sx q[0];
rz(-1.5082486) q[0];
sx q[0];
rz(3.0330171) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9134175) q[2];
sx q[2];
rz(-2.4552305) q[2];
sx q[2];
rz(-0.31191269) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.407302) q[1];
sx q[1];
rz(-1.764174) q[1];
sx q[1];
rz(-0.9949245) q[1];
rz(-2.5254123) q[3];
sx q[3];
rz(-2.1839704) q[3];
sx q[3];
rz(2.4626436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6370411) q[2];
sx q[2];
rz(-1.7402288) q[2];
sx q[2];
rz(-1.6999647) q[2];
rz(2.1291034) q[3];
sx q[3];
rz(-1.9952521) q[3];
sx q[3];
rz(-2.6700524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83990324) q[0];
sx q[0];
rz(-0.56115264) q[0];
sx q[0];
rz(1.0294234) q[0];
rz(1.3599716) q[1];
sx q[1];
rz(-1.1861035) q[1];
sx q[1];
rz(2.634868) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6678807) q[0];
sx q[0];
rz(-0.10598826) q[0];
sx q[0];
rz(1.3876268) q[0];
x q[1];
rz(1.1400181) q[2];
sx q[2];
rz(-0.74249858) q[2];
sx q[2];
rz(3.1378821) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6935061) q[1];
sx q[1];
rz(-2.5593981) q[1];
sx q[1];
rz(1.6901928) q[1];
rz(-pi) q[2];
rz(0.30015035) q[3];
sx q[3];
rz(-0.86816278) q[3];
sx q[3];
rz(0.62338487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3158675) q[2];
sx q[2];
rz(-0.59810144) q[2];
sx q[2];
rz(2.9244847) q[2];
rz(1.0202967) q[3];
sx q[3];
rz(-1.3291357) q[3];
sx q[3];
rz(0.98062688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8751136) q[0];
sx q[0];
rz(-1.3631692) q[0];
sx q[0];
rz(0.76817051) q[0];
rz(0.29113302) q[1];
sx q[1];
rz(-1.365064) q[1];
sx q[1];
rz(1.9545782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087045036) q[0];
sx q[0];
rz(-1.7816418) q[0];
sx q[0];
rz(1.6561942) q[0];
rz(0.044069604) q[2];
sx q[2];
rz(-1.6440834) q[2];
sx q[2];
rz(-1.2512887) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2666013) q[1];
sx q[1];
rz(-0.35540798) q[1];
sx q[1];
rz(-0.0796109) q[1];
rz(-pi) q[2];
rz(0.94908917) q[3];
sx q[3];
rz(-1.177985) q[3];
sx q[3];
rz(1.9167241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8698296) q[2];
sx q[2];
rz(-2.2972079) q[2];
sx q[2];
rz(2.1738906) q[2];
rz(-2.9076231) q[3];
sx q[3];
rz(-2.440019) q[3];
sx q[3];
rz(2.7847737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135603) q[0];
sx q[0];
rz(-1.5505294) q[0];
sx q[0];
rz(1.0796984) q[0];
rz(-2.8375541) q[1];
sx q[1];
rz(-1.9875151) q[1];
sx q[1];
rz(1.3160204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94155055) q[0];
sx q[0];
rz(-0.5467177) q[0];
sx q[0];
rz(-0.23970516) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6929246) q[2];
sx q[2];
rz(-0.59905648) q[2];
sx q[2];
rz(-2.0640822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5583937) q[1];
sx q[1];
rz(-2.2474562) q[1];
sx q[1];
rz(-0.79460245) q[1];
rz(-pi) q[2];
rz(1.1445266) q[3];
sx q[3];
rz(-1.5139587) q[3];
sx q[3];
rz(2.6714966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6198081) q[2];
sx q[2];
rz(-1.2336171) q[2];
sx q[2];
rz(-0.141315) q[2];
rz(2.5610793) q[3];
sx q[3];
rz(-1.4760735) q[3];
sx q[3];
rz(-2.9132402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572131) q[0];
sx q[0];
rz(-0.79638201) q[0];
sx q[0];
rz(-1.6656531) q[0];
rz(-0.93302226) q[1];
sx q[1];
rz(-0.7111744) q[1];
sx q[1];
rz(1.7500386) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6835502) q[0];
sx q[0];
rz(-1.453896) q[0];
sx q[0];
rz(-0.83302541) q[0];
rz(-1.4165164) q[2];
sx q[2];
rz(-1.1517467) q[2];
sx q[2];
rz(1.8928469) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4310225) q[1];
sx q[1];
rz(-2.6198513) q[1];
sx q[1];
rz(1.6304593) q[1];
x q[2];
rz(2.103533) q[3];
sx q[3];
rz(-1.8577855) q[3];
sx q[3];
rz(2.5764333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89620227) q[2];
sx q[2];
rz(-0.37084493) q[2];
sx q[2];
rz(0.24946269) q[2];
rz(-1.6438515) q[3];
sx q[3];
rz(-1.7353053) q[3];
sx q[3];
rz(-0.41516414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45667085) q[0];
sx q[0];
rz(-0.49600729) q[0];
sx q[0];
rz(-1.7556835) q[0];
rz(-1.2617525) q[1];
sx q[1];
rz(-1.2077786) q[1];
sx q[1];
rz(2.6127846) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518327) q[0];
sx q[0];
rz(-2.6963137) q[0];
sx q[0];
rz(1.5543429) q[0];
rz(-pi) q[1];
rz(2.5031935) q[2];
sx q[2];
rz(-1.0784454) q[2];
sx q[2];
rz(-1.3251208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
rz(-1.7011794) q[3];
sx q[3];
rz(-2.2516139) q[3];
sx q[3];
rz(-0.21475204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1888107) q[2];
sx q[2];
rz(-0.47241259) q[2];
sx q[2];
rz(-1.7702276) q[2];
rz(2.93907) q[3];
sx q[3];
rz(-1.4684418) q[3];
sx q[3];
rz(1.9523581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.33223575) q[0];
sx q[0];
rz(-1.7928596) q[0];
sx q[0];
rz(0.15383823) q[0];
rz(2.4065252) q[1];
sx q[1];
rz(-2.0497597) q[1];
sx q[1];
rz(0.14911266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6693838) q[0];
sx q[0];
rz(-0.65919224) q[0];
sx q[0];
rz(-1.4344331) q[0];
rz(-1.2066226) q[2];
sx q[2];
rz(-1.7330164) q[2];
sx q[2];
rz(-0.71031865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.393776) q[1];
sx q[1];
rz(-0.92499176) q[1];
sx q[1];
rz(0.95816905) q[1];
rz(-pi) q[2];
rz(-2.074769) q[3];
sx q[3];
rz(-1.4912581) q[3];
sx q[3];
rz(2.7635637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4601124) q[2];
sx q[2];
rz(-1.6296547) q[2];
sx q[2];
rz(0.45664772) q[2];
rz(-1.7234507) q[3];
sx q[3];
rz(-2.2599615) q[3];
sx q[3];
rz(0.68896967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6316471) q[0];
sx q[0];
rz(-2.8896285) q[0];
sx q[0];
rz(3.0889567) q[0];
rz(1.3098199) q[1];
sx q[1];
rz(-0.24886623) q[1];
sx q[1];
rz(-1.9356669) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9034957) q[0];
sx q[0];
rz(-0.66199979) q[0];
sx q[0];
rz(0.9107301) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5809999) q[2];
sx q[2];
rz(-1.3729323) q[2];
sx q[2];
rz(-2.3677111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0830906) q[1];
sx q[1];
rz(-1.7278226) q[1];
sx q[1];
rz(0.67879063) q[1];
x q[2];
rz(-0.94843978) q[3];
sx q[3];
rz(-1.3788169) q[3];
sx q[3];
rz(0.49302671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2974818) q[2];
sx q[2];
rz(-1.7606807) q[2];
sx q[2];
rz(0.55997509) q[2];
rz(-1.0567793) q[3];
sx q[3];
rz(-0.70928514) q[3];
sx q[3];
rz(0.81282508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.63447222) q[0];
sx q[0];
rz(-1.7422603) q[0];
sx q[0];
rz(-1.5768453) q[0];
rz(-1.5722081) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(2.3326468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841815) q[0];
sx q[0];
rz(-2.0657259) q[0];
sx q[0];
rz(1.5280426) q[0];
rz(-pi) q[1];
rz(0.39500631) q[2];
sx q[2];
rz(-0.85112903) q[2];
sx q[2];
rz(2.0499944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8524127) q[1];
sx q[1];
rz(-0.9599664) q[1];
sx q[1];
rz(3.1256366) q[1];
x q[2];
rz(-1.21502) q[3];
sx q[3];
rz(-0.18863931) q[3];
sx q[3];
rz(-1.1805746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6605777) q[2];
sx q[2];
rz(-2.1889841) q[2];
sx q[2];
rz(3.0926404) q[2];
rz(-1.281721) q[3];
sx q[3];
rz(-2.6661524) q[3];
sx q[3];
rz(1.4648645) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82185164) q[0];
sx q[0];
rz(-2.8806683) q[0];
sx q[0];
rz(1.2224181) q[0];
rz(1.4756731) q[1];
sx q[1];
rz(-1.2011352) q[1];
sx q[1];
rz(-0.45752057) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7675661) q[0];
sx q[0];
rz(-1.6609285) q[0];
sx q[0];
rz(0.73390168) q[0];
rz(-pi) q[1];
rz(-1.0636163) q[2];
sx q[2];
rz(-1.9263144) q[2];
sx q[2];
rz(1.4414444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19472518) q[1];
sx q[1];
rz(-1.3683142) q[1];
sx q[1];
rz(2.7848886) q[1];
rz(-1.1531257) q[3];
sx q[3];
rz(-2.2137797) q[3];
sx q[3];
rz(-2.2394799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7484625) q[2];
sx q[2];
rz(-0.25105432) q[2];
sx q[2];
rz(2.1012696) q[2];
rz(0.48804247) q[3];
sx q[3];
rz(-1.4405684) q[3];
sx q[3];
rz(0.53781167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77448612) q[0];
sx q[0];
rz(-2.0699061) q[0];
sx q[0];
rz(0.62181428) q[0];
rz(-1.294301) q[1];
sx q[1];
rz(-0.58302561) q[1];
sx q[1];
rz(2.2517712) q[1];
rz(2.842335) q[2];
sx q[2];
rz(-0.99874961) q[2];
sx q[2];
rz(-0.082538907) q[2];
rz(-2.4518012) q[3];
sx q[3];
rz(-1.2648911) q[3];
sx q[3];
rz(-2.779724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
