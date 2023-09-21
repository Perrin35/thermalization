OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(-0.76173705) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(-2.0643056) q[1];
sx q[1];
rz(0.74365562) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0775454) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(1.8403948) q[0];
x q[1];
rz(-0.6526297) q[2];
sx q[2];
rz(-1.666781) q[2];
sx q[2];
rz(-1.8362311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8929157) q[1];
sx q[1];
rz(-0.49202737) q[1];
sx q[1];
rz(0.9534652) q[1];
rz(-pi) q[2];
rz(-0.071804382) q[3];
sx q[3];
rz(-1.9339438) q[3];
sx q[3];
rz(2.8920409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3687392) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(-2.1570737) q[0];
rz(-pi) q[1];
rz(1.8789005) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(-1.845713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7296655) q[1];
sx q[1];
rz(-2.0265059) q[1];
sx q[1];
rz(-0.63974849) q[1];
rz(-pi) q[2];
rz(1.9713692) q[3];
sx q[3];
rz(-1.2188606) q[3];
sx q[3];
rz(0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-0.67726642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30401858) q[0];
sx q[0];
rz(-1.7390334) q[0];
sx q[0];
rz(-1.5181245) q[0];
x q[1];
rz(1.5718939) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(2.9874143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87264204) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(-0.93572576) q[1];
x q[2];
rz(-0.47795313) q[3];
sx q[3];
rz(-1.2198997) q[3];
sx q[3];
rz(1.4357476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(-0.81165195) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(-2.6507846) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.1725918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2617944) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(-1.7394702) q[0];
rz(-pi) q[1];
rz(-1.5863717) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(1.3941927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6231411) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(1.2007984) q[1];
x q[2];
rz(1.4106393) q[3];
sx q[3];
rz(-0.43678108) q[3];
sx q[3];
rz(-0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(2.6546997) q[2];
rz(-1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(2.4868734) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514873) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(-2.4672227) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7582558) q[2];
sx q[2];
rz(-1.9831295) q[2];
sx q[2];
rz(0.86134855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59963804) q[1];
sx q[1];
rz(-1.1918187) q[1];
sx q[1];
rz(-2.42948) q[1];
rz(-2.914364) q[3];
sx q[3];
rz(-0.39441808) q[3];
sx q[3];
rz(1.7928746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.6437795) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2346674) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(-0.90953565) q[0];
x q[1];
rz(2.3408893) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(0.18037361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3469997) q[1];
sx q[1];
rz(-2.0302677) q[1];
sx q[1];
rz(-0.32230349) q[1];
x q[2];
rz(0.27961126) q[3];
sx q[3];
rz(-1.3307443) q[3];
sx q[3];
rz(-2.1820501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(2.3251422) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291173) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-0.11418848) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(0.79089975) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.233477) q[0];
sx q[0];
rz(-2.18963) q[0];
sx q[0];
rz(-0.11615642) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7795243) q[2];
sx q[2];
rz(-2.3644591) q[2];
sx q[2];
rz(0.91285489) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0485718) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(3.0572592) q[1];
x q[2];
rz(0.88364717) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131322) q[0];
sx q[0];
rz(-1.5043133) q[0];
sx q[0];
rz(0.19595887) q[0];
rz(-pi) q[1];
rz(2.8587564) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(-1.1495513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64995631) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(1.7040571) q[1];
rz(-pi) q[2];
rz(1.4905606) q[3];
sx q[3];
rz(-0.9730162) q[3];
sx q[3];
rz(-2.762592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-0.33561486) q[2];
rz(0.32661682) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083387233) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(2.1597247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69753416) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(-0.27192893) q[0];
rz(3.0332546) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(0.94891753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6135785) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(1.8730875) q[1];
x q[2];
rz(0.710886) q[3];
sx q[3];
rz(-1.3804187) q[3];
sx q[3];
rz(-3.0829317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(-2.1441933) q[2];
rz(0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-0.26836747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8440588) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(-0.78155078) q[0];
x q[1];
rz(-2.1860113) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(-2.5126484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8908773) q[1];
sx q[1];
rz(-0.6579352) q[1];
sx q[1];
rz(-0.4472181) q[1];
rz(-pi) q[2];
rz(0.26465613) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(-2.2789126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(-2.5777585) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47678369) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(2.6396991) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
rz(-1.7715122) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
