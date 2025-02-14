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
rz(1.9873729) q[0];
sx q[0];
rz(4.6648751) q[0];
sx q[0];
rz(9.2821791) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(-2.4886517) q[1];
sx q[1];
rz(-1.0452193) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0526177) q[0];
sx q[0];
rz(-2.0988905) q[0];
sx q[0];
rz(-1.0333956) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93306834) q[2];
sx q[2];
rz(-1.6003813) q[2];
sx q[2];
rz(-2.3552908) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18148566) q[1];
sx q[1];
rz(-1.4131261) q[1];
sx q[1];
rz(0.10206359) q[1];
rz(-2.1395438) q[3];
sx q[3];
rz(-1.6709329) q[3];
sx q[3];
rz(-1.1796234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97592252) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(2.9284076) q[2];
rz(-2.2255157) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(-2.7540414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85584545) q[0];
sx q[0];
rz(-2.2540932) q[0];
sx q[0];
rz(2.6708653) q[0];
rz(1.1031021) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(-2.5618166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9410163) q[0];
sx q[0];
rz(-0.52768702) q[0];
sx q[0];
rz(-3.0224968) q[0];
rz(-pi) q[1];
rz(-2.6227399) q[2];
sx q[2];
rz(-0.60728549) q[2];
sx q[2];
rz(2.3407206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75232154) q[1];
sx q[1];
rz(-2.139341) q[1];
sx q[1];
rz(2.6804377) q[1];
rz(-pi) q[2];
rz(2.1768119) q[3];
sx q[3];
rz(-2.7005115) q[3];
sx q[3];
rz(-1.4853322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(3.0465872) q[2];
rz(-2.5659918) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(-1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59548241) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(-2.8693759) q[0];
rz(0.1046003) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(1.6786172) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3264831) q[0];
sx q[0];
rz(-0.72569752) q[0];
sx q[0];
rz(1.4952907) q[0];
x q[1];
rz(-1.6350503) q[2];
sx q[2];
rz(-0.78650606) q[2];
sx q[2];
rz(0.74555874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.453664) q[1];
sx q[1];
rz(-1.2896104) q[1];
sx q[1];
rz(2.4635876) q[1];
rz(1.5103064) q[3];
sx q[3];
rz(-1.20245) q[3];
sx q[3];
rz(2.9591536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0649123) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(0.40985516) q[2];
rz(0.05154933) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16875295) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(2.4421316) q[0];
rz(0.93060023) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(-1.0994937) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4874026) q[0];
sx q[0];
rz(-2.497597) q[0];
sx q[0];
rz(-2.0430768) q[0];
x q[1];
rz(-1.1874974) q[2];
sx q[2];
rz(-0.58247236) q[2];
sx q[2];
rz(0.7815643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2840957) q[1];
sx q[1];
rz(-1.2624143) q[1];
sx q[1];
rz(-0.33607884) q[1];
x q[2];
rz(0.52249281) q[3];
sx q[3];
rz(-0.71194369) q[3];
sx q[3];
rz(-2.5639736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76116556) q[2];
sx q[2];
rz(-2.0695504) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(-0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(0.57653069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.024427323) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(-1.0166136) q[0];
rz(-0.92574614) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(0.076676682) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163239) q[0];
sx q[0];
rz(-1.7364572) q[0];
sx q[0];
rz(-1.9693841) q[0];
rz(-1.3252649) q[2];
sx q[2];
rz(-1.3708198) q[2];
sx q[2];
rz(2.1918346) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5591084) q[1];
sx q[1];
rz(-2.541572) q[1];
sx q[1];
rz(-1.9514854) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.816964) q[3];
sx q[3];
rz(-1.7868407) q[3];
sx q[3];
rz(2.1271355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-0.7879476) q[2];
sx q[2];
rz(-0.12316556) q[2];
rz(0.67433107) q[3];
sx q[3];
rz(-2.6682523) q[3];
sx q[3];
rz(2.3136852) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8318361) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(-0.48102608) q[0];
rz(2.9006529) q[1];
sx q[1];
rz(-2.6048581) q[1];
sx q[1];
rz(-0.13490881) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99687679) q[0];
sx q[0];
rz(-0.92665387) q[0];
sx q[0];
rz(-2.5597325) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4675702) q[2];
sx q[2];
rz(-1.0380259) q[2];
sx q[2];
rz(2.4599883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0162867) q[1];
sx q[1];
rz(-1.8452438) q[1];
sx q[1];
rz(0.80933615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44996275) q[3];
sx q[3];
rz(-1.5057766) q[3];
sx q[3];
rz(-1.7700333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2065108) q[2];
sx q[2];
rz(-2.2057081) q[2];
sx q[2];
rz(-1.8180465) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226591) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(2.4830699) q[0];
rz(-1.5918484) q[1];
sx q[1];
rz(-1.0127944) q[1];
sx q[1];
rz(-2.6351567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69909912) q[0];
sx q[0];
rz(-2.7662219) q[0];
sx q[0];
rz(2.9604549) q[0];
rz(3.1108702) q[2];
sx q[2];
rz(-1.0892158) q[2];
sx q[2];
rz(-2.5831163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5750908) q[1];
sx q[1];
rz(-0.12037863) q[1];
sx q[1];
rz(-2.4574404) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5793367) q[3];
sx q[3];
rz(-1.6054573) q[3];
sx q[3];
rz(-0.12753419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(-0.55240101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34614554) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(-0.41918293) q[0];
rz(3.0508793) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(1.027164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61495435) q[0];
sx q[0];
rz(-1.1246343) q[0];
sx q[0];
rz(2.5102708) q[0];
rz(-pi) q[1];
rz(1.5787324) q[2];
sx q[2];
rz(-1.5936226) q[2];
sx q[2];
rz(0.16032444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0106869) q[1];
sx q[1];
rz(-1.3584777) q[1];
sx q[1];
rz(-3.0328018) q[1];
x q[2];
rz(1.5040565) q[3];
sx q[3];
rz(-0.68664521) q[3];
sx q[3];
rz(0.29717577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(-0.99009222) q[2];
rz(0.31791911) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7541499) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-2.5842174) q[0];
rz(-0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-0.16709669) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97826577) q[0];
sx q[0];
rz(-1.3045132) q[0];
sx q[0];
rz(2.8243218) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9252824) q[2];
sx q[2];
rz(-0.71291332) q[2];
sx q[2];
rz(2.0305433) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7726248) q[1];
sx q[1];
rz(-2.4030622) q[1];
sx q[1];
rz(-0.46486295) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4792419) q[3];
sx q[3];
rz(-2.4560438) q[3];
sx q[3];
rz(2.8753672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1780213) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(-0.42845217) q[2];
rz(2.4109449) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927602) q[0];
sx q[0];
rz(-1.6535783) q[0];
sx q[0];
rz(2.0220508) q[0];
rz(-0.66850942) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(-0.25984919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4762806) q[0];
sx q[0];
rz(-0.57054936) q[0];
sx q[0];
rz(-2.0670939) q[0];
x q[1];
rz(2.0374184) q[2];
sx q[2];
rz(-1.2248618) q[2];
sx q[2];
rz(-0.26901252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30505896) q[1];
sx q[1];
rz(-1.1593474) q[1];
sx q[1];
rz(1.9484503) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15774653) q[3];
sx q[3];
rz(-2.254527) q[3];
sx q[3];
rz(0.58303787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4708289) q[2];
sx q[2];
rz(-1.2329817) q[2];
sx q[2];
rz(2.0551576) q[2];
rz(-0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(2.4194748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5928741) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(-1.2552352) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(-1.7931425) q[2];
sx q[2];
rz(-2.7334474) q[2];
sx q[2];
rz(1.853142) q[2];
rz(-1.4896354) q[3];
sx q[3];
rz(-2.4774144) q[3];
sx q[3];
rz(-1.7437205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
