OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(2.0988965) q[0];
sx q[0];
rz(12.287696) q[0];
rz(4.2884045) q[1];
sx q[1];
rz(0.52462259) q[1];
sx q[1];
rz(8.1537032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575114) q[0];
sx q[0];
rz(-1.3102828) q[0];
sx q[0];
rz(1.5519014) q[0];
x q[1];
rz(2.991011) q[2];
sx q[2];
rz(-1.7007593) q[2];
sx q[2];
rz(0.58565631) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.179175) q[1];
sx q[1];
rz(-1.6500874) q[1];
sx q[1];
rz(2.2843664) q[1];
rz(-2.1500312) q[3];
sx q[3];
rz(-2.09822) q[3];
sx q[3];
rz(2.9468731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.110454) q[2];
sx q[2];
rz(-1.9240856) q[2];
sx q[2];
rz(0.26220599) q[2];
rz(2.0860705) q[3];
sx q[3];
rz(-1.8942098) q[3];
sx q[3];
rz(-0.087337703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5259041) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(2.9600034) q[0];
rz(-1.7407821) q[1];
sx q[1];
rz(-1.9488275) q[1];
sx q[1];
rz(0.84091944) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070905237) q[0];
sx q[0];
rz(-1.5884261) q[0];
sx q[0];
rz(-0.32004642) q[0];
x q[1];
rz(-2.2107) q[2];
sx q[2];
rz(-1.962365) q[2];
sx q[2];
rz(-1.9545133) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.052472981) q[1];
sx q[1];
rz(-1.6069672) q[1];
sx q[1];
rz(1.7391886) q[1];
rz(-pi) q[2];
rz(2.7737229) q[3];
sx q[3];
rz(-2.1506393) q[3];
sx q[3];
rz(2.232792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.759364) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(-1.0648897) q[2];
rz(-3.0801638) q[3];
sx q[3];
rz(-1.8062785) q[3];
sx q[3];
rz(-1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29838022) q[0];
sx q[0];
rz(-1.1792264) q[0];
sx q[0];
rz(-0.19905736) q[0];
rz(-0.82636034) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(-2.3645649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8733002) q[0];
sx q[0];
rz(-1.4140633) q[0];
sx q[0];
rz(-3.0746033) q[0];
rz(-pi) q[1];
rz(-2.1146745) q[2];
sx q[2];
rz(-1.9444398) q[2];
sx q[2];
rz(-2.0622562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9183082) q[1];
sx q[1];
rz(-1.3121843) q[1];
sx q[1];
rz(-2.1373939) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97717173) q[3];
sx q[3];
rz(-1.8026173) q[3];
sx q[3];
rz(2.1013732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0394773) q[2];
sx q[2];
rz(-1.1481043) q[2];
sx q[2];
rz(0.51304212) q[2];
rz(-0.78220621) q[3];
sx q[3];
rz(-1.9358044) q[3];
sx q[3];
rz(-1.3768844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9272598) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(1.4904892) q[0];
rz(0.53228846) q[1];
sx q[1];
rz(-2.5106301) q[1];
sx q[1];
rz(-1.5375563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89808116) q[0];
sx q[0];
rz(-1.4644721) q[0];
sx q[0];
rz(2.9778019) q[0];
rz(-pi) q[1];
rz(-0.093393383) q[2];
sx q[2];
rz(-1.0265145) q[2];
sx q[2];
rz(1.613008) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2080431) q[1];
sx q[1];
rz(-1.4999823) q[1];
sx q[1];
rz(0.93472247) q[1];
rz(-0.52738366) q[3];
sx q[3];
rz(-1.77447) q[3];
sx q[3];
rz(2.7319997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.535546) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(-3.1287076) q[2];
rz(1.2031201) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(-1.0217246) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0032229) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(-1.3997929) q[0];
rz(2.7895582) q[1];
sx q[1];
rz(-2.0875918) q[1];
sx q[1];
rz(-0.83784109) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1803446) q[0];
sx q[0];
rz(-0.77416468) q[0];
sx q[0];
rz(1.3526239) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1878686) q[2];
sx q[2];
rz(-1.0487818) q[2];
sx q[2];
rz(-1.2921289) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.074118) q[1];
sx q[1];
rz(-2.5606541) q[1];
sx q[1];
rz(0.23260637) q[1];
rz(0.16437634) q[3];
sx q[3];
rz(-1.7225186) q[3];
sx q[3];
rz(-0.80899948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1167404) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(1.9579197) q[2];
rz(-2.8708894) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(-1.5326327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026641332) q[0];
sx q[0];
rz(-0.6237492) q[0];
sx q[0];
rz(0.57914105) q[0];
rz(1.9193513) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(2.0319669) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752788) q[0];
sx q[0];
rz(-1.7357317) q[0];
sx q[0];
rz(-1.724302) q[0];
rz(0.30495651) q[2];
sx q[2];
rz(-1.5169608) q[2];
sx q[2];
rz(2.3862145) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5286104) q[1];
sx q[1];
rz(-0.41944606) q[1];
sx q[1];
rz(-2.4952013) q[1];
rz(-pi) q[2];
rz(-2.6251276) q[3];
sx q[3];
rz(-2.0544009) q[3];
sx q[3];
rz(2.2276239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1960725) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(-0.35745364) q[2];
rz(1.1899905) q[3];
sx q[3];
rz(-1.2001195) q[3];
sx q[3];
rz(-1.3397217) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59915197) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(2.3495112) q[0];
rz(-2.4374841) q[1];
sx q[1];
rz(-0.45583615) q[1];
sx q[1];
rz(1.5584996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38403758) q[0];
sx q[0];
rz(-0.4641986) q[0];
sx q[0];
rz(-3.1350101) q[0];
x q[1];
rz(-0.59410166) q[2];
sx q[2];
rz(-1.033448) q[2];
sx q[2];
rz(-2.0955476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5643107) q[1];
sx q[1];
rz(-0.4275215) q[1];
sx q[1];
rz(-0.083210817) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24105234) q[3];
sx q[3];
rz(-2.3755382) q[3];
sx q[3];
rz(-2.0360006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.030297) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(3.1374068) q[2];
rz(0.26675102) q[3];
sx q[3];
rz(-1.6175852) q[3];
sx q[3];
rz(0.72223103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58077145) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(0.15604493) q[0];
rz(-1.5566298) q[1];
sx q[1];
rz(-2.2792351) q[1];
sx q[1];
rz(0.59476605) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0069591) q[0];
sx q[0];
rz(-1.3158321) q[0];
sx q[0];
rz(0.50936551) q[0];
rz(-1.2807219) q[2];
sx q[2];
rz(-2.2422487) q[2];
sx q[2];
rz(-2.4460276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76164272) q[1];
sx q[1];
rz(-0.98230668) q[1];
sx q[1];
rz(-2.2698927) q[1];
x q[2];
rz(2.2355763) q[3];
sx q[3];
rz(-0.95905868) q[3];
sx q[3];
rz(0.25463984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56738371) q[2];
sx q[2];
rz(-1.5790066) q[2];
sx q[2];
rz(1.0279921) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-0.65516156) q[3];
sx q[3];
rz(1.3348234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26355711) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(-2.3433319) q[0];
rz(2.3410666) q[1];
sx q[1];
rz(-0.48897484) q[1];
sx q[1];
rz(-2.5953603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29833405) q[0];
sx q[0];
rz(-0.97163659) q[0];
sx q[0];
rz(-0.13183044) q[0];
rz(1.9990218) q[2];
sx q[2];
rz(-1.0248803) q[2];
sx q[2];
rz(-0.8949309) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45193397) q[1];
sx q[1];
rz(-2.0469189) q[1];
sx q[1];
rz(-1.0187638) q[1];
rz(-pi) q[2];
rz(2.915156) q[3];
sx q[3];
rz(-1.4996935) q[3];
sx q[3];
rz(-1.16217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74786782) q[2];
sx q[2];
rz(-2.0515029) q[2];
sx q[2];
rz(-3.0785576) q[2];
rz(2.4437599) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6460655) q[0];
sx q[0];
rz(-0.87764469) q[0];
sx q[0];
rz(-0.416042) q[0];
rz(-1.3620194) q[1];
sx q[1];
rz(-1.1241309) q[1];
sx q[1];
rz(1.8488041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36609367) q[0];
sx q[0];
rz(-1.3028569) q[0];
sx q[0];
rz(1.8084722) q[0];
rz(-pi) q[1];
rz(-2.6854158) q[2];
sx q[2];
rz(-1.8229228) q[2];
sx q[2];
rz(1.1751564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9465465) q[1];
sx q[1];
rz(-2.1113696) q[1];
sx q[1];
rz(-2.5953496) q[1];
rz(-pi) q[2];
rz(-2.8804194) q[3];
sx q[3];
rz(-1.9800485) q[3];
sx q[3];
rz(-0.60200426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4585939) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(1.2609743) q[3];
sx q[3];
rz(-0.57456273) q[3];
sx q[3];
rz(-1.449409) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4954729) q[0];
sx q[0];
rz(-1.6874122) q[0];
sx q[0];
rz(0.73943403) q[0];
rz(-2.6932035) q[1];
sx q[1];
rz(-1.8012451) q[1];
sx q[1];
rz(-1.9582122) q[1];
rz(2.420633) q[2];
sx q[2];
rz(-2.6868281) q[2];
sx q[2];
rz(1.6285298) q[2];
rz(-1.7200593) q[3];
sx q[3];
rz(-2.4076318) q[3];
sx q[3];
rz(-2.9416495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
