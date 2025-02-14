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
rz(-1.1542198) q[0];
sx q[0];
rz(-1.5232824) q[0];
sx q[0];
rz(0.14259882) q[0];
rz(-2.5481186) q[1];
sx q[1];
rz(-0.65294099) q[1];
sx q[1];
rz(-2.0963734) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3680688) q[0];
sx q[0];
rz(-1.1127825) q[0];
sx q[0];
rz(-2.5450443) q[0];
rz(-3.1047761) q[2];
sx q[2];
rz(-2.2082001) q[2];
sx q[2];
rz(0.76257818) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7683623) q[1];
sx q[1];
rz(-1.4700031) q[1];
sx q[1];
rz(1.7292777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3863643) q[3];
sx q[3];
rz(-2.5650555) q[3];
sx q[3];
rz(2.5952796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.97592252) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(0.21318501) q[2];
rz(2.2255157) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(-0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-2.6328937) q[1];
sx q[1];
rz(2.5618166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743917) q[0];
sx q[0];
rz(-1.5109332) q[0];
sx q[0];
rz(2.6169928) q[0];
rz(1.9026372) q[2];
sx q[2];
rz(-2.08925) q[2];
sx q[2];
rz(-1.4084852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75232154) q[1];
sx q[1];
rz(-2.139341) q[1];
sx q[1];
rz(-2.6804377) q[1];
rz(-0.96478077) q[3];
sx q[3];
rz(-2.7005115) q[3];
sx q[3];
rz(1.6562605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(0.095005438) q[2];
rz(-2.5659918) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5461102) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(-0.27221671) q[0];
rz(0.1046003) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(1.6786172) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.187785) q[0];
sx q[0];
rz(-1.5207131) q[0];
sx q[0];
rz(2.2950776) q[0];
rz(1.5065423) q[2];
sx q[2];
rz(-2.3550866) q[2];
sx q[2];
rz(2.3960339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.453664) q[1];
sx q[1];
rz(-1.2896104) q[1];
sx q[1];
rz(2.4635876) q[1];
rz(-pi) q[2];
rz(-2.9862254) q[3];
sx q[3];
rz(-2.7685363) q[3];
sx q[3];
rz(-3.1257926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0649123) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(-2.7317375) q[2];
rz(3.0900433) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(-2.4079017) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16875295) q[0];
sx q[0];
rz(-2.3641455) q[0];
sx q[0];
rz(-0.69946104) q[0];
rz(-0.93060023) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(-2.0420989) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65419009) q[0];
sx q[0];
rz(-2.497597) q[0];
sx q[0];
rz(2.0430768) q[0];
rz(-1.0223939) q[2];
sx q[2];
rz(-1.363593) q[2];
sx q[2];
rz(-1.1140547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18106724) q[1];
sx q[1];
rz(-1.8904422) q[1];
sx q[1];
rz(-1.8962217) q[1];
rz(2.4995125) q[3];
sx q[3];
rz(-1.2386981) q[3];
sx q[3];
rz(-1.737271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76116556) q[2];
sx q[2];
rz(-2.0695504) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(2.565062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.024427323) q[0];
sx q[0];
rz(-3.0553525) q[0];
sx q[0];
rz(1.0166136) q[0];
rz(-2.2158465) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(3.064916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163239) q[0];
sx q[0];
rz(-1.7364572) q[0];
sx q[0];
rz(1.1722086) q[0];
rz(-2.935604) q[2];
sx q[2];
rz(-1.811337) q[2];
sx q[2];
rz(-2.5702916) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5591084) q[1];
sx q[1];
rz(-0.6000207) q[1];
sx q[1];
rz(1.1901072) q[1];
x q[2];
rz(-1.2604976) q[3];
sx q[3];
rz(-0.7782794) q[3];
sx q[3];
rz(-0.33184127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3735247) q[2];
sx q[2];
rz(-0.7879476) q[2];
sx q[2];
rz(0.12316556) q[2];
rz(-0.67433107) q[3];
sx q[3];
rz(-2.6682523) q[3];
sx q[3];
rz(-2.3136852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3097565) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(-2.6605666) q[0];
rz(-2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(3.0066838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19770007) q[0];
sx q[0];
rz(-1.1158082) q[0];
sx q[0];
rz(2.3030192) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53510869) q[2];
sx q[2];
rz(-1.4819179) q[2];
sx q[2];
rz(0.8366226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9489158) q[1];
sx q[1];
rz(-0.8444311) q[1];
sx q[1];
rz(2.7706258) q[1];
x q[2];
rz(2.9929912) q[3];
sx q[3];
rz(-0.45431787) q[3];
sx q[3];
rz(0.065520614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2065108) q[2];
sx q[2];
rz(-2.2057081) q[2];
sx q[2];
rz(-1.3235462) q[2];
rz(2.8694618) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-0.14565295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018933522) q[0];
sx q[0];
rz(-2.3441362) q[0];
sx q[0];
rz(-2.4830699) q[0];
rz(-1.5497442) q[1];
sx q[1];
rz(-1.0127944) q[1];
sx q[1];
rz(-0.50643593) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424935) q[0];
sx q[0];
rz(-0.37537071) q[0];
sx q[0];
rz(-2.9604549) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1108702) q[2];
sx q[2];
rz(-2.0523768) q[2];
sx q[2];
rz(0.55847634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5750908) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(-0.68415227) q[1];
rz(-pi) q[2];
rz(1.6122088) q[3];
sx q[3];
rz(-0.99185252) q[3];
sx q[3];
rz(1.4659297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.2696421) q[2];
rz(1.7757724) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(-2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7954471) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(2.7224097) q[0];
rz(-0.090713352) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(-2.1144287) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8801422) q[0];
sx q[0];
rz(-2.1322726) q[0];
sx q[0];
rz(-2.1057157) q[0];
rz(1.5628603) q[2];
sx q[2];
rz(-1.5479701) q[2];
sx q[2];
rz(-2.9812682) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5324983) q[1];
sx q[1];
rz(-2.9033992) q[1];
sx q[1];
rz(2.0373809) q[1];
x q[2];
rz(-2.2563491) q[3];
sx q[3];
rz(-1.5285057) q[3];
sx q[3];
rz(1.919618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(-0.99009222) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(0.038343553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7541499) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-0.55737525) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-2.974496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1633269) q[0];
sx q[0];
rz(-1.8370795) q[0];
sx q[0];
rz(0.31727088) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4403247) q[2];
sx q[2];
rz(-1.4299562) q[2];
sx q[2];
rz(2.5171211) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7726248) q[1];
sx q[1];
rz(-2.4030622) q[1];
sx q[1];
rz(-0.46486295) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0669555) q[3];
sx q[3];
rz(-2.2529246) q[3];
sx q[3];
rz(-2.9934237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1780213) q[2];
sx q[2];
rz(-2.9489297) q[2];
sx q[2];
rz(-2.7131405) q[2];
rz(0.7306478) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(-0.60034269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44883248) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(2.0220508) q[0];
rz(-2.4730832) q[1];
sx q[1];
rz(-0.57066494) q[1];
sx q[1];
rz(-0.25984919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2370473) q[0];
sx q[0];
rz(-2.0656842) q[0];
sx q[0];
rz(2.8450235) q[0];
rz(2.2462559) q[2];
sx q[2];
rz(-0.5731715) q[2];
sx q[2];
rz(-1.2473388) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8365337) q[1];
sx q[1];
rz(-1.1593474) q[1];
sx q[1];
rz(-1.9484503) q[1];
x q[2];
rz(-1.7612475) q[3];
sx q[3];
rz(-0.69882876) q[3];
sx q[3];
rz(0.33635283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6707637) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(1.0864351) q[2];
rz(-0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(-0.72211784) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54871854) q[0];
sx q[0];
rz(-1.6328136) q[0];
sx q[0];
rz(1.6528224) q[0];
rz(-1.2552352) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(1.9699388) q[2];
sx q[2];
rz(-1.4831586) q[2];
sx q[2];
rz(-3.0638564) q[2];
rz(-0.063379824) q[3];
sx q[3];
rz(-2.2324003) q[3];
sx q[3];
rz(1.5008055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
