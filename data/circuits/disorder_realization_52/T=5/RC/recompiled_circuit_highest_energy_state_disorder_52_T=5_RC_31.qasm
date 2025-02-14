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
rz(1.7088543) q[0];
sx q[0];
rz(-2.813485) q[0];
sx q[0];
rz(-0.83972591) q[0];
rz(-1.7307164) q[1];
sx q[1];
rz(-1.6003992) q[1];
sx q[1];
rz(2.3720001) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037207863) q[0];
sx q[0];
rz(-1.5360556) q[0];
sx q[0];
rz(1.6503235) q[0];
x q[1];
rz(-3.1159471) q[2];
sx q[2];
rz(-0.8249976) q[2];
sx q[2];
rz(-0.96482575) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.855866) q[1];
sx q[1];
rz(-1.2190423) q[1];
sx q[1];
rz(0.00067623059) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0528209) q[3];
sx q[3];
rz(-2.149579) q[3];
sx q[3];
rz(2.5797648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8771693) q[2];
sx q[2];
rz(-1.8921655) q[2];
sx q[2];
rz(-1.9682311) q[2];
rz(2.7133283) q[3];
sx q[3];
rz(-0.56848017) q[3];
sx q[3];
rz(-2.4885524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0083171) q[0];
sx q[0];
rz(-0.44047099) q[0];
sx q[0];
rz(-2.0793197) q[0];
rz(-1.590689) q[1];
sx q[1];
rz(-2.5804602) q[1];
sx q[1];
rz(-1.0947469) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2424205) q[0];
sx q[0];
rz(-1.3160442) q[0];
sx q[0];
rz(3.0887449) q[0];
rz(-pi) q[1];
rz(0.44610141) q[2];
sx q[2];
rz(-2.1992536) q[2];
sx q[2];
rz(2.7035509) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9506663) q[1];
sx q[1];
rz(-2.6537173) q[1];
sx q[1];
rz(2.4937954) q[1];
rz(-pi) q[2];
rz(-2.3660819) q[3];
sx q[3];
rz(-1.417096) q[3];
sx q[3];
rz(2.6284144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4636479) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(2.2410683) q[2];
rz(1.0362222) q[3];
sx q[3];
rz(-3.0039054) q[3];
sx q[3];
rz(0.010995939) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4765332) q[0];
sx q[0];
rz(-1.1256951) q[0];
sx q[0];
rz(-1.4680468) q[0];
rz(2.2131069) q[1];
sx q[1];
rz(-1.0037582) q[1];
sx q[1];
rz(-1.7195864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4741459) q[0];
sx q[0];
rz(-2.8774539) q[0];
sx q[0];
rz(0.35281065) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68520968) q[2];
sx q[2];
rz(-2.1985719) q[2];
sx q[2];
rz(-0.9762828) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8633008) q[1];
sx q[1];
rz(-1.2108786) q[1];
sx q[1];
rz(2.5665096) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5676651) q[3];
sx q[3];
rz(-0.49387723) q[3];
sx q[3];
rz(-2.9536794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15098393) q[2];
sx q[2];
rz(-1.432212) q[2];
sx q[2];
rz(2.3217311) q[2];
rz(0.12623434) q[3];
sx q[3];
rz(-1.1773959) q[3];
sx q[3];
rz(1.4950776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4594629) q[0];
sx q[0];
rz(-1.7403025) q[0];
sx q[0];
rz(-1.6023741) q[0];
rz(1.0914717) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.1702671) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4226886) q[0];
sx q[0];
rz(-0.50245521) q[0];
sx q[0];
rz(-0.86645856) q[0];
x q[1];
rz(1.2762047) q[2];
sx q[2];
rz(-0.74199332) q[2];
sx q[2];
rz(-2.2592253) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5398769) q[1];
sx q[1];
rz(-1.4810432) q[1];
sx q[1];
rz(1.490277) q[1];
x q[2];
rz(2.7682414) q[3];
sx q[3];
rz(-1.3697094) q[3];
sx q[3];
rz(-1.1652077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.22152659) q[2];
sx q[2];
rz(-2.2130794) q[2];
sx q[2];
rz(-3.0583734) q[2];
rz(-0.65256882) q[3];
sx q[3];
rz(-0.021952732) q[3];
sx q[3];
rz(0.86161247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6457152) q[0];
sx q[0];
rz(-0.28384122) q[0];
sx q[0];
rz(-0.19677095) q[0];
rz(-1.1605877) q[1];
sx q[1];
rz(-0.44191688) q[1];
sx q[1];
rz(1.7534076) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3144238) q[0];
sx q[0];
rz(-1.4008796) q[0];
sx q[0];
rz(1.569237) q[0];
rz(0.6394425) q[2];
sx q[2];
rz(-2.8911984) q[2];
sx q[2];
rz(-2.3441803) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1886471) q[1];
sx q[1];
rz(-1.7428696) q[1];
sx q[1];
rz(-2.2285372) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84204414) q[3];
sx q[3];
rz(-1.2938515) q[3];
sx q[3];
rz(-2.7634948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56472003) q[2];
sx q[2];
rz(-1.9331965) q[2];
sx q[2];
rz(1.2811071) q[2];
rz(0.34230226) q[3];
sx q[3];
rz(-0.050693158) q[3];
sx q[3];
rz(2.2127693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36011919) q[0];
sx q[0];
rz(-2.0442648) q[0];
sx q[0];
rz(1.0327562) q[0];
rz(0.67689854) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(0.78183174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.483455) q[0];
sx q[0];
rz(-1.0080604) q[0];
sx q[0];
rz(-1.9938009) q[0];
rz(1.621945) q[2];
sx q[2];
rz(-1.5384073) q[2];
sx q[2];
rz(-2.6408104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1429174) q[1];
sx q[1];
rz(-1.8774722) q[1];
sx q[1];
rz(-0.59345133) q[1];
rz(-pi) q[2];
rz(-1.6138462) q[3];
sx q[3];
rz(-1.2331748) q[3];
sx q[3];
rz(1.4211602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3013762) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(2.2678579) q[2];
rz(2.3524763) q[3];
sx q[3];
rz(-1.5566166) q[3];
sx q[3];
rz(-2.6553787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9588722) q[0];
sx q[0];
rz(-3.0954269) q[0];
sx q[0];
rz(-0.090855457) q[0];
rz(-0.60984045) q[1];
sx q[1];
rz(-1.5900541) q[1];
sx q[1];
rz(2.5441433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4700902) q[0];
sx q[0];
rz(-1.4250942) q[0];
sx q[0];
rz(1.7007366) q[0];
x q[1];
rz(0.84709527) q[2];
sx q[2];
rz(-2.6604746) q[2];
sx q[2];
rz(-1.4916949) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56126334) q[1];
sx q[1];
rz(-2.4819907) q[1];
sx q[1];
rz(0.81239169) q[1];
rz(-1.6555637) q[3];
sx q[3];
rz(-2.7206507) q[3];
sx q[3];
rz(1.3011025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37683836) q[2];
sx q[2];
rz(-0.54543269) q[2];
sx q[2];
rz(3.0010014) q[2];
rz(-1.8600672) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(2.6144821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.918688) q[0];
sx q[0];
rz(-2.1253026) q[0];
sx q[0];
rz(1.5284398) q[0];
rz(2.3279922) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(2.334107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19014004) q[0];
sx q[0];
rz(-1.7684816) q[0];
sx q[0];
rz(0.19301599) q[0];
x q[1];
rz(-0.1296223) q[2];
sx q[2];
rz(-1.1309012) q[2];
sx q[2];
rz(2.9578046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85881104) q[1];
sx q[1];
rz(-0.75560299) q[1];
sx q[1];
rz(-2.1367461) q[1];
rz(2.7164338) q[3];
sx q[3];
rz(-1.0711373) q[3];
sx q[3];
rz(0.560958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7878824) q[2];
sx q[2];
rz(-1.3545802) q[2];
sx q[2];
rz(-0.14459571) q[2];
rz(1.1041798) q[3];
sx q[3];
rz(-0.2140597) q[3];
sx q[3];
rz(-1.0415227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60177326) q[0];
sx q[0];
rz(-2.0105392) q[0];
sx q[0];
rz(-2.5850776) q[0];
rz(-2.3717608) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(2.6803023) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78019136) q[0];
sx q[0];
rz(-2.5690329) q[0];
sx q[0];
rz(1.4828234) q[0];
rz(0.67619063) q[2];
sx q[2];
rz(-2.1874223) q[2];
sx q[2];
rz(-0.21756324) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3756421) q[1];
sx q[1];
rz(-1.8643171) q[1];
sx q[1];
rz(0.75325812) q[1];
rz(-pi) q[2];
rz(2.5784745) q[3];
sx q[3];
rz(-0.65059911) q[3];
sx q[3];
rz(1.7312778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7028659) q[2];
sx q[2];
rz(-0.30539572) q[2];
sx q[2];
rz(-0.76496441) q[2];
rz(1.9276098) q[3];
sx q[3];
rz(-0.58911222) q[3];
sx q[3];
rz(-2.7629619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42534378) q[0];
sx q[0];
rz(-0.049976293) q[0];
sx q[0];
rz(2.7783527) q[0];
rz(1.4092457) q[1];
sx q[1];
rz(-0.98594085) q[1];
sx q[1];
rz(-2.9149616) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55174819) q[0];
sx q[0];
rz(-1.4355816) q[0];
sx q[0];
rz(-1.5677794) q[0];
rz(-pi) q[1];
rz(2.4754074) q[2];
sx q[2];
rz(-2.2720784) q[2];
sx q[2];
rz(-1.9124857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5635101) q[1];
sx q[1];
rz(-2.3704766) q[1];
sx q[1];
rz(1.095849) q[1];
rz(-pi) q[2];
rz(-0.33785401) q[3];
sx q[3];
rz(-1.5607335) q[3];
sx q[3];
rz(-0.73051605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54214415) q[2];
sx q[2];
rz(-0.52079529) q[2];
sx q[2];
rz(-0.64860541) q[2];
rz(0.058622807) q[3];
sx q[3];
rz(-2.5128745) q[3];
sx q[3];
rz(0.64529836) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7764353) q[0];
sx q[0];
rz(-1.9226274) q[0];
sx q[0];
rz(2.6489039) q[0];
rz(-2.0877214) q[1];
sx q[1];
rz(-2.5851879) q[1];
sx q[1];
rz(-2.7256706) q[1];
rz(3.1171534) q[2];
sx q[2];
rz(-1.3346439) q[2];
sx q[2];
rz(1.1572184) q[2];
rz(0.875474) q[3];
sx q[3];
rz(-1.8322104) q[3];
sx q[3];
rz(-1.9598243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
