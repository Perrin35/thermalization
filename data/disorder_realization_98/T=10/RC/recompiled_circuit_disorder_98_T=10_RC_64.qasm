OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(2.9918848) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5877085) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(2.5250838) q[0];
rz(-pi) q[1];
rz(-2.1120464) q[2];
sx q[2];
rz(-2.0143348) q[2];
sx q[2];
rz(-2.8461547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7357199) q[1];
sx q[1];
rz(-2.957203) q[1];
sx q[1];
rz(1.7806446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19102328) q[3];
sx q[3];
rz(-2.1130307) q[3];
sx q[3];
rz(2.1477826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(-0.63252226) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1855542) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(1.5831328) q[0];
rz(-pi) q[1];
rz(-0.89276887) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(-0.54112753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3995041) q[1];
sx q[1];
rz(-1.794145) q[1];
sx q[1];
rz(-2.8469267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90918031) q[3];
sx q[3];
rz(-1.0717857) q[3];
sx q[3];
rz(-3.1391075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.9799505) q[2];
rz(-1.15796) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.3495548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7499381) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(1.8750989) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6367958) q[2];
sx q[2];
rz(-1.608035) q[2];
sx q[2];
rz(-0.50022349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6408491) q[1];
sx q[1];
rz(-1.3355681) q[1];
sx q[1];
rz(2.6231223) q[1];
rz(-pi) q[2];
rz(1.549167) q[3];
sx q[3];
rz(-0.4261741) q[3];
sx q[3];
rz(-2.222995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(2.5097805) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(-3.1052123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2934389) q[0];
sx q[0];
rz(-1.358195) q[0];
sx q[0];
rz(-0.92377499) q[0];
x q[1];
rz(-0.58156275) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(0.046422596) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42602793) q[1];
sx q[1];
rz(-1.4392816) q[1];
sx q[1];
rz(2.2711666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23457228) q[3];
sx q[3];
rz(-1.6664701) q[3];
sx q[3];
rz(-2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.9533763) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7017512) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(0.23434815) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4317961) q[0];
sx q[0];
rz(-1.0456107) q[0];
sx q[0];
rz(1.1774506) q[0];
x q[1];
rz(1.3946103) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(-1.2622152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4489331) q[1];
sx q[1];
rz(-1.4217136) q[1];
sx q[1];
rz(-0.74933185) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0626416) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(0.80432804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691539) q[0];
sx q[0];
rz(-1.5307431) q[0];
sx q[0];
rz(-0.37102951) q[0];
x q[1];
rz(1.8311062) q[2];
sx q[2];
rz(-0.91432768) q[2];
sx q[2];
rz(1.3442163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3559349) q[1];
sx q[1];
rz(-2.6895084) q[1];
sx q[1];
rz(-1.362734) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2915217) q[3];
sx q[3];
rz(-1.0696628) q[3];
sx q[3];
rz(0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(2.6521519) q[2];
rz(-0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5722826) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.2868767) q[0];
rz(0.6634179) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.2333262) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7147303) q[0];
sx q[0];
rz(-2.7573702) q[0];
sx q[0];
rz(-1.5412488) q[0];
x q[1];
rz(-2.2592505) q[2];
sx q[2];
rz(-0.8493648) q[2];
sx q[2];
rz(-1.7989858) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26112939) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(3.0304099) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4864278) q[3];
sx q[3];
rz(-1.8842116) q[3];
sx q[3];
rz(-3.0048971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90342605) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(-1.1358791) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7779185) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(-1.9922436) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82842365) q[1];
sx q[1];
rz(-2.1143882) q[1];
sx q[1];
rz(1.4491175) q[1];
rz(2.2860252) q[3];
sx q[3];
rz(-1.1547525) q[3];
sx q[3];
rz(0.95203979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(-2.941926) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(-0.02773157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32556191) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-0.82657878) q[0];
rz(-1.2559782) q[2];
sx q[2];
rz(-2.417649) q[2];
sx q[2];
rz(2.1997423) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30222826) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(-0.60217963) q[1];
rz(1.6373709) q[3];
sx q[3];
rz(-1.5449617) q[3];
sx q[3];
rz(2.0500101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(2.8840816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89850241) q[0];
sx q[0];
rz(-1.1872963) q[0];
sx q[0];
rz(2.6570508) q[0];
rz(2.4117878) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(0.44621106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7911437) q[1];
sx q[1];
rz(-2.2669683) q[1];
sx q[1];
rz(-1.552156) q[1];
rz(-pi) q[2];
rz(-1.1274687) q[3];
sx q[3];
rz(-2.6192198) q[3];
sx q[3];
rz(-0.70538196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-0.46802855) q[2];
sx q[2];
rz(-1.9149018) q[2];
sx q[2];
rz(0.55185774) q[2];
rz(-2.0541035) q[3];
sx q[3];
rz(-1.8046422) q[3];
sx q[3];
rz(1.3840152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];