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
rz(-0.14970782) q[1];
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
rz(-0.61650886) q[0];
x q[1];
rz(-2.1120464) q[2];
sx q[2];
rz(-2.0143348) q[2];
sx q[2];
rz(0.29543791) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.522361) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(0.038832263) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2655067) q[3];
sx q[3];
rz(-2.5698834) q[3];
sx q[3];
rz(-1.7892464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-0.70409888) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(0.44644341) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7783633) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(2.7345783) q[0];
rz(-pi) q[1];
rz(2.9035283) q[2];
sx q[2];
rz(-1.8556343) q[2];
sx q[2];
rz(-2.970398) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0374239) q[1];
sx q[1];
rz(-1.8579322) q[1];
sx q[1];
rz(-1.3377405) q[1];
rz(-pi) q[2];
rz(-0.84516256) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(-1.15796) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.3495548) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2230167) q[0];
sx q[0];
rz(-0.71764676) q[0];
sx q[0];
rz(-0.36803228) q[0];
rz(-1.0565874) q[2];
sx q[2];
rz(-0.075767013) q[2];
sx q[2];
rz(2.5839992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.459356) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(-0.45046803) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9968824) q[3];
sx q[3];
rz(-1.5618556) q[3];
sx q[3];
rz(0.63250354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(-0.036380336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.13658) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.9146634) q[0];
x q[1];
rz(0.47469791) q[2];
sx q[2];
rz(-0.63650741) q[2];
sx q[2];
rz(2.0091025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2991043) q[1];
sx q[1];
rz(-0.71055382) q[1];
sx q[1];
rz(1.3683661) q[1];
rz(-pi) q[2];
rz(2.7500238) q[3];
sx q[3];
rz(-2.8885926) q[3];
sx q[3];
rz(2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(-1.9533763) q[2];
rz(2.4711117) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-0.23434815) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74041022) q[0];
sx q[0];
rz(-0.64490841) q[0];
sx q[0];
rz(0.58437225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0767897) q[2];
sx q[2];
rz(-1.6550118) q[2];
sx q[2];
rz(-0.46352026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.25917945) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(1.773136) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3393199) q[3];
sx q[3];
rz(-2.8095062) q[3];
sx q[3];
rz(2.092923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-0.22053545) q[2];
rz(-2.7045414) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(1.1436499) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4139347) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(-1.5278221) q[0];
rz(-pi) q[1];
rz(1.3104865) q[2];
sx q[2];
rz(-0.91432768) q[2];
sx q[2];
rz(1.7973763) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5551344) q[1];
sx q[1];
rz(-2.0124334) q[1];
sx q[1];
rz(-0.099979062) q[1];
rz(1.2915217) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(2.9135381) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.2868767) q[0];
rz(-0.6634179) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(-1.9082665) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9702643) q[0];
sx q[0];
rz(-1.5597222) q[0];
sx q[0];
rz(1.1867255) q[0];
x q[1];
rz(-2.2912824) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(0.26956272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8334956) q[1];
sx q[1];
rz(-1.4596246) q[1];
sx q[1];
rz(1.5566467) q[1];
rz(-pi) q[2];
rz(2.535378) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-1.0160149) q[2];
rz(-0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.4455147) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.258237) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94851516) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(1.8817188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4969205) q[2];
sx q[2];
rz(-2.1053227) q[2];
sx q[2];
rz(1.0169741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80553493) q[1];
sx q[1];
rz(-1.6748669) q[1];
sx q[1];
rz(2.59471) q[1];
rz(-pi) q[2];
rz(2.6120841) q[3];
sx q[3];
rz(-0.92748517) q[3];
sx q[3];
rz(-0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68226472) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(-2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(-2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(0.2510221) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(3.1138611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160307) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-0.82657878) q[0];
x q[1];
rz(-2.8743923) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(-2.6097678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2003277) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(-1.6914678) q[1];
rz(-1.6373709) q[3];
sx q[3];
rz(-1.5449617) q[3];
sx q[3];
rz(1.0915826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-0.25751105) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8513545) q[0];
sx q[0];
rz(-0.60831735) q[0];
sx q[0];
rz(-0.71382199) q[0];
rz(0.72980482) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(2.6953816) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7911437) q[1];
sx q[1];
rz(-2.2669683) q[1];
sx q[1];
rz(-1.5894366) q[1];
rz(1.0912283) q[3];
sx q[3];
rz(-1.7864831) q[3];
sx q[3];
rz(1.8857764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(0.47484067) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-2.6735641) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(2.0444617) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
