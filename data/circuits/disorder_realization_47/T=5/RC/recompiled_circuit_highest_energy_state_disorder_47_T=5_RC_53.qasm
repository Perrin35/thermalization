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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(1.6562847) q[0];
rz(1.0506884) q[1];
sx q[1];
rz(-1.7424072) q[1];
sx q[1];
rz(2.4032226) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0169058) q[0];
sx q[0];
rz(-1.4198185) q[0];
sx q[0];
rz(-0.24589234) q[0];
x q[1];
rz(-1.1777056) q[2];
sx q[2];
rz(-1.4294659) q[2];
sx q[2];
rz(1.7766952) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0289606) q[1];
sx q[1];
rz(-0.46507177) q[1];
sx q[1];
rz(-2.8264224) q[1];
rz(-pi) q[2];
rz(-2.5220736) q[3];
sx q[3];
rz(-1.3196919) q[3];
sx q[3];
rz(2.8015603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(-2.7281249) q[2];
rz(0.56420285) q[3];
sx q[3];
rz(-1.4671624) q[3];
sx q[3];
rz(1.9570501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.41459945) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(-0.10398908) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-2.7242421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39825059) q[0];
sx q[0];
rz(-0.62055991) q[0];
sx q[0];
rz(-2.1055566) q[0];
rz(-3.1293489) q[2];
sx q[2];
rz(-1.2070388) q[2];
sx q[2];
rz(-3.0573483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88637832) q[1];
sx q[1];
rz(-0.30579771) q[1];
sx q[1];
rz(1.7096108) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52882282) q[3];
sx q[3];
rz(-1.0695056) q[3];
sx q[3];
rz(-2.5426585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.523681) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(-2.9376612) q[2];
rz(-1.6265053) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(1.6319857) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6025036) q[0];
sx q[0];
rz(-1.6747549) q[0];
sx q[0];
rz(2.7432826) q[0];
rz(-1.0771982) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(2.5881252) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.80478) q[0];
sx q[0];
rz(-0.83641988) q[0];
sx q[0];
rz(-1.5848716) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4301694) q[2];
sx q[2];
rz(-1.4069948) q[2];
sx q[2];
rz(-1.8970416) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7872214) q[1];
sx q[1];
rz(-2.1408892) q[1];
sx q[1];
rz(1.1524423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3970277) q[3];
sx q[3];
rz(-2.8855763) q[3];
sx q[3];
rz(2.721173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.074097721) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(1.2158166) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(1.2522987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85553402) q[0];
sx q[0];
rz(-1.0581886) q[0];
sx q[0];
rz(-1.4372987) q[0];
rz(-1.3358491) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(-0.45509532) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44667654) q[0];
sx q[0];
rz(-1.4908264) q[0];
sx q[0];
rz(0.6178426) q[0];
rz(-1.9432151) q[2];
sx q[2];
rz(-1.7657868) q[2];
sx q[2];
rz(-2.7253828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3434717) q[1];
sx q[1];
rz(-1.2373072) q[1];
sx q[1];
rz(-1.8300959) q[1];
rz(2.9330071) q[3];
sx q[3];
rz(-1.6634395) q[3];
sx q[3];
rz(-1.1535597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(0.8320128) q[2];
rz(0.653382) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117383) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(-2.4526556) q[0];
rz(2.9298933) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(0.45417085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9849761) q[0];
sx q[0];
rz(-2.1097221) q[0];
sx q[0];
rz(2.6205553) q[0];
rz(-pi) q[1];
rz(2.0549704) q[2];
sx q[2];
rz(-2.4459029) q[2];
sx q[2];
rz(-1.5351968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.28605697) q[1];
sx q[1];
rz(-2.2814676) q[1];
sx q[1];
rz(-0.90531207) q[1];
rz(-pi) q[2];
rz(1.0096142) q[3];
sx q[3];
rz(-2.8343763) q[3];
sx q[3];
rz(-3.0209783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5357369) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(1.2153252) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(-0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-3.0257492) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(0.51344839) q[0];
rz(2.2250941) q[1];
sx q[1];
rz(-1.1767574) q[1];
sx q[1];
rz(1.7880012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3037864) q[0];
sx q[0];
rz(-0.45188658) q[0];
sx q[0];
rz(-2.9511098) q[0];
rz(0.84557477) q[2];
sx q[2];
rz(-2.6161199) q[2];
sx q[2];
rz(0.89286823) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1648851) q[1];
sx q[1];
rz(-1.083255) q[1];
sx q[1];
rz(-1.6009925) q[1];
rz(-pi) q[2];
rz(-1.6597802) q[3];
sx q[3];
rz(-1.5060194) q[3];
sx q[3];
rz(-1.806206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6650271) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(3.108007) q[2];
rz(-2.636886) q[3];
sx q[3];
rz(-1.7028156) q[3];
sx q[3];
rz(2.3869042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018709239) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(1.1705742) q[0];
rz(-2.9508044) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(0.64291397) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87827728) q[0];
sx q[0];
rz(-1.3049135) q[0];
sx q[0];
rz(-0.10590597) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.229081) q[2];
sx q[2];
rz(-0.40900074) q[2];
sx q[2];
rz(-0.3106948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9289405) q[1];
sx q[1];
rz(-1.8096605) q[1];
sx q[1];
rz(-2.8293306) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6899818) q[3];
sx q[3];
rz(-1.168712) q[3];
sx q[3];
rz(-1.6329671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99870318) q[2];
sx q[2];
rz(-3.0316752) q[2];
sx q[2];
rz(1.9996803) q[2];
rz(0.029684639) q[3];
sx q[3];
rz(-1.6073062) q[3];
sx q[3];
rz(-0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2340045) q[0];
sx q[0];
rz(-1.1827844) q[0];
sx q[0];
rz(1.3080066) q[0];
rz(-1.1169149) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(0.21211472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.708688) q[0];
sx q[0];
rz(-1.8861265) q[0];
sx q[0];
rz(-0.84724119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4918152) q[2];
sx q[2];
rz(-0.37576518) q[2];
sx q[2];
rz(0.77451578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37112889) q[1];
sx q[1];
rz(-1.1327289) q[1];
sx q[1];
rz(0.99645946) q[1];
rz(-0.85101012) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(1.2860166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1954631) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(1.2809666) q[2];
rz(-1.7383176) q[3];
sx q[3];
rz(-2.2303228) q[3];
sx q[3];
rz(-2.4173071) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4487149) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(-2.2591059) q[0];
rz(-2.8485883) q[1];
sx q[1];
rz(-2.2182783) q[1];
sx q[1];
rz(2.015347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6544271) q[0];
sx q[0];
rz(-1.2058655) q[0];
sx q[0];
rz(-0.69661822) q[0];
rz(-pi) q[1];
rz(-0.0056929767) q[2];
sx q[2];
rz(-2.5416592) q[2];
sx q[2];
rz(-1.1545187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4418151) q[1];
sx q[1];
rz(-0.80301563) q[1];
sx q[1];
rz(0.074549874) q[1];
x q[2];
rz(2.7357941) q[3];
sx q[3];
rz(-1.278459) q[3];
sx q[3];
rz(0.97989156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5249411) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(-0.16723995) q[2];
rz(3.1104769) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(-0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50877082) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(-0.81047812) q[0];
rz(-1.8327911) q[1];
sx q[1];
rz(-2.2056613) q[1];
sx q[1];
rz(0.097361758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7435849) q[0];
sx q[0];
rz(-2.3438489) q[0];
sx q[0];
rz(-1.8037075) q[0];
x q[1];
rz(1.1896336) q[2];
sx q[2];
rz(-1.544017) q[2];
sx q[2];
rz(-1.3285411) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88579455) q[1];
sx q[1];
rz(-0.58572873) q[1];
sx q[1];
rz(0.08759193) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11730365) q[3];
sx q[3];
rz(-1.3051093) q[3];
sx q[3];
rz(2.9151288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8436766) q[2];
sx q[2];
rz(-2.127779) q[2];
sx q[2];
rz(-1.8664912) q[2];
rz(-0.44025931) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(0.27537235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-0.62200017) q[1];
sx q[1];
rz(-1.2983464) q[1];
sx q[1];
rz(1.3141528) q[1];
rz(-0.95876454) q[2];
sx q[2];
rz(-0.13560451) q[2];
sx q[2];
rz(2.3252631) q[2];
rz(-2.4980656) q[3];
sx q[3];
rz(-1.4326) q[3];
sx q[3];
rz(2.8958733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
