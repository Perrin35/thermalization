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
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(-2.1767148) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(-2.092195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8937738) q[0];
sx q[0];
rz(-0.6914833) q[0];
sx q[0];
rz(2.1201503) q[0];
x q[1];
rz(1.1683091) q[2];
sx q[2];
rz(-2.7912729) q[2];
sx q[2];
rz(2.868319) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8559554) q[1];
sx q[1];
rz(-0.95889839) q[1];
sx q[1];
rz(-0.4680856) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4780294) q[3];
sx q[3];
rz(-2.1481107) q[3];
sx q[3];
rz(-0.64914671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6875978) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(-3.122984) q[2];
rz(0.54443693) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(-0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3675156) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(0.53502214) q[0];
rz(0.53994838) q[1];
sx q[1];
rz(-0.63553634) q[1];
sx q[1];
rz(-1.3998869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0947123) q[0];
sx q[0];
rz(-2.0167354) q[0];
sx q[0];
rz(-2.711854) q[0];
rz(2.1715047) q[2];
sx q[2];
rz(-2.0510489) q[2];
sx q[2];
rz(2.9532972) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3532039) q[1];
sx q[1];
rz(-1.8192768) q[1];
sx q[1];
rz(2.7491991) q[1];
rz(-pi) q[2];
rz(2.4055491) q[3];
sx q[3];
rz(-0.76220817) q[3];
sx q[3];
rz(-2.8413642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(0.92607099) q[2];
rz(2.1615084) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(-2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3493018) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(1.6287623) q[0];
rz(0.28383645) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(-1.1121174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84217269) q[0];
sx q[0];
rz(-1.5779691) q[0];
sx q[0];
rz(-3.1345063) q[0];
rz(-pi) q[1];
rz(3.0653333) q[2];
sx q[2];
rz(-1.6789376) q[2];
sx q[2];
rz(-2.0881483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97858799) q[1];
sx q[1];
rz(-1.3842062) q[1];
sx q[1];
rz(0.48689894) q[1];
x q[2];
rz(-0.25896163) q[3];
sx q[3];
rz(-1.2185214) q[3];
sx q[3];
rz(-1.5332298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(-1.1019361) q[2];
rz(-2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(-2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717644) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(-2.4523822) q[0];
rz(0.31461942) q[1];
sx q[1];
rz(-2.2210821) q[1];
sx q[1];
rz(1.7291732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2570626) q[0];
sx q[0];
rz(-1.2572968) q[0];
sx q[0];
rz(-2.351981) q[0];
x q[1];
rz(0.41823776) q[2];
sx q[2];
rz(-0.23242885) q[2];
sx q[2];
rz(1.4168036) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3312348) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(0.15366252) q[1];
rz(-pi) q[2];
rz(-2.9474218) q[3];
sx q[3];
rz(-1.0946858) q[3];
sx q[3];
rz(0.74060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41669258) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(2.5856384) q[2];
rz(-1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3274662) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(-2.5471174) q[0];
rz(2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-2.1655703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9851889) q[0];
sx q[0];
rz(-1.9984584) q[0];
sx q[0];
rz(-2.7424863) q[0];
rz(-pi) q[1];
rz(2.3875434) q[2];
sx q[2];
rz(-1.5496786) q[2];
sx q[2];
rz(-1.6789951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12078602) q[1];
sx q[1];
rz(-0.73985064) q[1];
sx q[1];
rz(0.45005656) q[1];
rz(-pi) q[2];
rz(-0.61686744) q[3];
sx q[3];
rz(-0.36412334) q[3];
sx q[3];
rz(-1.1222276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6550265) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(2.5308934) q[2];
rz(0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(-1.8122199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82764757) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(1.4661283) q[0];
rz(-2.3620391) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-0.73878845) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8429564) q[0];
sx q[0];
rz(-1.8514575) q[0];
sx q[0];
rz(0.21224169) q[0];
rz(-pi) q[1];
rz(-2.8347978) q[2];
sx q[2];
rz(-1.5606797) q[2];
sx q[2];
rz(0.9900569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93529451) q[1];
sx q[1];
rz(-1.9351462) q[1];
sx q[1];
rz(0.97779001) q[1];
rz(-pi) q[2];
rz(-1.1973279) q[3];
sx q[3];
rz(-1.9091144) q[3];
sx q[3];
rz(-0.23972971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5256727) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(0.45977965) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(-1.1427243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90726844) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(-1.020485) q[0];
rz(3.0842969) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(-2.3540156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4754935) q[0];
sx q[0];
rz(-1.3331873) q[0];
sx q[0];
rz(-0.65639021) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26489218) q[2];
sx q[2];
rz(-2.8745626) q[2];
sx q[2];
rz(1.9445813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.139442) q[1];
sx q[1];
rz(-2.6137335) q[1];
sx q[1];
rz(1.9308286) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5191139) q[3];
sx q[3];
rz(-0.58850901) q[3];
sx q[3];
rz(-1.9543598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3202177) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(-2.5569432) q[2];
rz(0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(2.8133494) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.1513386) q[1];
sx q[1];
rz(-1.4788871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23173702) q[0];
sx q[0];
rz(-1.3074991) q[0];
sx q[0];
rz(-1.3366827) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4166862) q[2];
sx q[2];
rz(-0.6577684) q[2];
sx q[2];
rz(-0.36054128) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11301576) q[1];
sx q[1];
rz(-1.2927107) q[1];
sx q[1];
rz(0.098829513) q[1];
rz(-2.9313179) q[3];
sx q[3];
rz(-1.2089024) q[3];
sx q[3];
rz(2.0023605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(-0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64860827) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(-0.96187821) q[0];
rz(0.23421639) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-2.403517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80958145) q[0];
sx q[0];
rz(-1.1763078) q[0];
sx q[0];
rz(-1.8277192) q[0];
rz(1.1253276) q[2];
sx q[2];
rz(-2.1574852) q[2];
sx q[2];
rz(2.4054804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2710423) q[1];
sx q[1];
rz(-1.9890474) q[1];
sx q[1];
rz(-0.31914945) q[1];
rz(-pi) q[2];
rz(0.35154147) q[3];
sx q[3];
rz(-2.4443887) q[3];
sx q[3];
rz(-0.54920025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48478475) q[2];
sx q[2];
rz(-2.1190376) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(-1.8079181) q[3];
sx q[3];
rz(-2.1402054) q[3];
sx q[3];
rz(2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.242908) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(0.86724487) q[0];
rz(-1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(1.0702466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9673706) q[0];
sx q[0];
rz(-1.4486607) q[0];
sx q[0];
rz(-1.723126) q[0];
x q[1];
rz(0.82735586) q[2];
sx q[2];
rz(-2.1340886) q[2];
sx q[2];
rz(2.2652596) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6079191) q[1];
sx q[1];
rz(-1.9551139) q[1];
sx q[1];
rz(-0.55748765) q[1];
rz(-pi) q[2];
rz(3.002043) q[3];
sx q[3];
rz(-1.4613536) q[3];
sx q[3];
rz(-2.8982996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16537198) q[2];
sx q[2];
rz(-0.82087159) q[2];
sx q[2];
rz(1.2316068) q[2];
rz(1.6407137) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-2.4098082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-0.41481836) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(-2.6957569) q[2];
sx q[2];
rz(-1.6204964) q[2];
sx q[2];
rz(1.3514018) q[2];
rz(0.75178643) q[3];
sx q[3];
rz(-1.1204168) q[3];
sx q[3];
rz(0.71747019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
