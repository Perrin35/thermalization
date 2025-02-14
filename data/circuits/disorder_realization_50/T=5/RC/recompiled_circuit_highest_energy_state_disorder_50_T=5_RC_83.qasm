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
rz(0.69831508) q[0];
sx q[0];
rz(-0.3787711) q[0];
sx q[0];
rz(-1.9842499) q[0];
rz(-2.3565489) q[1];
sx q[1];
rz(-2.4554689) q[1];
sx q[1];
rz(0.52678984) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3279545) q[0];
sx q[0];
rz(-2.1305973) q[0];
sx q[0];
rz(0.35388057) q[0];
rz(-pi) q[1];
rz(-0.77157048) q[2];
sx q[2];
rz(-0.89961806) q[2];
sx q[2];
rz(-3.0085466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4545645) q[1];
sx q[1];
rz(-0.44378456) q[1];
sx q[1];
rz(3.0566178) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.029954682) q[3];
sx q[3];
rz(-2.1308793) q[3];
sx q[3];
rz(-2.4966979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36246768) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(2.989952) q[2];
rz(-0.14196299) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-1.1375455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66459429) q[0];
sx q[0];
rz(-0.73200309) q[0];
sx q[0];
rz(0.81749302) q[0];
rz(0.11257653) q[1];
sx q[1];
rz(-1.45603) q[1];
sx q[1];
rz(-0.89961019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0257638) q[0];
sx q[0];
rz(-1.2983822) q[0];
sx q[0];
rz(-1.8194356) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4075872) q[2];
sx q[2];
rz(-2.3235934) q[2];
sx q[2];
rz(-0.98246511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4790597) q[1];
sx q[1];
rz(-0.823349) q[1];
sx q[1];
rz(1.6475299) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95590653) q[3];
sx q[3];
rz(-2.3645176) q[3];
sx q[3];
rz(1.9333378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4704935) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(1.0279083) q[2];
rz(-2.4043064) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(-0.024287311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1360433) q[0];
sx q[0];
rz(-0.5558973) q[0];
sx q[0];
rz(-1.1935724) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.4793414) q[1];
sx q[1];
rz(1.4847635) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45423266) q[0];
sx q[0];
rz(-3.0890565) q[0];
sx q[0];
rz(2.0665069) q[0];
x q[1];
rz(-1.7822927) q[2];
sx q[2];
rz(-1.4471608) q[2];
sx q[2];
rz(2.5619626) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72276593) q[1];
sx q[1];
rz(-1.9120941) q[1];
sx q[1];
rz(-1.7143634) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3562702) q[3];
sx q[3];
rz(-0.37060043) q[3];
sx q[3];
rz(1.4166946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.16126157) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(-0.70183357) q[2];
rz(-0.069132239) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0858916) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(2.5162146) q[0];
rz(2.0471795) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(-1.8249003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0023481) q[0];
sx q[0];
rz(-1.7138094) q[0];
sx q[0];
rz(1.0225747) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6676297) q[2];
sx q[2];
rz(-1.2487354) q[2];
sx q[2];
rz(1.2351954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91176468) q[1];
sx q[1];
rz(-2.266464) q[1];
sx q[1];
rz(-1.8063481) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3415001) q[3];
sx q[3];
rz(-2.0043819) q[3];
sx q[3];
rz(0.80899245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(-0.030108062) q[2];
rz(2.7249469) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(-0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3629214) q[0];
sx q[0];
rz(-1.2782949) q[0];
sx q[0];
rz(1.9238506) q[0];
rz(-2.9171004) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(-0.40103689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.359085) q[0];
sx q[0];
rz(-1.1993202) q[0];
sx q[0];
rz(2.0187733) q[0];
x q[1];
rz(0.46868344) q[2];
sx q[2];
rz(-2.4309553) q[2];
sx q[2];
rz(2.8993894) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89304533) q[1];
sx q[1];
rz(-1.6181989) q[1];
sx q[1];
rz(-0.23706146) q[1];
x q[2];
rz(1.4698074) q[3];
sx q[3];
rz(-1.4091531) q[3];
sx q[3];
rz(-0.18835959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2991221) q[2];
sx q[2];
rz(-1.9186019) q[2];
sx q[2];
rz(-2.6341338) q[2];
rz(2.2211645) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(-0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64666635) q[0];
sx q[0];
rz(-0.05519069) q[0];
sx q[0];
rz(-2.1221509) q[0];
rz(1.9693718) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(-1.9196462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86496124) q[0];
sx q[0];
rz(-2.2233133) q[0];
sx q[0];
rz(-0.44256532) q[0];
x q[1];
rz(3.1062791) q[2];
sx q[2];
rz(-2.0853634) q[2];
sx q[2];
rz(0.071431486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6060562) q[1];
sx q[1];
rz(-2.5185985) q[1];
sx q[1];
rz(-2.6583903) q[1];
rz(1.1282519) q[3];
sx q[3];
rz(-0.62164111) q[3];
sx q[3];
rz(2.3592202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6601861) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(-1.0931724) q[2];
rz(-2.6045351) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401232) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(-1.681666) q[0];
rz(-1.6319252) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(-3.0622838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58188841) q[0];
sx q[0];
rz(-1.6515344) q[0];
sx q[0];
rz(1.333192) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7081291) q[2];
sx q[2];
rz(-1.7725362) q[2];
sx q[2];
rz(1.2936426) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39187688) q[1];
sx q[1];
rz(-1.8444711) q[1];
sx q[1];
rz(0.43507149) q[1];
rz(-1.7071869) q[3];
sx q[3];
rz(-2.2838784) q[3];
sx q[3];
rz(-0.73120414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1377533) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(2.7286781) q[2];
rz(-1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(2.7220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9846648) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(2.8562163) q[0];
rz(3.0889619) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(-2.9579128) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7516458) q[0];
sx q[0];
rz(-3.0245298) q[0];
sx q[0];
rz(-1.879891) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2854699) q[2];
sx q[2];
rz(-1.7186224) q[2];
sx q[2];
rz(-1.2814652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0252776) q[1];
sx q[1];
rz(-2.5417788) q[1];
sx q[1];
rz(0.13529899) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64818188) q[3];
sx q[3];
rz(-1.8055918) q[3];
sx q[3];
rz(-0.33644331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24667428) q[2];
sx q[2];
rz(-1.4072714) q[2];
sx q[2];
rz(1.0651917) q[2];
rz(1.6861606) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0378549) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(1.6023585) q[0];
rz(-0.94929758) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(1.7112214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1827138) q[0];
sx q[0];
rz(-2.4005167) q[0];
sx q[0];
rz(0.7345906) q[0];
rz(0.84202853) q[2];
sx q[2];
rz(-1.6269738) q[2];
sx q[2];
rz(-0.89955715) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9897223) q[1];
sx q[1];
rz(-0.29726004) q[1];
sx q[1];
rz(-0.10135915) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1422906) q[3];
sx q[3];
rz(-1.8290038) q[3];
sx q[3];
rz(2.9203897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1324233) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(-3.0150748) q[2];
rz(-0.33454076) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(-2.2693966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406521) q[0];
sx q[0];
rz(-0.88467389) q[0];
sx q[0];
rz(0.51710039) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.5093191) q[1];
sx q[1];
rz(-0.089769207) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1506648) q[0];
sx q[0];
rz(-2.1164843) q[0];
sx q[0];
rz(0.47713466) q[0];
rz(0.062762063) q[2];
sx q[2];
rz(-2.1307935) q[2];
sx q[2];
rz(2.7095209) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7081212) q[1];
sx q[1];
rz(-0.29954391) q[1];
sx q[1];
rz(-2.5764861) q[1];
rz(1.0221972) q[3];
sx q[3];
rz(-0.41893235) q[3];
sx q[3];
rz(1.054712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6101997) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(0.72719491) q[2];
rz(-0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(2.9288647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1995734) q[0];
sx q[0];
rz(-1.6006391) q[0];
sx q[0];
rz(-2.0508456) q[0];
rz(1.3923116) q[1];
sx q[1];
rz(-1.6284457) q[1];
sx q[1];
rz(1.5096691) q[1];
rz(3.0346995) q[2];
sx q[2];
rz(-0.86915599) q[2];
sx q[2];
rz(-0.55021777) q[2];
rz(0.58842622) q[3];
sx q[3];
rz(-1.0259494) q[3];
sx q[3];
rz(-1.1001724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
