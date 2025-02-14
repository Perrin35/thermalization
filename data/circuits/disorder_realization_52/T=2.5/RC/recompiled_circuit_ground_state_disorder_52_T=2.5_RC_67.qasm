OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(0.90149108) q[0];
rz(-1.787552) q[1];
sx q[1];
rz(-1.6156337) q[1];
sx q[1];
rz(1.3344596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1839361) q[0];
sx q[0];
rz(-1.5165197) q[0];
sx q[0];
rz(0.21845777) q[0];
x q[1];
rz(-0.15057474) q[2];
sx q[2];
rz(-1.5158487) q[2];
sx q[2];
rz(0.51412941) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.035369594) q[1];
sx q[1];
rz(-1.5179885) q[1];
sx q[1];
rz(1.04671) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9525312) q[3];
sx q[3];
rz(-0.24654085) q[3];
sx q[3];
rz(2.2061359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91288519) q[2];
sx q[2];
rz(-3.0444453) q[2];
sx q[2];
rz(2.482282) q[2];
rz(0.77624503) q[3];
sx q[3];
rz(-0.018748911) q[3];
sx q[3];
rz(-0.68500486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9126251) q[0];
sx q[0];
rz(-1.9321059) q[0];
sx q[0];
rz(-1.1932766) q[0];
rz(3.0972262) q[1];
sx q[1];
rz(-3.1293479) q[1];
sx q[1];
rz(-2.9150229) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8346342) q[0];
sx q[0];
rz(-0.0037841664) q[0];
sx q[0];
rz(-2.0223122) q[0];
rz(0.37906693) q[2];
sx q[2];
rz(-1.5866536) q[2];
sx q[2];
rz(-0.011034688) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5331796) q[1];
sx q[1];
rz(-0.4371818) q[1];
sx q[1];
rz(-0.1685779) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9162493) q[3];
sx q[3];
rz(-0.85593191) q[3];
sx q[3];
rz(-0.9957141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6988397) q[2];
sx q[2];
rz(-1.5696462) q[2];
sx q[2];
rz(1.6119831) q[2];
rz(-0.93305856) q[3];
sx q[3];
rz(-1.482684) q[3];
sx q[3];
rz(0.23811594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642692) q[0];
sx q[0];
rz(-0.015559109) q[0];
sx q[0];
rz(-2.9413057) q[0];
rz(-0.00042032584) q[1];
sx q[1];
rz(-0.93685189) q[1];
sx q[1];
rz(-0.013484152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1258365) q[0];
sx q[0];
rz(-1.5061597) q[0];
sx q[0];
rz(-1.9958853) q[0];
rz(-1.6137984) q[2];
sx q[2];
rz(-1.637405) q[2];
sx q[2];
rz(3.1284077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9428448) q[1];
sx q[1];
rz(-1.6403733) q[1];
sx q[1];
rz(1.7067616) q[1];
rz(-1.5192658) q[3];
sx q[3];
rz(-1.4084219) q[3];
sx q[3];
rz(-0.80013093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1991594) q[2];
sx q[2];
rz(-1.556267) q[2];
sx q[2];
rz(-1.5996492) q[2];
rz(-2.13983) q[3];
sx q[3];
rz(-2.9250513) q[3];
sx q[3];
rz(0.17222968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37429419) q[0];
sx q[0];
rz(-0.16574398) q[0];
sx q[0];
rz(-2.7941008) q[0];
rz(2.5047498) q[1];
sx q[1];
rz(-3.1359735) q[1];
sx q[1];
rz(1.9510795) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40909262) q[0];
sx q[0];
rz(-0.14483368) q[0];
sx q[0];
rz(-1.9961137) q[0];
x q[1];
rz(1.6830446) q[2];
sx q[2];
rz(-0.069321037) q[2];
sx q[2];
rz(1.6905418) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61100436) q[1];
sx q[1];
rz(-1.3827033) q[1];
sx q[1];
rz(-2.3744319) q[1];
rz(-pi) q[2];
rz(-2.7772831) q[3];
sx q[3];
rz(-1.3650044) q[3];
sx q[3];
rz(0.28401532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5570598) q[2];
sx q[2];
rz(-0.0396885) q[2];
sx q[2];
rz(-1.7812799) q[2];
rz(1.6654061) q[3];
sx q[3];
rz(-1.5675661) q[3];
sx q[3];
rz(-0.40569693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0535468) q[0];
sx q[0];
rz(-0.73943728) q[0];
sx q[0];
rz(-0.0060225688) q[0];
rz(-1.7209523) q[1];
sx q[1];
rz(-3.0907478) q[1];
sx q[1];
rz(-0.086070148) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0417418) q[0];
sx q[0];
rz(-1.6779416) q[0];
sx q[0];
rz(1.5876549) q[0];
rz(0.055069607) q[2];
sx q[2];
rz(-1.5447642) q[2];
sx q[2];
rz(-0.30773417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7988551) q[1];
sx q[1];
rz(-1.4794083) q[1];
sx q[1];
rz(1.6080086) q[1];
rz(-pi) q[2];
rz(0.9283916) q[3];
sx q[3];
rz(-0.07559055) q[3];
sx q[3];
rz(-0.20353157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.033279557) q[2];
sx q[2];
rz(-0.72282183) q[2];
sx q[2];
rz(1.3839728) q[2];
rz(2.7464187) q[3];
sx q[3];
rz(-0.049001781) q[3];
sx q[3];
rz(-1.9743732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.60642099) q[0];
sx q[0];
rz(-2.9223154) q[0];
sx q[0];
rz(-1.0004591) q[0];
rz(2.4011627) q[1];
sx q[1];
rz(-2.705997) q[1];
sx q[1];
rz(-2.7620517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4576806) q[0];
sx q[0];
rz(-3.0427986) q[0];
sx q[0];
rz(0.26215078) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8431115) q[2];
sx q[2];
rz(-2.0305995) q[2];
sx q[2];
rz(1.0155755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9866375) q[1];
sx q[1];
rz(-1.0670245) q[1];
sx q[1];
rz(-0.030435199) q[1];
rz(1.3336181) q[3];
sx q[3];
rz(-0.46018013) q[3];
sx q[3];
rz(1.4165914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1257816) q[2];
sx q[2];
rz(-2.8534079) q[2];
sx q[2];
rz(-3.06456) q[2];
rz(-0.020126255) q[3];
sx q[3];
rz(-0.052611668) q[3];
sx q[3];
rz(-2.3441815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1064442) q[0];
sx q[0];
rz(-0.012520944) q[0];
sx q[0];
rz(1.7092108) q[0];
rz(0.2969946) q[1];
sx q[1];
rz(-0.15428267) q[1];
sx q[1];
rz(0.061554734) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8650353) q[0];
sx q[0];
rz(-1.9752968) q[0];
sx q[0];
rz(2.1767031) q[0];
x q[1];
rz(-3.1157137) q[2];
sx q[2];
rz(-2.9270083) q[2];
sx q[2];
rz(1.9865004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6541432) q[1];
sx q[1];
rz(-1.2795078) q[1];
sx q[1];
rz(1.4994411) q[1];
x q[2];
rz(1.5331623) q[3];
sx q[3];
rz(-1.1363582) q[3];
sx q[3];
rz(1.561059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2986472) q[2];
sx q[2];
rz(-2.8892398) q[2];
sx q[2];
rz(1.7227777) q[2];
rz(-1.336054) q[3];
sx q[3];
rz(-3.1137443) q[3];
sx q[3];
rz(-1.4068039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1143188) q[0];
sx q[0];
rz(-2.424746) q[0];
sx q[0];
rz(-1.5007716) q[0];
rz(-2.0188792) q[1];
sx q[1];
rz(-0.32674679) q[1];
sx q[1];
rz(1.2935125) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5260391) q[0];
sx q[0];
rz(-0.74239555) q[0];
sx q[0];
rz(1.1008939) q[0];
x q[1];
rz(-0.40385623) q[2];
sx q[2];
rz(-2.4060898) q[2];
sx q[2];
rz(2.1063358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97630097) q[1];
sx q[1];
rz(-1.8433851) q[1];
sx q[1];
rz(-1.3957681) q[1];
x q[2];
rz(-2.4582793) q[3];
sx q[3];
rz(-0.97021996) q[3];
sx q[3];
rz(-1.7699458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5844172) q[2];
sx q[2];
rz(-0.36089218) q[2];
sx q[2];
rz(1.7822251) q[2];
rz(-0.44698295) q[3];
sx q[3];
rz(-0.042526571) q[3];
sx q[3];
rz(-0.16714787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8166703) q[0];
sx q[0];
rz(-1.336038) q[0];
sx q[0];
rz(-2.0409806) q[0];
rz(2.0657516) q[1];
sx q[1];
rz(-2.4918719) q[1];
sx q[1];
rz(-0.67695391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8487602) q[0];
sx q[0];
rz(-2.7037132) q[0];
sx q[0];
rz(-0.48012244) q[0];
x q[1];
rz(0.81440429) q[2];
sx q[2];
rz(-0.16781092) q[2];
sx q[2];
rz(1.4175159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9857603) q[1];
sx q[1];
rz(-1.5712156) q[1];
sx q[1];
rz(-3.1414933) q[1];
x q[2];
rz(-0.30218924) q[3];
sx q[3];
rz(-2.6679278) q[3];
sx q[3];
rz(0.74407265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0000275) q[2];
sx q[2];
rz(-0.0024777369) q[2];
sx q[2];
rz(2.4139717) q[2];
rz(-1.242312) q[3];
sx q[3];
rz(-0.036402313) q[3];
sx q[3];
rz(1.9696994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90945554) q[0];
sx q[0];
rz(-2.1196892) q[0];
sx q[0];
rz(-0.68956462) q[0];
rz(1.6652971) q[1];
sx q[1];
rz(-0.25054014) q[1];
sx q[1];
rz(-0.15588674) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4094226) q[0];
sx q[0];
rz(-0.66930938) q[0];
sx q[0];
rz(0.21440345) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8901398) q[2];
sx q[2];
rz(-1.4551468) q[2];
sx q[2];
rz(2.419099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5248651) q[1];
sx q[1];
rz(-1.5749802) q[1];
sx q[1];
rz(-1.5723438) q[1];
x q[2];
rz(3.0468349) q[3];
sx q[3];
rz(-1.4156431) q[3];
sx q[3];
rz(1.1984389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9547687) q[2];
sx q[2];
rz(-3.0195152) q[2];
sx q[2];
rz(-2.1464777) q[2];
rz(2.9156445) q[3];
sx q[3];
rz(-0.048361691) q[3];
sx q[3];
rz(-2.3154955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7786998) q[0];
sx q[0];
rz(-0.97465546) q[0];
sx q[0];
rz(1.4018651) q[0];
rz(-1.4317935) q[1];
sx q[1];
rz(-1.8451537) q[1];
sx q[1];
rz(0.61737212) q[1];
rz(1.7079034) q[2];
sx q[2];
rz(-0.89024407) q[2];
sx q[2];
rz(-2.2355516) q[2];
rz(-3.0052983) q[3];
sx q[3];
rz(-1.1850428) q[3];
sx q[3];
rz(-0.94120126) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
