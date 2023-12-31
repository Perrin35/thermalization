OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2479808) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(-3.0889838) q[0];
rz(1.5300418) q[2];
sx q[2];
rz(-1.1062804) q[2];
sx q[2];
rz(-1.0711311) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6005046) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(2.0946676) q[1];
rz(-pi) q[2];
rz(-1.7928042) q[3];
sx q[3];
rz(-1.7553925) q[3];
sx q[3];
rz(-2.8649462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4336006) q[0];
sx q[0];
rz(-0.029768243) q[0];
sx q[0];
rz(-3.004651) q[0];
x q[1];
rz(-2.8005373) q[2];
sx q[2];
rz(-0.58652069) q[2];
sx q[2];
rz(-1.7797433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95784159) q[1];
sx q[1];
rz(-1.9722003) q[1];
sx q[1];
rz(0.5387696) q[1];
x q[2];
rz(-0.12947793) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(-2.5667218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7188321) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(-2.7405222) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7558407) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(-1.2582448) q[0];
rz(0.24984078) q[2];
sx q[2];
rz(-0.65953883) q[2];
sx q[2];
rz(-1.8674873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2420826) q[1];
sx q[1];
rz(-0.84065719) q[1];
sx q[1];
rz(-1.5622557) q[1];
rz(2.9455455) q[3];
sx q[3];
rz(-1.9571597) q[3];
sx q[3];
rz(-2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.4439616) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(0.25340432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93638203) q[0];
sx q[0];
rz(-2.10996) q[0];
sx q[0];
rz(0.11986952) q[0];
rz(-pi) q[1];
rz(-1.0850111) q[2];
sx q[2];
rz(-1.5861142) q[2];
sx q[2];
rz(-2.3682396) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4761915) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(2.979216) q[1];
rz(-pi) q[2];
rz(-0.13417379) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(-1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(2.3542662) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-1.0095989) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-0.45809349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8100909) q[0];
sx q[0];
rz(-2.0134263) q[0];
sx q[0];
rz(-1.5682194) q[0];
rz(2.6890254) q[2];
sx q[2];
rz(-1.0107702) q[2];
sx q[2];
rz(-2.3134311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0702857) q[1];
sx q[1];
rz(-2.3100393) q[1];
sx q[1];
rz(0.15858312) q[1];
x q[2];
rz(-1.7430274) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(2.9120973) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8454682) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(-1.0494997) q[0];
rz(1.426258) q[2];
sx q[2];
rz(-0.83206165) q[2];
sx q[2];
rz(-1.9369672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8900745) q[1];
sx q[1];
rz(-1.2444082) q[1];
sx q[1];
rz(0.20784394) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.091286) q[3];
sx q[3];
rz(-1.3538085) q[3];
sx q[3];
rz(0.12245164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-1.0725347) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-2.0297208) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3277153) q[0];
sx q[0];
rz(-1.0881256) q[0];
sx q[0];
rz(-3.0359603) q[0];
x q[1];
rz(-2.578031) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(2.3129472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48549451) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(-0.2904201) q[1];
rz(-pi) q[2];
rz(-2.0734378) q[3];
sx q[3];
rz(-1.2253237) q[3];
sx q[3];
rz(-2.2964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(-2.4107966) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(2.3866167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2728111) q[0];
sx q[0];
rz(-2.4065354) q[0];
sx q[0];
rz(1.3520157) q[0];
rz(-pi) q[1];
rz(2.5078012) q[2];
sx q[2];
rz(-0.83665028) q[2];
sx q[2];
rz(-1.447669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8440486) q[1];
sx q[1];
rz(-2.5522759) q[1];
sx q[1];
rz(0.20357666) q[1];
rz(2.6492277) q[3];
sx q[3];
rz(-1.0400606) q[3];
sx q[3];
rz(-2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(-0.26646715) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-0.019502217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7050539) q[0];
sx q[0];
rz(-0.22163135) q[0];
sx q[0];
rz(1.1987232) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2150061) q[2];
sx q[2];
rz(-2.8928061) q[2];
sx q[2];
rz(-0.34005806) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49669493) q[1];
sx q[1];
rz(-0.95458191) q[1];
sx q[1];
rz(-2.9195021) q[1];
x q[2];
rz(0.1653413) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185547) q[0];
sx q[0];
rz(-1.7114637) q[0];
sx q[0];
rz(-1.592357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7933153) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(2.7811108) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1710098) q[1];
sx q[1];
rz(-0.90921558) q[1];
sx q[1];
rz(0.2172825) q[1];
x q[2];
rz(1.6230574) q[3];
sx q[3];
rz(-1.9196379) q[3];
sx q[3];
rz(2.7319752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(2.7764376) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(0.45583367) q[3];
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
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(-1.6297324) q[3];
sx q[3];
rz(-0.56093506) q[3];
sx q[3];
rz(-0.18118071) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
