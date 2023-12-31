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
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(1.9231208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89361184) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(-0.052608842) q[0];
rz(1.6115509) q[2];
sx q[2];
rz(-1.1062804) q[2];
sx q[2];
rz(-2.0704616) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6005046) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(-1.0469251) q[1];
rz(-pi) q[2];
rz(0.86745947) q[3];
sx q[3];
rz(-0.28775035) q[3];
sx q[3];
rz(0.61144637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(-1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(-2.2024343) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70799202) q[0];
sx q[0];
rz(-3.1118244) q[0];
sx q[0];
rz(-0.13694163) q[0];
rz(-pi) q[1];
rz(-0.34105532) q[2];
sx q[2];
rz(-0.58652069) q[2];
sx q[2];
rz(-1.3618493) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95784159) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(-0.5387696) q[1];
rz(-3.0121147) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(2.5667218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7188321) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(0.20673949) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553592) q[0];
sx q[0];
rz(-1.8952574) q[0];
sx q[0];
rz(-0.1089036) q[0];
x q[1];
rz(2.4972649) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(-2.6459141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.229278) q[1];
sx q[1];
rz(-2.4114128) q[1];
sx q[1];
rz(0.0095403949) q[1];
x q[2];
rz(1.9639575) q[3];
sx q[3];
rz(-1.3893681) q[3];
sx q[3];
rz(-1.7593311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-0.78021375) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.697631) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(2.8881883) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93638203) q[0];
sx q[0];
rz(-1.0316327) q[0];
sx q[0];
rz(0.11986952) q[0];
x q[1];
rz(1.5379982) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(0.82644586) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4761915) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(0.16237662) q[1];
rz(-pi) q[2];
rz(-3.0074189) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(2.3542662) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(0.062967904) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-0.45809349) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23819085) q[0];
sx q[0];
rz(-1.5684677) q[0];
sx q[0];
rz(-0.44263126) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96197084) q[2];
sx q[2];
rz(-1.95032) q[2];
sx q[2];
rz(2.6517207) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.749436) q[1];
sx q[1];
rz(-1.4538308) q[1];
sx q[1];
rz(-2.3163296) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7430274) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(-0.34190049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(2.9120973) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.5620419) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8454682) q[0];
sx q[0];
rz(-0.86588973) q[0];
sx q[0];
rz(1.0494997) q[0];
rz(-2.9847449) q[2];
sx q[2];
rz(-2.3914797) q[2];
sx q[2];
rz(1.7240806) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25176469) q[1];
sx q[1];
rz(-1.7675195) q[1];
sx q[1];
rz(-1.2377435) q[1];
rz(-pi) q[2];
rz(-0.24354981) q[3];
sx q[3];
rz(-1.1034414) q[3];
sx q[3];
rz(-1.5598284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(2.069058) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(-2.0297208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81387732) q[0];
sx q[0];
rz(-1.0881256) q[0];
sx q[0];
rz(-3.0359603) q[0];
x q[1];
rz(-1.9260336) q[2];
sx q[2];
rz(-1.0676427) q[2];
sx q[2];
rz(0.1728729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8060311) q[1];
sx q[1];
rz(-2.8415488) q[1];
sx q[1];
rz(-2.8801444) q[1];
rz(-2.0734378) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(2.2964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599931) q[0];
sx q[0];
rz(-0.85708517) q[0];
sx q[0];
rz(-2.9478361) q[0];
rz(-pi) q[1];
rz(-2.412699) q[2];
sx q[2];
rz(-2.0260099) q[2];
sx q[2];
rz(0.58065562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0874487) q[1];
sx q[1];
rz(-0.99522299) q[1];
sx q[1];
rz(-1.4364442) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6492277) q[3];
sx q[3];
rz(-2.1015321) q[3];
sx q[3];
rz(-2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.5554957) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-2.0027347) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(3.1220904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6396128) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(-1.7777068) q[0];
x q[1];
rz(2.2150061) q[2];
sx q[2];
rz(-0.24878657) q[2];
sx q[2];
rz(-0.34005806) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2724185) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(1.2692578) q[1];
x q[2];
rz(1.8832302) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(2.0563994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0120221) q[0];
sx q[0];
rz(-2.9992933) q[0];
sx q[0];
rz(0.15108959) q[0];
rz(-pi) q[1];
rz(-2.8414367) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(2.4094827) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3155568) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(1.300699) q[1];
x q[2];
rz(-2.9989472) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(0.25776097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2293573) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(-1.1744432) q[2];
sx q[2];
rz(-0.94559961) q[2];
sx q[2];
rz(-0.0018975817) q[2];
rz(0.036988463) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
