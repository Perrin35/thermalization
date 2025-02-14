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
rz(-2.9989938) q[0];
rz(-2.5481186) q[1];
sx q[1];
rz(-0.65294099) q[1];
sx q[1];
rz(1.0452193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3680688) q[0];
sx q[0];
rz(-2.0288101) q[0];
sx q[0];
rz(2.5450443) q[0];
x q[1];
rz(3.1047761) q[2];
sx q[2];
rz(-0.93339257) q[2];
sx q[2];
rz(-2.3790145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3821311) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(2.1407024) q[1];
rz(0.11868415) q[3];
sx q[3];
rz(-2.136345) q[3];
sx q[3];
rz(-2.8142336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97592252) q[2];
sx q[2];
rz(-1.0079577) q[2];
sx q[2];
rz(2.9284076) q[2];
rz(0.91607696) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(-0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857472) q[0];
sx q[0];
rz(-2.2540932) q[0];
sx q[0];
rz(2.6708653) q[0];
rz(2.0384906) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(2.5618166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.267201) q[0];
sx q[0];
rz(-1.6306595) q[0];
sx q[0];
rz(-2.6169928) q[0];
x q[1];
rz(2.5986541) q[2];
sx q[2];
rz(-1.2839054) q[2];
sx q[2];
rz(-0.33143383) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75232154) q[1];
sx q[1];
rz(-1.0022517) q[1];
sx q[1];
rz(2.6804377) q[1];
rz(2.1768119) q[3];
sx q[3];
rz(-2.7005115) q[3];
sx q[3];
rz(1.6562605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9819928) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(-1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5461102) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(-2.8693759) q[0];
rz(0.1046003) q[1];
sx q[1];
rz(-1.0403386) q[1];
sx q[1];
rz(1.4629755) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8151096) q[0];
sx q[0];
rz(-2.4158951) q[0];
sx q[0];
rz(1.4952907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3562695) q[2];
sx q[2];
rz(-1.6162655) q[2];
sx q[2];
rz(-0.77982219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.453664) q[1];
sx q[1];
rz(-1.2896104) q[1];
sx q[1];
rz(-2.4635876) q[1];
rz(-pi) q[2];
rz(-0.1553673) q[3];
sx q[3];
rz(-2.7685363) q[3];
sx q[3];
rz(-0.01580007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0649123) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(0.40985516) q[2];
rz(0.05154933) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(2.4079017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9728397) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(-2.4421316) q[0];
rz(2.2109924) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(-2.0420989) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52877766) q[0];
sx q[0];
rz(-1.2941501) q[0];
sx q[0];
rz(2.1601281) q[0];
rz(-pi) q[1];
rz(1.9540953) q[2];
sx q[2];
rz(-0.58247236) q[2];
sx q[2];
rz(0.7815643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1395542) q[1];
sx q[1];
rz(-2.6894662) q[1];
sx q[1];
rz(0.76805784) q[1];
rz(1.1641665) q[3];
sx q[3];
rz(-0.96897954) q[3];
sx q[3];
rz(-0.072657771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76116556) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(-0.94044828) q[2];
rz(-0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(0.57653069) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1171653) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(1.0166136) q[0];
rz(-2.2158465) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(-0.076676682) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8228418) q[0];
sx q[0];
rz(-0.42993507) q[0];
sx q[0];
rz(-1.1640401) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20598866) q[2];
sx q[2];
rz(-1.811337) q[2];
sx q[2];
rz(-2.5702916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0106153) q[1];
sx q[1];
rz(-2.1226127) q[1];
sx q[1];
rz(0.24893649) q[1];
x q[2];
rz(1.2604976) q[3];
sx q[3];
rz(-0.7782794) q[3];
sx q[3];
rz(-2.8097514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-0.7879476) q[2];
sx q[2];
rz(3.0184271) q[2];
rz(-2.4672616) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(0.82790747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8318361) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-0.48102608) q[0];
rz(0.24093974) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(3.0066838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3137993) q[0];
sx q[0];
rz(-2.302357) q[0];
sx q[0];
rz(-0.93904943) q[0];
x q[1];
rz(1.4675702) q[2];
sx q[2];
rz(-2.1035668) q[2];
sx q[2];
rz(-0.68160439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1926769) q[1];
sx q[1];
rz(-2.2971616) q[1];
sx q[1];
rz(0.37096687) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6429796) q[3];
sx q[3];
rz(-1.1218539) q[3];
sx q[3];
rz(-0.23060966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9350819) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(-1.3235462) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-2.2836253) q[3];
sx q[3];
rz(2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(2.6351567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424935) q[0];
sx q[0];
rz(-0.37537071) q[0];
sx q[0];
rz(-2.9604549) q[0];
rz(-pi) q[1];
rz(-2.0525706) q[2];
sx q[2];
rz(-1.543569) q[2];
sx q[2];
rz(-2.1435062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87878972) q[1];
sx q[1];
rz(-1.4775986) q[1];
sx q[1];
rz(1.4944939) q[1];
rz(-pi) q[2];
rz(3.0783404) q[3];
sx q[3];
rz(-0.58025415) q[3];
sx q[3];
rz(1.3903416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(-0.55240101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34614554) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(-0.41918293) q[0];
rz(-3.0508793) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(-1.027164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4226728) q[0];
sx q[0];
rz(-0.75506588) q[0];
sx q[0];
rz(-2.4605196) q[0];
x q[1];
rz(1.5628603) q[2];
sx q[2];
rz(-1.5936226) q[2];
sx q[2];
rz(-0.16032444) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5324983) q[1];
sx q[1];
rz(-0.23819345) q[1];
sx q[1];
rz(-1.1042117) q[1];
x q[2];
rz(0.88524359) q[3];
sx q[3];
rz(-1.5285057) q[3];
sx q[3];
rz(-1.2219747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1586228) q[2];
sx q[2];
rz(-1.0047487) q[2];
sx q[2];
rz(2.1515004) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(-3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38744277) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-2.5842174) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-2.974496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97826577) q[0];
sx q[0];
rz(-1.8370795) q[0];
sx q[0];
rz(2.8243218) q[0];
rz(1.3873151) q[2];
sx q[2];
rz(-2.263732) q[2];
sx q[2];
rz(-0.82833457) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22699478) q[1];
sx q[1];
rz(-0.92508864) q[1];
sx q[1];
rz(1.1832994) q[1];
x q[2];
rz(0.88730335) q[3];
sx q[3];
rz(-1.512882) q[3];
sx q[3];
rz(-1.7660727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1780213) q[2];
sx q[2];
rz(-2.9489297) q[2];
sx q[2];
rz(-2.7131405) q[2];
rz(-0.7306478) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(-2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927602) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(2.0220508) q[0];
rz(2.4730832) q[1];
sx q[1];
rz(-0.57066494) q[1];
sx q[1];
rz(0.25984919) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6194659) q[0];
sx q[0];
rz(-1.8308976) q[0];
sx q[0];
rz(1.0570265) q[0];
rz(-pi) q[1];
rz(-2.2462559) q[2];
sx q[2];
rz(-0.5731715) q[2];
sx q[2];
rz(-1.8942539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47626859) q[1];
sx q[1];
rz(-0.55108738) q[1];
sx q[1];
rz(0.70161201) q[1];
rz(-pi) q[2];
rz(0.15774653) q[3];
sx q[3];
rz(-0.88706568) q[3];
sx q[3];
rz(-2.5585548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4708289) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-2.0551576) q[2];
rz(-0.49232617) q[3];
sx q[3];
rz(-2.2958675) q[3];
sx q[3];
rz(0.72211784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928741) q[0];
sx q[0];
rz(-1.6328136) q[0];
sx q[0];
rz(1.6528224) q[0];
rz(1.8863574) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(-0.095070953) q[2];
sx q[2];
rz(-1.173272) q[2];
sx q[2];
rz(1.6116326) q[2];
rz(-1.6519573) q[3];
sx q[3];
rz(-0.66417821) q[3];
sx q[3];
rz(1.3978721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
