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
rz(0.084963381) q[0];
sx q[0];
rz(3.4439937) q[0];
sx q[0];
rz(6.2511282) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(-2.6377331) q[1];
sx q[1];
rz(-2.388968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0085474) q[0];
sx q[0];
rz(-0.38131443) q[0];
sx q[0];
rz(2.5757936) q[0];
rz(-0.36280323) q[2];
sx q[2];
rz(-1.9224836) q[2];
sx q[2];
rz(-1.2376518) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9960027) q[1];
sx q[1];
rz(-1.8217345) q[1];
sx q[1];
rz(0.35059021) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21423046) q[3];
sx q[3];
rz(-1.9211244) q[3];
sx q[3];
rz(3.0607579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(2.3215129) q[2];
rz(-2.7005633) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.012861982) q[0];
sx q[0];
rz(-0.91330376) q[0];
sx q[0];
rz(-2.7313857) q[0];
rz(2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.358323) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84285865) q[0];
sx q[0];
rz(-2.4141387) q[0];
sx q[0];
rz(-0.95746118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3097281) q[2];
sx q[2];
rz(-2.4565149) q[2];
sx q[2];
rz(0.58704228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1413026) q[1];
sx q[1];
rz(-1.4030255) q[1];
sx q[1];
rz(-3.0649158) q[1];
rz(-2.210564) q[3];
sx q[3];
rz(-1.6954386) q[3];
sx q[3];
rz(2.3384936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35839781) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(-0.50296339) q[2];
rz(-1.8424235) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(2.1657535) q[0];
rz(0.93636912) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(0.13793129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59527151) q[0];
sx q[0];
rz(-2.0450651) q[0];
sx q[0];
rz(2.7300937) q[0];
x q[1];
rz(1.4549655) q[2];
sx q[2];
rz(-1.1755014) q[2];
sx q[2];
rz(-0.17648736) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8616069) q[1];
sx q[1];
rz(-2.2642038) q[1];
sx q[1];
rz(-2.9704291) q[1];
rz(-1.8432328) q[3];
sx q[3];
rz(-1.1198634) q[3];
sx q[3];
rz(-2.65254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5522573) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(1.4790347) q[2];
rz(-0.79834437) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8524356) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(1.5929219) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(-0.41935316) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33808483) q[0];
sx q[0];
rz(-0.28059059) q[0];
sx q[0];
rz(2.1795033) q[0];
rz(-pi) q[1];
rz(2.6134562) q[2];
sx q[2];
rz(-1.4924876) q[2];
sx q[2];
rz(-2.9700043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6904628) q[1];
sx q[1];
rz(-1.9277384) q[1];
sx q[1];
rz(0.068515645) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0859503) q[3];
sx q[3];
rz(-1.8797698) q[3];
sx q[3];
rz(0.072871836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7833917) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(1.5979213) q[2];
rz(-1.915043) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2655547) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(0.75620404) q[0];
rz(-1.8887695) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(1.7255712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9276816) q[0];
sx q[0];
rz(-0.98386231) q[0];
sx q[0];
rz(-1.8326493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7977095) q[2];
sx q[2];
rz(-1.4522219) q[2];
sx q[2];
rz(-0.32394513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6452613) q[1];
sx q[1];
rz(-2.1095536) q[1];
sx q[1];
rz(-1.5625728) q[1];
x q[2];
rz(-2.6719499) q[3];
sx q[3];
rz(-0.40294632) q[3];
sx q[3];
rz(0.5136036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8004134) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(2.9608012) q[2];
rz(-1.6527269) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(-1.5159877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4329231) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(-2.3413626) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(0.12399331) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2469047) q[0];
sx q[0];
rz(-0.13291153) q[0];
sx q[0];
rz(-1.7869496) q[0];
rz(-pi) q[1];
rz(3.0818865) q[2];
sx q[2];
rz(-1.8054031) q[2];
sx q[2];
rz(1.4690746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1734487) q[1];
sx q[1];
rz(-0.3160797) q[1];
sx q[1];
rz(2.8760002) q[1];
rz(2.0923397) q[3];
sx q[3];
rz(-0.66415706) q[3];
sx q[3];
rz(2.5356027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(3.0653811) q[2];
rz(1.8501836) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(0.61856234) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16977075) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(0.75702697) q[0];
rz(-2.3454759) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-2.1597791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2404838) q[0];
sx q[0];
rz(-0.79190688) q[0];
sx q[0];
rz(1.9740941) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4635529) q[2];
sx q[2];
rz(-1.5011884) q[2];
sx q[2];
rz(-2.7977242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0714149) q[1];
sx q[1];
rz(-1.9767673) q[1];
sx q[1];
rz(2.0858048) q[1];
rz(0.14431868) q[3];
sx q[3];
rz(-1.8486406) q[3];
sx q[3];
rz(-1.8622189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41165274) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(0.5438424) q[2];
rz(-3.1206711) q[3];
sx q[3];
rz(-1.7048416) q[3];
sx q[3];
rz(0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92886096) q[0];
sx q[0];
rz(-2.2305771) q[0];
sx q[0];
rz(0.62335706) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-0.61505452) q[1];
sx q[1];
rz(0.26434937) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8835444) q[0];
sx q[0];
rz(-2.2730581) q[0];
sx q[0];
rz(2.8625898) q[0];
x q[1];
rz(3.1000443) q[2];
sx q[2];
rz(-1.2354406) q[2];
sx q[2];
rz(-1.2183508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0919288) q[1];
sx q[1];
rz(-2.290002) q[1];
sx q[1];
rz(-2.767241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13295138) q[3];
sx q[3];
rz(-1.2556517) q[3];
sx q[3];
rz(-0.12115762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.136772) q[2];
sx q[2];
rz(-0.68237582) q[2];
sx q[2];
rz(0.47478673) q[2];
rz(-2.7759806) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(2.9584208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356165) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(-2.6557652) q[0];
rz(2.0756508) q[1];
sx q[1];
rz(-1.424788) q[1];
sx q[1];
rz(-0.33445439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4518648) q[0];
sx q[0];
rz(-2.7948927) q[0];
sx q[0];
rz(2.6348216) q[0];
rz(-pi) q[1];
rz(-2.7031519) q[2];
sx q[2];
rz(-1.6980357) q[2];
sx q[2];
rz(1.3194989) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6600646) q[1];
sx q[1];
rz(-2.017952) q[1];
sx q[1];
rz(-1.7363973) q[1];
rz(0.79631898) q[3];
sx q[3];
rz(-1.2346754) q[3];
sx q[3];
rz(-1.4848061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7740384) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(2.1749707) q[2];
rz(-1.4878081) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(-0.070913471) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8129355) q[0];
sx q[0];
rz(-0.85449496) q[0];
sx q[0];
rz(0.27035126) q[0];
rz(-1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(0.62634748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1046902) q[0];
sx q[0];
rz(-2.2459163) q[0];
sx q[0];
rz(-2.7319607) q[0];
rz(-pi) q[1];
rz(0.93339351) q[2];
sx q[2];
rz(-1.9248171) q[2];
sx q[2];
rz(-0.4590946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8423564) q[1];
sx q[1];
rz(-1.6116214) q[1];
sx q[1];
rz(3.051332) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38624951) q[3];
sx q[3];
rz(-1.874141) q[3];
sx q[3];
rz(1.7084427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35499972) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(-1.5691441) q[2];
rz(0.7507503) q[3];
sx q[3];
rz(-2.7995977) q[3];
sx q[3];
rz(0.87740889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.5759721) q[0];
sx q[0];
rz(-1.8436057) q[0];
sx q[0];
rz(2.5634503) q[0];
rz(-1.3427973) q[1];
sx q[1];
rz(-0.31675757) q[1];
sx q[1];
rz(-2.996179) q[1];
rz(3.1411896) q[2];
sx q[2];
rz(-3.0192791) q[2];
sx q[2];
rz(3.064439) q[2];
rz(2.9407637) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
