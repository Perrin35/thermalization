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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(-0.65499175) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(-1.8185599) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060299035) q[0];
sx q[0];
rz(-0.64612389) q[0];
sx q[0];
rz(0.52559121) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.336786) q[2];
sx q[2];
rz(-1.26097) q[2];
sx q[2];
rz(1.3739623) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.536024) q[1];
sx q[1];
rz(-1.4711416) q[1];
sx q[1];
rz(-1.8479363) q[1];
rz(-pi) q[2];
rz(2.5526974) q[3];
sx q[3];
rz(-2.3213904) q[3];
sx q[3];
rz(-0.39862788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(2.4832671) q[2];
rz(-0.63411921) q[3];
sx q[3];
rz(-1.6987957) q[3];
sx q[3];
rz(-1.770795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089791678) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(3.064503) q[0];
rz(0.58473051) q[1];
sx q[1];
rz(-1.8010151) q[1];
sx q[1];
rz(-1.9416521) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8611885) q[0];
sx q[0];
rz(-1.1139835) q[0];
sx q[0];
rz(-0.21976075) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7351772) q[2];
sx q[2];
rz(-1.5417678) q[2];
sx q[2];
rz(0.2982225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89948326) q[1];
sx q[1];
rz(-2.1596908) q[1];
sx q[1];
rz(-0.45581006) q[1];
rz(-pi) q[2];
rz(0.097278519) q[3];
sx q[3];
rz(-2.5803498) q[3];
sx q[3];
rz(1.2429855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8889019) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(-2.4661031) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(-0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460687) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-2.38696) q[0];
rz(3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(-2.1307814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6957798) q[0];
sx q[0];
rz(-0.58923972) q[0];
sx q[0];
rz(1.0007829) q[0];
rz(-pi) q[1];
rz(1.5358244) q[2];
sx q[2];
rz(-2.4596697) q[2];
sx q[2];
rz(-0.65820314) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52151187) q[1];
sx q[1];
rz(-0.82640582) q[1];
sx q[1];
rz(1.6537731) q[1];
rz(-0.76356865) q[3];
sx q[3];
rz(-2.0393622) q[3];
sx q[3];
rz(-2.6811622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90091554) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(-0.50673103) q[2];
rz(1.2725376) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5469359) q[0];
sx q[0];
rz(-1.9173859) q[0];
sx q[0];
rz(-1.1154255) q[0];
rz(-1.2660654) q[1];
sx q[1];
rz(-1.5645809) q[1];
sx q[1];
rz(2.26684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47622555) q[0];
sx q[0];
rz(-1.3194808) q[0];
sx q[0];
rz(3.0500924) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1227448) q[2];
sx q[2];
rz(-0.91721469) q[2];
sx q[2];
rz(3.1342521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0656888) q[1];
sx q[1];
rz(-1.3049647) q[1];
sx q[1];
rz(0.42215438) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3618008) q[3];
sx q[3];
rz(-2.4242988) q[3];
sx q[3];
rz(1.6092666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(-1.215722) q[2];
rz(-1.3974961) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27595156) q[0];
sx q[0];
rz(-3.0396099) q[0];
sx q[0];
rz(1.4707461) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-0.78796092) q[1];
sx q[1];
rz(-1.4102304) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263162) q[0];
sx q[0];
rz(-1.363376) q[0];
sx q[0];
rz(2.0735969) q[0];
rz(0.67362154) q[2];
sx q[2];
rz(-1.7604019) q[2];
sx q[2];
rz(-3.0233235) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2454589) q[1];
sx q[1];
rz(-1.8832379) q[1];
sx q[1];
rz(-1.5724843) q[1];
rz(-0.35378176) q[3];
sx q[3];
rz(-1.1615586) q[3];
sx q[3];
rz(-1.080846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1269647) q[2];
sx q[2];
rz(-0.85863272) q[2];
sx q[2];
rz(-1.3237759) q[2];
rz(1.7823559) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(-3.1193745) q[0];
rz(-0.70023316) q[1];
sx q[1];
rz(-1.2513688) q[1];
sx q[1];
rz(2.1846695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6895804) q[0];
sx q[0];
rz(-2.4274024) q[0];
sx q[0];
rz(0.22269188) q[0];
rz(-pi) q[1];
rz(-1.7526723) q[2];
sx q[2];
rz(-1.7851014) q[2];
sx q[2];
rz(2.8150812) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8808525) q[1];
sx q[1];
rz(-1.8449952) q[1];
sx q[1];
rz(-1.1541648) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92659183) q[3];
sx q[3];
rz(-0.53172382) q[3];
sx q[3];
rz(2.7972972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2172829) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(-1.4420606) q[2];
rz(-2.0349272) q[3];
sx q[3];
rz(-2.536074) q[3];
sx q[3];
rz(-1.1520011) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42665136) q[0];
sx q[0];
rz(-1.277667) q[0];
sx q[0];
rz(-0.46698025) q[0];
rz(1.8309719) q[1];
sx q[1];
rz(-2.1881723) q[1];
sx q[1];
rz(0.39042815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434131) q[0];
sx q[0];
rz(-1.5382086) q[0];
sx q[0];
rz(-3.0010953) q[0];
rz(-pi) q[1];
rz(3.0340334) q[2];
sx q[2];
rz(-1.2778408) q[2];
sx q[2];
rz(0.40574542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8666208) q[1];
sx q[1];
rz(-2.3094147) q[1];
sx q[1];
rz(0.93865959) q[1];
x q[2];
rz(-1.7795227) q[3];
sx q[3];
rz(-0.62297076) q[3];
sx q[3];
rz(1.8546212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(1.0170271) q[2];
rz(0.93938604) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(2.2123607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(0.7153213) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(-1.9287047) q[0];
rz(0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(-2.718198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.903666) q[0];
sx q[0];
rz(-2.4366424) q[0];
sx q[0];
rz(-1.4506819) q[0];
x q[1];
rz(-1.7216484) q[2];
sx q[2];
rz(-1.6467588) q[2];
sx q[2];
rz(-1.9585889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5880206) q[1];
sx q[1];
rz(-1.7436073) q[1];
sx q[1];
rz(-0.2134448) q[1];
rz(-pi) q[2];
rz(0.043181852) q[3];
sx q[3];
rz(-2.9715952) q[3];
sx q[3];
rz(1.7640597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49070552) q[2];
sx q[2];
rz(-2.8287973) q[2];
sx q[2];
rz(0.97839626) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(-1.8470701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438743) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(3.0978715) q[0];
rz(1.9678736) q[1];
sx q[1];
rz(-2.3269188) q[1];
sx q[1];
rz(0.79183212) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4872015) q[0];
sx q[0];
rz(-0.7672317) q[0];
sx q[0];
rz(-0.61875288) q[0];
x q[1];
rz(0.9845244) q[2];
sx q[2];
rz(-0.35536623) q[2];
sx q[2];
rz(2.9815069) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2481195) q[1];
sx q[1];
rz(-1.9933102) q[1];
sx q[1];
rz(2.7759107) q[1];
rz(-pi) q[2];
rz(0.47747647) q[3];
sx q[3];
rz(-1.2315742) q[3];
sx q[3];
rz(-2.4216258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.23058471) q[2];
sx q[2];
rz(-0.30969301) q[2];
sx q[2];
rz(-1.9045551) q[2];
rz(-0.48219314) q[3];
sx q[3];
rz(-0.92415205) q[3];
sx q[3];
rz(0.48875109) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.024254) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(-1.8035969) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(1.9911912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8088732) q[0];
sx q[0];
rz(-2.3925663) q[0];
sx q[0];
rz(-2.3179884) q[0];
rz(1.8878292) q[2];
sx q[2];
rz(-1.46508) q[2];
sx q[2];
rz(-0.52717613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4070421) q[1];
sx q[1];
rz(-1.7104539) q[1];
sx q[1];
rz(1.9662844) q[1];
rz(1.0329594) q[3];
sx q[3];
rz(-0.67194429) q[3];
sx q[3];
rz(0.6368466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.368025) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(0.31039882) q[2];
rz(-2.5144905) q[3];
sx q[3];
rz(-2.1576364) q[3];
sx q[3];
rz(-1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-1.2720168) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(-1.2177474) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(-1.0711674) q[2];
sx q[2];
rz(-1.3860605) q[2];
sx q[2];
rz(-0.99169127) q[2];
rz(0.18219215) q[3];
sx q[3];
rz(-0.32705083) q[3];
sx q[3];
rz(0.22352945) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
