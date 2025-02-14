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
rz(1.3230327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060299035) q[0];
sx q[0];
rz(-0.64612389) q[0];
sx q[0];
rz(-2.6160014) q[0];
x q[1];
rz(1.8048067) q[2];
sx q[2];
rz(-1.26097) q[2];
sx q[2];
rz(-1.7676304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3016794) q[1];
sx q[1];
rz(-2.8475146) q[1];
sx q[1];
rz(-1.2204351) q[1];
rz(-pi) q[2];
rz(1.0336797) q[3];
sx q[3];
rz(-2.2244646) q[3];
sx q[3];
rz(-1.9680227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(-2.4832671) q[2];
rz(0.63411921) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(-1.770795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.051801) q[0];
sx q[0];
rz(-0.3312411) q[0];
sx q[0];
rz(-3.064503) q[0];
rz(-0.58473051) q[1];
sx q[1];
rz(-1.8010151) q[1];
sx q[1];
rz(-1.1999406) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18835078) q[0];
sx q[0];
rz(-0.50354275) q[0];
sx q[0];
rz(1.9882697) q[0];
x q[1];
rz(1.7351772) q[2];
sx q[2];
rz(-1.5998249) q[2];
sx q[2];
rz(-2.8433702) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2421094) q[1];
sx q[1];
rz(-0.98190183) q[1];
sx q[1];
rz(-0.45581006) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0443141) q[3];
sx q[3];
rz(-0.56124282) q[3];
sx q[3];
rz(-1.2429855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8889019) q[2];
sx q[2];
rz(-2.1609047) q[2];
sx q[2];
rz(-0.67548951) q[2];
rz(-0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49552396) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(1.0108112) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7769612) q[0];
sx q[0];
rz(-1.2662132) q[0];
sx q[0];
rz(2.0833894) q[0];
x q[1];
rz(3.1132142) q[2];
sx q[2];
rz(-2.252223) q[2];
sx q[2];
rz(-2.4383557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1055731) q[1];
sx q[1];
rz(-1.509799) q[1];
sx q[1];
rz(2.3954844) q[1];
rz(-pi) q[2];
rz(-0.95960404) q[3];
sx q[3];
rz(-2.2356847) q[3];
sx q[3];
rz(-1.5184107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90091554) q[2];
sx q[2];
rz(-2.7132576) q[2];
sx q[2];
rz(-0.50673103) q[2];
rz(1.2725376) q[3];
sx q[3];
rz(-1.3576018) q[3];
sx q[3];
rz(-0.27819628) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5946567) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(-1.1154255) q[0];
rz(-1.2660654) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-2.26684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0717569) q[0];
sx q[0];
rz(-1.6594145) q[0];
sx q[0];
rz(-1.3184692) q[0];
x q[1];
rz(-1.546193) q[2];
sx q[2];
rz(-0.65381351) q[2];
sx q[2];
rz(0.023651274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0656888) q[1];
sx q[1];
rz(-1.836628) q[1];
sx q[1];
rz(0.42215438) q[1];
rz(-0.86438365) q[3];
sx q[3];
rz(-1.4339851) q[3];
sx q[3];
rz(-0.19696008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-2.6648882) q[2];
sx q[2];
rz(1.215722) q[2];
rz(-1.7440965) q[3];
sx q[3];
rz(-1.6179061) q[3];
sx q[3];
rz(1.7106748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8656411) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.4707461) q[0];
rz(-1.7631081) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.4102304) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2988457) q[0];
sx q[0];
rz(-2.0618467) q[0];
sx q[0];
rz(2.9058855) q[0];
rz(-1.8115787) q[2];
sx q[2];
rz(-0.91139873) q[2];
sx q[2];
rz(-1.8383775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3258563) q[1];
sx q[1];
rz(-1.5691901) q[1];
sx q[1];
rz(-0.31244197) q[1];
rz(-pi) q[2];
rz(-2.0038807) q[3];
sx q[3];
rz(-1.8942465) q[3];
sx q[3];
rz(-2.5057305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0146279) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(1.8178168) q[2];
rz(-1.3592367) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(2.7544379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112261) q[0];
sx q[0];
rz(-1.68196) q[0];
sx q[0];
rz(0.02221814) q[0];
rz(-0.70023316) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(-2.1846695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8533405) q[0];
sx q[0];
rz(-1.4256251) q[0];
sx q[0];
rz(-2.4397544) q[0];
rz(0.69338696) q[2];
sx q[2];
rz(-0.28017226) q[2];
sx q[2];
rz(-2.1020773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80920389) q[1];
sx q[1];
rz(-1.1706377) q[1];
sx q[1];
rz(-2.8431812) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33958667) q[3];
sx q[3];
rz(-1.9882147) q[3];
sx q[3];
rz(-1.0610895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2172829) q[2];
sx q[2];
rz(-1.1223015) q[2];
sx q[2];
rz(-1.4420606) q[2];
rz(-2.0349272) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(-1.9895915) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149413) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(-2.6746124) q[0];
rz(-1.3106208) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(2.7511645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4425499) q[0];
sx q[0];
rz(-2.9973898) q[0];
sx q[0];
rz(-0.22871916) q[0];
rz(-pi) q[1];
rz(3.0340334) q[2];
sx q[2];
rz(-1.2778408) q[2];
sx q[2];
rz(0.40574542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5891287) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(2.5659849) q[1];
rz(-0.14777811) q[3];
sx q[3];
rz(-0.96333233) q[3];
sx q[3];
rz(-2.1097418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19162576) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(-2.1245655) q[2];
rz(-0.93938604) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(-0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7153213) q[0];
sx q[0];
rz(-2.4763698) q[0];
sx q[0];
rz(-1.9287047) q[0];
rz(2.5384278) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903666) q[0];
sx q[0];
rz(-0.70495021) q[0];
sx q[0];
rz(1.6909107) q[0];
rz(-pi) q[1];
rz(-1.7216484) q[2];
sx q[2];
rz(-1.6467588) q[2];
sx q[2];
rz(1.1830038) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.054476995) q[1];
sx q[1];
rz(-1.3605788) q[1];
sx q[1];
rz(-1.7475378) q[1];
rz(0.043181852) q[3];
sx q[3];
rz(-2.9715952) q[3];
sx q[3];
rz(-1.377533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49070552) q[2];
sx q[2];
rz(-2.8287973) q[2];
sx q[2];
rz(0.97839626) q[2];
rz(-0.69495106) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(-1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438743) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(3.0978715) q[0];
rz(-1.1737191) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(-0.79183212) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7515562) q[0];
sx q[0];
rz(-1.9851713) q[0];
sx q[0];
rz(2.4757371) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2709684) q[2];
sx q[2];
rz(-1.764503) q[2];
sx q[2];
rz(1.1739588) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14230141) q[1];
sx q[1];
rz(-0.55146927) q[1];
sx q[1];
rz(-2.2427008) q[1];
x q[2];
rz(-2.6641162) q[3];
sx q[3];
rz(-1.2315742) q[3];
sx q[3];
rz(-2.4216258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9110079) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(1.9045551) q[2];
rz(2.6593995) q[3];
sx q[3];
rz(-2.2174406) q[3];
sx q[3];
rz(2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.024254) q[0];
sx q[0];
rz(-2.2145705) q[0];
sx q[0];
rz(-1.8035969) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(-1.1504014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3327194) q[0];
sx q[0];
rz(-2.3925663) q[0];
sx q[0];
rz(-2.3179884) q[0];
rz(-pi) q[1];
rz(1.2427203) q[2];
sx q[2];
rz(-0.33362937) q[2];
sx q[2];
rz(-1.3547813) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89429606) q[1];
sx q[1];
rz(-1.1793696) q[1];
sx q[1];
rz(2.9904234) q[1];
rz(-pi) q[2];
rz(-2.1701067) q[3];
sx q[3];
rz(-1.8953634) q[3];
sx q[3];
rz(-2.6443986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77356768) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(-2.8311938) q[2];
rz(-0.6271022) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(-1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8695759) q[0];
sx q[0];
rz(-1.5675114) q[0];
sx q[0];
rz(-1.5935224) q[0];
rz(1.2177474) q[1];
sx q[1];
rz(-1.6883313) q[1];
sx q[1];
rz(2.2399166) q[1];
rz(-1.0711674) q[2];
sx q[2];
rz(-1.3860605) q[2];
sx q[2];
rz(-0.99169127) q[2];
rz(2.8195856) q[3];
sx q[3];
rz(-1.6290355) q[3];
sx q[3];
rz(1.6215948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
