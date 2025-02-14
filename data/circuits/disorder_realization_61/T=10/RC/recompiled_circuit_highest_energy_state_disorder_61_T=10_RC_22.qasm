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
rz(1.6545777) q[1];
sx q[1];
rz(4.5643212) q[1];
sx q[1];
rz(8.1017452) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9441863) q[0];
sx q[0];
rz(-1.8776769) q[0];
sx q[0];
rz(0.57800622) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5147314) q[2];
sx q[2];
rz(-2.755609) q[2];
sx q[2];
rz(-2.4311993) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3016794) q[1];
sx q[1];
rz(-0.29407802) q[1];
sx q[1];
rz(1.9211576) q[1];
rz(-pi) q[2];
rz(-2.5526974) q[3];
sx q[3];
rz(-2.3213904) q[3];
sx q[3];
rz(0.39862788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089791678) q[0];
sx q[0];
rz(-0.3312411) q[0];
sx q[0];
rz(-0.077089699) q[0];
rz(2.5568621) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.9416521) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18835078) q[0];
sx q[0];
rz(-0.50354275) q[0];
sx q[0];
rz(1.9882697) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7351772) q[2];
sx q[2];
rz(-1.5417678) q[2];
sx q[2];
rz(2.8433702) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5189831) q[1];
sx q[1];
rz(-2.4138192) q[1];
sx q[1];
rz(-2.1534797) q[1];
rz(-0.097278519) q[3];
sx q[3];
rz(-0.56124282) q[3];
sx q[3];
rz(-1.8986072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.2526907) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(2.4661031) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(3.1233136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.49552396) q[0];
sx q[0];
rz(-1.6697474) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(0.010628788) q[1];
sx q[1];
rz(-1.7517215) q[1];
sx q[1];
rz(-2.1307814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0389688) q[0];
sx q[0];
rz(-2.0576697) q[0];
sx q[0];
rz(0.34619934) q[0];
x q[1];
rz(-3.1132142) q[2];
sx q[2];
rz(-0.88936964) q[2];
sx q[2];
rz(0.703237) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0360196) q[1];
sx q[1];
rz(-1.509799) q[1];
sx q[1];
rz(2.3954844) q[1];
rz(-pi) q[2];
rz(-0.63186462) q[3];
sx q[3];
rz(-0.87040983) q[3];
sx q[3];
rz(-2.4720342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90091554) q[2];
sx q[2];
rz(-2.7132576) q[2];
sx q[2];
rz(-2.6348616) q[2];
rz(1.8690551) q[3];
sx q[3];
rz(-1.7839909) q[3];
sx q[3];
rz(2.8633964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5469359) q[0];
sx q[0];
rz(-1.9173859) q[0];
sx q[0];
rz(1.1154255) q[0];
rz(1.2660654) q[1];
sx q[1];
rz(-1.5770117) q[1];
sx q[1];
rz(-0.87475264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6653671) q[0];
sx q[0];
rz(-1.8221119) q[0];
sx q[0];
rz(-0.091500207) q[0];
x q[1];
rz(1.5953996) q[2];
sx q[2];
rz(-2.4877791) q[2];
sx q[2];
rz(3.1179414) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0759039) q[1];
sx q[1];
rz(-1.3049647) q[1];
sx q[1];
rz(2.7194383) q[1];
x q[2];
rz(-1.3618008) q[3];
sx q[3];
rz(-0.71729388) q[3];
sx q[3];
rz(-1.6092666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0227585) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(-1.215722) q[2];
rz(1.3974961) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(1.4309179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27595156) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(-1.4707461) q[0];
rz(-1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(-1.4102304) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84274693) q[0];
sx q[0];
rz(-2.0618467) q[0];
sx q[0];
rz(-0.23570717) q[0];
rz(-pi) q[1];
rz(1.330014) q[2];
sx q[2];
rz(-2.2301939) q[2];
sx q[2];
rz(1.8383775) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3258563) q[1];
sx q[1];
rz(-1.5724025) q[1];
sx q[1];
rz(-2.8291507) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2447884) q[3];
sx q[3];
rz(-2.6072579) q[3];
sx q[3];
rz(-0.33269445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1269647) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(-1.3237759) q[2];
rz(-1.3592367) q[3];
sx q[3];
rz(-1.8424572) q[3];
sx q[3];
rz(0.38715473) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(-0.02221814) q[0];
rz(-0.70023316) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(-2.1846695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2882521) q[0];
sx q[0];
rz(-1.4256251) q[0];
sx q[0];
rz(2.4397544) q[0];
rz(-pi) q[1];
rz(-0.21778743) q[2];
sx q[2];
rz(-1.3931257) q[2];
sx q[2];
rz(-1.9363994) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80920389) q[1];
sx q[1];
rz(-1.1706377) q[1];
sx q[1];
rz(2.8431812) q[1];
rz(2.2150008) q[3];
sx q[3];
rz(-2.6098688) q[3];
sx q[3];
rz(0.34429541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9243098) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(-1.4420606) q[2];
rz(2.0349272) q[3];
sx q[3];
rz(-0.60551867) q[3];
sx q[3];
rz(-1.1520011) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149413) q[0];
sx q[0];
rz(-1.8639257) q[0];
sx q[0];
rz(0.46698025) q[0];
rz(1.3106208) q[1];
sx q[1];
rz(-2.1881723) q[1];
sx q[1];
rz(-0.39042815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0434131) q[0];
sx q[0];
rz(-1.5382086) q[0];
sx q[0];
rz(0.14049732) q[0];
rz(-pi) q[1];
rz(0.10755928) q[2];
sx q[2];
rz(-1.2778408) q[2];
sx q[2];
rz(2.7358472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5891287) q[1];
sx q[1];
rz(-2.2099582) q[1];
sx q[1];
rz(0.57560779) q[1];
rz(-pi) q[2];
rz(1.3620699) q[3];
sx q[3];
rz(-0.62297076) q[3];
sx q[3];
rz(1.8546212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9499669) q[2];
sx q[2];
rz(-1.2167296) q[2];
sx q[2];
rz(2.1245655) q[2];
rz(-2.2022066) q[3];
sx q[3];
rz(-0.87299577) q[3];
sx q[3];
rz(0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7153213) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(-1.2128879) q[0];
rz(0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(0.42339465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903666) q[0];
sx q[0];
rz(-2.4366424) q[0];
sx q[0];
rz(-1.6909107) q[0];
rz(-pi) q[1];
rz(-1.4199443) q[2];
sx q[2];
rz(-1.6467588) q[2];
sx q[2];
rz(-1.1830038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4881542) q[1];
sx q[1];
rz(-2.8677925) q[1];
sx q[1];
rz(-0.68922148) q[1];
rz(0.16984197) q[3];
sx q[3];
rz(-1.563493) q[3];
sx q[3];
rz(-0.23582349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49070552) q[2];
sx q[2];
rz(-2.8287973) q[2];
sx q[2];
rz(0.97839626) q[2];
rz(-2.4466416) q[3];
sx q[3];
rz(-1.9415104) q[3];
sx q[3];
rz(1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4872015) q[0];
sx q[0];
rz(-2.374361) q[0];
sx q[0];
rz(2.5228398) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9390807) q[2];
sx q[2];
rz(-1.2767451) q[2];
sx q[2];
rz(2.685315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14230141) q[1];
sx q[1];
rz(-0.55146927) q[1];
sx q[1];
rz(0.89889185) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6641162) q[3];
sx q[3];
rz(-1.9100185) q[3];
sx q[3];
rz(-0.71996688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9110079) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(-1.9045551) q[2];
rz(-2.6593995) q[3];
sx q[3];
rz(-0.92415205) q[3];
sx q[3];
rz(-0.48875109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1173387) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(1.8035969) q[0];
rz(2.2894739) q[1];
sx q[1];
rz(-1.2282649) q[1];
sx q[1];
rz(-1.9911912) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.710708) q[0];
sx q[0];
rz(-2.0938494) q[0];
sx q[0];
rz(-2.5780748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2427203) q[2];
sx q[2];
rz(-0.33362937) q[2];
sx q[2];
rz(1.7868113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2472966) q[1];
sx q[1];
rz(-1.9622231) q[1];
sx q[1];
rz(0.15116926) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97148599) q[3];
sx q[3];
rz(-1.8953634) q[3];
sx q[3];
rz(0.49719405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77356768) q[2];
sx q[2];
rz(-2.1358392) q[2];
sx q[2];
rz(2.8311938) q[2];
rz(-0.6271022) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(-1.9703126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2720168) q[0];
sx q[0];
rz(-1.5740812) q[0];
sx q[0];
rz(1.5480702) q[0];
rz(-1.2177474) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(1.942684) q[2];
sx q[2];
rz(-0.52996427) q[2];
sx q[2];
rz(0.90373273) q[2];
rz(1.5094093) q[3];
sx q[3];
rz(-1.2493549) q[3];
sx q[3];
rz(-3.1102104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
