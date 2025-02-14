OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(-1.2195725) q[0];
rz(3.976877) q[1];
sx q[1];
rz(4.4983954) q[1];
sx q[1];
rz(8.5228336) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25489461) q[0];
sx q[0];
rz(-2.5097846) q[0];
sx q[0];
rz(2.8796701) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3839533) q[2];
sx q[2];
rz(-1.9451491) q[2];
sx q[2];
rz(2.6461305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28022596) q[1];
sx q[1];
rz(-1.8609443) q[1];
sx q[1];
rz(-0.18158438) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91633032) q[3];
sx q[3];
rz(-2.188465) q[3];
sx q[3];
rz(0.82181069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1787662) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(-0.82628769) q[2];
rz(-1.4055584) q[3];
sx q[3];
rz(-0.30705753) q[3];
sx q[3];
rz(-0.74339408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58251441) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(-1.4003117) q[0];
rz(2.0179613) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(-2.0434911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208291) q[0];
sx q[0];
rz(-1.6801103) q[0];
sx q[0];
rz(2.9846322) q[0];
rz(1.9015354) q[2];
sx q[2];
rz(-1.8914701) q[2];
sx q[2];
rz(1.641524) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73126331) q[1];
sx q[1];
rz(-1.6917233) q[1];
sx q[1];
rz(1.4951453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9722749) q[3];
sx q[3];
rz(-0.40630925) q[3];
sx q[3];
rz(-0.10029785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0565679) q[2];
sx q[2];
rz(-0.3173863) q[2];
sx q[2];
rz(-1.8918096) q[2];
rz(-1.7791087) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(-1.484163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054841) q[0];
sx q[0];
rz(-0.69378575) q[0];
sx q[0];
rz(2.8662477) q[0];
rz(-2.6823726) q[1];
sx q[1];
rz(-2.1870435) q[1];
sx q[1];
rz(3.088248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33063525) q[0];
sx q[0];
rz(-2.914408) q[0];
sx q[0];
rz(0.35463984) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6079536) q[2];
sx q[2];
rz(-2.6803663) q[2];
sx q[2];
rz(-0.30468582) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5819323) q[1];
sx q[1];
rz(-1.6114863) q[1];
sx q[1];
rz(-2.5896124) q[1];
rz(-pi) q[2];
rz(0.20823222) q[3];
sx q[3];
rz(-1.2484115) q[3];
sx q[3];
rz(0.65402189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3205545) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(-1.5041171) q[2];
rz(1.5004246) q[3];
sx q[3];
rz(-1.0502366) q[3];
sx q[3];
rz(0.27568451) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42030537) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(0.48459184) q[0];
rz(1.8991607) q[1];
sx q[1];
rz(-0.74598765) q[1];
sx q[1];
rz(-2.7925083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32767998) q[0];
sx q[0];
rz(-0.55024863) q[0];
sx q[0];
rz(1.6525066) q[0];
x q[1];
rz(-0.075320638) q[2];
sx q[2];
rz(-1.9131129) q[2];
sx q[2];
rz(1.3665733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4306077) q[1];
sx q[1];
rz(-0.29919708) q[1];
sx q[1];
rz(-0.69327142) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7233299) q[3];
sx q[3];
rz(-2.1304479) q[3];
sx q[3];
rz(2.5734316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7829973) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(2.6915468) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(0.1990327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49804509) q[0];
sx q[0];
rz(-1.4692551) q[0];
sx q[0];
rz(0.099844649) q[0];
rz(-1.8210583) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(0.60755306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4293885) q[0];
sx q[0];
rz(-0.11183248) q[0];
sx q[0];
rz(0.73356103) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25600146) q[2];
sx q[2];
rz(-2.5814272) q[2];
sx q[2];
rz(-1.7988811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2098908) q[1];
sx q[1];
rz(-1.7875515) q[1];
sx q[1];
rz(-2.5561129) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9266379) q[3];
sx q[3];
rz(-0.27493335) q[3];
sx q[3];
rz(2.7235018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9472092) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-2.6109931) q[2];
rz(0.44140205) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(-0.94943625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7049103) q[0];
sx q[0];
rz(-1.7969776) q[0];
sx q[0];
rz(-0.5740903) q[0];
rz(2.2615945) q[1];
sx q[1];
rz(-2.0206385) q[1];
sx q[1];
rz(-1.3425739) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18595565) q[0];
sx q[0];
rz(-1.370472) q[0];
sx q[0];
rz(-2.0485417) q[0];
x q[1];
rz(1.9832575) q[2];
sx q[2];
rz(-1.165373) q[2];
sx q[2];
rz(2.6087922) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7751392) q[1];
sx q[1];
rz(-2.7045858) q[1];
sx q[1];
rz(-0.98458146) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4785512) q[3];
sx q[3];
rz(-1.6812857) q[3];
sx q[3];
rz(-0.35943951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67419702) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(0.44710844) q[2];
rz(2.1111264) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36987385) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(-0.41279992) q[0];
rz(-2.0665456) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(-0.80361754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5301054) q[0];
sx q[0];
rz(-2.3854333) q[0];
sx q[0];
rz(-2.3458781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56411479) q[2];
sx q[2];
rz(-1.0061641) q[2];
sx q[2];
rz(-2.5052793) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10104837) q[1];
sx q[1];
rz(-0.78436995) q[1];
sx q[1];
rz(-3.1215828) q[1];
rz(-0.17937029) q[3];
sx q[3];
rz(-1.6458047) q[3];
sx q[3];
rz(-1.6841494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2794118) q[2];
sx q[2];
rz(-1.923424) q[2];
sx q[2];
rz(1.5009521) q[2];
rz(-0.62266478) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(-1.6525432) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72614661) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(2.5836482) q[0];
rz(0.56124148) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(1.7254613) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5197882) q[0];
sx q[0];
rz(-2.0178231) q[0];
sx q[0];
rz(2.8267953) q[0];
x q[1];
rz(0.29406644) q[2];
sx q[2];
rz(-1.9379788) q[2];
sx q[2];
rz(2.4096556) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3741403) q[1];
sx q[1];
rz(-1.6686286) q[1];
sx q[1];
rz(1.2409452) q[1];
x q[2];
rz(-1.775557) q[3];
sx q[3];
rz(-0.8258709) q[3];
sx q[3];
rz(1.5958169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2531565) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-2.1412795) q[2];
rz(-0.12217626) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(2.9848671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67824739) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(-0.004322411) q[0];
rz(2.9777572) q[1];
sx q[1];
rz(-0.42525735) q[1];
sx q[1];
rz(1.77553) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23200792) q[0];
sx q[0];
rz(-0.9238832) q[0];
sx q[0];
rz(0.034897403) q[0];
rz(-2.7694584) q[2];
sx q[2];
rz(-1.2805802) q[2];
sx q[2];
rz(-1.2843101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.26936618) q[1];
sx q[1];
rz(-0.25372836) q[1];
sx q[1];
rz(0.14464186) q[1];
x q[2];
rz(2.4716464) q[3];
sx q[3];
rz(-0.94196999) q[3];
sx q[3];
rz(-0.82312246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(2.5223993) q[2];
rz(-2.5891417) q[3];
sx q[3];
rz(-1.3055472) q[3];
sx q[3];
rz(2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7904952) q[0];
sx q[0];
rz(-0.18432291) q[0];
sx q[0];
rz(2.6628394) q[0];
rz(-1.9888196) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(2.1943888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21175948) q[0];
sx q[0];
rz(-0.18778983) q[0];
sx q[0];
rz(-1.7540356) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9302773) q[2];
sx q[2];
rz(-0.61865846) q[2];
sx q[2];
rz(1.2115237) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45631623) q[1];
sx q[1];
rz(-0.79689491) q[1];
sx q[1];
rz(-2.0048098) q[1];
x q[2];
rz(-0.60444215) q[3];
sx q[3];
rz(-0.42057188) q[3];
sx q[3];
rz(2.4673691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(-2.0742778) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0530479) q[0];
sx q[0];
rz(-1.6041258) q[0];
sx q[0];
rz(-1.5969101) q[0];
rz(2.0883941) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(1.5736035) q[2];
sx q[2];
rz(-1.8127828) q[2];
sx q[2];
rz(0.91290963) q[2];
rz(1.347318) q[3];
sx q[3];
rz(-2.1383994) q[3];
sx q[3];
rz(2.4169328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
