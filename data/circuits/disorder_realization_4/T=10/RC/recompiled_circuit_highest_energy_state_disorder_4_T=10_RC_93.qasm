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
rz(-1.1308489) q[0];
sx q[0];
rz(-1.3698438) q[0];
sx q[0];
rz(-0.1432336) q[0];
rz(-0.020996006) q[1];
sx q[1];
rz(-0.58240533) q[1];
sx q[1];
rz(-2.6666759) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015939322) q[0];
sx q[0];
rz(-1.6467735) q[0];
sx q[0];
rz(-1.6267908) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036587759) q[2];
sx q[2];
rz(-0.6152336) q[2];
sx q[2];
rz(1.7177291) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2432952) q[1];
sx q[1];
rz(-1.6800632) q[1];
sx q[1];
rz(2.3508049) q[1];
rz(2.0359614) q[3];
sx q[3];
rz(-1.1316) q[3];
sx q[3];
rz(0.10094563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.016606) q[2];
sx q[2];
rz(-1.6281444) q[2];
sx q[2];
rz(0.58546698) q[2];
rz(-2.3440907) q[3];
sx q[3];
rz(-1.5944642) q[3];
sx q[3];
rz(-0.40369478) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1095235) q[0];
sx q[0];
rz(-1.7970947) q[0];
sx q[0];
rz(2.4865785) q[0];
rz(-3.1380999) q[1];
sx q[1];
rz(-0.73492903) q[1];
sx q[1];
rz(3.1266812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61136041) q[0];
sx q[0];
rz(-1.6515628) q[0];
sx q[0];
rz(1.5967036) q[0];
rz(-pi) q[1];
rz(-2.0582366) q[2];
sx q[2];
rz(-2.3771046) q[2];
sx q[2];
rz(-2.9202094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0449734) q[1];
sx q[1];
rz(-2.3311021) q[1];
sx q[1];
rz(-1.0721447) q[1];
x q[2];
rz(-0.064830975) q[3];
sx q[3];
rz(-2.3309543) q[3];
sx q[3];
rz(-0.027934542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1917176) q[2];
sx q[2];
rz(-1.8052552) q[2];
sx q[2];
rz(-0.28676644) q[2];
rz(3.1112572) q[3];
sx q[3];
rz(-1.5807296) q[3];
sx q[3];
rz(-1.102977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86998087) q[0];
sx q[0];
rz(-2.5219707) q[0];
sx q[0];
rz(1.1645338) q[0];
rz(-0.33186913) q[1];
sx q[1];
rz(-1.7121366) q[1];
sx q[1];
rz(-1.9934995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6004922) q[0];
sx q[0];
rz(-2.9229188) q[0];
sx q[0];
rz(-1.4906916) q[0];
rz(1.2562468) q[2];
sx q[2];
rz(-0.49731055) q[2];
sx q[2];
rz(1.0337551) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61473523) q[1];
sx q[1];
rz(-1.8737443) q[1];
sx q[1];
rz(-0.2016068) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3018191) q[3];
sx q[3];
rz(-2.9570974) q[3];
sx q[3];
rz(0.46637812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0120734) q[2];
sx q[2];
rz(-1.4612863) q[2];
sx q[2];
rz(0.47895437) q[2];
rz(0.78914133) q[3];
sx q[3];
rz(-2.456587) q[3];
sx q[3];
rz(1.7710549) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1016462) q[0];
sx q[0];
rz(-0.78831035) q[0];
sx q[0];
rz(-2.3729861) q[0];
rz(2.6885314) q[1];
sx q[1];
rz(-2.1078347) q[1];
sx q[1];
rz(-0.57910848) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3192037) q[0];
sx q[0];
rz(-1.3487848) q[0];
sx q[0];
rz(-0.25698203) q[0];
rz(-pi) q[1];
rz(-1.0592346) q[2];
sx q[2];
rz(-1.2280012) q[2];
sx q[2];
rz(-2.2276218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4873638) q[1];
sx q[1];
rz(-2.0647065) q[1];
sx q[1];
rz(3.0904675) q[1];
rz(1.7594537) q[3];
sx q[3];
rz(-1.9629729) q[3];
sx q[3];
rz(-2.3463263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1733178) q[2];
sx q[2];
rz(-1.7705132) q[2];
sx q[2];
rz(-1.0260065) q[2];
rz(-0.46938986) q[3];
sx q[3];
rz(-2.1211076) q[3];
sx q[3];
rz(2.6494086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39407179) q[0];
sx q[0];
rz(-1.6265656) q[0];
sx q[0];
rz(-2.1409905) q[0];
rz(-1.8383149) q[1];
sx q[1];
rz(-2.3293827) q[1];
sx q[1];
rz(-2.2607048) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7422646) q[0];
sx q[0];
rz(-2.3652746) q[0];
sx q[0];
rz(-0.72350435) q[0];
rz(-2.4816031) q[2];
sx q[2];
rz(-1.7777598) q[2];
sx q[2];
rz(2.3530838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.31277675) q[1];
sx q[1];
rz(-2.34561) q[1];
sx q[1];
rz(-0.46013855) q[1];
x q[2];
rz(-0.88663574) q[3];
sx q[3];
rz(-0.58707844) q[3];
sx q[3];
rz(2.1140703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96054968) q[2];
sx q[2];
rz(-2.0281823) q[2];
sx q[2];
rz(1.95365) q[2];
rz(-0.12399593) q[3];
sx q[3];
rz(-1.32722) q[3];
sx q[3];
rz(2.7441062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7749629) q[0];
sx q[0];
rz(-0.22659817) q[0];
sx q[0];
rz(-2.9887548) q[0];
rz(-2.9951908) q[1];
sx q[1];
rz(-1.791879) q[1];
sx q[1];
rz(0.70405594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8237127) q[0];
sx q[0];
rz(-1.6005662) q[0];
sx q[0];
rz(-1.5882306) q[0];
rz(-0.91135773) q[2];
sx q[2];
rz(-1.8375085) q[2];
sx q[2];
rz(0.12202036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4761281) q[1];
sx q[1];
rz(-1.4515467) q[1];
sx q[1];
rz(0.075983451) q[1];
rz(-pi) q[2];
rz(2.0364159) q[3];
sx q[3];
rz(-2.4676968) q[3];
sx q[3];
rz(-0.020581882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0600837) q[2];
sx q[2];
rz(-2.1174105) q[2];
sx q[2];
rz(-2.7626959) q[2];
rz(1.0259519) q[3];
sx q[3];
rz(-2.43695) q[3];
sx q[3];
rz(1.1161233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41554552) q[0];
sx q[0];
rz(-1.3274095) q[0];
sx q[0];
rz(-0.56336796) q[0];
rz(-2.7401961) q[1];
sx q[1];
rz(-1.0029663) q[1];
sx q[1];
rz(2.4664403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7679325) q[0];
sx q[0];
rz(-1.5392173) q[0];
sx q[0];
rz(-3.1302559) q[0];
rz(-pi) q[1];
rz(-1.364351) q[2];
sx q[2];
rz(-0.9725625) q[2];
sx q[2];
rz(0.75144671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4298378) q[1];
sx q[1];
rz(-2.7378509) q[1];
sx q[1];
rz(-2.212119) q[1];
rz(0.83480699) q[3];
sx q[3];
rz(-1.2865598) q[3];
sx q[3];
rz(0.72407297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99951619) q[2];
sx q[2];
rz(-0.88366214) q[2];
sx q[2];
rz(-2.7189972) q[2];
rz(-2.8790867) q[3];
sx q[3];
rz(-1.0650977) q[3];
sx q[3];
rz(1.0814103) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074908681) q[0];
sx q[0];
rz(-2.1796362) q[0];
sx q[0];
rz(0.99895507) q[0];
rz(-1.2133489) q[1];
sx q[1];
rz(-1.9478925) q[1];
sx q[1];
rz(0.34128571) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0564108) q[0];
sx q[0];
rz(-0.23427948) q[0];
sx q[0];
rz(-2.5368669) q[0];
rz(-2.2407124) q[2];
sx q[2];
rz(-1.3645126) q[2];
sx q[2];
rz(-1.8697949) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7012257) q[1];
sx q[1];
rz(-1.6858674) q[1];
sx q[1];
rz(-0.23903317) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058505614) q[3];
sx q[3];
rz(-1.2924177) q[3];
sx q[3];
rz(-0.52245058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21218941) q[2];
sx q[2];
rz(-1.5985649) q[2];
sx q[2];
rz(-2.8821778) q[2];
rz(-0.2855531) q[3];
sx q[3];
rz(-2.4694337) q[3];
sx q[3];
rz(-2.0937505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015273) q[0];
sx q[0];
rz(-2.6492388) q[0];
sx q[0];
rz(1.0858076) q[0];
rz(-1.0049413) q[1];
sx q[1];
rz(-1.137038) q[1];
sx q[1];
rz(2.5850632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426612) q[0];
sx q[0];
rz(-3.001881) q[0];
sx q[0];
rz(1.8675787) q[0];
x q[1];
rz(-0.6350216) q[2];
sx q[2];
rz(-2.1974034) q[2];
sx q[2];
rz(-2.602586) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5270015) q[1];
sx q[1];
rz(-1.6629991) q[1];
sx q[1];
rz(2.8648977) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79995221) q[3];
sx q[3];
rz(-2.4728007) q[3];
sx q[3];
rz(0.92000414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0050547) q[2];
sx q[2];
rz(-1.7359066) q[2];
sx q[2];
rz(1.271099) q[2];
rz(2.2942885) q[3];
sx q[3];
rz(-1.3264791) q[3];
sx q[3];
rz(-0.29720753) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6274881) q[0];
sx q[0];
rz(-1.4493554) q[0];
sx q[0];
rz(0.76747146) q[0];
rz(-0.96495676) q[1];
sx q[1];
rz(-1.283353) q[1];
sx q[1];
rz(-2.8585785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3393769) q[0];
sx q[0];
rz(-2.2664323) q[0];
sx q[0];
rz(0.2494158) q[0];
rz(-pi) q[1];
rz(0.04472132) q[2];
sx q[2];
rz(-0.43202094) q[2];
sx q[2];
rz(0.17145874) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.231608) q[1];
sx q[1];
rz(-2.5043813) q[1];
sx q[1];
rz(-0.63225327) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55693407) q[3];
sx q[3];
rz(-0.86268988) q[3];
sx q[3];
rz(-2.1754337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1538126) q[2];
sx q[2];
rz(-0.38975468) q[2];
sx q[2];
rz(-1.7392996) q[2];
rz(0.66579372) q[3];
sx q[3];
rz(-1.2716764) q[3];
sx q[3];
rz(-1.3027035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5182198) q[0];
sx q[0];
rz(-1.8898531) q[0];
sx q[0];
rz(1.3367597) q[0];
rz(-1.549859) q[1];
sx q[1];
rz(-1.5668329) q[1];
sx q[1];
rz(-0.11650539) q[1];
rz(0.58407513) q[2];
sx q[2];
rz(-1.948195) q[2];
sx q[2];
rz(0.52577166) q[2];
rz(0.20768349) q[3];
sx q[3];
rz(-0.80266914) q[3];
sx q[3];
rz(-1.707984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
