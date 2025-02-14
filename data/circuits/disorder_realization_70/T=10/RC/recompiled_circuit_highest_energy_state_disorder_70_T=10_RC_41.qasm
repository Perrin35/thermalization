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
rz(1.7477859) q[0];
sx q[0];
rz(4.0417606) q[0];
sx q[0];
rz(10.834122) q[0];
rz(0.70977587) q[1];
sx q[1];
rz(3.4748454) q[1];
sx q[1];
rz(9.8069053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724498) q[0];
sx q[0];
rz(-0.74688321) q[0];
sx q[0];
rz(-0.19799671) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13488954) q[2];
sx q[2];
rz(-1.3185878) q[2];
sx q[2];
rz(1.8602961) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.00013) q[1];
sx q[1];
rz(-1.9381818) q[1];
sx q[1];
rz(-0.39333435) q[1];
rz(-pi) q[2];
rz(0.11973937) q[3];
sx q[3];
rz(-2.2594707) q[3];
sx q[3];
rz(0.39018351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0673151) q[2];
sx q[2];
rz(-1.6396319) q[2];
sx q[2];
rz(3.066646) q[2];
rz(-1.8986757) q[3];
sx q[3];
rz(-2.8800745) q[3];
sx q[3];
rz(0.3956795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0089834) q[0];
sx q[0];
rz(-0.41985303) q[0];
sx q[0];
rz(-2.5123151) q[0];
rz(-2.3043326) q[1];
sx q[1];
rz(-2.3910797) q[1];
sx q[1];
rz(-2.5619521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.048174) q[0];
sx q[0];
rz(-2.3701128) q[0];
sx q[0];
rz(-0.29399612) q[0];
rz(-pi) q[1];
rz(1.3084035) q[2];
sx q[2];
rz(-1.9540522) q[2];
sx q[2];
rz(2.8079536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15687277) q[1];
sx q[1];
rz(-2.619209) q[1];
sx q[1];
rz(-0.70981437) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0218815) q[3];
sx q[3];
rz(-0.16489794) q[3];
sx q[3];
rz(-0.69530693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89977449) q[2];
sx q[2];
rz(-0.16336975) q[2];
sx q[2];
rz(0.74431288) q[2];
rz(2.7742079) q[3];
sx q[3];
rz(-1.8495879) q[3];
sx q[3];
rz(-1.415409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431417) q[0];
sx q[0];
rz(-2.1911868) q[0];
sx q[0];
rz(-2.1344192) q[0];
rz(0.011064359) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(0.99753582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38356203) q[0];
sx q[0];
rz(-0.13995384) q[0];
sx q[0];
rz(-1.3849694) q[0];
x q[1];
rz(-1.4513955) q[2];
sx q[2];
rz(-0.7660256) q[2];
sx q[2];
rz(1.8325782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1147194) q[1];
sx q[1];
rz(-2.623154) q[1];
sx q[1];
rz(-3.1264831) q[1];
rz(-pi) q[2];
rz(3.1155048) q[3];
sx q[3];
rz(-2.6296277) q[3];
sx q[3];
rz(-1.0796368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1204388) q[2];
sx q[2];
rz(-0.28205559) q[2];
sx q[2];
rz(2.2998478) q[2];
rz(-0.65819955) q[3];
sx q[3];
rz(-0.87116146) q[3];
sx q[3];
rz(-0.021520821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4206674) q[0];
sx q[0];
rz(-3.0166716) q[0];
sx q[0];
rz(-0.7290054) q[0];
rz(2.3782102) q[1];
sx q[1];
rz(-0.63012505) q[1];
sx q[1];
rz(-0.36300945) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1141407) q[0];
sx q[0];
rz(-1.8233366) q[0];
sx q[0];
rz(2.8414498) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55018376) q[2];
sx q[2];
rz(-0.90298684) q[2];
sx q[2];
rz(0.097152348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69907585) q[1];
sx q[1];
rz(-1.867379) q[1];
sx q[1];
rz(-1.7188598) q[1];
rz(2.7458526) q[3];
sx q[3];
rz(-3.0043663) q[3];
sx q[3];
rz(-1.9326222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18569854) q[2];
sx q[2];
rz(-2.318435) q[2];
sx q[2];
rz(-1.4996747) q[2];
rz(0.58756346) q[3];
sx q[3];
rz(-0.9891808) q[3];
sx q[3];
rz(2.6232918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3041621) q[0];
sx q[0];
rz(-2.4346133) q[0];
sx q[0];
rz(0.28513232) q[0];
rz(-2.888491) q[1];
sx q[1];
rz(-1.0292091) q[1];
sx q[1];
rz(1.0458127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8305215) q[0];
sx q[0];
rz(-1.9609299) q[0];
sx q[0];
rz(-2.5907787) q[0];
rz(-pi) q[1];
rz(-0.4407119) q[2];
sx q[2];
rz(-1.1870702) q[2];
sx q[2];
rz(-1.1659932) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48165136) q[1];
sx q[1];
rz(-1.5143947) q[1];
sx q[1];
rz(-2.7475471) q[1];
x q[2];
rz(-2.8941514) q[3];
sx q[3];
rz(-2.4486922) q[3];
sx q[3];
rz(-0.019236658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4046341) q[2];
sx q[2];
rz(-2.0817231) q[2];
sx q[2];
rz(-0.57445478) q[2];
rz(0.61070853) q[3];
sx q[3];
rz(-2.6210531) q[3];
sx q[3];
rz(-1.0569388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.33559281) q[0];
sx q[0];
rz(-1.7570423) q[0];
sx q[0];
rz(0.35032508) q[0];
rz(2.8490745) q[1];
sx q[1];
rz(-0.12648335) q[1];
sx q[1];
rz(-0.77936053) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.821601) q[0];
sx q[0];
rz(-1.6896833) q[0];
sx q[0];
rz(-3.0827808) q[0];
rz(0.47996491) q[2];
sx q[2];
rz(-1.4128608) q[2];
sx q[2];
rz(1.8851033) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8694589) q[1];
sx q[1];
rz(-1.773352) q[1];
sx q[1];
rz(0.32702469) q[1];
rz(-pi) q[2];
rz(2.8661714) q[3];
sx q[3];
rz(-2.6044327) q[3];
sx q[3];
rz(-0.28322434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23955762) q[2];
sx q[2];
rz(-1.5152479) q[2];
sx q[2];
rz(0.96417344) q[2];
rz(2.9686109) q[3];
sx q[3];
rz(-1.0123342) q[3];
sx q[3];
rz(3.0103736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5445589) q[0];
sx q[0];
rz(-1.9439161) q[0];
sx q[0];
rz(-0.069393754) q[0];
rz(1.767905) q[1];
sx q[1];
rz(-1.8927788) q[1];
sx q[1];
rz(-0.51838851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4189506) q[0];
sx q[0];
rz(-1.4012453) q[0];
sx q[0];
rz(0.8954312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80194998) q[2];
sx q[2];
rz(-0.63849245) q[2];
sx q[2];
rz(-1.0494159) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9897108) q[1];
sx q[1];
rz(-0.1406142) q[1];
sx q[1];
rz(2.9563044) q[1];
rz(-1.5126918) q[3];
sx q[3];
rz(-1.4062506) q[3];
sx q[3];
rz(-1.8376415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9531276) q[2];
sx q[2];
rz(-1.943482) q[2];
sx q[2];
rz(-0.24448621) q[2];
rz(1.2178347) q[3];
sx q[3];
rz(-0.094450258) q[3];
sx q[3];
rz(2.794246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0511047) q[0];
sx q[0];
rz(-0.29260391) q[0];
sx q[0];
rz(-0.69945139) q[0];
rz(-2.4199016) q[1];
sx q[1];
rz(-1.3612008) q[1];
sx q[1];
rz(-1.9500505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724567) q[0];
sx q[0];
rz(-1.4459608) q[0];
sx q[0];
rz(-0.21376507) q[0];
rz(-1.4100685) q[2];
sx q[2];
rz(-1.3597915) q[2];
sx q[2];
rz(-0.69549045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9758975) q[1];
sx q[1];
rz(-1.6553406) q[1];
sx q[1];
rz(1.0313453) q[1];
x q[2];
rz(-1.9189214) q[3];
sx q[3];
rz(-1.1468107) q[3];
sx q[3];
rz(-2.4676187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2193853) q[2];
sx q[2];
rz(-0.81874138) q[2];
sx q[2];
rz(-1.7259664) q[2];
rz(-2.7042232) q[3];
sx q[3];
rz(-0.24961095) q[3];
sx q[3];
rz(2.6387446) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1943844) q[0];
sx q[0];
rz(-1.9358862) q[0];
sx q[0];
rz(-2.6956287) q[0];
rz(-1.3392316) q[1];
sx q[1];
rz(-0.65473348) q[1];
sx q[1];
rz(-1.9816678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98731536) q[0];
sx q[0];
rz(-1.2270708) q[0];
sx q[0];
rz(2.4489342) q[0];
rz(-1.3092988) q[2];
sx q[2];
rz(-1.6435677) q[2];
sx q[2];
rz(-2.6852599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36148188) q[1];
sx q[1];
rz(-1.0358397) q[1];
sx q[1];
rz(-1.8423716) q[1];
rz(0.91083093) q[3];
sx q[3];
rz(-2.0487222) q[3];
sx q[3];
rz(2.7895989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40598536) q[2];
sx q[2];
rz(-1.0286101) q[2];
sx q[2];
rz(2.410991) q[2];
rz(-0.90100151) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(-3.1072531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63448298) q[0];
sx q[0];
rz(-0.49840885) q[0];
sx q[0];
rz(-2.679017) q[0];
rz(1.0628465) q[1];
sx q[1];
rz(-1.3674419) q[1];
sx q[1];
rz(3.0752693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6505867) q[0];
sx q[0];
rz(-2.1438476) q[0];
sx q[0];
rz(-1.9377524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0424329) q[2];
sx q[2];
rz(-2.7688469) q[2];
sx q[2];
rz(-0.78102222) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5196701) q[1];
sx q[1];
rz(-2.1035476) q[1];
sx q[1];
rz(0.68343917) q[1];
rz(-pi) q[2];
rz(-1.0995501) q[3];
sx q[3];
rz(-0.46898919) q[3];
sx q[3];
rz(2.2820306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1041717) q[2];
sx q[2];
rz(-2.6080242) q[2];
sx q[2];
rz(-1.2332234) q[2];
rz(2.6836266) q[3];
sx q[3];
rz(-0.27126867) q[3];
sx q[3];
rz(0.39877322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11252277) q[0];
sx q[0];
rz(-1.4334913) q[0];
sx q[0];
rz(1.7423472) q[0];
rz(-2.9872672) q[1];
sx q[1];
rz(-1.7185153) q[1];
sx q[1];
rz(-1.1850866) q[1];
rz(0.90171705) q[2];
sx q[2];
rz(-0.48090413) q[2];
sx q[2];
rz(3.0851875) q[2];
rz(-0.88944351) q[3];
sx q[3];
rz(-0.69557299) q[3];
sx q[3];
rz(2.099149) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
