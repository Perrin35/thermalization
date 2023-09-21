OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5320839) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-2.8340333) q[0];
x q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(-0.27054271) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38274256) q[1];
sx q[1];
rz(-1.1945801) q[1];
sx q[1];
rz(1.7099027) q[1];
x q[2];
rz(-2.9620142) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(-1.4048525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75630674) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(0.73861648) q[0];
rz(2.0915394) q[2];
sx q[2];
rz(-0.69486952) q[2];
sx q[2];
rz(-0.75887242) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0305811) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(0.76534033) q[1];
x q[2];
rz(-2.9566544) q[3];
sx q[3];
rz(-1.8293081) q[3];
sx q[3];
rz(-1.0227433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.5154243) q[0];
rz(0.60107636) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9515581) q[0];
sx q[0];
rz(-2.1265731) q[0];
sx q[0];
rz(0.65727289) q[0];
rz(2.2268779) q[2];
sx q[2];
rz(-1.9064184) q[2];
sx q[2];
rz(0.4682954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6669238) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(1.0679507) q[1];
rz(-pi) q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.2325866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6757641) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(-2.7990544) q[0];
rz(0.68508673) q[2];
sx q[2];
rz(-1.6685467) q[2];
sx q[2];
rz(0.098066559) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1102317) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(-1.5178174) q[1];
x q[2];
rz(-1.2918606) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(-3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.189165) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2244959) q[0];
sx q[0];
rz(-1.3670237) q[0];
sx q[0];
rz(-2.8208371) q[0];
rz(-pi) q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(0.86415926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89228499) q[1];
sx q[1];
rz(-1.882949) q[1];
sx q[1];
rz(2.409163) q[1];
rz(-pi) q[2];
rz(1.0088483) q[3];
sx q[3];
rz(-1.7147439) q[3];
sx q[3];
rz(-0.89509237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(0.67374054) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(1.649883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30116044) q[0];
sx q[0];
rz(-1.4727117) q[0];
sx q[0];
rz(0.43099404) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6368124) q[2];
sx q[2];
rz(-0.74540388) q[2];
sx q[2];
rz(2.8842852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3867492) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(2.2350603) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10156472) q[3];
sx q[3];
rz(-1.933681) q[3];
sx q[3];
rz(-0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(0.86589083) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(3.1076028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86173979) q[0];
sx q[0];
rz(-1.732677) q[0];
sx q[0];
rz(-0.78761657) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(-2.5069782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1696724) q[1];
sx q[1];
rz(-1.6213413) q[1];
sx q[1];
rz(0.89269841) q[1];
rz(2.5266685) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(1.8254335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.4454909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60364265) q[0];
sx q[0];
rz(-2.0041487) q[0];
sx q[0];
rz(-0.005714697) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5567899) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(-0.98758299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1576924) q[1];
sx q[1];
rz(-1.0219814) q[1];
sx q[1];
rz(-1.3859315) q[1];
x q[2];
rz(0.66283488) q[3];
sx q[3];
rz(-2.088306) q[3];
sx q[3];
rz(-2.585632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(2.8919162) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99684925) q[0];
sx q[0];
rz(-1.0486756) q[0];
sx q[0];
rz(1.7168619) q[0];
x q[1];
rz(0.0094527761) q[2];
sx q[2];
rz(-1.7477034) q[2];
sx q[2];
rz(-1.1284019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8621091) q[1];
sx q[1];
rz(-2.9276507) q[1];
sx q[1];
rz(-1.2941542) q[1];
x q[2];
rz(0.59659776) q[3];
sx q[3];
rz(-1.3571686) q[3];
sx q[3];
rz(-1.9006157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(1.6419798) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23696454) q[0];
sx q[0];
rz(-1.4693854) q[0];
sx q[0];
rz(-1.2786091) q[0];
x q[1];
rz(2.3775616) q[2];
sx q[2];
rz(-1.6198297) q[2];
sx q[2];
rz(-1.3438091) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2227576) q[1];
sx q[1];
rz(-1.7822052) q[1];
sx q[1];
rz(1.097015) q[1];
x q[2];
rz(-0.0025000574) q[3];
sx q[3];
rz(-1.7435939) q[3];
sx q[3];
rz(2.1713184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-2.1949027) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(0.46438607) q[2];
sx q[2];
rz(-1.1190363) q[2];
sx q[2];
rz(-2.6418532) q[2];
rz(-1.0740888) q[3];
sx q[3];
rz(-0.91377331) q[3];
sx q[3];
rz(2.5627315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];