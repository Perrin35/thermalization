OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28461449) q[0];
sx q[0];
rz(-2.1224788) q[0];
sx q[0];
rz(0.3509699) q[0];
rz(-pi) q[1];
rz(-1.0640261) q[2];
sx q[2];
rz(-2.1760586) q[2];
sx q[2];
rz(1.93047) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0776699) q[1];
sx q[1];
rz(-0.7709255) q[1];
sx q[1];
rz(1.116712) q[1];
x q[2];
rz(1.6139908) q[3];
sx q[3];
rz(-2.3213263) q[3];
sx q[3];
rz(-0.84212069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2261752) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.0579257) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-0.57587409) q[0];
rz(-1.2469762) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(1.974568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.5123899) q[0];
sx q[0];
rz(-1.8286684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0300693) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(-0.18268798) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0547304) q[1];
sx q[1];
rz(-0.80838258) q[1];
sx q[1];
rz(3.0562835) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0608276) q[3];
sx q[3];
rz(-2.1809275) q[3];
sx q[3];
rz(-0.78435635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4784933) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944511) q[0];
sx q[0];
rz(-1.6269636) q[0];
sx q[0];
rz(2.4980314) q[0];
x q[1];
rz(2.0929167) q[2];
sx q[2];
rz(-0.45959696) q[2];
sx q[2];
rz(2.2082579) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6940569) q[1];
sx q[1];
rz(-2.5390365) q[1];
sx q[1];
rz(-1.7130997) q[1];
rz(-1.2546) q[3];
sx q[3];
rz(-0.71101515) q[3];
sx q[3];
rz(2.7544114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-2.2423559) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(-0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(2.7096601) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(0.63582173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83946562) q[0];
sx q[0];
rz(-1.3913904) q[0];
sx q[0];
rz(-1.0666215) q[0];
rz(-1.3859149) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(-0.76279574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3969877) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(2.8810487) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5084247) q[3];
sx q[3];
rz(-0.35053262) q[3];
sx q[3];
rz(1.7115953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(0.27854663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771839) q[0];
sx q[0];
rz(-1.9946949) q[0];
sx q[0];
rz(2.3197078) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6804382) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(1.4102175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3814195) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(2.6890523) q[1];
rz(-pi) q[2];
rz(-2.3260818) q[3];
sx q[3];
rz(-1.9933356) q[3];
sx q[3];
rz(0.58327196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0329523) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(-3.0335398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5393339) q[0];
sx q[0];
rz(-0.56319153) q[0];
sx q[0];
rz(-2.5186033) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4322386) q[2];
sx q[2];
rz(-1.9493305) q[2];
sx q[2];
rz(-2.0356503) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5080155) q[1];
sx q[1];
rz(-2.0532236) q[1];
sx q[1];
rz(-1.1303933) q[1];
rz(-pi) q[2];
rz(0.26515682) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(1.4491855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(1.8704869) q[2];
rz(0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69041396) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(0.50509392) q[0];
x q[1];
rz(-1.0929843) q[2];
sx q[2];
rz(-1.3078948) q[2];
sx q[2];
rz(2.773657) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.212008) q[1];
sx q[1];
rz(-1.7716904) q[1];
sx q[1];
rz(-2.4561873) q[1];
rz(1.4736389) q[3];
sx q[3];
rz(-0.81340862) q[3];
sx q[3];
rz(-2.017445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(2.612109) q[2];
rz(-0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(-2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.9876678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740828) q[0];
sx q[0];
rz(-0.65069288) q[0];
sx q[0];
rz(-3.0853737) q[0];
rz(-2.9032193) q[2];
sx q[2];
rz(-0.56406883) q[2];
sx q[2];
rz(0.43018815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0358419) q[1];
sx q[1];
rz(-0.81139794) q[1];
sx q[1];
rz(1.6559385) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.047394) q[3];
sx q[3];
rz(-2.4123441) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5899137) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8001051) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-2.5019116) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1146678) q[0];
sx q[0];
rz(-1.9244734) q[0];
sx q[0];
rz(2.7915519) q[0];
rz(-pi) q[1];
rz(0.77634546) q[2];
sx q[2];
rz(-1.7661957) q[2];
sx q[2];
rz(2.5459144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6744086) q[1];
sx q[1];
rz(-1.7519752) q[1];
sx q[1];
rz(-0.49023899) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0036003) q[3];
sx q[3];
rz(-0.90413168) q[3];
sx q[3];
rz(0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-2.6616667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878764) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(2.2399708) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9058814) q[2];
sx q[2];
rz(-1.7582298) q[2];
sx q[2];
rz(2.4311709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(-0.6154284) q[1];
rz(-pi) q[2];
rz(-2.2009833) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(1.9711232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(-1.3170362) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(-0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(0.97933979) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(-0.32817763) q[2];
sx q[2];
rz(-2.6371418) q[2];
sx q[2];
rz(-0.20124659) q[2];
rz(-0.76673037) q[3];
sx q[3];
rz(-0.37692108) q[3];
sx q[3];
rz(2.2212096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
