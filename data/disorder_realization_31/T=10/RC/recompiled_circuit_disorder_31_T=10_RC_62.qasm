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
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9561477) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(-2.3953715) q[0];
x q[1];
rz(-1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(0.27054271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.018651389) q[1];
sx q[1];
rz(-2.7416347) q[1];
sx q[1];
rz(0.33756983) q[1];
x q[2];
rz(1.6928715) q[3];
sx q[3];
rz(-0.97994084) q[3];
sx q[3];
rz(1.6216244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(0.17949417) q[2];
rz(1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(1.1546086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1240631) q[0];
sx q[0];
rz(-1.5895827) q[0];
sx q[0];
rz(-1.5879052) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0500533) q[2];
sx q[2];
rz(-0.69486952) q[2];
sx q[2];
rz(2.3827202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6007538) q[1];
sx q[1];
rz(-0.80889091) q[1];
sx q[1];
rz(1.4536588) q[1];
x q[2];
rz(-1.8335908) q[3];
sx q[3];
rz(-1.7495219) q[3];
sx q[3];
rz(-2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(-2.6702821) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(2.0498958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(-0.65727289) q[0];
x q[1];
rz(0.41468427) q[2];
sx q[2];
rz(-2.1846002) q[2];
sx q[2];
rz(-2.2874122) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97869067) q[1];
sx q[1];
rz(-0.86135094) q[1];
sx q[1];
rz(2.6497926) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5063498) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(0.32999048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096075637) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(-1.7128574) q[0];
rz(-2.9878346) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(1.549987) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.51018184) q[1];
sx q[1];
rz(-1.5522172) q[1];
sx q[1];
rz(-1.2127962) q[1];
rz(-pi) q[2];
x q[2];
rz(1.849732) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0115396) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(-1.0774353) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(-2.2713984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91709671) q[0];
sx q[0];
rz(-1.774569) q[0];
sx q[0];
rz(-0.32075551) q[0];
rz(1.7587897) q[2];
sx q[2];
rz(-2.2586939) q[2];
sx q[2];
rz(-2.5685513) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1340027) q[1];
sx q[1];
rz(-0.78467272) q[1];
sx q[1];
rz(2.6919634) q[1];
rz(-pi) q[2];
rz(1.0088483) q[3];
sx q[3];
rz(-1.7147439) q[3];
sx q[3];
rz(-0.89509237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-2.4678521) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(1.4917096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8404322) q[0];
sx q[0];
rz(-1.668881) q[0];
sx q[0];
rz(0.43099404) q[0];
rz(-pi) q[1];
rz(-0.67955534) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(-0.92781767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41960934) q[1];
sx q[1];
rz(-1.0227385) q[1];
sx q[1];
rz(-2.4649058) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0400279) q[3];
sx q[3];
rz(-1.933681) q[3];
sx q[3];
rz(-0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-2.2475524) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-3.1076028) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995178) q[0];
sx q[0];
rz(-2.3410428) q[0];
sx q[0];
rz(0.22649015) q[0];
rz(-2.2529644) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(1.6737446) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46376343) q[1];
sx q[1];
rz(-2.4619108) q[1];
sx q[1];
rz(-1.6512647) q[1];
rz(-pi) q[2];
rz(-2.9339318) q[3];
sx q[3];
rz(-0.62519473) q[3];
sx q[3];
rz(3.0561662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(2.2650488) q[2];
rz(0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.4454909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59003419) q[0];
sx q[0];
rz(-2.708205) q[0];
sx q[0];
rz(1.5584459) q[0];
x q[1];
rz(2.5567899) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(0.98758299) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.48965028) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(-2.5850992) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94533841) q[3];
sx q[3];
rz(-1.006554) q[3];
sx q[3];
rz(2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.8898213) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315902) q[0];
sx q[0];
rz(-2.6012523) q[0];
sx q[0];
rz(-0.24775981) q[0];
x q[1];
rz(-1.6236213) q[2];
sx q[2];
rz(-0.17715684) q[2];
sx q[2];
rz(1.0747386) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.579638) q[1];
sx q[1];
rz(-1.5127752) q[1];
sx q[1];
rz(1.7768363) q[1];
rz(1.8272607) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(2.6687711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(-2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050215125) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6584872) q[0];
sx q[0];
rz(-0.30880901) q[0];
sx q[0];
rz(-1.2312066) q[0];
rz(2.3775616) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(-1.7977835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26350281) q[1];
sx q[1];
rz(-2.6260758) q[1];
sx q[1];
rz(-1.131119) q[1];
x q[2];
rz(1.5851192) q[3];
sx q[3];
rz(-0.17281547) q[3];
sx q[3];
rz(-2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88400921) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35836999) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(2.6772066) q[2];
sx q[2];
rz(-2.0225564) q[2];
sx q[2];
rz(0.49973942) q[2];
rz(-0.72017097) q[3];
sx q[3];
rz(-1.1838893) q[3];
sx q[3];
rz(0.67223687) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];