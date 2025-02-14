OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7810516) q[0];
sx q[0];
rz(-2.5268351) q[0];
sx q[0];
rz(1.0714666) q[0];
rz(2.7331424) q[1];
sx q[1];
rz(-0.93745679) q[1];
sx q[1];
rz(-1.6931005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1030185) q[0];
sx q[0];
rz(-1.9644613) q[0];
sx q[0];
rz(2.066924) q[0];
rz(1.229847) q[2];
sx q[2];
rz(-1.023479) q[2];
sx q[2];
rz(0.42423074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.422157) q[1];
sx q[1];
rz(-1.4220752) q[1];
sx q[1];
rz(1.6610816) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43464147) q[3];
sx q[3];
rz(-0.9081525) q[3];
sx q[3];
rz(1.545411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.19746) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(2.3485363) q[2];
rz(2.872725) q[3];
sx q[3];
rz(-1.3258508) q[3];
sx q[3];
rz(2.1813006) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307513) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(-2.261396) q[0];
rz(-0.63086069) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(-1.0989443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99115935) q[0];
sx q[0];
rz(-2.9637103) q[0];
sx q[0];
rz(-0.028602914) q[0];
rz(1.3395202) q[2];
sx q[2];
rz(-0.63098365) q[2];
sx q[2];
rz(-2.4914666) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6827755) q[1];
sx q[1];
rz(-1.7331682) q[1];
sx q[1];
rz(-1.1332676) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97084011) q[3];
sx q[3];
rz(-2.4106541) q[3];
sx q[3];
rz(2.746664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38689303) q[2];
sx q[2];
rz(-1.4852445) q[2];
sx q[2];
rz(-2.5999542) q[2];
rz(-0.70332876) q[3];
sx q[3];
rz(-2.7024305) q[3];
sx q[3];
rz(-0.37731236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2331053) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(-3.044686) q[0];
rz(2.5540409) q[1];
sx q[1];
rz(-0.89043003) q[1];
sx q[1];
rz(0.60120916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9438697) q[0];
sx q[0];
rz(-1.7515148) q[0];
sx q[0];
rz(-0.69036071) q[0];
x q[1];
rz(2.7695724) q[2];
sx q[2];
rz(-0.17401234) q[2];
sx q[2];
rz(3.0559828) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86089462) q[1];
sx q[1];
rz(-1.1105683) q[1];
sx q[1];
rz(-1.5621508) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1264301) q[3];
sx q[3];
rz(-2.4345102) q[3];
sx q[3];
rz(2.9418687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84843695) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(0.71451521) q[2];
rz(1.6589288) q[3];
sx q[3];
rz(-0.27479333) q[3];
sx q[3];
rz(-2.7627435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191384) q[0];
sx q[0];
rz(-1.2353354) q[0];
sx q[0];
rz(1.5041014) q[0];
rz(1.0927041) q[1];
sx q[1];
rz(-1.8605109) q[1];
sx q[1];
rz(-0.5853931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8849759) q[0];
sx q[0];
rz(-2.4451849) q[0];
sx q[0];
rz(0.43088669) q[0];
rz(-pi) q[1];
rz(-3.131065) q[2];
sx q[2];
rz(-1.6750355) q[2];
sx q[2];
rz(-0.57360211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8429381) q[1];
sx q[1];
rz(-2.2982592) q[1];
sx q[1];
rz(2.2266989) q[1];
rz(-pi) q[2];
rz(-2.0616509) q[3];
sx q[3];
rz(-0.98831359) q[3];
sx q[3];
rz(2.7459308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66712159) q[2];
sx q[2];
rz(-2.4690364) q[2];
sx q[2];
rz(2.7453864) q[2];
rz(-2.951156) q[3];
sx q[3];
rz(-2.2021144) q[3];
sx q[3];
rz(2.6207391) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025295479) q[0];
sx q[0];
rz(-0.46023661) q[0];
sx q[0];
rz(-0.39475557) q[0];
rz(-2.1341628) q[1];
sx q[1];
rz(-1.9837244) q[1];
sx q[1];
rz(-1.6068858) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0368113) q[0];
sx q[0];
rz(-0.80950538) q[0];
sx q[0];
rz(-2.7466989) q[0];
x q[1];
rz(1.2053648) q[2];
sx q[2];
rz(-0.88813587) q[2];
sx q[2];
rz(-0.78675573) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64749657) q[1];
sx q[1];
rz(-2.7733884) q[1];
sx q[1];
rz(1.0009264) q[1];
x q[2];
rz(2.6768489) q[3];
sx q[3];
rz(-1.9806974) q[3];
sx q[3];
rz(-1.003554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69636238) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(2.9635079) q[2];
rz(0.94545025) q[3];
sx q[3];
rz(-0.33792308) q[3];
sx q[3];
rz(-2.7054355) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6984542) q[0];
sx q[0];
rz(-2.4789424) q[0];
sx q[0];
rz(-2.9887001) q[0];
rz(2.1235509) q[1];
sx q[1];
rz(-2.1986304) q[1];
sx q[1];
rz(2.5315888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40326443) q[0];
sx q[0];
rz(-2.3685507) q[0];
sx q[0];
rz(-0.18456642) q[0];
rz(-pi) q[1];
rz(-1.583408) q[2];
sx q[2];
rz(-0.57390814) q[2];
sx q[2];
rz(1.2645666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2392366) q[1];
sx q[1];
rz(-0.78241759) q[1];
sx q[1];
rz(-0.33496015) q[1];
rz(-3.043706) q[3];
sx q[3];
rz(-1.4892508) q[3];
sx q[3];
rz(1.3697325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4623922) q[2];
sx q[2];
rz(-1.8588763) q[2];
sx q[2];
rz(-2.882615) q[2];
rz(0.14608832) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(2.5025388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6530957) q[0];
sx q[0];
rz(-0.86576068) q[0];
sx q[0];
rz(2.3897032) q[0];
rz(2.409626) q[1];
sx q[1];
rz(-1.7611793) q[1];
sx q[1];
rz(-0.52925777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37664304) q[0];
sx q[0];
rz(-1.8762454) q[0];
sx q[0];
rz(0.16965387) q[0];
rz(-pi) q[1];
x q[1];
rz(1.115832) q[2];
sx q[2];
rz(-1.6935529) q[2];
sx q[2];
rz(1.794874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0433474) q[1];
sx q[1];
rz(-1.6520513) q[1];
sx q[1];
rz(0.01290613) q[1];
rz(-pi) q[2];
rz(0.31792171) q[3];
sx q[3];
rz(-2.0057851) q[3];
sx q[3];
rz(3.0299713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25972128) q[2];
sx q[2];
rz(-1.0292116) q[2];
sx q[2];
rz(-0.79505801) q[2];
rz(-1.9184387) q[3];
sx q[3];
rz(-1.3431679) q[3];
sx q[3];
rz(-0.77643001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20599468) q[0];
sx q[0];
rz(-1.5933651) q[0];
sx q[0];
rz(-3.026631) q[0];
rz(3.0523172) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(2.9866536) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0099338) q[0];
sx q[0];
rz(-1.4408875) q[0];
sx q[0];
rz(0.1013426) q[0];
rz(-pi) q[1];
rz(2.2346588) q[2];
sx q[2];
rz(-1.1198545) q[2];
sx q[2];
rz(-2.5189357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7748002) q[1];
sx q[1];
rz(-1.8272562) q[1];
sx q[1];
rz(0.40089295) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84568923) q[3];
sx q[3];
rz(-1.6250623) q[3];
sx q[3];
rz(2.3783663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(-0.6883626) q[2];
rz(-2.5734731) q[3];
sx q[3];
rz(-1.0024242) q[3];
sx q[3];
rz(2.4585371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4532918) q[0];
sx q[0];
rz(-0.81128565) q[0];
sx q[0];
rz(1.9360833) q[0];
rz(-2.0506809) q[1];
sx q[1];
rz(-2.5542407) q[1];
sx q[1];
rz(1.6039414) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086370416) q[0];
sx q[0];
rz(-2.0049262) q[0];
sx q[0];
rz(0.40121292) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1411016) q[2];
sx q[2];
rz(-1.0028936) q[2];
sx q[2];
rz(2.9141324) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22828211) q[1];
sx q[1];
rz(-2.3582715) q[1];
sx q[1];
rz(1.7287152) q[1];
rz(0.64301771) q[3];
sx q[3];
rz(-2.4367691) q[3];
sx q[3];
rz(-0.091136668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9566112) q[2];
sx q[2];
rz(-2.0272171) q[2];
sx q[2];
rz(1.1448917) q[2];
rz(-1.6635118) q[3];
sx q[3];
rz(-0.76668113) q[3];
sx q[3];
rz(2.0293106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0803273) q[0];
sx q[0];
rz(-2.5054131) q[0];
sx q[0];
rz(-1.7058477) q[0];
rz(1.1371293) q[1];
sx q[1];
rz(-2.4500193) q[1];
sx q[1];
rz(0.93708509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3242151) q[0];
sx q[0];
rz(-2.827001) q[0];
sx q[0];
rz(-0.28633519) q[0];
rz(-pi) q[1];
rz(-1.996973) q[2];
sx q[2];
rz(-1.1010896) q[2];
sx q[2];
rz(-1.6324279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91793375) q[1];
sx q[1];
rz(-1.5964701) q[1];
sx q[1];
rz(2.3707546) q[1];
rz(-pi) q[2];
rz(2.5255975) q[3];
sx q[3];
rz(-1.8503891) q[3];
sx q[3];
rz(-0.071690138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3252141) q[2];
sx q[2];
rz(-2.4983695) q[2];
sx q[2];
rz(-0.31036672) q[2];
rz(-2.5988633) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(2.824596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9515789) q[0];
sx q[0];
rz(-1.0931451) q[0];
sx q[0];
rz(-2.6059294) q[0];
rz(2.426892) q[1];
sx q[1];
rz(-1.8938046) q[1];
sx q[1];
rz(-1.6751777) q[1];
rz(-2.428029) q[2];
sx q[2];
rz(-2.5893671) q[2];
sx q[2];
rz(1.6895369) q[2];
rz(2.4126903) q[3];
sx q[3];
rz(-1.6160252) q[3];
sx q[3];
rz(3.0795815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
