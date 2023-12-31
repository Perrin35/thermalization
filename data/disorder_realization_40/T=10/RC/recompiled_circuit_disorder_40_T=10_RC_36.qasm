OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(0.11178804) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(-1.4844126) q[1];
sx q[1];
rz(2.9878374) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68054799) q[0];
sx q[0];
rz(-1.8616315) q[0];
sx q[0];
rz(-2.6936147) q[0];
rz(-pi) q[1];
rz(-0.076924952) q[2];
sx q[2];
rz(-1.3300606) q[2];
sx q[2];
rz(-1.7488637) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1195214) q[1];
sx q[1];
rz(-1.6011366) q[1];
sx q[1];
rz(-1.8896709) q[1];
rz(2.2114803) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(-1.9213284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.443632) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(-1.704818) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-2.2170128) q[0];
rz(2.1444767) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.3234214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01100563) q[0];
sx q[0];
rz(-0.6479833) q[0];
sx q[0];
rz(2.381071) q[0];
rz(-pi) q[1];
rz(-3.0201868) q[2];
sx q[2];
rz(-2.7896023) q[2];
sx q[2];
rz(-0.44167232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18113187) q[1];
sx q[1];
rz(-0.84003969) q[1];
sx q[1];
rz(1.9366656) q[1];
rz(-pi) q[2];
rz(1.0074535) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(-0.37936488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20415846) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-0.75817529) q[2];
rz(2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(-1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4085098) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-0.68840233) q[0];
rz(-3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(-2.6115131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2358658) q[0];
sx q[0];
rz(-1.4639963) q[0];
sx q[0];
rz(-1.8624767) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0968139) q[2];
sx q[2];
rz(-1.194343) q[2];
sx q[2];
rz(-1.9386171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59280076) q[1];
sx q[1];
rz(-2.2081516) q[1];
sx q[1];
rz(2.7194276) q[1];
rz(-pi) q[2];
rz(-1.8102874) q[3];
sx q[3];
rz(-1.4756087) q[3];
sx q[3];
rz(-0.59323192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(-0.83550134) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(-1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-0.13664666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9204606) q[0];
sx q[0];
rz(-1.4401299) q[0];
sx q[0];
rz(-2.0544102) q[0];
rz(-pi) q[1];
rz(2.6776671) q[2];
sx q[2];
rz(-1.8564965) q[2];
sx q[2];
rz(0.49109101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5977051) q[1];
sx q[1];
rz(-1.6391616) q[1];
sx q[1];
rz(0.023936546) q[1];
rz(1.7131545) q[3];
sx q[3];
rz(-0.58066237) q[3];
sx q[3];
rz(2.2992087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(0.56048918) q[2];
rz(3.1292606) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43301582) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.4076153) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27027425) q[0];
sx q[0];
rz(-2.6042013) q[0];
sx q[0];
rz(-0.73211615) q[0];
rz(-0.3048004) q[2];
sx q[2];
rz(-2.1314203) q[2];
sx q[2];
rz(2.0914818) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69406063) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.8156169) q[1];
x q[2];
rz(1.0293343) q[3];
sx q[3];
rz(-0.41356219) q[3];
sx q[3];
rz(3.001861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(0.042479854) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(-0.4831627) q[0];
rz(2.0893611) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(0.56484708) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72950596) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(-1.8338404) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0658423) q[2];
sx q[2];
rz(-1.1510013) q[2];
sx q[2];
rz(-0.23725739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5833252) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(-1.7130909) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7259898) q[3];
sx q[3];
rz(-1.7994013) q[3];
sx q[3];
rz(2.5582563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5715282) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(1.3048874) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(-0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17334443) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(-0.54779732) q[0];
rz(2.3563747) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(0.26842591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82314202) q[0];
sx q[0];
rz(-2.8822691) q[0];
sx q[0];
rz(0.32487049) q[0];
rz(-pi) q[1];
rz(2.354291) q[2];
sx q[2];
rz(-0.95677081) q[2];
sx q[2];
rz(-1.3348483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0722326) q[1];
sx q[1];
rz(-1.4727122) q[1];
sx q[1];
rz(2.1588438) q[1];
rz(-pi) q[2];
rz(-1.0915756) q[3];
sx q[3];
rz(-2.5476533) q[3];
sx q[3];
rz(-1.493243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(2.3274373) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58586621) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-2.1222173) q[0];
rz(2.2881919) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(2.6928435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1837346) q[0];
sx q[0];
rz(-2.0439889) q[0];
sx q[0];
rz(-2.8790022) q[0];
rz(2.2737695) q[2];
sx q[2];
rz(-2.7542369) q[2];
sx q[2];
rz(2.5572436) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1467421) q[1];
sx q[1];
rz(-0.59616201) q[1];
sx q[1];
rz(-0.46390987) q[1];
rz(-pi) q[2];
rz(2.429871) q[3];
sx q[3];
rz(-0.91152836) q[3];
sx q[3];
rz(-0.014051138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2723508) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(2.8273919) q[2];
rz(2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(1.2605793) q[0];
rz(-2.4977327) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(-3.1226645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277055) q[0];
sx q[0];
rz(-0.2642309) q[0];
sx q[0];
rz(3.0731191) q[0];
rz(-pi) q[1];
rz(-1.47255) q[2];
sx q[2];
rz(-1.4155404) q[2];
sx q[2];
rz(-1.9948024) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70871204) q[1];
sx q[1];
rz(-1.6212665) q[1];
sx q[1];
rz(-2.9570079) q[1];
rz(-pi) q[2];
rz(-0.67846672) q[3];
sx q[3];
rz(-2.7760091) q[3];
sx q[3];
rz(1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5921322) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(-0.67031676) q[2];
rz(-2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-2.0589028) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95613854) q[0];
sx q[0];
rz(-2.1573665) q[0];
sx q[0];
rz(-1.6295208) q[0];
x q[1];
rz(0.033109025) q[2];
sx q[2];
rz(-2.3220255) q[2];
sx q[2];
rz(2.5493252) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5233366) q[1];
sx q[1];
rz(-0.54392951) q[1];
sx q[1];
rz(3.1052599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38090221) q[3];
sx q[3];
rz(-1.8758095) q[3];
sx q[3];
rz(1.1327133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.164914) q[2];
sx q[2];
rz(-0.51188525) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(2.0846562) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95505161) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(0.74116771) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(-1.0319866) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(1.742733) q[3];
sx q[3];
rz(-2.3098683) q[3];
sx q[3];
rz(1.5666425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
