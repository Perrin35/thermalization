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
rz(1.3266069) q[0];
sx q[0];
rz(3.366037) q[0];
sx q[0];
rz(11.233575) q[0];
rz(2.3092071) q[1];
sx q[1];
rz(-1.7008984) q[1];
sx q[1];
rz(2.031215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719037) q[0];
sx q[0];
rz(-2.9326322) q[0];
sx q[0];
rz(0.67061575) q[0];
rz(-pi) q[1];
rz(1.9259077) q[2];
sx q[2];
rz(-0.74374226) q[2];
sx q[2];
rz(-1.6031934) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.346525) q[1];
sx q[1];
rz(-1.5988886) q[1];
sx q[1];
rz(-3.0884519) q[1];
x q[2];
rz(1.2625804) q[3];
sx q[3];
rz(-2.1860414) q[3];
sx q[3];
rz(-2.8160994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9139468) q[2];
sx q[2];
rz(-3.1326742) q[2];
sx q[2];
rz(1.7019567) q[2];
rz(-1.7281744) q[3];
sx q[3];
rz(-3.1296802) q[3];
sx q[3];
rz(-0.2695151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.444376) q[0];
sx q[0];
rz(-1.5378636) q[0];
sx q[0];
rz(-0.61475301) q[0];
rz(2.6055824) q[1];
sx q[1];
rz(-3.1159846) q[1];
sx q[1];
rz(-0.33686179) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0885411) q[0];
sx q[0];
rz(-1.7053002) q[0];
sx q[0];
rz(-0.10945871) q[0];
rz(3.1104149) q[2];
sx q[2];
rz(-0.47652361) q[2];
sx q[2];
rz(2.2533803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.24241703) q[1];
sx q[1];
rz(-1.6149033) q[1];
sx q[1];
rz(-3.1257948) q[1];
x q[2];
rz(-2.9514489) q[3];
sx q[3];
rz(-1.2570253) q[3];
sx q[3];
rz(0.47167512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88379318) q[2];
sx q[2];
rz(-3.1287441) q[2];
sx q[2];
rz(2.4419355) q[2];
rz(-0.16957016) q[3];
sx q[3];
rz(-3.1278059) q[3];
sx q[3];
rz(-0.56210303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70497847) q[0];
sx q[0];
rz(-0.545937) q[0];
sx q[0];
rz(0.28526947) q[0];
rz(-0.26372313) q[1];
sx q[1];
rz(-0.00058760651) q[1];
sx q[1];
rz(-2.347351) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7274311) q[0];
sx q[0];
rz(-2.2530956) q[0];
sx q[0];
rz(-1.4201565) q[0];
rz(-pi) q[1];
rz(-1.2474485) q[2];
sx q[2];
rz(-1.4686606) q[2];
sx q[2];
rz(0.43854005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8869739) q[1];
sx q[1];
rz(-1.6052263) q[1];
sx q[1];
rz(-1.5643584) q[1];
rz(0.77190001) q[3];
sx q[3];
rz(-0.83166122) q[3];
sx q[3];
rz(-1.2496304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0245725) q[2];
sx q[2];
rz(-3.0954376) q[2];
sx q[2];
rz(-0.864492) q[2];
rz(-2.3965059) q[3];
sx q[3];
rz(-0.86747187) q[3];
sx q[3];
rz(0.043070506) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2660148) q[0];
sx q[0];
rz(-0.046253007) q[0];
sx q[0];
rz(0.86638802) q[0];
rz(-0.59151793) q[1];
sx q[1];
rz(-0.94647995) q[1];
sx q[1];
rz(1.0513069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49120228) q[0];
sx q[0];
rz(-1.5687546) q[0];
sx q[0];
rz(-1.5709086) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5685097) q[2];
sx q[2];
rz(-1.5703452) q[2];
sx q[2];
rz(1.566103) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1278197) q[1];
sx q[1];
rz(-0.88370815) q[1];
sx q[1];
rz(-0.52718157) q[1];
rz(2.4254708) q[3];
sx q[3];
rz(-2.4181626) q[3];
sx q[3];
rz(1.0264068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6136578) q[2];
sx q[2];
rz(-0.064585678) q[2];
sx q[2];
rz(2.2876372) q[2];
rz(2.4438786) q[3];
sx q[3];
rz(-1.3359952) q[3];
sx q[3];
rz(0.31049389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39603221) q[0];
sx q[0];
rz(-0.10265352) q[0];
sx q[0];
rz(2.7295617) q[0];
rz(1.283006) q[1];
sx q[1];
rz(-2.367815) q[1];
sx q[1];
rz(-0.7512908) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4921732) q[0];
sx q[0];
rz(-2.87531) q[0];
sx q[0];
rz(-2.6800781) q[0];
x q[1];
rz(2.2504248) q[2];
sx q[2];
rz(-1.749265) q[2];
sx q[2];
rz(2.8986069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48482286) q[1];
sx q[1];
rz(-0.98451603) q[1];
sx q[1];
rz(0.9723248) q[1];
rz(-pi) q[2];
rz(-2.2734518) q[3];
sx q[3];
rz(-2.9538547) q[3];
sx q[3];
rz(2.7526698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6777307) q[2];
sx q[2];
rz(-3.1158713) q[2];
sx q[2];
rz(0.99979293) q[2];
rz(-0.60100466) q[3];
sx q[3];
rz(-3.0809564) q[3];
sx q[3];
rz(0.84990466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26337013) q[0];
sx q[0];
rz(-3.0510986) q[0];
sx q[0];
rz(-0.40060842) q[0];
rz(1.7349617) q[1];
sx q[1];
rz(-0.35146439) q[1];
sx q[1];
rz(1.7469143) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6215537) q[0];
sx q[0];
rz(-1.5661608) q[0];
sx q[0];
rz(1.5511447) q[0];
rz(-pi) q[1];
rz(1.2591061) q[2];
sx q[2];
rz(-2.3505728) q[2];
sx q[2];
rz(-2.493801) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7813599) q[1];
sx q[1];
rz(-0.1905687) q[1];
sx q[1];
rz(-3.1182454) q[1];
rz(-0.73598273) q[3];
sx q[3];
rz(-1.6998359) q[3];
sx q[3];
rz(-2.1202212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3957735) q[2];
sx q[2];
rz(-2.5534111) q[2];
sx q[2];
rz(-1.0751209) q[2];
rz(-2.6221258) q[3];
sx q[3];
rz(-0.1778917) q[3];
sx q[3];
rz(-1.5943257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290685) q[0];
sx q[0];
rz(-1.2146177) q[0];
sx q[0];
rz(1.6125096) q[0];
rz(0.56135881) q[1];
sx q[1];
rz(-0.012265597) q[1];
sx q[1];
rz(0.57048172) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5382696) q[0];
sx q[0];
rz(-0.13165671) q[0];
sx q[0];
rz(-1.1696474) q[0];
rz(-1.8884185) q[2];
sx q[2];
rz(-1.3537242) q[2];
sx q[2];
rz(-2.0906015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.884932) q[1];
sx q[1];
rz(-3.1260371) q[1];
sx q[1];
rz(-1.6473098) q[1];
x q[2];
rz(2.3164767) q[3];
sx q[3];
rz(-1.9843915) q[3];
sx q[3];
rz(0.90892413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.059171112) q[2];
sx q[2];
rz(-2.8736281) q[2];
sx q[2];
rz(-1.9783665) q[2];
rz(-1.6022812) q[3];
sx q[3];
rz(-0.10207615) q[3];
sx q[3];
rz(2.1479837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38692835) q[0];
sx q[0];
rz(-2.8480777) q[0];
sx q[0];
rz(2.5013404) q[0];
rz(-0.86296588) q[1];
sx q[1];
rz(-0.14185618) q[1];
sx q[1];
rz(2.364025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3024738) q[0];
sx q[0];
rz(-0.57838744) q[0];
sx q[0];
rz(0.79721398) q[0];
rz(-pi) q[1];
x q[1];
rz(1.338358) q[2];
sx q[2];
rz(-0.82045499) q[2];
sx q[2];
rz(-2.6636843) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9915087) q[1];
sx q[1];
rz(-1.5822295) q[1];
sx q[1];
rz(1.6437794) q[1];
rz(-pi) q[2];
rz(2.0361774) q[3];
sx q[3];
rz(-1.1148561) q[3];
sx q[3];
rz(-2.3400644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.61206478) q[2];
sx q[2];
rz(-2.7175588) q[2];
sx q[2];
rz(-0.5303793) q[2];
rz(0.027675962) q[3];
sx q[3];
rz(-3.0961302) q[3];
sx q[3];
rz(-1.8224705) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767122) q[0];
sx q[0];
rz(-2.9806529) q[0];
sx q[0];
rz(-2.2093534) q[0];
rz(-1.5981916) q[1];
sx q[1];
rz(-0.80359572) q[1];
sx q[1];
rz(2.7995321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8810976) q[0];
sx q[0];
rz(-1.6552345) q[0];
sx q[0];
rz(1.5859833) q[0];
rz(-pi) q[1];
x q[1];
rz(0.069829536) q[2];
sx q[2];
rz(-0.73903767) q[2];
sx q[2];
rz(0.42934092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6530214) q[1];
sx q[1];
rz(-2.584051) q[1];
sx q[1];
rz(-2.7776021) q[1];
rz(-pi) q[2];
rz(-1.6430278) q[3];
sx q[3];
rz(-2.9119303) q[3];
sx q[3];
rz(-3.0185452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.082569294) q[2];
sx q[2];
rz(-3.1408568) q[2];
sx q[2];
rz(-1.2081344) q[2];
rz(2.1900603) q[3];
sx q[3];
rz(-3.1335242) q[3];
sx q[3];
rz(2.4212196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35219881) q[0];
sx q[0];
rz(-2.3645526) q[0];
sx q[0];
rz(3.0598031) q[0];
rz(-2.8261322) q[1];
sx q[1];
rz(-0.052736484) q[1];
sx q[1];
rz(-1.9116521) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30222505) q[0];
sx q[0];
rz(-1.4124083) q[0];
sx q[0];
rz(-1.7296289) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17773262) q[2];
sx q[2];
rz(-1.7077674) q[2];
sx q[2];
rz(0.21868071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3620748) q[1];
sx q[1];
rz(-1.6884585) q[1];
sx q[1];
rz(-0.090434342) q[1];
x q[2];
rz(-2.1020402) q[3];
sx q[3];
rz(-1.4508411) q[3];
sx q[3];
rz(2.1841072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54084593) q[2];
sx q[2];
rz(-3.1223065) q[2];
sx q[2];
rz(-2.02796) q[2];
rz(-3.1019548) q[3];
sx q[3];
rz(-0.0097291917) q[3];
sx q[3];
rz(-2.5447194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780846) q[0];
sx q[0];
rz(-1.9263374) q[0];
sx q[0];
rz(1.3081464) q[0];
rz(2.5254163) q[1];
sx q[1];
rz(-2.0694852) q[1];
sx q[1];
rz(-2.9328666) q[1];
rz(-2.2439416) q[2];
sx q[2];
rz(-2.8406526) q[2];
sx q[2];
rz(0.1546897) q[2];
rz(1.3499089) q[3];
sx q[3];
rz(-1.532423) q[3];
sx q[3];
rz(1.5929089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
