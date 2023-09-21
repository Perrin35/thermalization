OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1238414) q[0];
sx q[0];
rz(-3.0525644) q[0];
sx q[0];
rz(0.71319367) q[0];
x q[1];
rz(0.97857742) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(0.65939553) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.19385399) q[1];
sx q[1];
rz(-2.1734936) q[1];
sx q[1];
rz(2.1147453) q[1];
rz(2.1300415) q[3];
sx q[3];
rz(-0.64239255) q[3];
sx q[3];
rz(1.1005644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1258939) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-0.57693276) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.4651441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0365021) q[0];
sx q[0];
rz(-0.72421342) q[0];
sx q[0];
rz(3.1378531) q[0];
x q[1];
rz(-1.0826153) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(2.8607334) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.022545594) q[1];
sx q[1];
rz(-2.3733768) q[1];
sx q[1];
rz(0.70336282) q[1];
rz(-pi) q[2];
rz(-0.37104718) q[3];
sx q[3];
rz(-0.78073946) q[3];
sx q[3];
rz(0.77663976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-0.1097651) q[2];
rz(0.62260735) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96317545) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-0.80054545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313989) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(-1.2244768) q[0];
rz(-pi) q[1];
rz(2.2764552) q[2];
sx q[2];
rz(-2.1102998) q[2];
sx q[2];
rz(1.3441966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86299455) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(1.5763271) q[1];
rz(-2.839746) q[3];
sx q[3];
rz(-0.89248025) q[3];
sx q[3];
rz(-0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.7819972) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490705) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(-0.46491369) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-2.0565313) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59180075) q[0];
sx q[0];
rz(-1.8013957) q[0];
sx q[0];
rz(-2.6016298) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98767878) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(2.0476066) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.81740582) q[1];
sx q[1];
rz(-1.7122014) q[1];
sx q[1];
rz(-1.4694587) q[1];
rz(0.61120175) q[3];
sx q[3];
rz(-1.2308321) q[3];
sx q[3];
rz(1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.00502914) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3571346) q[0];
sx q[0];
rz(-2.031209) q[0];
sx q[0];
rz(-2.5255192) q[0];
rz(-pi) q[1];
rz(2.9833897) q[2];
sx q[2];
rz(-1.8100097) q[2];
sx q[2];
rz(1.4520793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7179467) q[1];
sx q[1];
rz(-1.4232993) q[1];
sx q[1];
rz(0.46848483) q[1];
rz(1.4622299) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(1.1443646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-2.7691675) q[0];
rz(-1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(-2.9763124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.107347) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(-0.29883595) q[0];
rz(-pi) q[1];
rz(1.6406306) q[2];
sx q[2];
rz(-0.88431057) q[2];
sx q[2];
rz(-1.7718466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.936441) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(-2.5564297) q[1];
rz(-pi) q[2];
rz(-2.9061756) q[3];
sx q[3];
rz(-2.0093007) q[3];
sx q[3];
rz(2.8188843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(-1.6983263) q[2];
rz(0.36744395) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-0.34926397) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(2.4051037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6155646) q[0];
sx q[0];
rz(-3.0897339) q[0];
sx q[0];
rz(1.2349013) q[0];
x q[1];
rz(-2.8924106) q[2];
sx q[2];
rz(-1.8998002) q[2];
sx q[2];
rz(-2.5788384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2939824) q[1];
sx q[1];
rz(-0.81969205) q[1];
sx q[1];
rz(-1.3241029) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1398846) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(-1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(-2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.2954856) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(2.334306) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(-2.2198026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1944273) q[0];
sx q[0];
rz(-0.88473407) q[0];
sx q[0];
rz(1.7454733) q[0];
x q[1];
rz(1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.8956172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31729749) q[1];
sx q[1];
rz(-1.824914) q[1];
sx q[1];
rz(1.3586678) q[1];
rz(-pi) q[2];
rz(2.962035) q[3];
sx q[3];
rz(-2.2710544) q[3];
sx q[3];
rz(-1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(0.67970651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8863227) q[0];
sx q[0];
rz(-0.5605883) q[0];
sx q[0];
rz(-2.7532996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87551261) q[2];
sx q[2];
rz(-1.4317703) q[2];
sx q[2];
rz(1.306844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6913773) q[1];
sx q[1];
rz(-0.9141578) q[1];
sx q[1];
rz(-1.1378098) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9548095) q[3];
sx q[3];
rz(-0.98815742) q[3];
sx q[3];
rz(2.5481176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(-2.4662468) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.9649327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.048621) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(-0.049157814) q[0];
rz(-2.2628225) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(-2.9825485) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1706108) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(1.7112205) q[1];
rz(-2.6305466) q[3];
sx q[3];
rz(-1.647445) q[3];
sx q[3];
rz(0.49728909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-2.4049315) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(0.35076326) q[3];
sx q[3];
rz(-1.5630994) q[3];
sx q[3];
rz(0.85170436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
