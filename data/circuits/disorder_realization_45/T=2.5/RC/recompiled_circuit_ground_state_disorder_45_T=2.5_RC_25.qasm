OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.328182) q[0];
sx q[0];
rz(5.2433104) q[0];
sx q[0];
rz(5.1007895) q[0];
rz(-2.8662968) q[1];
sx q[1];
rz(-1.0010012) q[1];
sx q[1];
rz(-1.6943975) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793854) q[0];
sx q[0];
rz(-1.1027063) q[0];
sx q[0];
rz(1.9759959) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2252407) q[2];
sx q[2];
rz(-2.3961903) q[2];
sx q[2];
rz(0.8236664) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.027982546) q[1];
sx q[1];
rz(-1.1259127) q[1];
sx q[1];
rz(1.0335733) q[1];
rz(-1.0436023) q[3];
sx q[3];
rz(-1.2498244) q[3];
sx q[3];
rz(2.2466602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64624661) q[2];
sx q[2];
rz(-2.4301811) q[2];
sx q[2];
rz(1.4900253) q[2];
rz(-1.1680158) q[3];
sx q[3];
rz(-1.7707337) q[3];
sx q[3];
rz(-1.542154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1646351) q[0];
sx q[0];
rz(-2.5049306) q[0];
sx q[0];
rz(-2.013999) q[0];
rz(-0.98059869) q[1];
sx q[1];
rz(-0.53739986) q[1];
sx q[1];
rz(1.5375536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011833103) q[0];
sx q[0];
rz(-0.3922222) q[0];
sx q[0];
rz(-1.5935728) q[0];
rz(-pi) q[1];
rz(1.6631378) q[2];
sx q[2];
rz(-1.5870816) q[2];
sx q[2];
rz(0.0013601842) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95629803) q[1];
sx q[1];
rz(-1.0377656) q[1];
sx q[1];
rz(-3.0009439) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62894024) q[3];
sx q[3];
rz(-2.6835614) q[3];
sx q[3];
rz(-2.0322134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.299122) q[2];
sx q[2];
rz(-0.8490347) q[2];
sx q[2];
rz(-2.9910198) q[2];
rz(1.8922837) q[3];
sx q[3];
rz(-2.1842897) q[3];
sx q[3];
rz(-1.5540436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097178) q[0];
sx q[0];
rz(-1.4282325) q[0];
sx q[0];
rz(-0.10794434) q[0];
rz(2.3297294) q[1];
sx q[1];
rz(-0.64473647) q[1];
sx q[1];
rz(0.20400861) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0013179) q[0];
sx q[0];
rz(-1.5840127) q[0];
sx q[0];
rz(-0.93854882) q[0];
rz(2.7008867) q[2];
sx q[2];
rz(-1.9153898) q[2];
sx q[2];
rz(-1.454892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9274689) q[1];
sx q[1];
rz(-1.0368616) q[1];
sx q[1];
rz(0.1410036) q[1];
x q[2];
rz(-1.9526945) q[3];
sx q[3];
rz(-2.3586267) q[3];
sx q[3];
rz(3.0375089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.418499) q[2];
sx q[2];
rz(-0.30947026) q[2];
sx q[2];
rz(-1.9401559) q[2];
rz(-0.57032436) q[3];
sx q[3];
rz(-1.0446965) q[3];
sx q[3];
rz(-2.9111011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-3.0975321) q[0];
sx q[0];
rz(-2.1756797) q[0];
sx q[0];
rz(0.048138054) q[0];
rz(0.33787456) q[1];
sx q[1];
rz(-2.477555) q[1];
sx q[1];
rz(-2.7545676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88833798) q[0];
sx q[0];
rz(-1.3562886) q[0];
sx q[0];
rz(-2.2127445) q[0];
rz(1.1311283) q[2];
sx q[2];
rz(-2.2626468) q[2];
sx q[2];
rz(-0.95702167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.59874) q[1];
sx q[1];
rz(-2.4495995) q[1];
sx q[1];
rz(2.180468) q[1];
rz(-pi) q[2];
rz(0.13366661) q[3];
sx q[3];
rz(-1.8456621) q[3];
sx q[3];
rz(0.2256338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43992299) q[2];
sx q[2];
rz(-2.5950626) q[2];
sx q[2];
rz(0.92865357) q[2];
rz(2.6141911) q[3];
sx q[3];
rz(-1.4832067) q[3];
sx q[3];
rz(1.5206913) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034390282) q[0];
sx q[0];
rz(-0.90029383) q[0];
sx q[0];
rz(0.1874371) q[0];
rz(1.772359) q[1];
sx q[1];
rz(-2.1010459) q[1];
sx q[1];
rz(1.7231411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013896) q[0];
sx q[0];
rz(-1.6121776) q[0];
sx q[0];
rz(0.94766728) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4985053) q[2];
sx q[2];
rz(-1.4355112) q[2];
sx q[2];
rz(0.39150086) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7781648) q[1];
sx q[1];
rz(-1.3194934) q[1];
sx q[1];
rz(2.9775724) q[1];
x q[2];
rz(0.62331919) q[3];
sx q[3];
rz(-1.5006362) q[3];
sx q[3];
rz(-0.96621418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9886542) q[2];
sx q[2];
rz(-2.8230437) q[2];
sx q[2];
rz(-2.5195276) q[2];
rz(-2.7234744) q[3];
sx q[3];
rz(-1.5112292) q[3];
sx q[3];
rz(2.6869584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76573265) q[0];
sx q[0];
rz(-1.3272165) q[0];
sx q[0];
rz(0.8263998) q[0];
rz(1.3580258) q[1];
sx q[1];
rz(-2.2260901) q[1];
sx q[1];
rz(-0.079039097) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86980692) q[0];
sx q[0];
rz(-3.0328864) q[0];
sx q[0];
rz(-2.4034017) q[0];
rz(-pi) q[1];
rz(-2.6174829) q[2];
sx q[2];
rz(-0.75499207) q[2];
sx q[2];
rz(-2.4885677) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0896645) q[1];
sx q[1];
rz(-0.66510495) q[1];
sx q[1];
rz(1.0699349) q[1];
x q[2];
rz(-1.7675381) q[3];
sx q[3];
rz(-1.1181346) q[3];
sx q[3];
rz(2.6552185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0283811) q[2];
sx q[2];
rz(-1.6190642) q[2];
sx q[2];
rz(-1.1979423) q[2];
rz(1.0269264) q[3];
sx q[3];
rz(-0.17470655) q[3];
sx q[3];
rz(3.0026156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2748264) q[0];
sx q[0];
rz(-0.7678031) q[0];
sx q[0];
rz(-0.36198947) q[0];
rz(-0.29701862) q[1];
sx q[1];
rz(-1.2535008) q[1];
sx q[1];
rz(2.6716935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74772111) q[0];
sx q[0];
rz(-2.031051) q[0];
sx q[0];
rz(2.3097765) q[0];
rz(1.5748575) q[2];
sx q[2];
rz(-0.60535678) q[2];
sx q[2];
rz(-0.53077215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5928157) q[1];
sx q[1];
rz(-1.5632742) q[1];
sx q[1];
rz(-1.8846017) q[1];
rz(-pi) q[2];
rz(3.0918202) q[3];
sx q[3];
rz(-0.94954606) q[3];
sx q[3];
rz(-1.6175818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.517259) q[2];
sx q[2];
rz(-0.91906157) q[2];
sx q[2];
rz(-0.84226912) q[2];
rz(-1.2256578) q[3];
sx q[3];
rz(-1.813348) q[3];
sx q[3];
rz(-1.8153502) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4351742) q[0];
sx q[0];
rz(-1.59732) q[0];
sx q[0];
rz(0.73961863) q[0];
rz(2.836272) q[1];
sx q[1];
rz(-1.0605597) q[1];
sx q[1];
rz(1.7230497) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.212767) q[0];
sx q[0];
rz(-0.97090844) q[0];
sx q[0];
rz(1.8119353) q[0];
rz(-2.9786359) q[2];
sx q[2];
rz(-0.96471404) q[2];
sx q[2];
rz(-1.9418756) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67603079) q[1];
sx q[1];
rz(-1.9622161) q[1];
sx q[1];
rz(-1.2483601) q[1];
rz(-3.0036054) q[3];
sx q[3];
rz(-1.7847698) q[3];
sx q[3];
rz(0.54420769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50489396) q[2];
sx q[2];
rz(-2.2765474) q[2];
sx q[2];
rz(-2.0828341) q[2];
rz(2.5508358) q[3];
sx q[3];
rz(-2.1754913) q[3];
sx q[3];
rz(3.0920715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1307986) q[0];
sx q[0];
rz(-2.8191691) q[0];
sx q[0];
rz(1.7266493) q[0];
rz(-2.2825799) q[1];
sx q[1];
rz(-2.450727) q[1];
sx q[1];
rz(2.1376999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52951661) q[0];
sx q[0];
rz(-2.8359543) q[0];
sx q[0];
rz(-2.5661657) q[0];
x q[1];
rz(0.42385139) q[2];
sx q[2];
rz(-0.7025439) q[2];
sx q[2];
rz(-2.4830816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2768663) q[1];
sx q[1];
rz(-1.4819615) q[1];
sx q[1];
rz(0.72901841) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3968302) q[3];
sx q[3];
rz(-1.6339374) q[3];
sx q[3];
rz(2.5936332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88177195) q[2];
sx q[2];
rz(-1.5544145) q[2];
sx q[2];
rz(0.47015321) q[2];
rz(1.638089) q[3];
sx q[3];
rz(-0.94109002) q[3];
sx q[3];
rz(0.21465429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44899061) q[0];
sx q[0];
rz(-2.7916491) q[0];
sx q[0];
rz(-1.0119447) q[0];
rz(0.28255209) q[1];
sx q[1];
rz(-2.2239182) q[1];
sx q[1];
rz(0.73175398) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7912023) q[0];
sx q[0];
rz(-1.6178432) q[0];
sx q[0];
rz(1.5806517) q[0];
rz(-pi) q[1];
rz(0.75586478) q[2];
sx q[2];
rz(-2.8939092) q[2];
sx q[2];
rz(-2.8510954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2866757) q[1];
sx q[1];
rz(-2.1647506) q[1];
sx q[1];
rz(2.0718365) q[1];
x q[2];
rz(1.326272) q[3];
sx q[3];
rz(-2.7774307) q[3];
sx q[3];
rz(-2.9579666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67873502) q[2];
sx q[2];
rz(-1.5588458) q[2];
sx q[2];
rz(-0.25657121) q[2];
rz(-2.3738142) q[3];
sx q[3];
rz(-1.0662181) q[3];
sx q[3];
rz(-2.0496875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11138042) q[0];
sx q[0];
rz(-2.31507) q[0];
sx q[0];
rz(2.4872381) q[0];
rz(2.4334473) q[1];
sx q[1];
rz(-1.000052) q[1];
sx q[1];
rz(2.4194385) q[1];
rz(0.80864334) q[2];
sx q[2];
rz(-1.8075347) q[2];
sx q[2];
rz(-1.6703477) q[2];
rz(1.9740029) q[3];
sx q[3];
rz(-1.5807865) q[3];
sx q[3];
rz(-2.577232) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
