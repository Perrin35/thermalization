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
rz(-1.8104115) q[0];
sx q[0];
rz(-1.6348732) q[0];
sx q[0];
rz(2.9659086) q[0];
rz(-2.076258) q[1];
sx q[1];
rz(-1.330436) q[1];
sx q[1];
rz(-2.4863844) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6106843) q[0];
sx q[0];
rz(-1.887868) q[0];
sx q[0];
rz(-1.0560587) q[0];
rz(-pi) q[1];
rz(-3.0426435) q[2];
sx q[2];
rz(-1.5155158) q[2];
sx q[2];
rz(1.7951622) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53319327) q[1];
sx q[1];
rz(-0.13020588) q[1];
sx q[1];
rz(-0.59334468) q[1];
x q[2];
rz(-0.96027042) q[3];
sx q[3];
rz(-1.8132079) q[3];
sx q[3];
rz(-1.9376955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4483999) q[2];
sx q[2];
rz(-0.11036631) q[2];
sx q[2];
rz(-2.8006862) q[2];
rz(-2.36813) q[3];
sx q[3];
rz(-2.1229459) q[3];
sx q[3];
rz(0.10507467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648436) q[0];
sx q[0];
rz(-0.28443795) q[0];
sx q[0];
rz(0.25962096) q[0];
rz(-2.3530841) q[1];
sx q[1];
rz(-0.85616833) q[1];
sx q[1];
rz(1.8964918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9246747) q[0];
sx q[0];
rz(-0.62546906) q[0];
sx q[0];
rz(-1.3979925) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8444031) q[2];
sx q[2];
rz(-2.8126394) q[2];
sx q[2];
rz(-2.7847814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4739405) q[1];
sx q[1];
rz(-0.86939916) q[1];
sx q[1];
rz(2.6768973) q[1];
x q[2];
rz(-0.50778265) q[3];
sx q[3];
rz(-1.4243817) q[3];
sx q[3];
rz(3.0203569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.239324) q[2];
sx q[2];
rz(-2.6332492) q[2];
sx q[2];
rz(0.90947214) q[2];
rz(-0.66696143) q[3];
sx q[3];
rz(-1.704155) q[3];
sx q[3];
rz(-0.10233574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7143836) q[0];
sx q[0];
rz(-0.41218555) q[0];
sx q[0];
rz(-1.3877731) q[0];
rz(-2.1448403) q[1];
sx q[1];
rz(-2.4501188) q[1];
sx q[1];
rz(-0.48666418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.023611) q[0];
sx q[0];
rz(-1.9754462) q[0];
sx q[0];
rz(0.16594529) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66101199) q[2];
sx q[2];
rz(-2.1901988) q[2];
sx q[2];
rz(0.41601478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4506766) q[1];
sx q[1];
rz(-1.6918536) q[1];
sx q[1];
rz(0.52074213) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1975945) q[3];
sx q[3];
rz(-1.6907915) q[3];
sx q[3];
rz(-1.7606408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3643058) q[2];
sx q[2];
rz(-1.2591668) q[2];
sx q[2];
rz(-2.6510009) q[2];
rz(-0.81654882) q[3];
sx q[3];
rz(-0.74211636) q[3];
sx q[3];
rz(-2.1013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.33646026) q[0];
sx q[0];
rz(-1.8959624) q[0];
sx q[0];
rz(-2.1601987) q[0];
rz(-0.40311748) q[1];
sx q[1];
rz(-2.6336481) q[1];
sx q[1];
rz(0.53781646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11975461) q[0];
sx q[0];
rz(-1.7203385) q[0];
sx q[0];
rz(-1.8200021) q[0];
x q[1];
rz(2.2914406) q[2];
sx q[2];
rz(-1.4701029) q[2];
sx q[2];
rz(-0.062391524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.855558) q[1];
sx q[1];
rz(-1.0893309) q[1];
sx q[1];
rz(0.43586327) q[1];
rz(-2.245491) q[3];
sx q[3];
rz(-1.3235914) q[3];
sx q[3];
rz(0.02442115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2526523) q[2];
sx q[2];
rz(-2.6344968) q[2];
sx q[2];
rz(2.0895152) q[2];
rz(-2.5893411) q[3];
sx q[3];
rz(-1.735732) q[3];
sx q[3];
rz(1.9210057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1110765) q[0];
sx q[0];
rz(-1.3330326) q[0];
sx q[0];
rz(-0.037574969) q[0];
rz(-2.3887718) q[1];
sx q[1];
rz(-1.5767187) q[1];
sx q[1];
rz(-3.0573696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0287474) q[0];
sx q[0];
rz(-2.0984834) q[0];
sx q[0];
rz(-0.28808388) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8657612) q[2];
sx q[2];
rz(-1.2216481) q[2];
sx q[2];
rz(-2.5087506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7561235) q[1];
sx q[1];
rz(-1.9944571) q[1];
sx q[1];
rz(3.0961354) q[1];
rz(1.7315167) q[3];
sx q[3];
rz(-1.6165716) q[3];
sx q[3];
rz(1.20884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5273744) q[2];
sx q[2];
rz(-0.97443333) q[2];
sx q[2];
rz(0.43240377) q[2];
rz(1.4934941) q[3];
sx q[3];
rz(-1.7688388) q[3];
sx q[3];
rz(2.4026332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768141) q[0];
sx q[0];
rz(-1.5274763) q[0];
sx q[0];
rz(-2.3798808) q[0];
rz(2.1181882) q[1];
sx q[1];
rz(-2.6506212) q[1];
sx q[1];
rz(2.8010119) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42804952) q[0];
sx q[0];
rz(-2.3556956) q[0];
sx q[0];
rz(-3.0160496) q[0];
rz(2.5804065) q[2];
sx q[2];
rz(-1.4380939) q[2];
sx q[2];
rz(-0.36280879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2557166) q[1];
sx q[1];
rz(-1.2310561) q[1];
sx q[1];
rz(-2.8541616) q[1];
rz(-1.7238925) q[3];
sx q[3];
rz(-0.82624895) q[3];
sx q[3];
rz(-1.561867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1731825) q[2];
sx q[2];
rz(-2.0081655) q[2];
sx q[2];
rz(-0.75562149) q[2];
rz(-0.38718811) q[3];
sx q[3];
rz(-1.3914934) q[3];
sx q[3];
rz(2.3473327) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8887535) q[0];
sx q[0];
rz(-3.0536953) q[0];
sx q[0];
rz(-1.2095691) q[0];
rz(-1.6190489) q[1];
sx q[1];
rz(-2.3698898) q[1];
sx q[1];
rz(1.3873772) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4693869) q[0];
sx q[0];
rz(-1.2325365) q[0];
sx q[0];
rz(2.8126008) q[0];
rz(-pi) q[1];
rz(0.31529324) q[2];
sx q[2];
rz(-2.1209426) q[2];
sx q[2];
rz(-1.2619051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3183115) q[1];
sx q[1];
rz(-2.2522587) q[1];
sx q[1];
rz(1.0266515) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8903743) q[3];
sx q[3];
rz(-2.0154023) q[3];
sx q[3];
rz(0.79892677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9428955) q[2];
sx q[2];
rz(-0.2936475) q[2];
sx q[2];
rz(2.7315308) q[2];
rz(-2.6041218) q[3];
sx q[3];
rz(-2.0987857) q[3];
sx q[3];
rz(1.0679831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033567) q[0];
sx q[0];
rz(-1.9060059) q[0];
sx q[0];
rz(2.1349452) q[0];
rz(-1.135723) q[1];
sx q[1];
rz(-2.6508811) q[1];
sx q[1];
rz(-2.8699285) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1170809) q[0];
sx q[0];
rz(-0.70444626) q[0];
sx q[0];
rz(-2.4570877) q[0];
rz(-pi) q[1];
rz(0.92046787) q[2];
sx q[2];
rz(-1.0050541) q[2];
sx q[2];
rz(2.5793902) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0783196) q[1];
sx q[1];
rz(-2.1542179) q[1];
sx q[1];
rz(2.2607445) q[1];
rz(-pi) q[2];
rz(-2.397698) q[3];
sx q[3];
rz(-1.7366323) q[3];
sx q[3];
rz(-1.3506571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9378308) q[2];
sx q[2];
rz(-2.2088642) q[2];
sx q[2];
rz(-2.0358098) q[2];
rz(1.3984937) q[3];
sx q[3];
rz(-2.1574557) q[3];
sx q[3];
rz(-0.85247803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49649134) q[0];
sx q[0];
rz(-0.95106769) q[0];
sx q[0];
rz(-2.7880461) q[0];
rz(-1.722466) q[1];
sx q[1];
rz(-1.5357693) q[1];
sx q[1];
rz(0.062573418) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1997027) q[0];
sx q[0];
rz(-2.1158005) q[0];
sx q[0];
rz(-2.1502994) q[0];
rz(-pi) q[1];
x q[1];
rz(2.563971) q[2];
sx q[2];
rz(-0.84813877) q[2];
sx q[2];
rz(-1.8236782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8282916) q[1];
sx q[1];
rz(-0.57692674) q[1];
sx q[1];
rz(-0.58842701) q[1];
rz(-1.6639641) q[3];
sx q[3];
rz(-0.9103295) q[3];
sx q[3];
rz(-0.70445337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.320437) q[2];
sx q[2];
rz(-1.374036) q[2];
sx q[2];
rz(-2.721411) q[2];
rz(0.34475103) q[3];
sx q[3];
rz(-0.77778608) q[3];
sx q[3];
rz(2.2043118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6688113) q[0];
sx q[0];
rz(-2.9720699) q[0];
sx q[0];
rz(-1.8602759) q[0];
rz(-1.5877113) q[1];
sx q[1];
rz(-2.3322767) q[1];
sx q[1];
rz(-2.4344427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2408694) q[0];
sx q[0];
rz(-1.5828195) q[0];
sx q[0];
rz(-2.2012086) q[0];
rz(-pi) q[1];
rz(2.1538582) q[2];
sx q[2];
rz(-1.1942004) q[2];
sx q[2];
rz(-1.131112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36701074) q[1];
sx q[1];
rz(-0.38644192) q[1];
sx q[1];
rz(-2.5143753) q[1];
rz(-pi) q[2];
rz(0.78451802) q[3];
sx q[3];
rz(-0.79567474) q[3];
sx q[3];
rz(-0.28233847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3215434) q[2];
sx q[2];
rz(-2.0115435) q[2];
sx q[2];
rz(-2.6450487) q[2];
rz(-1.8494362) q[3];
sx q[3];
rz(-2.8437331) q[3];
sx q[3];
rz(0.37798247) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15688607) q[0];
sx q[0];
rz(-0.26837415) q[0];
sx q[0];
rz(-1.0347086) q[0];
rz(2.5490419) q[1];
sx q[1];
rz(-0.68065803) q[1];
sx q[1];
rz(-1.5417644) q[1];
rz(0.040475673) q[2];
sx q[2];
rz(-1.802609) q[2];
sx q[2];
rz(2.6229057) q[2];
rz(-0.32101202) q[3];
sx q[3];
rz(-0.62855084) q[3];
sx q[3];
rz(1.4756065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
