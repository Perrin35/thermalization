OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642759) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(-3.0741865) q[0];
rz(-1.2201266) q[2];
sx q[2];
rz(-1.7979243) q[2];
sx q[2];
rz(-0.35958689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19385399) q[1];
sx q[1];
rz(-2.1734936) q[1];
sx q[1];
rz(-2.1147453) q[1];
x q[2];
rz(-2.7636823) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(-0.4370673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.4651441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11008308) q[0];
sx q[0];
rz(-0.84658909) q[0];
sx q[0];
rz(-1.5674885) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0589774) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(-0.28085923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0959024) q[1];
sx q[1];
rz(-2.0369139) q[1];
sx q[1];
rz(-2.5065266) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.915669) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(-1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(-2.3410472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313989) q[0];
sx q[0];
rz(-2.7079765) q[0];
sx q[0];
rz(1.9171159) q[0];
rz(-pi) q[1];
rz(-2.2764552) q[2];
sx q[2];
rz(-1.0312928) q[2];
sx q[2];
rz(-1.7973961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2785981) q[1];
sx q[1];
rz(-1.1578214) q[1];
sx q[1];
rz(-1.5652656) q[1];
rz(-0.86977264) q[3];
sx q[3];
rz(-1.8043892) q[3];
sx q[3];
rz(-1.9529395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67733726) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(1.7819972) q[2];
rz(-0.96757403) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(-1.7165855) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(-2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(1.0850614) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61468716) q[0];
sx q[0];
rz(-0.58262107) q[0];
sx q[0];
rz(0.42838642) q[0];
rz(-pi) q[1];
rz(-2.1539139) q[2];
sx q[2];
rz(-2.7104125) q[2];
sx q[2];
rz(1.093986) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76772056) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(-0.14212455) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(0.51182169) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-0.98168215) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22589332) q[0];
sx q[0];
rz(-0.75076538) q[0];
sx q[0];
rz(0.70930003) q[0];
x q[1];
rz(-1.328674) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(2.9850933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2770734) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(-0.31788748) q[1];
rz(3.0675689) q[3];
sx q[3];
rz(-2.166966) q[3];
sx q[3];
rz(2.1285469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-2.9763124) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76089345) q[0];
sx q[0];
rz(-1.6000415) q[0];
sx q[0];
rz(0.0951035) q[0];
x q[1];
rz(-3.0566453) q[2];
sx q[2];
rz(-2.4521378) q[2];
sx q[2];
rz(1.2598318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11152553) q[1];
sx q[1];
rz(-1.0377874) q[1];
sx q[1];
rz(2.0433321) q[1];
rz(-pi) q[2];
rz(-2.0201488) q[3];
sx q[3];
rz(-1.7835788) q[3];
sx q[3];
rz(-1.7920115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(0.73648891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7613477) q[0];
sx q[0];
rz(-1.5537098) q[0];
sx q[0];
rz(-1.5218309) q[0];
rz(-pi) q[1];
rz(2.1963504) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(-1.2298825) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0349717) q[1];
sx q[1];
rz(-1.391341) q[1];
sx q[1];
rz(-2.3751395) q[1];
rz(0.69865366) q[3];
sx q[3];
rz(-0.62056382) q[3];
sx q[3];
rz(-0.26649775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.2954856) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078995973) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(2.2198026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066276) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(0.28618257) q[2];
sx q[2];
rz(-2.2879062) q[2];
sx q[2];
rz(-1.66695) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39107716) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(-0.68117546) q[1];
x q[2];
rz(-0.17955762) q[3];
sx q[3];
rz(-2.2710544) q[3];
sx q[3];
rz(-1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(-1.3575859) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(-0.67970651) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25526991) q[0];
sx q[0];
rz(-2.5810044) q[0];
sx q[0];
rz(-0.38829304) q[0];
rz(-pi) q[1];
x q[1];
rz(2.26608) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(-1.306844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19787994) q[1];
sx q[1];
rz(-2.3730952) q[1];
sx q[1];
rz(0.49853034) q[1];
rz(-pi) q[2];
rz(0.67919517) q[3];
sx q[3];
rz(-2.0742356) q[3];
sx q[3];
rz(-2.5356843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(2.2311907) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.048621) q[0];
sx q[0];
rz(-1.7341359) q[0];
sx q[0];
rz(-0.049157814) q[0];
x q[1];
rz(0.54729692) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(-1.0123569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6224222) q[1];
sx q[1];
rz(-1.4960438) q[1];
sx q[1];
rz(-1.012411) q[1];
rz(-pi) q[2];
rz(0.51104607) q[3];
sx q[3];
rz(-1.647445) q[3];
sx q[3];
rz(0.49728909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(1.8489807) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
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
