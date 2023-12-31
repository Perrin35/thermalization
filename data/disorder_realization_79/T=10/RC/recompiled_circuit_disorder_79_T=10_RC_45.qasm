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
rz(7.1516501) q[0];
sx q[0];
rz(9.2317543) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773168) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(3.0741865) q[0];
rz(0.97857742) q[2];
sx q[2];
rz(-2.7263612) q[2];
sx q[2];
rz(2.4821971) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19385399) q[1];
sx q[1];
rz(-0.96809909) q[1];
sx q[1];
rz(-2.1147453) q[1];
rz(-pi) q[2];
rz(1.0055553) q[3];
sx q[3];
rz(-1.894265) q[3];
sx q[3];
rz(2.2068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(-2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-2.0400955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9803479) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.4651441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315096) q[0];
sx q[0];
rz(-2.2950036) q[0];
sx q[0];
rz(1.5741041) q[0];
x q[1];
rz(-2.0589774) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(-2.8607334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.84491731) q[1];
sx q[1];
rz(-1.0122609) q[1];
sx q[1];
rz(1.0122453) q[1];
x q[2];
rz(1.915669) q[3];
sx q[3];
rz(-2.2861835) q[3];
sx q[3];
rz(1.8638924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(1.4573147) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96317545) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(0.80054545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0101937) q[0];
sx q[0];
rz(-2.7079765) q[0];
sx q[0];
rz(-1.2244768) q[0];
x q[1];
rz(2.3163047) q[2];
sx q[2];
rz(-0.85916677) q[2];
sx q[2];
rz(-2.3724144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4360106) q[1];
sx q[1];
rz(-1.5758621) q[1];
sx q[1];
rz(-2.7286121) q[1];
rz(-1.9242026) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(1.0850614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0264498) q[0];
sx q[0];
rz(-1.0466252) q[0];
sx q[0];
rz(-1.8379704) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9374574) q[2];
sx q[2];
rz(-1.3385834) q[2];
sx q[2];
rz(-1.0166849) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81740582) q[1];
sx q[1];
rz(-1.7122014) q[1];
sx q[1];
rz(-1.6721339) q[1];
x q[2];
rz(1.163108) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(-2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(2.1599105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9156993) q[0];
sx q[0];
rz(-2.3908273) q[0];
sx q[0];
rz(0.70930003) q[0];
rz(-1.328674) q[2];
sx q[2];
rz(-1.4171346) q[2];
sx q[2];
rz(0.15649934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2213858) q[1];
sx q[1];
rz(-2.033794) q[1];
sx q[1];
rz(-1.7358001) q[1];
rz(-pi) q[2];
rz(-0.97335191) q[3];
sx q[3];
rz(-1.6320328) q[3];
sx q[3];
rz(0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-2.7691675) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-2.9763124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344791) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(-1.5414184) q[0];
rz(-pi) q[1];
rz(0.084947305) q[2];
sx q[2];
rz(-0.68945486) q[2];
sx q[2];
rz(1.8817608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.936441) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(-0.58516296) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0201488) q[3];
sx q[3];
rz(-1.3580139) q[3];
sx q[3];
rz(-1.3495812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(2.929556) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38024494) q[0];
sx q[0];
rz(-1.5537098) q[0];
sx q[0];
rz(1.6197617) q[0];
x q[1];
rz(2.8924106) q[2];
sx q[2];
rz(-1.8998002) q[2];
sx q[2];
rz(2.5788384) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0349717) q[1];
sx q[1];
rz(-1.391341) q[1];
sx q[1];
rz(-0.76645318) q[1];
rz(-pi) q[2];
rz(-1.1398846) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0050469) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.078995973) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.2980365) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-2.2198026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1944273) q[0];
sx q[0];
rz(-2.2568586) q[0];
sx q[0];
rz(1.7454733) q[0];
rz(2.3085262) q[2];
sx q[2];
rz(-1.3563915) q[2];
sx q[2];
rz(0.094878541) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3075879) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(-0.25968857) q[1];
x q[2];
rz(-1.7796302) q[3];
sx q[3];
rz(-2.4224671) q[3];
sx q[3];
rz(-1.5987087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(2.4618861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4363791) q[0];
sx q[0];
rz(-1.0562911) q[0];
sx q[0];
rz(-1.3374469) q[0];
rz(2.9613413) q[2];
sx q[2];
rz(-2.2580574) q[2];
sx q[2];
rz(-2.7625411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.745984) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(-0.70396522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71984843) q[3];
sx q[3];
rz(-0.82092972) q[3];
sx q[3];
rz(-1.6380701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1264964) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(0.85062406) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991853) q[0];
sx q[0];
rz(-2.9710794) q[0];
sx q[0];
rz(1.8605581) q[0];
x q[1];
rz(-2.5942957) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(-1.0123569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1365876) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(0.088076061) q[1];
x q[2];
rz(-1.6586187) q[3];
sx q[3];
rz(-1.0613958) q[3];
sx q[3];
rz(-2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(-0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-2.7453616) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(-2.4049315) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(-1.5626004) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
