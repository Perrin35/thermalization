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
rz(0.95712823) q[0];
sx q[0];
rz(4.794802) q[0];
sx q[0];
rz(10.573536) q[0];
rz(0.59612885) q[1];
sx q[1];
rz(-0.40607536) q[1];
sx q[1];
rz(-1.4546855) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6172778) q[0];
sx q[0];
rz(-1.2067267) q[0];
sx q[0];
rz(0.84366701) q[0];
rz(1.6802915) q[2];
sx q[2];
rz(-2.3286162) q[2];
sx q[2];
rz(2.2535498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0815989) q[1];
sx q[1];
rz(-2.3761056) q[1];
sx q[1];
rz(2.9503534) q[1];
x q[2];
rz(-1.4740397) q[3];
sx q[3];
rz(-2.2265937) q[3];
sx q[3];
rz(-0.024160926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.181432) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(1.9350447) q[2];
rz(-1.1613965) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777765) q[0];
sx q[0];
rz(-1.125536) q[0];
sx q[0];
rz(2.1695082) q[0];
rz(-2.8731335) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(-0.45480248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7642915) q[0];
sx q[0];
rz(-1.8669085) q[0];
sx q[0];
rz(0.13237093) q[0];
x q[1];
rz(2.9846086) q[2];
sx q[2];
rz(-1.9548226) q[2];
sx q[2];
rz(2.2228754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.763096) q[1];
sx q[1];
rz(-0.46751577) q[1];
sx q[1];
rz(-3.0276831) q[1];
rz(0.76105705) q[3];
sx q[3];
rz(-0.71963862) q[3];
sx q[3];
rz(-2.3005405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0316169) q[2];
sx q[2];
rz(-0.72828186) q[2];
sx q[2];
rz(-0.99962437) q[2];
rz(0.085974606) q[3];
sx q[3];
rz(-1.5857668) q[3];
sx q[3];
rz(0.011367817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407532) q[0];
sx q[0];
rz(-1.4680306) q[0];
sx q[0];
rz(0.22415796) q[0];
rz(1.5466746) q[1];
sx q[1];
rz(-1.5432576) q[1];
sx q[1];
rz(-1.3067783) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9973361) q[0];
sx q[0];
rz(-1.4486509) q[0];
sx q[0];
rz(2.5886378) q[0];
rz(-pi) q[1];
rz(-2.3036912) q[2];
sx q[2];
rz(-2.1616677) q[2];
sx q[2];
rz(1.6333579) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3114709) q[1];
sx q[1];
rz(-1.1563627) q[1];
sx q[1];
rz(0.45763514) q[1];
x q[2];
rz(0.31587028) q[3];
sx q[3];
rz(-0.69413737) q[3];
sx q[3];
rz(-2.4296076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.811502) q[2];
sx q[2];
rz(-1.5218488) q[2];
sx q[2];
rz(-2.1211993) q[2];
rz(-1.329782) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(-1.5248732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40600768) q[0];
sx q[0];
rz(-2.0841053) q[0];
sx q[0];
rz(-1.1757346) q[0];
rz(2.8658087) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(0.63794678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9383134) q[0];
sx q[0];
rz(-1.5601451) q[0];
sx q[0];
rz(-0.032708688) q[0];
rz(-2.5204896) q[2];
sx q[2];
rz(-1.3912462) q[2];
sx q[2];
rz(-0.68815069) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6762661) q[1];
sx q[1];
rz(-0.96083896) q[1];
sx q[1];
rz(1.9708939) q[1];
x q[2];
rz(-0.31973394) q[3];
sx q[3];
rz(-1.481062) q[3];
sx q[3];
rz(-3.0872726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63567579) q[2];
sx q[2];
rz(-1.5636779) q[2];
sx q[2];
rz(1.4471311) q[2];
rz(-2.9315089) q[3];
sx q[3];
rz(-0.45322067) q[3];
sx q[3];
rz(0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2337445) q[0];
sx q[0];
rz(-0.15234983) q[0];
sx q[0];
rz(1.0688758) q[0];
rz(-1.2999889) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(-1.2058421) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.962709) q[0];
sx q[0];
rz(-2.7307352) q[0];
sx q[0];
rz(-1.7715447) q[0];
rz(3.1320747) q[2];
sx q[2];
rz(-2.2235907) q[2];
sx q[2];
rz(1.4981086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0413549) q[1];
sx q[1];
rz(-0.85548988) q[1];
sx q[1];
rz(-3.0042786) q[1];
rz(-2.6603094) q[3];
sx q[3];
rz(-1.7168772) q[3];
sx q[3];
rz(1.8446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.892889) q[2];
sx q[2];
rz(-0.37962571) q[2];
sx q[2];
rz(-3.0730263) q[2];
rz(0.93962234) q[3];
sx q[3];
rz(-1.7328123) q[3];
sx q[3];
rz(1.9127964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3786316) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(2.6201541) q[0];
rz(1.6261082) q[1];
sx q[1];
rz(-1.7518967) q[1];
sx q[1];
rz(0.62561402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8396225) q[0];
sx q[0];
rz(-1.6018943) q[0];
sx q[0];
rz(-2.1108225) q[0];
rz(-pi) q[1];
rz(-2.7623903) q[2];
sx q[2];
rz(-1.450945) q[2];
sx q[2];
rz(2.3764936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.292899) q[1];
sx q[1];
rz(-1.1176511) q[1];
sx q[1];
rz(0.21993665) q[1];
x q[2];
rz(3.0066177) q[3];
sx q[3];
rz(-0.42536456) q[3];
sx q[3];
rz(-2.7736026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0079415) q[2];
sx q[2];
rz(-2.5456754) q[2];
sx q[2];
rz(0.1055183) q[2];
rz(2.1775235) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(-0.9849557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069020011) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(1.7286812) q[0];
rz(2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(0.55317318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1346388) q[0];
sx q[0];
rz(-2.3592822) q[0];
sx q[0];
rz(-0.99690848) q[0];
rz(-2.5752284) q[2];
sx q[2];
rz(-2.9743122) q[2];
sx q[2];
rz(-2.9110514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87482086) q[1];
sx q[1];
rz(-1.8541341) q[1];
sx q[1];
rz(-1.5288407) q[1];
x q[2];
rz(-0.91316789) q[3];
sx q[3];
rz(-0.32507691) q[3];
sx q[3];
rz(-0.43077786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46099123) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(-2.7294532) q[2];
rz(0.66017094) q[3];
sx q[3];
rz(-2.6001402) q[3];
sx q[3];
rz(-0.49303833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2039345) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(0.20189051) q[0];
rz(-1.0578602) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(-0.26022628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79151151) q[0];
sx q[0];
rz(-2.4790194) q[0];
sx q[0];
rz(2.4988079) q[0];
rz(-pi) q[1];
rz(0.8809907) q[2];
sx q[2];
rz(-0.63369232) q[2];
sx q[2];
rz(-1.6735759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34781814) q[1];
sx q[1];
rz(-1.8500237) q[1];
sx q[1];
rz(-0.3600959) q[1];
rz(-2.9256938) q[3];
sx q[3];
rz(-0.67189081) q[3];
sx q[3];
rz(-3.1300621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1883833) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(-0.58237135) q[2];
rz(-2.6992056) q[3];
sx q[3];
rz(-1.235639) q[3];
sx q[3];
rz(0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9198832) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(-1.7145994) q[0];
rz(-1.7859979) q[1];
sx q[1];
rz(-1.6537063) q[1];
sx q[1];
rz(-1.4367163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5858108) q[0];
sx q[0];
rz(-1.2511051) q[0];
sx q[0];
rz(3.1204819) q[0];
rz(0.53578761) q[2];
sx q[2];
rz(-2.5965207) q[2];
sx q[2];
rz(1.4591914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9475262) q[1];
sx q[1];
rz(-1.8119651) q[1];
sx q[1];
rz(-0.014726676) q[1];
x q[2];
rz(-2.0974227) q[3];
sx q[3];
rz(-1.4344525) q[3];
sx q[3];
rz(-2.6201893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9068678) q[2];
sx q[2];
rz(-1.8610672) q[2];
sx q[2];
rz(1.5391763) q[2];
rz(2.5336044) q[3];
sx q[3];
rz(-1.4092849) q[3];
sx q[3];
rz(-0.68979818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36668396) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(-2.3181424) q[0];
rz(-2.2460294) q[1];
sx q[1];
rz(-1.7915553) q[1];
sx q[1];
rz(-0.38111883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6247647) q[0];
sx q[0];
rz(-0.73113686) q[0];
sx q[0];
rz(0.69964377) q[0];
x q[1];
rz(-1.4032768) q[2];
sx q[2];
rz(-1.1741271) q[2];
sx q[2];
rz(-1.828891) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1824519) q[1];
sx q[1];
rz(-2.509626) q[1];
sx q[1];
rz(2.0672634) q[1];
x q[2];
rz(-0.70638871) q[3];
sx q[3];
rz(-0.76157842) q[3];
sx q[3];
rz(1.0672081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8054008) q[2];
sx q[2];
rz(-0.41173428) q[2];
sx q[2];
rz(2.1545048) q[2];
rz(-1.5385212) q[3];
sx q[3];
rz(-0.58731949) q[3];
sx q[3];
rz(0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22454746) q[0];
sx q[0];
rz(-1.2397091) q[0];
sx q[0];
rz(1.5201257) q[0];
rz(-0.70980258) q[1];
sx q[1];
rz(-0.55955049) q[1];
sx q[1];
rz(1.4581663) q[1];
rz(-1.456121) q[2];
sx q[2];
rz(-1.2085452) q[2];
sx q[2];
rz(2.1160407) q[2];
rz(-2.4003528) q[3];
sx q[3];
rz(-1.1164573) q[3];
sx q[3];
rz(-2.1637049) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
