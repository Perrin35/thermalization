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
rz(-0.91322652) q[0];
sx q[0];
rz(-0.93728596) q[0];
sx q[0];
rz(2.8159141) q[0];
rz(3.3477793) q[1];
sx q[1];
rz(3.4833796) q[1];
sx q[1];
rz(5.7835328) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9099115) q[0];
sx q[0];
rz(-1.2308111) q[0];
sx q[0];
rz(2.489734) q[0];
rz(-pi) q[1];
rz(-3.1392873) q[2];
sx q[2];
rz(-2.1426149) q[2];
sx q[2];
rz(-0.44849685) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.4573483) q[1];
sx q[1];
rz(-2.6169371) q[1];
sx q[1];
rz(-2.2150127) q[1];
rz(-pi) q[2];
rz(1.1127171) q[3];
sx q[3];
rz(-1.6290974) q[3];
sx q[3];
rz(1.3227303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8832522) q[2];
sx q[2];
rz(-0.84104937) q[2];
sx q[2];
rz(-2.945245) q[2];
rz(1.9718862) q[3];
sx q[3];
rz(-1.615808) q[3];
sx q[3];
rz(-1.4811146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.779988) q[0];
sx q[0];
rz(-1.9444281) q[0];
sx q[0];
rz(-2.8874604) q[0];
rz(-0.77468553) q[1];
sx q[1];
rz(-1.7278371) q[1];
sx q[1];
rz(2.4515801) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440702) q[0];
sx q[0];
rz(-1.43514) q[0];
sx q[0];
rz(-2.7619477) q[0];
rz(-pi) q[1];
rz(0.60890095) q[2];
sx q[2];
rz(-1.5594859) q[2];
sx q[2];
rz(-3.0593556) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0010052) q[1];
sx q[1];
rz(-2.6125357) q[1];
sx q[1];
rz(2.463813) q[1];
rz(2.0919741) q[3];
sx q[3];
rz(-1.3173977) q[3];
sx q[3];
rz(-2.1508642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0320354) q[2];
sx q[2];
rz(-3.1041225) q[2];
sx q[2];
rz(-2.1407342) q[2];
rz(-1.9199269) q[3];
sx q[3];
rz(-1.6358717) q[3];
sx q[3];
rz(-3.1083623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0669393) q[0];
sx q[0];
rz(-0.017843857) q[0];
sx q[0];
rz(-2.6113206) q[0];
rz(-0.52337921) q[1];
sx q[1];
rz(-1.7019848) q[1];
sx q[1];
rz(-1.9022) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059748404) q[0];
sx q[0];
rz(-0.91343266) q[0];
sx q[0];
rz(1.2685724) q[0];
x q[1];
rz(1.7943514) q[2];
sx q[2];
rz(-1.8164779) q[2];
sx q[2];
rz(-1.1576736) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5927913) q[1];
sx q[1];
rz(-0.32262412) q[1];
sx q[1];
rz(1.1760787) q[1];
rz(1.1887458) q[3];
sx q[3];
rz(-1.9460367) q[3];
sx q[3];
rz(0.98179152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8221028) q[2];
sx q[2];
rz(-2.4097061) q[2];
sx q[2];
rz(-3.0243995) q[2];
rz(-0.19290899) q[3];
sx q[3];
rz(-1.984805) q[3];
sx q[3];
rz(2.9624654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8837226) q[0];
sx q[0];
rz(-1.9921046) q[0];
sx q[0];
rz(0.024118751) q[0];
rz(-2.3564677) q[1];
sx q[1];
rz(-1.5631915) q[1];
sx q[1];
rz(-1.9804573) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1021082) q[0];
sx q[0];
rz(-2.2359936) q[0];
sx q[0];
rz(-2.7519144) q[0];
rz(-pi) q[1];
rz(-0.1382702) q[2];
sx q[2];
rz(-1.5492596) q[2];
sx q[2];
rz(-2.1295631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9930193) q[1];
sx q[1];
rz(-1.2482185) q[1];
sx q[1];
rz(0.60068074) q[1];
rz(2.9579453) q[3];
sx q[3];
rz(-2.1602294) q[3];
sx q[3];
rz(2.4494684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.44196096) q[2];
sx q[2];
rz(-0.72579759) q[2];
sx q[2];
rz(-0.3886784) q[2];
rz(-1.9937438) q[3];
sx q[3];
rz(-1.1286705) q[3];
sx q[3];
rz(1.8584937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1175784) q[0];
sx q[0];
rz(-2.3021181) q[0];
sx q[0];
rz(-1.041254) q[0];
rz(-2.4196692) q[1];
sx q[1];
rz(-1.8903774) q[1];
sx q[1];
rz(1.3232683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6206317) q[0];
sx q[0];
rz(-2.8537321) q[0];
sx q[0];
rz(-0.55572148) q[0];
rz(-pi) q[1];
rz(2.213654) q[2];
sx q[2];
rz(-1.7747806) q[2];
sx q[2];
rz(1.217733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3419575) q[1];
sx q[1];
rz(-0.87552908) q[1];
sx q[1];
rz(-1.2065843) q[1];
rz(-1.3106662) q[3];
sx q[3];
rz(-1.8046265) q[3];
sx q[3];
rz(-1.1375858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2920275) q[2];
sx q[2];
rz(-2.5490856) q[2];
sx q[2];
rz(-2.3805857) q[2];
rz(-1.5946397) q[3];
sx q[3];
rz(-1.2706815) q[3];
sx q[3];
rz(-2.7775619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4407235) q[0];
sx q[0];
rz(-2.4969164) q[0];
sx q[0];
rz(-2.0194637) q[0];
rz(2.2286277) q[1];
sx q[1];
rz(-2.2367621) q[1];
sx q[1];
rz(-0.95776552) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65993138) q[0];
sx q[0];
rz(-0.587981) q[0];
sx q[0];
rz(1.5065145) q[0];
rz(-1.945044) q[2];
sx q[2];
rz(-0.23434445) q[2];
sx q[2];
rz(-2.1017016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44555741) q[1];
sx q[1];
rz(-0.51338351) q[1];
sx q[1];
rz(-2.6460893) q[1];
rz(-pi) q[2];
rz(-1.692633) q[3];
sx q[3];
rz(-1.7434396) q[3];
sx q[3];
rz(-0.71820503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21445788) q[2];
sx q[2];
rz(-2.1467777) q[2];
sx q[2];
rz(2.4289995) q[2];
rz(-1.7566977) q[3];
sx q[3];
rz(-1.8637928) q[3];
sx q[3];
rz(-2.8539997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6433542) q[0];
sx q[0];
rz(-2.9297628) q[0];
sx q[0];
rz(-2.8969452) q[0];
rz(3.0554166) q[1];
sx q[1];
rz(-2.0241604) q[1];
sx q[1];
rz(0.92775956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13138765) q[0];
sx q[0];
rz(-2.3213815) q[0];
sx q[0];
rz(1.4641692) q[0];
rz(-pi) q[1];
rz(2.5439569) q[2];
sx q[2];
rz(-2.8888787) q[2];
sx q[2];
rz(-1.739335) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5416762) q[1];
sx q[1];
rz(-2.031026) q[1];
sx q[1];
rz(-0.61111889) q[1];
x q[2];
rz(2.4403839) q[3];
sx q[3];
rz(-1.0558914) q[3];
sx q[3];
rz(-1.9659404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9693532) q[2];
sx q[2];
rz(-0.66102207) q[2];
sx q[2];
rz(-1.4491402) q[2];
rz(-1.0714072) q[3];
sx q[3];
rz(-1.7457733) q[3];
sx q[3];
rz(1.9381819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0485359) q[0];
sx q[0];
rz(-3.0003248) q[0];
sx q[0];
rz(2.5988044) q[0];
rz(-0.12611783) q[1];
sx q[1];
rz(-1.3821955) q[1];
sx q[1];
rz(0.19475591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31349438) q[0];
sx q[0];
rz(-2.6107105) q[0];
sx q[0];
rz(-2.8110225) q[0];
rz(2.3210377) q[2];
sx q[2];
rz(-0.61167292) q[2];
sx q[2];
rz(0.95407971) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.69079) q[1];
sx q[1];
rz(-2.6581785) q[1];
sx q[1];
rz(0.40287896) q[1];
rz(-pi) q[2];
rz(0.1925999) q[3];
sx q[3];
rz(-2.3705774) q[3];
sx q[3];
rz(0.23480496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9572218) q[2];
sx q[2];
rz(-1.7879282) q[2];
sx q[2];
rz(-1.4380737) q[2];
rz(2.461869) q[3];
sx q[3];
rz(-0.34346911) q[3];
sx q[3];
rz(-0.22469416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1723802) q[0];
sx q[0];
rz(-0.42177105) q[0];
sx q[0];
rz(-1.0129741) q[0];
rz(-0.29533932) q[1];
sx q[1];
rz(-1.7177918) q[1];
sx q[1];
rz(1.0257592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00088770168) q[0];
sx q[0];
rz(-1.4460588) q[0];
sx q[0];
rz(-0.38356218) q[0];
rz(-pi) q[1];
rz(0.59495937) q[2];
sx q[2];
rz(-1.042871) q[2];
sx q[2];
rz(-2.4558112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.900454) q[1];
sx q[1];
rz(-1.608888) q[1];
sx q[1];
rz(2.898192) q[1];
rz(-pi) q[2];
rz(2.4186198) q[3];
sx q[3];
rz(-1.3275258) q[3];
sx q[3];
rz(-1.9169501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50461304) q[2];
sx q[2];
rz(-0.29256088) q[2];
sx q[2];
rz(-2.2620849) q[2];
rz(-0.28359908) q[3];
sx q[3];
rz(-1.2181543) q[3];
sx q[3];
rz(-1.846419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.857169) q[0];
sx q[0];
rz(-1.5623982) q[0];
sx q[0];
rz(-0.6012342) q[0];
rz(-0.55016905) q[1];
sx q[1];
rz(-2.1255122) q[1];
sx q[1];
rz(-0.30204958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4539393) q[0];
sx q[0];
rz(-2.2317356) q[0];
sx q[0];
rz(-2.3783408) q[0];
rz(2.1303802) q[2];
sx q[2];
rz(-2.4200984) q[2];
sx q[2];
rz(-1.341429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7833134) q[1];
sx q[1];
rz(-2.9960592) q[1];
sx q[1];
rz(-0.029749327) q[1];
rz(-pi) q[2];
rz(2.1393023) q[3];
sx q[3];
rz(-1.3703111) q[3];
sx q[3];
rz(-2.025163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.067513116) q[2];
sx q[2];
rz(-0.16177598) q[2];
sx q[2];
rz(1.3241241) q[2];
rz(0.4829123) q[3];
sx q[3];
rz(-0.78871471) q[3];
sx q[3];
rz(0.27831349) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0778238) q[0];
sx q[0];
rz(-1.1652729) q[0];
sx q[0];
rz(-1.2399974) q[0];
rz(-2.3700312) q[1];
sx q[1];
rz(-1.1411219) q[1];
sx q[1];
rz(0.065446767) q[1];
rz(-0.83648079) q[2];
sx q[2];
rz(-1.4118839) q[2];
sx q[2];
rz(-1.6250162) q[2];
rz(0.92486935) q[3];
sx q[3];
rz(-2.625941) q[3];
sx q[3];
rz(0.51437638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
