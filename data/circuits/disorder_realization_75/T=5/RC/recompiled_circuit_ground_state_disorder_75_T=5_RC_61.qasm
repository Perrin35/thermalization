OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45159856) q[0];
sx q[0];
rz(-0.30878433) q[0];
sx q[0];
rz(-0.2398332) q[0];
rz(2.580515) q[1];
sx q[1];
rz(-0.74220389) q[1];
sx q[1];
rz(-1.6286558) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5133249) q[0];
sx q[0];
rz(-1.4428992) q[0];
sx q[0];
rz(1.9403752) q[0];
x q[1];
rz(-0.26145491) q[2];
sx q[2];
rz(-2.0294218) q[2];
sx q[2];
rz(-2.0482091) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1354243) q[1];
sx q[1];
rz(-1.0841771) q[1];
sx q[1];
rz(1.7371337) q[1];
rz(-1.4024847) q[3];
sx q[3];
rz(-1.5136079) q[3];
sx q[3];
rz(-2.7528499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98770398) q[2];
sx q[2];
rz(-0.99606267) q[2];
sx q[2];
rz(-1.055701) q[2];
rz(-0.40134564) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(-1.6373985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5098679) q[0];
sx q[0];
rz(-0.26917502) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(2.3954605) q[1];
sx q[1];
rz(-1.1765307) q[1];
sx q[1];
rz(3.0768118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21532741) q[0];
sx q[0];
rz(-0.61867061) q[0];
sx q[0];
rz(-2.4149405) q[0];
rz(-2.4617875) q[2];
sx q[2];
rz(-1.7035988) q[2];
sx q[2];
rz(2.6390569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.244811) q[1];
sx q[1];
rz(-2.9771388) q[1];
sx q[1];
rz(2.0816148) q[1];
rz(-pi) q[2];
rz(-1.359932) q[3];
sx q[3];
rz(-1.5393889) q[3];
sx q[3];
rz(1.5109911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88227162) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-3.1003013) q[2];
rz(-1.1193554) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-0.4275221) q[0];
sx q[0];
rz(0.69931716) q[0];
rz(2.518867) q[1];
sx q[1];
rz(-2.3284349) q[1];
sx q[1];
rz(-0.25845382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759443) q[0];
sx q[0];
rz(-1.5830402) q[0];
sx q[0];
rz(-0.6338288) q[0];
rz(-pi) q[1];
rz(-3.0270789) q[2];
sx q[2];
rz(-2.2199759) q[2];
sx q[2];
rz(-2.442292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.041641673) q[1];
sx q[1];
rz(-2.0196242) q[1];
sx q[1];
rz(-1.5484518) q[1];
rz(-pi) q[2];
rz(1.4914037) q[3];
sx q[3];
rz(-0.43118048) q[3];
sx q[3];
rz(3.0357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2059325) q[2];
sx q[2];
rz(-1.6814597) q[2];
sx q[2];
rz(2.4350731) q[2];
rz(2.5126854) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(-1.1227054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0141456) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(-0.97417796) q[0];
rz(1.4840508) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(1.0008224) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580622) q[0];
sx q[0];
rz(-0.62072004) q[0];
sx q[0];
rz(-1.3002943) q[0];
rz(-0.77032178) q[2];
sx q[2];
rz(-2.1674967) q[2];
sx q[2];
rz(1.6092827) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82881935) q[1];
sx q[1];
rz(-2.3760894) q[1];
sx q[1];
rz(-0.70154066) q[1];
rz(1.6203247) q[3];
sx q[3];
rz(-1.2504745) q[3];
sx q[3];
rz(1.3763826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4546844) q[2];
sx q[2];
rz(-1.800622) q[2];
sx q[2];
rz(-0.81721133) q[2];
rz(2.5456083) q[3];
sx q[3];
rz(-1.7242566) q[3];
sx q[3];
rz(1.7433085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3572094) q[0];
sx q[0];
rz(-1.7359808) q[0];
sx q[0];
rz(-2.0696409) q[0];
rz(2.0042073) q[1];
sx q[1];
rz(-1.6268566) q[1];
sx q[1];
rz(0.17328182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.539285) q[0];
sx q[0];
rz(-1.1207016) q[0];
sx q[0];
rz(-2.0099239) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91867515) q[2];
sx q[2];
rz(-0.5917509) q[2];
sx q[2];
rz(2.7482035) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4049591) q[1];
sx q[1];
rz(-1.6425321) q[1];
sx q[1];
rz(0.81291764) q[1];
rz(2.625218) q[3];
sx q[3];
rz(-2.9547915) q[3];
sx q[3];
rz(-1.1774225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6413573) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(0.87453169) q[2];
rz(2.4747961) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(-0.82505208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2905529) q[0];
sx q[0];
rz(-2.9930826) q[0];
sx q[0];
rz(-0.17459757) q[0];
rz(2.5945276) q[1];
sx q[1];
rz(-2.6262296) q[1];
sx q[1];
rz(1.863716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1741076) q[0];
sx q[0];
rz(-0.293403) q[0];
sx q[0];
rz(1.3551329) q[0];
rz(-pi) q[1];
rz(0.49680423) q[2];
sx q[2];
rz(-2.3251109) q[2];
sx q[2];
rz(-1.4085438) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6656832) q[1];
sx q[1];
rz(-1.2194677) q[1];
sx q[1];
rz(-2.372588) q[1];
rz(-0.90323351) q[3];
sx q[3];
rz(-2.0732911) q[3];
sx q[3];
rz(-0.63695217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3829019) q[2];
sx q[2];
rz(-2.8768657) q[2];
sx q[2];
rz(2.7721789) q[2];
rz(0.82768011) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5675885) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(-0.23707238) q[0];
rz(1.7489307) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(2.1377835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35492453) q[0];
sx q[0];
rz(-1.1369497) q[0];
sx q[0];
rz(1.7545486) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0324442) q[2];
sx q[2];
rz(-0.91141111) q[2];
sx q[2];
rz(-0.75379363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79738338) q[1];
sx q[1];
rz(-0.64192574) q[1];
sx q[1];
rz(-2.5688237) q[1];
rz(-2.5828894) q[3];
sx q[3];
rz(-1.1538343) q[3];
sx q[3];
rz(-2.1897763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.5281333) q[2];
sx q[2];
rz(-1.5626835) q[2];
sx q[2];
rz(2.6252739) q[2];
rz(-0.83827072) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(0.18394884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(2.2273492) q[0];
rz(2.5573348) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(0.90726888) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4582918) q[0];
sx q[0];
rz(-2.2193546) q[0];
sx q[0];
rz(0.022819937) q[0];
x q[1];
rz(1.3970988) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(-0.71719682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71675352) q[1];
sx q[1];
rz(-0.83152926) q[1];
sx q[1];
rz(-0.6296954) q[1];
rz(-1.9531293) q[3];
sx q[3];
rz(-2.5075932) q[3];
sx q[3];
rz(1.7651737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3153136) q[2];
sx q[2];
rz(-1.4297012) q[2];
sx q[2];
rz(-0.86177525) q[2];
rz(1.9050542) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(-1.2110565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(1.030141) q[0];
rz(2.8252699) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(0.28824678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22209805) q[0];
sx q[0];
rz(-1.4256501) q[0];
sx q[0];
rz(0.40444379) q[0];
rz(-pi) q[1];
rz(2.4924303) q[2];
sx q[2];
rz(-1.9948975) q[2];
sx q[2];
rz(2.8555388) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8178351) q[1];
sx q[1];
rz(-2.164748) q[1];
sx q[1];
rz(-0.51333921) q[1];
x q[2];
rz(2.4699798) q[3];
sx q[3];
rz(-1.3928431) q[3];
sx q[3];
rz(-0.88615299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2890702) q[2];
sx q[2];
rz(-1.4129637) q[2];
sx q[2];
rz(-2.9150325) q[2];
rz(1.6864927) q[3];
sx q[3];
rz(-0.89613599) q[3];
sx q[3];
rz(-2.4494825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15014547) q[0];
sx q[0];
rz(-1.8658072) q[0];
sx q[0];
rz(-2.5164497) q[0];
rz(1.217968) q[1];
sx q[1];
rz(-2.4324799) q[1];
sx q[1];
rz(1.2120754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8344515) q[0];
sx q[0];
rz(-1.092415) q[0];
sx q[0];
rz(2.4605453) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4702414) q[2];
sx q[2];
rz(-1.7749447) q[2];
sx q[2];
rz(2.8189557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21704504) q[1];
sx q[1];
rz(-2.6967561) q[1];
sx q[1];
rz(-1.43631) q[1];
rz(-pi) q[2];
rz(2.8404866) q[3];
sx q[3];
rz(-1.8499377) q[3];
sx q[3];
rz(-3.0413586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43442279) q[2];
sx q[2];
rz(-1.7662798) q[2];
sx q[2];
rz(-2.0909069) q[2];
rz(0.55189842) q[3];
sx q[3];
rz(-2.4514908) q[3];
sx q[3];
rz(-1.2845854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5156749) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(2.3552409) q[1];
sx q[1];
rz(-0.66023371) q[1];
sx q[1];
rz(0.22088851) q[1];
rz(3.1037504) q[2];
sx q[2];
rz(-1.2076245) q[2];
sx q[2];
rz(-2.5943499) q[2];
rz(0.68610739) q[3];
sx q[3];
rz(-1.5024363) q[3];
sx q[3];
rz(-0.92492044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
