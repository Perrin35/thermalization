OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.759999) q[0];
sx q[0];
rz(-0.17188369) q[0];
sx q[0];
rz(2.5266393) q[0];
rz(1.3568658) q[1];
sx q[1];
rz(-1.5306127) q[1];
sx q[1];
rz(0.15712486) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9120795) q[0];
sx q[0];
rz(-2.3621971) q[0];
sx q[0];
rz(0.81672658) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31034361) q[2];
sx q[2];
rz(-2.9039798) q[2];
sx q[2];
rz(1.0050736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8871044) q[1];
sx q[1];
rz(-0.65782065) q[1];
sx q[1];
rz(-1.595934) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4218153) q[3];
sx q[3];
rz(-1.7474637) q[3];
sx q[3];
rz(-2.9797785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60363808) q[2];
sx q[2];
rz(-2.8261638) q[2];
sx q[2];
rz(1.8559378) q[2];
rz(-1.8042608) q[3];
sx q[3];
rz(-0.95271102) q[3];
sx q[3];
rz(0.11346909) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86785698) q[0];
sx q[0];
rz(-1.8524167) q[0];
sx q[0];
rz(-2.5751172) q[0];
rz(-1.6784338) q[1];
sx q[1];
rz(-0.17371829) q[1];
sx q[1];
rz(0.64243752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89680371) q[0];
sx q[0];
rz(-1.519108) q[0];
sx q[0];
rz(-2.0336853) q[0];
rz(-pi) q[1];
rz(1.240483) q[2];
sx q[2];
rz(-2.1871302) q[2];
sx q[2];
rz(-1.358145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21684619) q[1];
sx q[1];
rz(-1.5202151) q[1];
sx q[1];
rz(2.2994726) q[1];
rz(2.4107433) q[3];
sx q[3];
rz(-1.0789144) q[3];
sx q[3];
rz(3.0036894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88175875) q[2];
sx q[2];
rz(-0.79462785) q[2];
sx q[2];
rz(0.83414042) q[2];
rz(-0.086890876) q[3];
sx q[3];
rz(-2.0565242) q[3];
sx q[3];
rz(-2.0104008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6777545) q[0];
sx q[0];
rz(-1.0967655) q[0];
sx q[0];
rz(2.0820397) q[0];
rz(-0.73486596) q[1];
sx q[1];
rz(-0.015345416) q[1];
sx q[1];
rz(-0.14618348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37426078) q[0];
sx q[0];
rz(-2.0442937) q[0];
sx q[0];
rz(-2.9625508) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6017125) q[2];
sx q[2];
rz(-0.1242736) q[2];
sx q[2];
rz(-3.0310563) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7173269) q[1];
sx q[1];
rz(-1.4098412) q[1];
sx q[1];
rz(0.15086727) q[1];
rz(2.9133818) q[3];
sx q[3];
rz(-1.8442698) q[3];
sx q[3];
rz(0.21450689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77201468) q[2];
sx q[2];
rz(-0.39083189) q[2];
sx q[2];
rz(-2.3005627) q[2];
rz(3.0817025) q[3];
sx q[3];
rz(-1.6570897) q[3];
sx q[3];
rz(2.4976775) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.197072) q[0];
sx q[0];
rz(-0.50839013) q[0];
sx q[0];
rz(3.0149241) q[0];
rz(1.4177167) q[1];
sx q[1];
rz(-3.1210777) q[1];
sx q[1];
rz(-2.1518478) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5602172) q[0];
sx q[0];
rz(-0.84218431) q[0];
sx q[0];
rz(-1.5993662) q[0];
rz(-pi) q[1];
rz(0.085096882) q[2];
sx q[2];
rz(-1.9665045) q[2];
sx q[2];
rz(2.1425785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43112484) q[1];
sx q[1];
rz(-1.3763274) q[1];
sx q[1];
rz(-1.1072913) q[1];
rz(0.56405385) q[3];
sx q[3];
rz(-2.8587935) q[3];
sx q[3];
rz(-2.0334311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4385472) q[2];
sx q[2];
rz(-0.35312167) q[2];
sx q[2];
rz(-2.9959196) q[2];
rz(0.4438256) q[3];
sx q[3];
rz(-1.5747109) q[3];
sx q[3];
rz(2.023196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038427453) q[0];
sx q[0];
rz(-1.9157836) q[0];
sx q[0];
rz(-2.3577754) q[0];
rz(-1.8812284) q[1];
sx q[1];
rz(-0.0030219373) q[1];
sx q[1];
rz(0.22617117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40091915) q[0];
sx q[0];
rz(-0.52477389) q[0];
sx q[0];
rz(2.0361542) q[0];
x q[1];
rz(2.8402244) q[2];
sx q[2];
rz(-1.971481) q[2];
sx q[2];
rz(3.1189975) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5345787) q[1];
sx q[1];
rz(-1.4958989) q[1];
sx q[1];
rz(-1.0981047) q[1];
rz(-pi) q[2];
rz(-0.67490863) q[3];
sx q[3];
rz(-2.2332472) q[3];
sx q[3];
rz(-1.1250594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0164531) q[2];
sx q[2];
rz(-2.9235702) q[2];
sx q[2];
rz(1.3583292) q[2];
rz(-0.45768091) q[3];
sx q[3];
rz(-1.7031368) q[3];
sx q[3];
rz(2.5911736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0871171) q[0];
sx q[0];
rz(-3.1396907) q[0];
sx q[0];
rz(3.0892293) q[0];
rz(0.26854435) q[1];
sx q[1];
rz(-1.9895376) q[1];
sx q[1];
rz(0.23121887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58329952) q[0];
sx q[0];
rz(-0.31210408) q[0];
sx q[0];
rz(0.60404481) q[0];
rz(-pi) q[1];
rz(1.8753239) q[2];
sx q[2];
rz(-2.2923959) q[2];
sx q[2];
rz(0.90291427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0655122) q[1];
sx q[1];
rz(-1.339318) q[1];
sx q[1];
rz(-2.8006018) q[1];
x q[2];
rz(-2.1778132) q[3];
sx q[3];
rz(-1.4207109) q[3];
sx q[3];
rz(-0.039473783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3117567) q[2];
sx q[2];
rz(-0.51218963) q[2];
sx q[2];
rz(-0.83489418) q[2];
rz(3.0541776) q[3];
sx q[3];
rz(-1.089047) q[3];
sx q[3];
rz(-1.0138698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3111303) q[0];
sx q[0];
rz(-0.99263793) q[0];
sx q[0];
rz(-2.2845238) q[0];
rz(-2.3509707) q[1];
sx q[1];
rz(-3.1309083) q[1];
sx q[1];
rz(2.0649921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83504377) q[0];
sx q[0];
rz(-1.2694512) q[0];
sx q[0];
rz(-0.27923601) q[0];
rz(-pi) q[1];
rz(-2.1250379) q[2];
sx q[2];
rz(-0.86807251) q[2];
sx q[2];
rz(-1.6464485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8995114) q[1];
sx q[1];
rz(-0.82846009) q[1];
sx q[1];
rz(0.39023413) q[1];
rz(-pi) q[2];
rz(1.4348599) q[3];
sx q[3];
rz(-1.2316717) q[3];
sx q[3];
rz(-1.9423758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2662346) q[2];
sx q[2];
rz(-0.42576063) q[2];
sx q[2];
rz(-0.88179624) q[2];
rz(1.734123) q[3];
sx q[3];
rz(-2.2018645) q[3];
sx q[3];
rz(2.5772429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.8565916) q[0];
sx q[0];
rz(-0.0017310062) q[0];
sx q[0];
rz(0.28975394) q[0];
rz(-1.2334791) q[1];
sx q[1];
rz(-1.8461485) q[1];
sx q[1];
rz(2.8626056) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5279319) q[0];
sx q[0];
rz(-0.17027357) q[0];
sx q[0];
rz(-1.7199675) q[0];
rz(-2.5395168) q[2];
sx q[2];
rz(-0.47594949) q[2];
sx q[2];
rz(-2.0747831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.853914) q[1];
sx q[1];
rz(-2.8077237) q[1];
sx q[1];
rz(1.7823269) q[1];
rz(-pi) q[2];
rz(1.0881422) q[3];
sx q[3];
rz(-0.78018809) q[3];
sx q[3];
rz(1.7232945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4616619) q[2];
sx q[2];
rz(-1.8671904) q[2];
sx q[2];
rz(2.8604841) q[2];
rz(-0.57349652) q[3];
sx q[3];
rz(-2.1888013) q[3];
sx q[3];
rz(-0.41962418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24116521) q[0];
sx q[0];
rz(-0.92210162) q[0];
sx q[0];
rz(1.7036555) q[0];
rz(2.9519713) q[1];
sx q[1];
rz(-3.1260243) q[1];
sx q[1];
rz(1.9059034) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3251117) q[0];
sx q[0];
rz(-0.50704038) q[0];
sx q[0];
rz(1.1757686) q[0];
x q[1];
rz(0.81053712) q[2];
sx q[2];
rz(-0.97522465) q[2];
sx q[2];
rz(0.11071225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37285629) q[1];
sx q[1];
rz(-2.8017178) q[1];
sx q[1];
rz(2.0716785) q[1];
rz(-pi) q[2];
rz(-3.1240033) q[3];
sx q[3];
rz(-1.5776199) q[3];
sx q[3];
rz(1.1428892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29888657) q[2];
sx q[2];
rz(-1.3624531) q[2];
sx q[2];
rz(-2.5617808) q[2];
rz(-1.3696356) q[3];
sx q[3];
rz(-0.42391351) q[3];
sx q[3];
rz(0.4637318) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.46691) q[0];
sx q[0];
rz(-1.5002102) q[0];
sx q[0];
rz(1.7036194) q[0];
rz(-1.2314318) q[1];
sx q[1];
rz(-2.2302901) q[1];
sx q[1];
rz(-1.6447172) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2220331) q[0];
sx q[0];
rz(-0.78503937) q[0];
sx q[0];
rz(1.9505984) q[0];
rz(2.4980372) q[2];
sx q[2];
rz(-0.6236667) q[2];
sx q[2];
rz(0.31094956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0383427) q[1];
sx q[1];
rz(-1.5885067) q[1];
sx q[1];
rz(-1.5736386) q[1];
rz(-pi) q[2];
rz(-2.0541463) q[3];
sx q[3];
rz(-1.8385244) q[3];
sx q[3];
rz(-0.30743956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9897495) q[2];
sx q[2];
rz(-1.5630378) q[2];
sx q[2];
rz(2.6452276) q[2];
rz(2.1967891) q[3];
sx q[3];
rz(-3.1365518) q[3];
sx q[3];
rz(0.70281023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6481358) q[0];
sx q[0];
rz(-1.7127767) q[0];
sx q[0];
rz(2.0145804) q[0];
rz(1.5827178) q[1];
sx q[1];
rz(-2.8235148) q[1];
sx q[1];
rz(-2.943218) q[1];
rz(0.025540431) q[2];
sx q[2];
rz(-1.3141704) q[2];
sx q[2];
rz(-2.8258509) q[2];
rz(-0.55806969) q[3];
sx q[3];
rz(-1.2574099) q[3];
sx q[3];
rz(-1.3539061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
