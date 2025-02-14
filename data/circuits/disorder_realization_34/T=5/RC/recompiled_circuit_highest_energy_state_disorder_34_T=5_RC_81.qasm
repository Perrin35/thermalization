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
rz(-1.7582769) q[0];
sx q[0];
rz(-1.7051237) q[0];
sx q[0];
rz(-0.96487784) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(-0.42958346) q[1];
sx q[1];
rz(2.092195) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222264) q[0];
sx q[0];
rz(-0.99577409) q[0];
sx q[0];
rz(-0.40798835) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8951178) q[2];
sx q[2];
rz(-1.7056381) q[2];
sx q[2];
rz(-2.2244649) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2856372) q[1];
sx q[1];
rz(-0.95889839) q[1];
sx q[1];
rz(0.4680856) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4780294) q[3];
sx q[3];
rz(-0.99348196) q[3];
sx q[3];
rz(2.4924459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6875978) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(-0.018608658) q[2];
rz(0.54443693) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740771) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(-0.53502214) q[0];
rz(2.6016443) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.7417057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2311418) q[0];
sx q[0];
rz(-0.60908356) q[0];
sx q[0];
rz(-2.2876491) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1715047) q[2];
sx q[2];
rz(-2.0510489) q[2];
sx q[2];
rz(2.9532972) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4606159) q[1];
sx q[1];
rz(-1.1910805) q[1];
sx q[1];
rz(-1.3028076) q[1];
rz(-1.0008282) q[3];
sx q[3];
rz(-2.1080351) q[3];
sx q[3];
rz(1.9443823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0672368) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(2.2155217) q[2];
rz(-2.1615084) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3493018) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(1.6287623) q[0];
rz(-0.28383645) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(-1.1121174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4130198) q[0];
sx q[0];
rz(-1.5637102) q[0];
sx q[0];
rz(1.5779693) q[0];
x q[1];
rz(3.0653333) q[2];
sx q[2];
rz(-1.462655) q[2];
sx q[2];
rz(2.0881483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2123472) q[1];
sx q[1];
rz(-0.51873365) q[1];
sx q[1];
rz(0.38350819) q[1];
rz(-pi) q[2];
rz(-0.25896163) q[3];
sx q[3];
rz(-1.9230712) q[3];
sx q[3];
rz(-1.6083628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6185559) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(-0.88895041) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717644) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(2.4523822) q[0];
rz(2.8269732) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.4124195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1293761) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(1.1393121) q[0];
x q[1];
rz(2.7233549) q[2];
sx q[2];
rz(-0.23242885) q[2];
sx q[2];
rz(1.7247891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1656944) q[1];
sx q[1];
rz(-1.9503352) q[1];
sx q[1];
rz(1.6326153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9474218) q[3];
sx q[3];
rz(-1.0946858) q[3];
sx q[3];
rz(0.74060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7249001) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(-2.5856384) q[2];
rz(-1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(-2.8506193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274662) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(0.59447527) q[0];
rz(-0.56198436) q[1];
sx q[1];
rz(-2.2832506) q[1];
sx q[1];
rz(2.1655703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291479) q[0];
sx q[0];
rz(-1.2093822) q[0];
sx q[0];
rz(1.1114208) q[0];
rz(-pi) q[1];
x q[1];
rz(1.541829) q[2];
sx q[2];
rz(-0.81695643) q[2];
sx q[2];
rz(-3.0532233) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1072717) q[1];
sx q[1];
rz(-1.2731421) q[1];
sx q[1];
rz(-0.68796449) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30140437) q[3];
sx q[3];
rz(-1.7782974) q[3];
sx q[3];
rz(-1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48656616) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(-2.5308934) q[2];
rz(-0.30580172) q[3];
sx q[3];
rz(-2.1761201) q[3];
sx q[3];
rz(1.3293728) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82764757) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(1.4661283) q[0];
rz(-2.3620391) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-0.73878845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29863623) q[0];
sx q[0];
rz(-1.8514575) q[0];
sx q[0];
rz(-2.929351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5601842) q[2];
sx q[2];
rz(-1.2640177) q[2];
sx q[2];
rz(0.58394428) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.14911095) q[1];
sx q[1];
rz(-2.457239) q[1];
sx q[1];
rz(0.97196399) q[1];
rz(2.3380558) q[3];
sx q[3];
rz(-2.6430686) q[3];
sx q[3];
rz(0.62832181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5256727) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(0.8626779) q[2];
rz(-2.681813) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.9988683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90726844) q[0];
sx q[0];
rz(-2.7643804) q[0];
sx q[0];
rz(1.020485) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(-2.3540156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0840698) q[0];
sx q[0];
rz(-0.935874) q[0];
sx q[0];
rz(-1.8674839) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6422945) q[2];
sx q[2];
rz(-1.8282969) q[2];
sx q[2];
rz(1.471164) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1169402) q[1];
sx q[1];
rz(-1.7491915) q[1];
sx q[1];
rz(-2.0703038) q[1];
rz(-0.98290409) q[3];
sx q[3];
rz(-1.5994775) q[3];
sx q[3];
rz(2.8010288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82137498) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(2.5569432) q[2];
rz(-2.4892877) q[3];
sx q[3];
rz(-1.9125331) q[3];
sx q[3];
rz(-1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24340165) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(2.8133494) q[0];
rz(-2.7499061) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(-1.4788871) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680964) q[0];
sx q[0];
rz(-0.35051051) q[0];
sx q[0];
rz(-0.71061937) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52430341) q[2];
sx q[2];
rz(-1.1534165) q[2];
sx q[2];
rz(0.5988754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6819671) q[1];
sx q[1];
rz(-2.8468968) q[1];
sx q[1];
rz(-1.2380283) q[1];
rz(-pi) q[2];
rz(1.0669054) q[3];
sx q[3];
rz(-2.7253838) q[3];
sx q[3];
rz(-1.4598626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0216003) q[2];
sx q[2];
rz(-1.768521) q[2];
sx q[2];
rz(-0.78869406) q[2];
rz(-0.24104077) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(-2.327976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64860827) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(2.1797144) q[0];
rz(2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-0.73807565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3320112) q[0];
sx q[0];
rz(-1.9652848) q[0];
sx q[0];
rz(1.3138735) q[0];
x q[1];
rz(1.1253276) q[2];
sx q[2];
rz(-2.1574852) q[2];
sx q[2];
rz(2.4054804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9540087) q[1];
sx q[1];
rz(-2.6212647) q[1];
sx q[1];
rz(2.1854758) q[1];
rz(-pi) q[2];
rz(2.4752615) q[3];
sx q[3];
rz(-1.7937346) q[3];
sx q[3];
rz(1.8458888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6568079) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(-1.3336746) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(-0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8986847) q[0];
sx q[0];
rz(-0.81166357) q[0];
sx q[0];
rz(2.2743478) q[0];
rz(2.1278837) q[1];
sx q[1];
rz(-1.3288682) q[1];
sx q[1];
rz(2.0713461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9673706) q[0];
sx q[0];
rz(-1.692932) q[0];
sx q[0];
rz(-1.4184667) q[0];
x q[1];
rz(-2.4325244) q[2];
sx q[2];
rz(-0.96154172) q[2];
sx q[2];
rz(2.9035062) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5336736) q[1];
sx q[1];
rz(-1.1864788) q[1];
sx q[1];
rz(-2.584105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66859122) q[3];
sx q[3];
rz(-2.9644659) q[3];
sx q[3];
rz(1.1531545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16537198) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(1.9099859) q[2];
rz(-1.5008789) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(2.4098082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-0.41481836) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(0.11484405) q[2];
sx q[2];
rz(-2.6931802) q[2];
sx q[2];
rz(-0.32297011) q[2];
rz(0.75178643) q[3];
sx q[3];
rz(-1.1204168) q[3];
sx q[3];
rz(0.71747019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
