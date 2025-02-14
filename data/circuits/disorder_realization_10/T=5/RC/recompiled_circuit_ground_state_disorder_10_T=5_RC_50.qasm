OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6634231) q[0];
sx q[0];
rz(-2.2278251) q[0];
sx q[0];
rz(-1.0650286) q[0];
rz(-1.5364667) q[1];
sx q[1];
rz(4.2470266) q[1];
sx q[1];
rz(12.469835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6632412) q[0];
sx q[0];
rz(-0.4718726) q[0];
sx q[0];
rz(2.0357516) q[0];
rz(2.7235254) q[2];
sx q[2];
rz(-2.5804248) q[2];
sx q[2];
rz(-2.9594066) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7207011) q[1];
sx q[1];
rz(-1.8681453) q[1];
sx q[1];
rz(2.4908476) q[1];
rz(-pi) q[2];
rz(-0.99081087) q[3];
sx q[3];
rz(-1.4926877) q[3];
sx q[3];
rz(-3.1171321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1221293) q[2];
sx q[2];
rz(-0.18522842) q[2];
sx q[2];
rz(-1.2629925) q[2];
rz(2.7428135) q[3];
sx q[3];
rz(-1.140927) q[3];
sx q[3];
rz(-2.1122475) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1102092) q[0];
sx q[0];
rz(-0.6321913) q[0];
sx q[0];
rz(1.6803918) q[0];
rz(3.0714495) q[1];
sx q[1];
rz(-1.5488011) q[1];
sx q[1];
rz(2.7102914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.778743) q[0];
sx q[0];
rz(-1.6322789) q[0];
sx q[0];
rz(0.009470609) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8441418) q[2];
sx q[2];
rz(-1.1701726) q[2];
sx q[2];
rz(-1.144415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0491534) q[1];
sx q[1];
rz(-0.91040694) q[1];
sx q[1];
rz(-2.4699442) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5502843) q[3];
sx q[3];
rz(-1.6153803) q[3];
sx q[3];
rz(-1.4449028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18616072) q[2];
sx q[2];
rz(-1.4107076) q[2];
sx q[2];
rz(0.53039941) q[2];
rz(-1.0839328) q[3];
sx q[3];
rz(-1.089774) q[3];
sx q[3];
rz(-1.2497905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0386061) q[0];
sx q[0];
rz(-0.57612053) q[0];
sx q[0];
rz(1.9392133) q[0];
rz(-2.1309958) q[1];
sx q[1];
rz(-1.3253515) q[1];
sx q[1];
rz(-3.1114263) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8081431) q[0];
sx q[0];
rz(-1.6177869) q[0];
sx q[0];
rz(0.065568173) q[0];
rz(-pi) q[1];
rz(-1.6232477) q[2];
sx q[2];
rz(-1.5546335) q[2];
sx q[2];
rz(-0.22833041) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29459965) q[1];
sx q[1];
rz(-1.2275808) q[1];
sx q[1];
rz(-1.7061451) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7272337) q[3];
sx q[3];
rz(-2.5632901) q[3];
sx q[3];
rz(2.3012637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5304337) q[2];
sx q[2];
rz(-0.23808372) q[2];
sx q[2];
rz(-2.7670624) q[2];
rz(1.1206867) q[3];
sx q[3];
rz(-1.4846385) q[3];
sx q[3];
rz(2.4499031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.197914) q[0];
sx q[0];
rz(-1.2926084) q[0];
sx q[0];
rz(-1.9669272) q[0];
rz(-1.3150175) q[1];
sx q[1];
rz(-2.0349793) q[1];
sx q[1];
rz(2.9768501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.574802) q[0];
sx q[0];
rz(-1.0079449) q[0];
sx q[0];
rz(0.99250162) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4061178) q[2];
sx q[2];
rz(-2.0782172) q[2];
sx q[2];
rz(-0.31854445) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10197266) q[1];
sx q[1];
rz(-2.0594739) q[1];
sx q[1];
rz(-2.7571329) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.03892558) q[3];
sx q[3];
rz(-1.4994326) q[3];
sx q[3];
rz(-1.4369278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4508007) q[2];
sx q[2];
rz(-1.5361293) q[2];
sx q[2];
rz(2.6353321) q[2];
rz(0.51044232) q[3];
sx q[3];
rz(-2.3531395) q[3];
sx q[3];
rz(0.81152624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4617758) q[0];
sx q[0];
rz(-2.8051069) q[0];
sx q[0];
rz(0.80999723) q[0];
rz(-3.0432155) q[1];
sx q[1];
rz(-2.7521303) q[1];
sx q[1];
rz(0.32442763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3379441) q[0];
sx q[0];
rz(-1.6932825) q[0];
sx q[0];
rz(3.1396559) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7067028) q[2];
sx q[2];
rz(-1.7909484) q[2];
sx q[2];
rz(-2.1957601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4091332) q[1];
sx q[1];
rz(-2.1420825) q[1];
sx q[1];
rz(-1.7833227) q[1];
rz(-0.39119519) q[3];
sx q[3];
rz(-2.6325691) q[3];
sx q[3];
rz(2.5768821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7847932) q[2];
sx q[2];
rz(-1.5609799) q[2];
sx q[2];
rz(-0.55850935) q[2];
rz(-1.8465346) q[3];
sx q[3];
rz(-1.4174771) q[3];
sx q[3];
rz(-2.7891125) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2297939) q[0];
sx q[0];
rz(-2.7138382) q[0];
sx q[0];
rz(1.3707004) q[0];
rz(1.3937344) q[1];
sx q[1];
rz(-1.0153898) q[1];
sx q[1];
rz(-0.14399354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8407366) q[0];
sx q[0];
rz(-1.6413781) q[0];
sx q[0];
rz(1.4665804) q[0];
rz(-pi) q[1];
rz(-1.1156111) q[2];
sx q[2];
rz(-1.3521638) q[2];
sx q[2];
rz(1.3022547) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3578971) q[1];
sx q[1];
rz(-1.6829957) q[1];
sx q[1];
rz(0.96331255) q[1];
rz(-pi) q[2];
rz(2.6089588) q[3];
sx q[3];
rz(-2.4264779) q[3];
sx q[3];
rz(1.2330513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7250942) q[2];
sx q[2];
rz(-2.7855253) q[2];
sx q[2];
rz(0.98102513) q[2];
rz(0.52870005) q[3];
sx q[3];
rz(-1.7873535) q[3];
sx q[3];
rz(-0.56021148) q[3];
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
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33787802) q[0];
sx q[0];
rz(-0.80334544) q[0];
sx q[0];
rz(1.2038318) q[0];
rz(3.1375258) q[1];
sx q[1];
rz(-1.7151567) q[1];
sx q[1];
rz(-0.34122658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26449152) q[0];
sx q[0];
rz(-0.84853338) q[0];
sx q[0];
rz(2.5654656) q[0];
x q[1];
rz(-1.5690825) q[2];
sx q[2];
rz(-0.69134635) q[2];
sx q[2];
rz(1.6328263) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7497037) q[1];
sx q[1];
rz(-0.14924696) q[1];
sx q[1];
rz(3.0796389) q[1];
rz(-pi) q[2];
rz(-1.4248203) q[3];
sx q[3];
rz(-1.2314704) q[3];
sx q[3];
rz(1.4682352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5228806) q[2];
sx q[2];
rz(-1.2423542) q[2];
sx q[2];
rz(0.27113554) q[2];
rz(1.6297657) q[3];
sx q[3];
rz(-2.1745493) q[3];
sx q[3];
rz(0.40201521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4790344) q[0];
sx q[0];
rz(-2.5738578) q[0];
sx q[0];
rz(0.16831368) q[0];
rz(1.0100826) q[1];
sx q[1];
rz(-1.1770266) q[1];
sx q[1];
rz(3.041306) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71609488) q[0];
sx q[0];
rz(-1.8779813) q[0];
sx q[0];
rz(-1.2745122) q[0];
rz(-pi) q[1];
rz(-1.3919207) q[2];
sx q[2];
rz(-2.2958307) q[2];
sx q[2];
rz(-2.7357312) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0834962) q[1];
sx q[1];
rz(-1.5952949) q[1];
sx q[1];
rz(2.9398838) q[1];
rz(-pi) q[2];
rz(0.69117213) q[3];
sx q[3];
rz(-0.6414957) q[3];
sx q[3];
rz(2.3506899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68891305) q[2];
sx q[2];
rz(-0.49892384) q[2];
sx q[2];
rz(-1.720724) q[2];
rz(-2.0776599) q[3];
sx q[3];
rz(-1.7199793) q[3];
sx q[3];
rz(-0.091778278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38561472) q[0];
sx q[0];
rz(-1.6843963) q[0];
sx q[0];
rz(0.036855999) q[0];
rz(-1.8846177) q[1];
sx q[1];
rz(-1.6391862) q[1];
sx q[1];
rz(-1.3704376) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12808558) q[0];
sx q[0];
rz(-0.31541809) q[0];
sx q[0];
rz(-2.0924139) q[0];
x q[1];
rz(-2.7671934) q[2];
sx q[2];
rz(-2.7485195) q[2];
sx q[2];
rz(-0.094660195) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6681666) q[1];
sx q[1];
rz(-1.2915732) q[1];
sx q[1];
rz(2.2609026) q[1];
x q[2];
rz(0.044058756) q[3];
sx q[3];
rz(-1.9875196) q[3];
sx q[3];
rz(-2.8686881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86159414) q[2];
sx q[2];
rz(-1.7656606) q[2];
sx q[2];
rz(2.9691248) q[2];
rz(2.5837768) q[3];
sx q[3];
rz(-0.54111257) q[3];
sx q[3];
rz(-1.4403758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7114792) q[0];
sx q[0];
rz(-2.5028296) q[0];
sx q[0];
rz(-0.3399671) q[0];
rz(-2.4760447) q[1];
sx q[1];
rz(-1.7219209) q[1];
sx q[1];
rz(-1.8284214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198044) q[0];
sx q[0];
rz(-1.7504842) q[0];
sx q[0];
rz(-1.6186886) q[0];
rz(-pi) q[1];
rz(-0.27965172) q[2];
sx q[2];
rz(-0.2708098) q[2];
sx q[2];
rz(0.93284399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39496702) q[1];
sx q[1];
rz(-1.2685742) q[1];
sx q[1];
rz(2.8923467) q[1];
rz(-pi) q[2];
rz(-0.35452133) q[3];
sx q[3];
rz(-2.7149622) q[3];
sx q[3];
rz(1.420295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36249915) q[2];
sx q[2];
rz(-2.7069147) q[2];
sx q[2];
rz(0.67081007) q[2];
rz(0.32495156) q[3];
sx q[3];
rz(-0.81348014) q[3];
sx q[3];
rz(-1.0130829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.308607) q[0];
sx q[0];
rz(-1.2815463) q[0];
sx q[0];
rz(1.4776342) q[0];
rz(1.0302522) q[1];
sx q[1];
rz(-1.6696842) q[1];
sx q[1];
rz(-1.3546863) q[1];
rz(2.3823635) q[2];
sx q[2];
rz(-0.6741796) q[2];
sx q[2];
rz(0.37311935) q[2];
rz(-1.9454789) q[3];
sx q[3];
rz(-2.119172) q[3];
sx q[3];
rz(-0.86772002) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
