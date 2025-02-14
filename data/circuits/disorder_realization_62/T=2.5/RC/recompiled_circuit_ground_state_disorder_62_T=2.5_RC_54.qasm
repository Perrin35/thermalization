OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5476721) q[0];
sx q[0];
rz(-0.5991109) q[0];
sx q[0];
rz(-0.17679086) q[0];
rz(2.0556567) q[1];
sx q[1];
rz(-2.63201) q[1];
sx q[1];
rz(-0.086960763) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3172042) q[0];
sx q[0];
rz(-1.9476711) q[0];
sx q[0];
rz(0.61610846) q[0];
rz(-0.29524191) q[2];
sx q[2];
rz(-0.59165162) q[2];
sx q[2];
rz(2.0830784) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1198169) q[1];
sx q[1];
rz(-1.4939185) q[1];
sx q[1];
rz(1.4426115) q[1];
x q[2];
rz(-1.8636892) q[3];
sx q[3];
rz(-0.97145069) q[3];
sx q[3];
rz(-1.0224316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95734346) q[2];
sx q[2];
rz(-2.6441296) q[2];
sx q[2];
rz(-2.5498665) q[2];
rz(2.7033778) q[3];
sx q[3];
rz(-0.44132909) q[3];
sx q[3];
rz(-2.2653968) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1205207) q[0];
sx q[0];
rz(-0.35728917) q[0];
sx q[0];
rz(-1.0907809) q[0];
rz(2.4202994) q[1];
sx q[1];
rz(-1.483016) q[1];
sx q[1];
rz(0.24478197) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1274968) q[0];
sx q[0];
rz(-1.2359706) q[0];
sx q[0];
rz(2.0398519) q[0];
rz(-0.82142395) q[2];
sx q[2];
rz(-1.9000179) q[2];
sx q[2];
rz(-0.034857817) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2657425) q[1];
sx q[1];
rz(-0.85803723) q[1];
sx q[1];
rz(-0.73709388) q[1];
rz(-pi) q[2];
rz(-0.29676389) q[3];
sx q[3];
rz(-1.7113842) q[3];
sx q[3];
rz(0.97749099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.81405866) q[2];
sx q[2];
rz(-2.5746097) q[2];
sx q[2];
rz(0.84612334) q[2];
rz(0.36492473) q[3];
sx q[3];
rz(-0.42607421) q[3];
sx q[3];
rz(2.164771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.1035128) q[0];
sx q[0];
rz(-0.80779034) q[0];
sx q[0];
rz(-0.14347759) q[0];
rz(1.4682651) q[1];
sx q[1];
rz(-1.1500618) q[1];
sx q[1];
rz(-0.066468261) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22035711) q[0];
sx q[0];
rz(-2.4548303) q[0];
sx q[0];
rz(-0.14089091) q[0];
rz(-pi) q[1];
rz(2.8861553) q[2];
sx q[2];
rz(-1.0424926) q[2];
sx q[2];
rz(-2.262399) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79083645) q[1];
sx q[1];
rz(-0.15201026) q[1];
sx q[1];
rz(-0.69614567) q[1];
rz(-0.18174882) q[3];
sx q[3];
rz(-2.3929993) q[3];
sx q[3];
rz(-0.68457097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1116144) q[2];
sx q[2];
rz(-0.51631236) q[2];
sx q[2];
rz(-2.8288793) q[2];
rz(-0.2615658) q[3];
sx q[3];
rz(-1.4489737) q[3];
sx q[3];
rz(0.79791445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0161491) q[0];
sx q[0];
rz(-0.29368547) q[0];
sx q[0];
rz(3.0545767) q[0];
rz(-1.7866987) q[1];
sx q[1];
rz(-1.5630629) q[1];
sx q[1];
rz(0.048390128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66358405) q[0];
sx q[0];
rz(-0.80246325) q[0];
sx q[0];
rz(-1.3924696) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87746967) q[2];
sx q[2];
rz(-2.9089768) q[2];
sx q[2];
rz(-0.86384976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9632513) q[1];
sx q[1];
rz(-1.6251441) q[1];
sx q[1];
rz(-2.4456294) q[1];
x q[2];
rz(2.2807924) q[3];
sx q[3];
rz(-1.080092) q[3];
sx q[3];
rz(1.2482289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4309569) q[2];
sx q[2];
rz(-2.841876) q[2];
sx q[2];
rz(2.7002913) q[2];
rz(-0.86853164) q[3];
sx q[3];
rz(-1.7732311) q[3];
sx q[3];
rz(1.8080447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61442536) q[0];
sx q[0];
rz(-1.704957) q[0];
sx q[0];
rz(-0.72702485) q[0];
rz(1.6983039) q[1];
sx q[1];
rz(-1.2579505) q[1];
sx q[1];
rz(-1.7194933) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073339065) q[0];
sx q[0];
rz(-1.710482) q[0];
sx q[0];
rz(0.54648593) q[0];
rz(2.431972) q[2];
sx q[2];
rz(-0.31236744) q[2];
sx q[2];
rz(1.720495) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4053012) q[1];
sx q[1];
rz(-1.4854238) q[1];
sx q[1];
rz(-3.0194204) q[1];
rz(-pi) q[2];
rz(-0.90713769) q[3];
sx q[3];
rz(-2.0454413) q[3];
sx q[3];
rz(1.8812219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21165851) q[2];
sx q[2];
rz(-2.5514166) q[2];
sx q[2];
rz(1.5606073) q[2];
rz(1.3382781) q[3];
sx q[3];
rz(-2.9542597) q[3];
sx q[3];
rz(-0.54792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.018589858) q[0];
sx q[0];
rz(-2.328861) q[0];
sx q[0];
rz(2.4216968) q[0];
rz(-2.9010991) q[1];
sx q[1];
rz(-1.1613107) q[1];
sx q[1];
rz(0.24615157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1937597) q[0];
sx q[0];
rz(-1.5475377) q[0];
sx q[0];
rz(0.2808397) q[0];
x q[1];
rz(-2.7163137) q[2];
sx q[2];
rz(-2.3375247) q[2];
sx q[2];
rz(-0.67415392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0726021) q[1];
sx q[1];
rz(-1.700987) q[1];
sx q[1];
rz(2.4058002) q[1];
x q[2];
rz(-1.9897377) q[3];
sx q[3];
rz(-2.2139858) q[3];
sx q[3];
rz(0.96574984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64421946) q[2];
sx q[2];
rz(-2.5770598) q[2];
sx q[2];
rz(-0.79968828) q[2];
rz(2.6175446) q[3];
sx q[3];
rz(-0.3862114) q[3];
sx q[3];
rz(-0.016949765) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2972357) q[0];
sx q[0];
rz(-3.0581664) q[0];
sx q[0];
rz(-2.2204087) q[0];
rz(-1.4211897) q[1];
sx q[1];
rz(-2.4589296) q[1];
sx q[1];
rz(-0.99501077) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1608022) q[0];
sx q[0];
rz(-1.1931927) q[0];
sx q[0];
rz(-0.23997216) q[0];
x q[1];
rz(-1.7172377) q[2];
sx q[2];
rz(-2.7173923) q[2];
sx q[2];
rz(-2.7643124) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9239569) q[1];
sx q[1];
rz(-1.991365) q[1];
sx q[1];
rz(2.5687508) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1345484) q[3];
sx q[3];
rz(-1.9552257) q[3];
sx q[3];
rz(-2.4151797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2034188) q[2];
sx q[2];
rz(-0.19762453) q[2];
sx q[2];
rz(-2.7040238) q[2];
rz(-2.3214052) q[3];
sx q[3];
rz(-1.5486251) q[3];
sx q[3];
rz(-0.13535132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1831128) q[0];
sx q[0];
rz(-3.0218229) q[0];
sx q[0];
rz(-2.130765) q[0];
rz(-3.0567567) q[1];
sx q[1];
rz(-1.1654221) q[1];
sx q[1];
rz(2.5929677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6609782) q[0];
sx q[0];
rz(-1.7465495) q[0];
sx q[0];
rz(0.032175933) q[0];
rz(-pi) q[1];
rz(0.63603129) q[2];
sx q[2];
rz(-2.6205728) q[2];
sx q[2];
rz(-0.65083671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1774166) q[1];
sx q[1];
rz(-0.89277041) q[1];
sx q[1];
rz(-1.3864338) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0579257) q[3];
sx q[3];
rz(-1.4691969) q[3];
sx q[3];
rz(-1.6179832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34714547) q[2];
sx q[2];
rz(-0.5793137) q[2];
sx q[2];
rz(0.53317201) q[2];
rz(-1.0779856) q[3];
sx q[3];
rz(-0.92107934) q[3];
sx q[3];
rz(-2.6326411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0549523) q[0];
sx q[0];
rz(-0.3594048) q[0];
sx q[0];
rz(2.3593498) q[0];
rz(3.0714463) q[1];
sx q[1];
rz(-0.47824305) q[1];
sx q[1];
rz(0.22629647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5180888) q[0];
sx q[0];
rz(-1.5020348) q[0];
sx q[0];
rz(-3.070773) q[0];
rz(-0.18780577) q[2];
sx q[2];
rz(-0.90131288) q[2];
sx q[2];
rz(1.737843) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8574684) q[1];
sx q[1];
rz(-2.5924666) q[1];
sx q[1];
rz(0.21128564) q[1];
rz(-pi) q[2];
x q[2];
rz(0.022458301) q[3];
sx q[3];
rz(-1.8216672) q[3];
sx q[3];
rz(-0.24652265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1066863) q[2];
sx q[2];
rz(-2.8539113) q[2];
sx q[2];
rz(1.7314343) q[2];
rz(0.26257026) q[3];
sx q[3];
rz(-1.5574484) q[3];
sx q[3];
rz(0.13172758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054963741) q[0];
sx q[0];
rz(-2.9294736) q[0];
sx q[0];
rz(0.18375272) q[0];
rz(0.19206583) q[1];
sx q[1];
rz(-1.4552677) q[1];
sx q[1];
rz(2.5180838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79270169) q[0];
sx q[0];
rz(-1.1109795) q[0];
sx q[0];
rz(1.9135273) q[0];
x q[1];
rz(-2.170606) q[2];
sx q[2];
rz(-1.0375334) q[2];
sx q[2];
rz(-0.87109921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19621721) q[1];
sx q[1];
rz(-1.5935285) q[1];
sx q[1];
rz(-3.0356221) q[1];
rz(-pi) q[2];
rz(3.1106408) q[3];
sx q[3];
rz(-1.1963468) q[3];
sx q[3];
rz(-2.6778145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3925675) q[2];
sx q[2];
rz(-0.6124658) q[2];
sx q[2];
rz(-0.037671063) q[2];
rz(0.41845775) q[3];
sx q[3];
rz(-0.28067121) q[3];
sx q[3];
rz(2.9627964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2035718) q[0];
sx q[0];
rz(-0.9104712) q[0];
sx q[0];
rz(-0.18785432) q[0];
rz(2.6899295) q[1];
sx q[1];
rz(-1.3126806) q[1];
sx q[1];
rz(-1.5246593) q[1];
rz(2.6653566) q[2];
sx q[2];
rz(-0.58544896) q[2];
sx q[2];
rz(0.61665012) q[2];
rz(-2.8658297) q[3];
sx q[3];
rz(-1.388474) q[3];
sx q[3];
rz(0.17920517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
