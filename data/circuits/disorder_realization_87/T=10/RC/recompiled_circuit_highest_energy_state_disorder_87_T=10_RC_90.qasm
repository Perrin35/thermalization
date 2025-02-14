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
rz(1.4662161) q[0];
sx q[0];
rz(5.3386547) q[0];
sx q[0];
rz(9.6261779) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(-1.720517) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11036377) q[0];
sx q[0];
rz(-2.5284934) q[0];
sx q[0];
rz(-1.4754063) q[0];
x q[1];
rz(-2.6236963) q[2];
sx q[2];
rz(-2.1764048) q[2];
sx q[2];
rz(-0.54010682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.046467543) q[1];
sx q[1];
rz(-0.78739843) q[1];
sx q[1];
rz(-0.56708972) q[1];
rz(-pi) q[2];
rz(2.8086923) q[3];
sx q[3];
rz(-1.1984636) q[3];
sx q[3];
rz(1.8474735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.034915514) q[2];
sx q[2];
rz(-0.80433977) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(0.15160027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5308373) q[0];
sx q[0];
rz(-2.294367) q[0];
sx q[0];
rz(0.24638677) q[0];
rz(1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.5066719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6794066) q[0];
sx q[0];
rz(-0.52216086) q[0];
sx q[0];
rz(-1.0617816) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7072524) q[2];
sx q[2];
rz(-2.8359957) q[2];
sx q[2];
rz(-1.9836677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40827258) q[1];
sx q[1];
rz(-0.36767861) q[1];
sx q[1];
rz(1.4094844) q[1];
x q[2];
rz(-2.706932) q[3];
sx q[3];
rz(-2.0512329) q[3];
sx q[3];
rz(0.91033376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(1.9096036) q[2];
rz(2.039382) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(-1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6983637) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(2.2337636) q[0];
rz(-2.5746131) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(3.0388015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.028057) q[0];
sx q[0];
rz(-0.60657078) q[0];
sx q[0];
rz(-0.3740749) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91102131) q[2];
sx q[2];
rz(-3.0098923) q[2];
sx q[2];
rz(0.98051276) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6425848) q[1];
sx q[1];
rz(-2.4053768) q[1];
sx q[1];
rz(2.5166488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1249173) q[3];
sx q[3];
rz(-1.8401056) q[3];
sx q[3];
rz(-0.96162063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28321442) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(-1.4374479) q[2];
rz(-2.7490859) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(-0.2984305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0482386) q[0];
sx q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(1.5123222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8522569) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(2.1136978) q[0];
rz(-pi) q[1];
rz(-1.9061162) q[2];
sx q[2];
rz(-2.309707) q[2];
sx q[2];
rz(-0.12140935) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61002953) q[1];
sx q[1];
rz(-1.3852784) q[1];
sx q[1];
rz(-0.97037913) q[1];
rz(1.063594) q[3];
sx q[3];
rz(-1.385687) q[3];
sx q[3];
rz(2.0349658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7901223) q[2];
sx q[2];
rz(-1.0834379) q[2];
sx q[2];
rz(-2.4646087) q[2];
rz(-1.8949159) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(-1.6251132) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93382728) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(3.1166792) q[0];
rz(1.0554396) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-2.8581462) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0830529) q[0];
sx q[0];
rz(-1.827936) q[0];
sx q[0];
rz(-1.8399946) q[0];
x q[1];
rz(1.7666498) q[2];
sx q[2];
rz(-1.4417267) q[2];
sx q[2];
rz(1.376525) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1920212) q[1];
sx q[1];
rz(-1.5402435) q[1];
sx q[1];
rz(1.0699238) q[1];
rz(-0.33971649) q[3];
sx q[3];
rz(-0.92631683) q[3];
sx q[3];
rz(-3.0655053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31263605) q[2];
sx q[2];
rz(-2.7241311) q[2];
sx q[2];
rz(-2.8157595) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7376937) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(2.7666336) q[0];
rz(-2.6629579) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(-2.7714444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86287303) q[0];
sx q[0];
rz(-1.5939043) q[0];
sx q[0];
rz(-1.5824759) q[0];
rz(-pi) q[1];
rz(-0.13932087) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(0.56688165) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83997516) q[1];
sx q[1];
rz(-2.2711922) q[1];
sx q[1];
rz(2.5216334) q[1];
x q[2];
rz(-2.7137854) q[3];
sx q[3];
rz(-1.5753581) q[3];
sx q[3];
rz(3.004194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39773539) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.7605304) q[2];
rz(-2.0928275) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(2.5019808) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25310707) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-2.2174368) q[0];
rz(-2.2249075) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(-1.4013269) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2271254) q[0];
sx q[0];
rz(-1.0773398) q[0];
sx q[0];
rz(-1.5052133) q[0];
rz(-0.6551459) q[2];
sx q[2];
rz(-1.4339951) q[2];
sx q[2];
rz(3.0828147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37937134) q[1];
sx q[1];
rz(-1.7306384) q[1];
sx q[1];
rz(-2.400378) q[1];
x q[2];
rz(2.9244366) q[3];
sx q[3];
rz(-1.0921061) q[3];
sx q[3];
rz(-1.1622805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.4618123) q[2];
sx q[2];
rz(-1.224996) q[2];
rz(3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5834354) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(0.48026568) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(0.86437782) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40369331) q[0];
sx q[0];
rz(-1.4072352) q[0];
sx q[0];
rz(-3.0944194) q[0];
rz(2.7988465) q[2];
sx q[2];
rz(-0.10048332) q[2];
sx q[2];
rz(-3.036694) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56039366) q[1];
sx q[1];
rz(-0.59948363) q[1];
sx q[1];
rz(-2.0066714) q[1];
rz(-pi) q[2];
rz(0.79473991) q[3];
sx q[3];
rz(-2.7594372) q[3];
sx q[3];
rz(2.4433745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(2.5065191) q[2];
rz(-2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(-2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315345) q[0];
sx q[0];
rz(-2.3681971) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(0.36525137) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(-1.7100547) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99125049) q[0];
sx q[0];
rz(-1.3936773) q[0];
sx q[0];
rz(2.9464673) q[0];
rz(-2.4875239) q[2];
sx q[2];
rz(-2.0383012) q[2];
sx q[2];
rz(1.2142177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4540951) q[1];
sx q[1];
rz(-2.3907067) q[1];
sx q[1];
rz(2.6862957) q[1];
rz(-pi) q[2];
rz(-2.6142653) q[3];
sx q[3];
rz(-1.9550712) q[3];
sx q[3];
rz(1.7155855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26620904) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(1.1727775) q[2];
rz(-1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(-1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239546) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(0.59463516) q[0];
rz(-1.0058962) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(-0.88919052) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62458663) q[0];
sx q[0];
rz(-1.3505624) q[0];
sx q[0];
rz(-0.051924719) q[0];
x q[1];
rz(-1.4290733) q[2];
sx q[2];
rz(-0.46438875) q[2];
sx q[2];
rz(1.022867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2500221) q[1];
sx q[1];
rz(-1.6305939) q[1];
sx q[1];
rz(-1.476261) q[1];
x q[2];
rz(-0.51639207) q[3];
sx q[3];
rz(-1.603873) q[3];
sx q[3];
rz(2.3865478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(-0.54538837) q[2];
rz(-0.96327463) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49190285) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(0.84359618) q[1];
sx q[1];
rz(-2.0633162) q[1];
sx q[1];
rz(1.8819527) q[1];
rz(-1.9764523) q[2];
sx q[2];
rz(-0.47274796) q[2];
sx q[2];
rz(-1.8190073) q[2];
rz(-0.58335494) q[3];
sx q[3];
rz(-1.1370549) q[3];
sx q[3];
rz(2.1255253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
