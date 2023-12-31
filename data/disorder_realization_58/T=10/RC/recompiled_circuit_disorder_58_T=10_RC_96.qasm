OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2031189) q[0];
sx q[0];
rz(-1.9460558) q[0];
sx q[0];
rz(1.24036) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(-2.8505461) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0349866) q[1];
sx q[1];
rz(-2.258746) q[1];
sx q[1];
rz(-2.7290542) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72676267) q[3];
sx q[3];
rz(-2.5708377) q[3];
sx q[3];
rz(1.4445514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-0.1201771) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90601901) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(2.2554382) q[0];
rz(-pi) q[1];
rz(1.1268483) q[2];
sx q[2];
rz(-1.4269281) q[2];
sx q[2];
rz(-0.680188) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3152299) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(0.96674322) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88300206) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(-1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(1.057391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085416) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(3.1197381) q[0];
rz(-pi) q[1];
rz(-0.34122841) q[2];
sx q[2];
rz(-1.382302) q[2];
sx q[2];
rz(-0.54268062) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1765269) q[1];
sx q[1];
rz(-0.9170734) q[1];
sx q[1];
rz(0.42018004) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4226513) q[3];
sx q[3];
rz(-1.395441) q[3];
sx q[3];
rz(2.5734176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.3624297) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-2.9761956) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20277682) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(-1.3846272) q[0];
x q[1];
rz(2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(-1.6844815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0927825) q[1];
sx q[1];
rz(-2.4895992) q[1];
sx q[1];
rz(2.584143) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(-1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(0.35476312) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(-2.9096471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9493128) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(-0.85286661) q[0];
rz(-pi) q[1];
rz(-0.69552341) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(2.7930789) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9310589) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(-1.348043) q[1];
x q[2];
rz(-0.050458126) q[3];
sx q[3];
rz(-2.7726463) q[3];
sx q[3];
rz(-1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0971138) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4076685) q[0];
sx q[0];
rz(-1.4364388) q[0];
sx q[0];
rz(0.863048) q[0];
x q[1];
rz(2.7315797) q[2];
sx q[2];
rz(-0.7747246) q[2];
sx q[2];
rz(-1.1020401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4221229) q[1];
sx q[1];
rz(-0.49233961) q[1];
sx q[1];
rz(-1.0455529) q[1];
x q[2];
rz(-0.44316767) q[3];
sx q[3];
rz(-2.245095) q[3];
sx q[3];
rz(-0.5815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-0.61378941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44871556) q[0];
sx q[0];
rz(-0.84445092) q[0];
sx q[0];
rz(-2.8918299) q[0];
rz(-1.6348398) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(-1.3751021) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1544513) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(0.51388545) q[1];
x q[2];
rz(0.11717637) q[3];
sx q[3];
rz(-2.0332608) q[3];
sx q[3];
rz(1.5069435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293837) q[0];
sx q[0];
rz(-1.0249656) q[0];
sx q[0];
rz(1.5270385) q[0];
rz(2.7266399) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(0.50843898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24482803) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(-pi) q[2];
rz(-0.55926178) q[3];
sx q[3];
rz(-2.5945633) q[3];
sx q[3];
rz(-0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(2.7507239) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(0.92591441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1122702) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(0.14365833) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7943031) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(1.6910451) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3013819) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(-1.128003) q[1];
rz(0.52243201) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(-2.0731376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9459076) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(2.6055028) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(2.7465903) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(2.1283456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.917406) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(-1.805472) q[0];
x q[1];
rz(1.6147862) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(-2.8994438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80876795) q[1];
sx q[1];
rz(-1.1507532) q[1];
sx q[1];
rz(1.4648449) q[1];
rz(-pi) q[2];
rz(-2.3771044) q[3];
sx q[3];
rz(-1.0478684) q[3];
sx q[3];
rz(-0.58512277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41291819) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(0.57327523) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(2.2686601) q[3];
sx q[3];
rz(-1.6193661) q[3];
sx q[3];
rz(-1.247874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
