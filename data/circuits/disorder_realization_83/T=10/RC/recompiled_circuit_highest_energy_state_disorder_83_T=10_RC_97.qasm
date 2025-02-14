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
rz(-0.42521617) q[0];
sx q[0];
rz(-1.8073616) q[0];
sx q[0];
rz(-0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(0.92460728) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3596852) q[0];
sx q[0];
rz(-1.230002) q[0];
sx q[0];
rz(-2.4016892) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0546646) q[2];
sx q[2];
rz(-1.4251815) q[2];
sx q[2];
rz(2.9235385) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.315359) q[1];
sx q[1];
rz(-0.72972882) q[1];
sx q[1];
rz(1.5118096) q[1];
rz(-pi) q[2];
rz(-0.12792952) q[3];
sx q[3];
rz(-1.5618213) q[3];
sx q[3];
rz(-0.25024807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1251462) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(-0.98801405) q[2];
rz(0.2068578) q[3];
sx q[3];
rz(-1.7346953) q[3];
sx q[3];
rz(-0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67966953) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(2.8894506) q[0];
rz(3.120046) q[1];
sx q[1];
rz(-2.4485059) q[1];
sx q[1];
rz(-0.85539877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148672) q[0];
sx q[0];
rz(-1.5692737) q[0];
sx q[0];
rz(-1.5437838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62521817) q[2];
sx q[2];
rz(-0.65186912) q[2];
sx q[2];
rz(2.9232581) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8586848) q[1];
sx q[1];
rz(-0.58382422) q[1];
sx q[1];
rz(-1.1633384) q[1];
rz(-pi) q[2];
rz(2.7383201) q[3];
sx q[3];
rz(-0.72411116) q[3];
sx q[3];
rz(1.4693174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.97027332) q[2];
sx q[2];
rz(-0.0042985175) q[2];
sx q[2];
rz(-0.26101905) q[2];
rz(-0.039693443) q[3];
sx q[3];
rz(-1.7429765) q[3];
sx q[3];
rz(2.5980914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4259341) q[0];
sx q[0];
rz(-1.6682699) q[0];
sx q[0];
rz(0.69450992) q[0];
rz(-1.1272686) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(-0.99004254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309826) q[0];
sx q[0];
rz(-1.6525804) q[0];
sx q[0];
rz(-2.7922996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8672529) q[2];
sx q[2];
rz(-0.60288376) q[2];
sx q[2];
rz(-2.4566513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1676869) q[1];
sx q[1];
rz(-1.1264633) q[1];
sx q[1];
rz(-1.9698039) q[1];
x q[2];
rz(-2.6022791) q[3];
sx q[3];
rz(-0.94180543) q[3];
sx q[3];
rz(-1.1719538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8433044) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(-2.7847024) q[3];
sx q[3];
rz(-2.7481952) q[3];
sx q[3];
rz(1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056034293) q[0];
sx q[0];
rz(-1.9558676) q[0];
sx q[0];
rz(2.2991142) q[0];
rz(-0.8017686) q[1];
sx q[1];
rz(-0.27920488) q[1];
sx q[1];
rz(0.78559771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39488897) q[0];
sx q[0];
rz(-2.3869042) q[0];
sx q[0];
rz(-2.7695023) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7217721) q[2];
sx q[2];
rz(-0.85938007) q[2];
sx q[2];
rz(-0.47455088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29348865) q[1];
sx q[1];
rz(-1.2624718) q[1];
sx q[1];
rz(-0.56667324) q[1];
rz(-pi) q[2];
rz(-2.9369041) q[3];
sx q[3];
rz(-2.31143) q[3];
sx q[3];
rz(1.1799174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12513146) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(-2.1330736) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(-1.0803224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652902) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(-1.0705795) q[0];
rz(-0.77752441) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(-1.0106962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3883481) q[0];
sx q[0];
rz(-2.0439612) q[0];
sx q[0];
rz(-1.9761258) q[0];
rz(-pi) q[1];
rz(-3.0762482) q[2];
sx q[2];
rz(-1.3962708) q[2];
sx q[2];
rz(-3.0734398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7091211) q[1];
sx q[1];
rz(-0.84408954) q[1];
sx q[1];
rz(-2.5120887) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56212496) q[3];
sx q[3];
rz(-1.486612) q[3];
sx q[3];
rz(0.69417324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(-1.5555752) q[2];
rz(-1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(-1.4122081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3509336) q[0];
sx q[0];
rz(-0.18336329) q[0];
sx q[0];
rz(2.3405128) q[0];
rz(0.13310295) q[1];
sx q[1];
rz(-1.3214654) q[1];
sx q[1];
rz(-2.7395693) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4426081) q[0];
sx q[0];
rz(-0.8359209) q[0];
sx q[0];
rz(2.0280272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0021462321) q[2];
sx q[2];
rz(-1.8615926) q[2];
sx q[2];
rz(2.4418497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0401329) q[1];
sx q[1];
rz(-1.8454843) q[1];
sx q[1];
rz(2.0707612) q[1];
x q[2];
rz(1.6958713) q[3];
sx q[3];
rz(-0.91642028) q[3];
sx q[3];
rz(-0.47548166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.63783995) q[2];
sx q[2];
rz(-1.4053586) q[2];
sx q[2];
rz(-0.77581882) q[2];
rz(-1.6860298) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(-1.0078526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044947226) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(-0.67895472) q[0];
rz(1.2848162) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(-1.3791893) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5011936) q[0];
sx q[0];
rz(-1.7441347) q[0];
sx q[0];
rz(-0.79357432) q[0];
x q[1];
rz(1.2913338) q[2];
sx q[2];
rz(-1.6510291) q[2];
sx q[2];
rz(-0.81278518) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0028982) q[1];
sx q[1];
rz(-2.0636673) q[1];
sx q[1];
rz(-1.0652131) q[1];
rz(-pi) q[2];
rz(1.3706743) q[3];
sx q[3];
rz(-1.6255696) q[3];
sx q[3];
rz(1.2204264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8046367) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(1.2389368) q[2];
rz(1.3321446) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-0.24470394) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7234583) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(3.0031257) q[0];
rz(-2.9578517) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(0.26652452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375232) q[0];
sx q[0];
rz(-0.36629391) q[0];
sx q[0];
rz(2.0384602) q[0];
x q[1];
rz(1.1600003) q[2];
sx q[2];
rz(-2.5345384) q[2];
sx q[2];
rz(2.0701054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5268685) q[1];
sx q[1];
rz(-2.7182332) q[1];
sx q[1];
rz(1.7080659) q[1];
rz(-pi) q[2];
rz(2.4146904) q[3];
sx q[3];
rz(-1.5642526) q[3];
sx q[3];
rz(0.18101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0900241) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(-0.23452342) q[2];
rz(1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(-0.75638151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804857) q[0];
sx q[0];
rz(-0.8661626) q[0];
sx q[0];
rz(-0.79972237) q[0];
rz(0.58894482) q[1];
sx q[1];
rz(-0.84732333) q[1];
sx q[1];
rz(-0.2074997) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8556535) q[0];
sx q[0];
rz(-2.7253599) q[0];
sx q[0];
rz(-2.8134384) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8325808) q[2];
sx q[2];
rz(-0.88243077) q[2];
sx q[2];
rz(-2.7613044) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5821857) q[1];
sx q[1];
rz(-1.1731802) q[1];
sx q[1];
rz(0.84487265) q[1];
rz(-0.028771632) q[3];
sx q[3];
rz(-1.1011916) q[3];
sx q[3];
rz(-2.4503277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3346682) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(2.7743288) q[2];
rz(-1.4372829) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(-1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.87840286) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(0.81018418) q[0];
rz(1.1148249) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(2.5883163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5779752) q[0];
sx q[0];
rz(-1.5153236) q[0];
sx q[0];
rz(2.2091921) q[0];
rz(-pi) q[1];
rz(-2.192334) q[2];
sx q[2];
rz(-0.87998968) q[2];
sx q[2];
rz(2.031428) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3685963) q[1];
sx q[1];
rz(-0.34075156) q[1];
sx q[1];
rz(-2.1590538) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3555967) q[3];
sx q[3];
rz(-1.82535) q[3];
sx q[3];
rz(2.8543775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1222003) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(0.23966399) q[2];
rz(1.3278809) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(-2.6740668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(-2.486034) q[1];
sx q[1];
rz(-1.1759023) q[1];
sx q[1];
rz(-2.1942153) q[1];
rz(-2.7262676) q[2];
sx q[2];
rz(-1.2653744) q[2];
sx q[2];
rz(1.1933586) q[2];
rz(0.58047337) q[3];
sx q[3];
rz(-1.3571285) q[3];
sx q[3];
rz(0.50783689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
