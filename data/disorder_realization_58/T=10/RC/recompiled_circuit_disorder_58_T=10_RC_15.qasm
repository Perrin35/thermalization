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
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18618628) q[0];
sx q[0];
rz(-0.49476981) q[0];
sx q[0];
rz(0.68899378) q[0];
x q[1];
rz(0.51672207) q[2];
sx q[2];
rz(-2.3331254) q[2];
sx q[2];
rz(2.7441623) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0349866) q[1];
sx q[1];
rz(-0.88284661) q[1];
sx q[1];
rz(0.41253849) q[1];
x q[2];
rz(-0.44736638) q[3];
sx q[3];
rz(-1.2036185) q[3];
sx q[3];
rz(0.51608738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(0.91180116) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(0.3266913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625898) q[0];
sx q[0];
rz(-1.0674745) q[0];
sx q[0];
rz(-0.83053629) q[0];
rz(-1.2454883) q[2];
sx q[2];
rz(-0.46519687) q[2];
sx q[2];
rz(-0.59782366) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5527199) q[1];
sx q[1];
rz(-2.1761804) q[1];
sx q[1];
rz(-0.49792109) q[1];
rz(-3.0558415) q[3];
sx q[3];
rz(-0.88480703) q[3];
sx q[3];
rz(-3.0119197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111423) q[0];
sx q[0];
rz(-1.550195) q[0];
sx q[0];
rz(1.9111454) q[0];
rz(-1.7705275) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(2.0470326) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87103802) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(-2.2689181) q[1];
x q[2];
rz(0.40773817) q[3];
sx q[3];
rz(-2.6860415) q[3];
sx q[3];
rz(-1.3726485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388158) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(1.7569655) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(0.93726678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0619229) q[1];
sx q[1];
rz(-1.8975782) q[1];
sx q[1];
rz(-2.5667739) q[1];
rz(-pi) q[2];
rz(-1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(-1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(0.35476312) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(2.9096471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7888768) q[0];
sx q[0];
rz(-1.0838325) q[0];
sx q[0];
rz(-0.48077521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4300214) q[2];
sx q[2];
rz(-0.88015926) q[2];
sx q[2];
rz(1.3123133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1272808) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(2.36253) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4076685) q[0];
sx q[0];
rz(-1.7051538) q[0];
sx q[0];
rz(2.2785447) q[0];
x q[1];
rz(0.41001292) q[2];
sx q[2];
rz(-0.7747246) q[2];
sx q[2];
rz(-2.0395525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.84038489) q[1];
sx q[1];
rz(-1.1494698) q[1];
sx q[1];
rz(0.2627443) q[1];
rz(0.44316767) q[3];
sx q[3];
rz(-2.245095) q[3];
sx q[3];
rz(0.5815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(2.732892) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(1.2020948) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(2.5278032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873338) q[0];
sx q[0];
rz(-1.7566534) q[0];
sx q[0];
rz(2.3128187) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28366144) q[2];
sx q[2];
rz(-1.6322871) q[2];
sx q[2];
rz(-0.21360699) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63749667) q[1];
sx q[1];
rz(-1.4667965) q[1];
sx q[1];
rz(-1.5118447) q[1];
rz(-pi) q[2];
rz(-3.0244163) q[3];
sx q[3];
rz(-2.0332608) q[3];
sx q[3];
rz(1.5069435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7229066) q[0];
sx q[0];
rz(-1.5333999) q[0];
sx q[0];
rz(2.5953369) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41495277) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(2.6331537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22673785) q[1];
sx q[1];
rz(-1.9361155) q[1];
sx q[1];
rz(1.564333) q[1];
rz(1.2582614) q[3];
sx q[3];
rz(-1.1141889) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(-0.92591441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6322964) q[0];
sx q[0];
rz(-0.25521454) q[0];
sx q[0];
rz(-2.1584828) q[0];
rz(2.3472896) q[2];
sx q[2];
rz(-2.1199193) q[2];
sx q[2];
rz(1.4505475) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33854252) q[1];
sx q[1];
rz(-1.1330789) q[1];
sx q[1];
rz(2.9796757) q[1];
rz(-2.9374398) q[3];
sx q[3];
rz(-1.4545868) q[3];
sx q[3];
rz(-1.0128563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(2.7465903) q[0];
rz(-1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(1.013247) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-0.14302111) q[0];
rz(-pi) q[1];
rz(-3.0229438) q[2];
sx q[2];
rz(-0.35602202) q[2];
sx q[2];
rz(0.3686541) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4229065) q[1];
sx q[1];
rz(-1.6675073) q[1];
sx q[1];
rz(2.7194517) q[1];
x q[2];
rz(0.89684422) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(-0.57327523) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-0.076889597) q[2];
sx q[2];
rz(-1.1270317) q[2];
sx q[2];
rz(1.8539853) q[2];
rz(1.6462973) q[3];
sx q[3];
rz(-2.4423238) q[3];
sx q[3];
rz(-2.8764976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
