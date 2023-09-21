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
rz(0.83067218) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75738534) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(-0.39461179) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0933502) q[2];
sx q[2];
rz(-2.2508143) q[2];
sx q[2];
rz(-2.8505461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9484529) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(-0.83955168) q[1];
x q[2];
rz(1.1675695) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(0.91180116) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(2.8149014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.806843) q[0];
sx q[0];
rz(-0.93870367) q[0];
sx q[0];
rz(-2.500781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(0.95869267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82636278) q[1];
sx q[1];
rz(-2.3781207) q[1];
sx q[1];
rz(-2.1748494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88300206) q[3];
sx q[3];
rz(-1.6371173) q[3];
sx q[3];
rz(1.4955213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(2.8116995) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(-2.0842016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49845093) q[0];
sx q[0];
rz(-2.8006449) q[0];
sx q[0];
rz(1.6324415) q[0];
x q[1];
rz(0.34122841) q[2];
sx q[2];
rz(-1.7592906) q[2];
sx q[2];
rz(-0.54268062) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96506572) q[1];
sx q[1];
rz(-2.2245193) q[1];
sx q[1];
rz(-2.7214126) q[1];
rz(-pi) q[2];
rz(-1.3789165) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3588336) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20277682) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(-1.3846272) q[0];
rz(-pi) q[1];
rz(-2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(1.6844815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0619229) q[1];
sx q[1];
rz(-1.8975782) q[1];
sx q[1];
rz(-0.57481874) q[1];
rz(-pi) q[2];
rz(-2.1256251) q[3];
sx q[3];
rz(-1.3362243) q[3];
sx q[3];
rz(0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7455204) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(0.23194557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35271586) q[0];
sx q[0];
rz(-1.0838325) q[0];
sx q[0];
rz(0.48077521) q[0];
rz(0.69552341) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-0.34851375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0143118) q[1];
sx q[1];
rz(-2.8297272) q[1];
sx q[1];
rz(-2.36253) q[1];
rz(-1.5902953) q[3];
sx q[3];
rz(-1.9392506) q[3];
sx q[3];
rz(1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0971138) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(2.1441377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0925804) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(2.9655365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73156725) q[2];
sx q[2];
rz(-1.8533857) q[2];
sx q[2];
rz(2.3716795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62085405) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(-1.1362856) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2951179) q[3];
sx q[3];
rz(-1.9123565) q[3];
sx q[3];
rz(2.4404756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(0.61378941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6928771) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(0.24976272) q[0];
rz(-1.6348398) q[2];
sx q[2];
rz(-1.8539068) q[2];
sx q[2];
rz(-1.7664906) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.504096) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(0.4883782) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293837) q[0];
sx q[0];
rz(-1.0249656) q[0];
sx q[0];
rz(-1.6145541) q[0];
rz(2.7266399) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(0.50843898) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8967646) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(3.1246964) q[1];
rz(-pi) q[2];
rz(-1.2582614) q[3];
sx q[3];
rz(-2.0274037) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-0.92591441) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5110916) q[0];
sx q[0];
rz(-1.4303659) q[0];
sx q[0];
rz(1.7846084) q[0];
x q[1];
rz(-2.2885867) q[2];
sx q[2];
rz(-2.2249613) q[2];
sx q[2];
rz(0.6086364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.029286413) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(-1.9025365) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9374398) q[3];
sx q[3];
rz(-1.4545868) q[3];
sx q[3];
rz(1.0128563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-1.013247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2241867) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(-1.805472) q[0];
rz(-pi) q[1];
rz(-1.5268065) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(0.24214889) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71868616) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(-0.42214091) q[1];
x q[2];
rz(2.2447484) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-2.3821793) q[2];
rz(-1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(-2.5683174) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
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