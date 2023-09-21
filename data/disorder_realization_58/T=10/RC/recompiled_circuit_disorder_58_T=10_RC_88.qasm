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
rz(2.7352754) q[0];
sx q[0];
rz(8.6046594) q[0];
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
rz(2.3842073) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(0.39461179) q[0];
rz(-pi) q[1];
rz(-2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(-0.29104656) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9484529) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(2.302041) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44736638) q[3];
sx q[3];
rz(-1.2036185) q[3];
sx q[3];
rz(-0.51608738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-2.2226298) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(0.3266913) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2355736) q[0];
sx q[0];
rz(-2.2740907) q[0];
sx q[0];
rz(2.2554382) q[0];
x q[1];
rz(-0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(-0.95869267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28194004) q[1];
sx q[1];
rz(-1.167206) q[1];
sx q[1];
rz(-0.90359009) q[1];
rz(-1.4665524) q[3];
sx q[3];
rz(-0.6904656) q[3];
sx q[3];
rz(-3.1363827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
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
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-2.0842016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111423) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(1.2304473) q[0];
rz(-pi) q[1];
rz(2.8003642) q[2];
sx q[2];
rz(-1.7592906) q[2];
sx q[2];
rz(0.54268062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2705546) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(-0.87267455) q[1];
x q[2];
rz(-1.3789165) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-2.2172745) q[3];
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
rz(1.3624297) q[2];
rz(-2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(0.16539703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388158) q[0];
sx q[0];
rz(-1.6949777) q[0];
sx q[0];
rz(1.7569655) q[0];
rz(1.7521162) q[2];
sx q[2];
rz(-1.4099979) q[2];
sx q[2];
rz(-0.60418512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0619229) q[1];
sx q[1];
rz(-1.8975782) q[1];
sx q[1];
rz(2.5667739) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27407077) q[3];
sx q[3];
rz(-2.1087397) q[3];
sx q[3];
rz(-2.2171973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(2.9361434) q[2];
rz(-2.0139587) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(1.9873437) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35271586) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(-0.48077521) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69552341) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(0.34851375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1272808) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(0.77906268) q[1];
x q[2];
rz(1.5512974) q[3];
sx q[3];
rz(-1.9392506) q[3];
sx q[3];
rz(1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(1.3034472) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(-2.1525106) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(2.1441377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31842445) q[0];
sx q[0];
rz(-0.71821763) q[0];
sx q[0];
rz(-1.7757925) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9428271) q[2];
sx q[2];
rz(-0.87429201) q[2];
sx q[2];
rz(0.55559413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3012078) q[1];
sx q[1];
rz(-1.1494698) q[1];
sx q[1];
rz(2.8788484) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84647471) q[3];
sx q[3];
rz(-1.2292362) q[3];
sx q[3];
rz(0.70111707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(-0.4513936) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-0.62414449) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-0.61378941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44871556) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(0.24976272) q[0];
x q[1];
rz(-2.9250547) q[2];
sx q[2];
rz(-2.8515186) q[2];
sx q[2];
rz(-1.9921583) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92717273) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(-3.0374132) q[1];
rz(-pi) q[2];
rz(1.3404487) q[3];
sx q[3];
rz(-0.47603546) q[3];
sx q[3];
rz(1.3766833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.9667352) q[2];
rz(0.60837778) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293837) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(-1.5270385) q[0];
x q[1];
rz(0.77163561) q[2];
sx q[2];
rz(-1.8688335) q[2];
sx q[2];
rz(2.3723797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22673785) q[1];
sx q[1];
rz(-1.2054772) q[1];
sx q[1];
rz(-1.5772596) q[1];
x q[2];
rz(-1.2582614) q[3];
sx q[3];
rz(-2.0274037) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-0.92591441) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6322964) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(-0.98310982) q[0];
rz(-0.7943031) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(1.6910451) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1123062) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(-1.2390562) q[1];
rz(-1.4521452) q[3];
sx q[3];
rz(-1.7735529) q[3];
sx q[3];
rz(-2.6076536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(0.53608981) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(2.1283456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625044) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(2.1129235) q[0];
rz(0.35372325) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(-1.8281787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80876795) q[1];
sx q[1];
rz(-1.9908394) q[1];
sx q[1];
rz(-1.4648449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4471531) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(-1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(-1.7761207) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(-0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41291819) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(0.076889597) q[2];
sx q[2];
rz(-2.0145609) q[2];
sx q[2];
rz(-1.2876074) q[2];
rz(0.063354062) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];