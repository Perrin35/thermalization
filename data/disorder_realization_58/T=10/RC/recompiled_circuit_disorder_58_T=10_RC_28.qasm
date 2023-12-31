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
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(pi/2) q[2];
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
rz(-0.51672207) q[2];
sx q[2];
rz(-0.80846723) q[2];
sx q[2];
rz(2.7441623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1066061) q[1];
sx q[1];
rz(-0.88284661) q[1];
sx q[1];
rz(2.7290542) q[1];
x q[2];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-0.1201771) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(0.78951019) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90601901) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(2.2554382) q[0];
x q[1];
rz(0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(-2.1829) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28194004) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(2.2380026) q[1];
rz(0.085751199) q[3];
sx q[3];
rz(-2.2567856) q[3];
sx q[3];
rz(3.0119197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(2.0842016) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1304504) q[0];
sx q[0];
rz(-1.550195) q[0];
sx q[0];
rz(-1.9111454) q[0];
x q[1];
rz(-2.8003642) q[2];
sx q[2];
rz(-1.7592906) q[2];
sx q[2];
rz(-0.54268062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8101465) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(-2.0600832) q[1];
rz(-pi) q[2];
rz(1.7626761) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3588336) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(1.3624297) q[2];
rz(-2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-2.9761956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1918614) q[0];
sx q[0];
rz(-0.22338578) q[0];
sx q[0];
rz(-0.97747691) q[0];
rz(-pi) q[1];
rz(-2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(0.93726678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0488102) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(0.55744967) q[1];
x q[2];
rz(0.27407077) q[3];
sx q[3];
rz(-1.032853) q[3];
sx q[3];
rz(-0.9243954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(-1.127634) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.9873437) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97868279) q[0];
sx q[0];
rz(-1.9918348) q[0];
sx q[0];
rz(-1.0324423) q[0];
rz(-pi) q[1];
rz(1.7115713) q[2];
sx q[2];
rz(-2.2614334) q[2];
sx q[2];
rz(1.8292793) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9310589) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(1.7935497) q[1];
x q[2];
rz(-0.050458126) q[3];
sx q[3];
rz(-0.36894635) q[3];
sx q[3];
rz(1.7794533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.8381455) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(2.1441377) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8231682) q[0];
sx q[0];
rz(-2.423375) q[0];
sx q[0];
rz(1.7757925) q[0];
rz(-pi) q[1];
rz(-1.1987655) q[2];
sx q[2];
rz(-2.2673006) q[2];
sx q[2];
rz(0.55559413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4221229) q[1];
sx q[1];
rz(-0.49233961) q[1];
sx q[1];
rz(-2.0960397) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(-1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9138907) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873338) q[0];
sx q[0];
rz(-1.3849392) q[0];
sx q[0];
rz(0.82877393) q[0];
rz(-0.21653793) q[2];
sx q[2];
rz(-2.8515186) q[2];
sx q[2];
rz(-1.1494344) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92717273) q[1];
sx q[1];
rz(-1.5121636) q[1];
sx q[1];
rz(-0.10417948) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11717637) q[3];
sx q[3];
rz(-2.0332608) q[3];
sx q[3];
rz(1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(-2.5332149) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.2414918) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2135293) q[0];
sx q[0];
rz(-2.5941879) q[0];
sx q[0];
rz(-0.071896032) q[0];
rz(-2.369957) q[2];
sx q[2];
rz(-1.8688335) q[2];
sx q[2];
rz(2.3723797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3417495) q[1];
sx q[1];
rz(-1.5647596) q[1];
sx q[1];
rz(0.36532613) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8833313) q[3];
sx q[3];
rz(-1.1141889) q[3];
sx q[3];
rz(1.4708335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6305011) q[0];
sx q[0];
rz(-1.4303659) q[0];
sx q[0];
rz(1.3569843) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8530059) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(2.5329563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1123062) q[1];
sx q[1];
rz(-0.46488133) q[1];
sx q[1];
rz(1.9025365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6894475) q[3];
sx q[3];
rz(-1.3680397) q[3];
sx q[3];
rz(2.6076536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(0.53608981) q[2];
rz(-0.4195956) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855899) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-1.013247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625044) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(1.0286691) q[0];
x q[1];
rz(-0.35372325) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(-1.8281787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4229065) q[1];
sx q[1];
rz(-1.6675073) q[1];
sx q[1];
rz(-2.7194517) q[1];
x q[2];
rz(-2.3771044) q[3];
sx q[3];
rz(-2.0937243) q[3];
sx q[3];
rz(-2.5564699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(-1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(-3.0782386) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
