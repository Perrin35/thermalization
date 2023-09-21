OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3842073) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(2.7469809) q[0];
x q[1];
rz(-2.0482424) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(2.8505461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7100916) q[1];
sx q[1];
rz(-2.3570865) q[1];
sx q[1];
rz(2.0246519) q[1];
rz(1.1675695) q[3];
sx q[3];
rz(-1.1551757) q[3];
sx q[3];
rz(0.88413873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43710199) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(2.2297915) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-0.3266913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90601901) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(-0.88615449) q[0];
rz(-pi) q[1];
rz(0.15906449) q[2];
sx q[2];
rz(-1.1317562) q[2];
sx q[2];
rz(0.95869267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8596526) q[1];
sx q[1];
rz(-1.167206) q[1];
sx q[1];
rz(-2.2380026) q[1];
rz(0.085751199) q[3];
sx q[3];
rz(-2.2567856) q[3];
sx q[3];
rz(3.0119197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111423) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(-1.2304473) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7705275) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(1.0945601) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2705546) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(0.87267455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7626761) q[3];
sx q[3];
rz(-1.1550316) q[3];
sx q[3];
rz(-2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.3624297) q[2];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-2.6065361) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(2.9761956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3446942) q[0];
sx q[0];
rz(-1.7555153) q[0];
sx q[0];
rz(3.0152507) q[0];
x q[1];
rz(-2.303316) q[2];
sx q[2];
rz(-0.24176134) q[2];
sx q[2];
rz(1.4571112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.71388984) q[1];
sx q[1];
rz(-2.1117003) q[1];
sx q[1];
rz(-1.9546024) q[1];
rz(-pi) q[2];
rz(-1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-2.9361434) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(2.9096471) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1629099) q[0];
sx q[0];
rz(-1.1497578) q[0];
sx q[0];
rz(2.1091503) q[0];
x q[1];
rz(-2.4460692) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(-2.7930789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8307454) q[1];
sx q[1];
rz(-1.7880882) q[1];
sx q[1];
rz(2.9160935) q[1];
rz(-2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-0.99745497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(-2.9655365) q[0];
rz(-pi) q[1];
rz(-1.9428271) q[2];
sx q[2];
rz(-2.2673006) q[2];
sx q[2];
rz(2.5859985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7194697) q[1];
sx q[1];
rz(-0.49233961) q[1];
sx q[1];
rz(-1.0455529) q[1];
rz(2.0632083) q[3];
sx q[3];
rz(-2.3541854) q[3];
sx q[3];
rz(1.9099964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(0.4513936) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873338) q[0];
sx q[0];
rz(-1.3849392) q[0];
sx q[0];
rz(-0.82877393) q[0];
rz(-pi) q[1];
rz(-1.5067528) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(-1.7664906) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1544513) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-0.51388545) q[1];
x q[2];
rz(1.1055787) q[3];
sx q[3];
rz(-1.4659766) q[3];
sx q[3];
rz(0.011381586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.9667352) q[2];
rz(0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229066) q[0];
sx q[0];
rz(-1.5333999) q[0];
sx q[0];
rz(0.54625578) q[0];
rz(-pi) q[1];
rz(-0.77163561) q[2];
sx q[2];
rz(-1.2727591) q[2];
sx q[2];
rz(2.3723797) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7998432) q[1];
sx q[1];
rz(-1.5647596) q[1];
sx q[1];
rz(-2.7762665) q[1];
rz(-1.8833313) q[3];
sx q[3];
rz(-1.1141889) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4349334) q[2];
sx q[2];
rz(-1.443053) q[2];
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
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-0.6859268) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5092963) q[0];
sx q[0];
rz(-2.8863781) q[0];
sx q[0];
rz(2.1584828) q[0];
x q[1];
rz(0.70897734) q[2];
sx q[2];
rz(-0.93009863) q[2];
sx q[2];
rz(0.35352732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1123062) q[1];
sx q[1];
rz(-0.46488133) q[1];
sx q[1];
rz(-1.2390562) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6191606) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(-2.0731376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9459076) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(0.39500239) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(1.013247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7618338) q[0];
sx q[0];
rz(-1.338431) q[0];
sx q[0];
rz(-0.14302111) q[0];
rz(-pi) q[1];
rz(-1.5268065) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(0.24214889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80876795) q[1];
sx q[1];
rz(-1.1507532) q[1];
sx q[1];
rz(-1.6767477) q[1];
rz(2.2447484) q[3];
sx q[3];
rz(-2.2138811) q[3];
sx q[3];
rz(-2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(0.56994462) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.7309932) q[2];
sx q[2];
rz(-2.6916531) q[2];
sx q[2];
rz(2.0315363) q[2];
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