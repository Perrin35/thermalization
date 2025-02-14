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
rz(2.4616315) q[0];
sx q[0];
rz(7.1890561) q[0];
sx q[0];
rz(11.472975) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(1.9688164) q[1];
sx q[1];
rz(9.8493675) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6494006) q[0];
sx q[0];
rz(-1.7723506) q[0];
sx q[0];
rz(-3.1026476) q[0];
rz(-pi) q[1];
rz(-1.8867887) q[2];
sx q[2];
rz(-2.9796763) q[2];
sx q[2];
rz(-2.201061) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72480768) q[1];
sx q[1];
rz(-2.8832286) q[1];
sx q[1];
rz(-1.2577673) q[1];
rz(2.2810117) q[3];
sx q[3];
rz(-1.4690377) q[3];
sx q[3];
rz(-1.4836131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37497416) q[2];
sx q[2];
rz(-1.5593636) q[2];
sx q[2];
rz(-1.9826822) q[2];
rz(-2.9186987) q[3];
sx q[3];
rz(-1.736172) q[3];
sx q[3];
rz(-0.29002732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952154) q[0];
sx q[0];
rz(-1.5257436) q[0];
sx q[0];
rz(-1.079153) q[0];
rz(-1.7201299) q[1];
sx q[1];
rz(-0.68669569) q[1];
sx q[1];
rz(-0.064528331) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20708974) q[0];
sx q[0];
rz(-0.10182589) q[0];
sx q[0];
rz(-2.1125211) q[0];
x q[1];
rz(2.6653397) q[2];
sx q[2];
rz(-1.5956399) q[2];
sx q[2];
rz(2.3551539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5777305) q[1];
sx q[1];
rz(-0.81474333) q[1];
sx q[1];
rz(1.9063063) q[1];
rz(-1.7713624) q[3];
sx q[3];
rz(-0.75955694) q[3];
sx q[3];
rz(-2.650039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9885538) q[2];
sx q[2];
rz(-1.7861853) q[2];
sx q[2];
rz(1.1174196) q[2];
rz(-0.85876632) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(2.7045238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7607464) q[0];
sx q[0];
rz(-2.8570211) q[0];
sx q[0];
rz(2.9685156) q[0];
rz(0.95678798) q[1];
sx q[1];
rz(-2.1134977) q[1];
sx q[1];
rz(2.079336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1011216) q[0];
sx q[0];
rz(-0.72325828) q[0];
sx q[0];
rz(-1.5612119) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4281279) q[2];
sx q[2];
rz(-0.77859813) q[2];
sx q[2];
rz(2.8215591) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77327496) q[1];
sx q[1];
rz(-2.2466806) q[1];
sx q[1];
rz(2.2818517) q[1];
x q[2];
rz(0.84433664) q[3];
sx q[3];
rz(-0.82714836) q[3];
sx q[3];
rz(-1.3177804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.411285) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(0.74472204) q[2];
rz(2.5005285) q[3];
sx q[3];
rz(-1.3499667) q[3];
sx q[3];
rz(0.49066576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915801) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(-2.2963754) q[0];
rz(-2.7382964) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(0.32803112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16763891) q[0];
sx q[0];
rz(-2.9120647) q[0];
sx q[0];
rz(-2.2295802) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9035527) q[2];
sx q[2];
rz(-1.5063707) q[2];
sx q[2];
rz(0.062459613) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8106336) q[1];
sx q[1];
rz(-1.0077268) q[1];
sx q[1];
rz(1.710239) q[1];
rz(-2.115504) q[3];
sx q[3];
rz(-1.8063365) q[3];
sx q[3];
rz(2.1229471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8959117) q[2];
sx q[2];
rz(-0.74378219) q[2];
sx q[2];
rz(-2.9823859) q[2];
rz(-3.1253452) q[3];
sx q[3];
rz(-1.0280321) q[3];
sx q[3];
rz(0.73133674) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943587) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(2.0514945) q[0];
rz(3.0116426) q[1];
sx q[1];
rz(-2.3268685) q[1];
sx q[1];
rz(1.5257588) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65924683) q[0];
sx q[0];
rz(-2.1123943) q[0];
sx q[0];
rz(2.0603176) q[0];
x q[1];
rz(-1.9984869) q[2];
sx q[2];
rz(-0.18015716) q[2];
sx q[2];
rz(-3.0017972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.255521) q[1];
sx q[1];
rz(-0.93734765) q[1];
sx q[1];
rz(-1.6048163) q[1];
x q[2];
rz(-0.035338621) q[3];
sx q[3];
rz(-1.3900847) q[3];
sx q[3];
rz(0.23739761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7502363) q[2];
sx q[2];
rz(-0.82166925) q[2];
sx q[2];
rz(2.6643122) q[2];
rz(0.20398772) q[3];
sx q[3];
rz(-2.9450649) q[3];
sx q[3];
rz(2.5175214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1962965) q[0];
sx q[0];
rz(-0.81552234) q[0];
sx q[0];
rz(0.77769172) q[0];
rz(2.5241191) q[1];
sx q[1];
rz(-1.4910881) q[1];
sx q[1];
rz(2.2106574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5701911) q[0];
sx q[0];
rz(-0.63136111) q[0];
sx q[0];
rz(-2.3124113) q[0];
x q[1];
rz(0.41895509) q[2];
sx q[2];
rz(-1.1280931) q[2];
sx q[2];
rz(2.1328164) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8691784) q[1];
sx q[1];
rz(-0.34124103) q[1];
sx q[1];
rz(-1.7029087) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0369426) q[3];
sx q[3];
rz(-0.6932689) q[3];
sx q[3];
rz(1.1806837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65333873) q[2];
sx q[2];
rz(-1.3335756) q[2];
sx q[2];
rz(-2.874157) q[2];
rz(-2.2087162) q[3];
sx q[3];
rz(-1.2837912) q[3];
sx q[3];
rz(2.647184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6701508) q[0];
sx q[0];
rz(-1.9981367) q[0];
sx q[0];
rz(2.4757521) q[0];
rz(-2.937607) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(-1.9542255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0788308) q[0];
sx q[0];
rz(-1.4785066) q[0];
sx q[0];
rz(-0.4706443) q[0];
rz(-pi) q[1];
rz(0.97839956) q[2];
sx q[2];
rz(-0.80139388) q[2];
sx q[2];
rz(-1.8155542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0625035) q[1];
sx q[1];
rz(-2.3983208) q[1];
sx q[1];
rz(0.39689245) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8245876) q[3];
sx q[3];
rz(-1.6204066) q[3];
sx q[3];
rz(1.3146159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5953956) q[2];
sx q[2];
rz(-2.6232145) q[2];
sx q[2];
rz(-2.4714244) q[2];
rz(-2.8403122) q[3];
sx q[3];
rz(-1.2944841) q[3];
sx q[3];
rz(-2.7231351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30348521) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(1.0656892) q[0];
rz(-2.4844555) q[1];
sx q[1];
rz(-0.70825759) q[1];
sx q[1];
rz(-1.7117975) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2524714) q[0];
sx q[0];
rz(-1.234113) q[0];
sx q[0];
rz(0.97674979) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33139511) q[2];
sx q[2];
rz(-1.3989395) q[2];
sx q[2];
rz(2.9239025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.011466) q[1];
sx q[1];
rz(-1.7892171) q[1];
sx q[1];
rz(0.41832025) q[1];
rz(-pi) q[2];
rz(-0.18971209) q[3];
sx q[3];
rz(-1.7125152) q[3];
sx q[3];
rz(-0.82822463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26057217) q[2];
sx q[2];
rz(-1.842097) q[2];
sx q[2];
rz(2.1232429) q[2];
rz(-1.5923422) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(2.3840267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26579478) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(2.1739668) q[0];
rz(-0.34010092) q[1];
sx q[1];
rz(-0.83931559) q[1];
sx q[1];
rz(-2.1077154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0268896) q[0];
sx q[0];
rz(-1.7091161) q[0];
sx q[0];
rz(-1.0113495) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56615717) q[2];
sx q[2];
rz(-0.70933178) q[2];
sx q[2];
rz(-0.69401238) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47017322) q[1];
sx q[1];
rz(-0.88693383) q[1];
sx q[1];
rz(-2.629423) q[1];
x q[2];
rz(2.9176209) q[3];
sx q[3];
rz(-1.4272318) q[3];
sx q[3];
rz(-0.75907133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0391482) q[2];
sx q[2];
rz(-1.4463964) q[2];
sx q[2];
rz(0.038912494) q[2];
rz(0.14032042) q[3];
sx q[3];
rz(-0.084429927) q[3];
sx q[3];
rz(1.8261568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4843531) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(1.8835618) q[0];
rz(-2.925442) q[1];
sx q[1];
rz(-0.70544568) q[1];
sx q[1];
rz(-0.36453882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.011019) q[0];
sx q[0];
rz(-2.6219061) q[0];
sx q[0];
rz(2.0849865) q[0];
x q[1];
rz(-2.6993518) q[2];
sx q[2];
rz(-0.85390515) q[2];
sx q[2];
rz(0.88767641) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8382524) q[1];
sx q[1];
rz(-1.6263744) q[1];
sx q[1];
rz(-1.5838632) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50915368) q[3];
sx q[3];
rz(-3.0753071) q[3];
sx q[3];
rz(1.2627511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2898499) q[2];
sx q[2];
rz(-1.936329) q[2];
sx q[2];
rz(2.5002948) q[2];
rz(1.4801721) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(0.67354584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73364532) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(-0.98450487) q[1];
sx q[1];
rz(-1.4823109) q[1];
sx q[1];
rz(-1.0479814) q[1];
rz(2.099809) q[2];
sx q[2];
rz(-0.22780252) q[2];
sx q[2];
rz(-3.0115423) q[2];
rz(1.0897286) q[3];
sx q[3];
rz(-1.0952598) q[3];
sx q[3];
rz(2.918603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
