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
rz(-1.9788111) q[0];
sx q[0];
rz(-2.2012308) q[0];
sx q[0];
rz(-0.19402394) q[0];
rz(2.5972875) q[1];
sx q[1];
rz(-1.5736009) q[1];
sx q[1];
rz(0.1967217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0531874) q[0];
sx q[0];
rz(-2.1413099) q[0];
sx q[0];
rz(-0.14630099) q[0];
rz(2.9241184) q[2];
sx q[2];
rz(-1.5494692) q[2];
sx q[2];
rz(2.3623717) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87952033) q[1];
sx q[1];
rz(-1.889887) q[1];
sx q[1];
rz(-1.8533005) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9858304) q[3];
sx q[3];
rz(-2.1418396) q[3];
sx q[3];
rz(-0.076249853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0804312) q[2];
sx q[2];
rz(-1.72074) q[2];
sx q[2];
rz(-3.1283992) q[2];
rz(-2.0773928) q[3];
sx q[3];
rz(-2.1049757) q[3];
sx q[3];
rz(-2.9707151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212293) q[0];
sx q[0];
rz(-2.5907488) q[0];
sx q[0];
rz(-0.4826104) q[0];
rz(-2.6699325) q[1];
sx q[1];
rz(-0.99716798) q[1];
sx q[1];
rz(0.68798033) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15287887) q[0];
sx q[0];
rz(-1.3804624) q[0];
sx q[0];
rz(2.8318263) q[0];
x q[1];
rz(-1.2763763) q[2];
sx q[2];
rz(-2.3465119) q[2];
sx q[2];
rz(-2.9816422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5045484) q[1];
sx q[1];
rz(-2.6682862) q[1];
sx q[1];
rz(0.096122336) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3673931) q[3];
sx q[3];
rz(-1.4869833) q[3];
sx q[3];
rz(2.0544685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.017435) q[2];
sx q[2];
rz(-2.0686801) q[2];
sx q[2];
rz(2.9166481) q[2];
rz(2.0624835) q[3];
sx q[3];
rz(-0.93828097) q[3];
sx q[3];
rz(2.6911531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17063046) q[0];
sx q[0];
rz(-0.58929515) q[0];
sx q[0];
rz(-1.2901837) q[0];
rz(-1.2782485) q[1];
sx q[1];
rz(-1.1561013) q[1];
sx q[1];
rz(1.2826756) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034866355) q[0];
sx q[0];
rz(-1.9283656) q[0];
sx q[0];
rz(0.94199543) q[0];
x q[1];
rz(3.01802) q[2];
sx q[2];
rz(-1.4915492) q[2];
sx q[2];
rz(1.0956956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8508529) q[1];
sx q[1];
rz(-0.026844414) q[1];
sx q[1];
rz(2.31183) q[1];
rz(-0.81961378) q[3];
sx q[3];
rz(-2.1791239) q[3];
sx q[3];
rz(-1.5386594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67905417) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(2.6666717) q[2];
rz(-1.9654407) q[3];
sx q[3];
rz(-2.2852496) q[3];
sx q[3];
rz(-0.28158751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3731641) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(0.1678634) q[0];
rz(-0.31790512) q[1];
sx q[1];
rz(-1.8524421) q[1];
sx q[1];
rz(-2.1071404) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4108334) q[0];
sx q[0];
rz(-2.877088) q[0];
sx q[0];
rz(1.4368617) q[0];
x q[1];
rz(1.9475161) q[2];
sx q[2];
rz(-2.2206306) q[2];
sx q[2];
rz(-1.5776933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75986162) q[1];
sx q[1];
rz(-2.029756) q[1];
sx q[1];
rz(-0.81736395) q[1];
rz(-pi) q[2];
rz(-1.4698704) q[3];
sx q[3];
rz(-0.96138326) q[3];
sx q[3];
rz(2.963629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1075403) q[2];
sx q[2];
rz(-2.5120021) q[2];
sx q[2];
rz(-0.26629392) q[2];
rz(-0.1693503) q[3];
sx q[3];
rz(-1.063187) q[3];
sx q[3];
rz(-0.17899409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.2696445) q[0];
sx q[0];
rz(-2.9565411) q[0];
sx q[0];
rz(1.8484775) q[0];
rz(-3.0710908) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(-0.806113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7370395) q[0];
sx q[0];
rz(-1.4463639) q[0];
sx q[0];
rz(-1.3621773) q[0];
x q[1];
rz(1.2048804) q[2];
sx q[2];
rz(-2.9822449) q[2];
sx q[2];
rz(-1.6312576) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15547968) q[1];
sx q[1];
rz(-2.1633185) q[1];
sx q[1];
rz(1.3377124) q[1];
rz(-pi) q[2];
rz(2.6025196) q[3];
sx q[3];
rz(-1.7601437) q[3];
sx q[3];
rz(-2.2089778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21141323) q[2];
sx q[2];
rz(-1.6665919) q[2];
sx q[2];
rz(0.30964568) q[2];
rz(-0.51112255) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(2.3206594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822134) q[0];
sx q[0];
rz(-0.62726227) q[0];
sx q[0];
rz(-2.7045767) q[0];
rz(2.6868195) q[1];
sx q[1];
rz(-0.79650703) q[1];
sx q[1];
rz(0.88266596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32188181) q[0];
sx q[0];
rz(-1.30153) q[0];
sx q[0];
rz(0.092094941) q[0];
rz(1.9493616) q[2];
sx q[2];
rz(-1.6003216) q[2];
sx q[2];
rz(1.5020509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5957634) q[1];
sx q[1];
rz(-0.59844164) q[1];
sx q[1];
rz(-0.95495895) q[1];
rz(2.1737133) q[3];
sx q[3];
rz(-1.9390188) q[3];
sx q[3];
rz(-0.24711025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5060045) q[2];
sx q[2];
rz(-1.4281861) q[2];
sx q[2];
rz(2.7006855) q[2];
rz(2.0697557) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(-2.5409839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275682) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(-1.2822275) q[0];
rz(2.0558689) q[1];
sx q[1];
rz(-1.6174569) q[1];
sx q[1];
rz(-1.6283183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075628932) q[0];
sx q[0];
rz(-2.3018078) q[0];
sx q[0];
rz(-0.15878598) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37570132) q[2];
sx q[2];
rz(-0.73601572) q[2];
sx q[2];
rz(-1.0003566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3320184) q[1];
sx q[1];
rz(-0.93832131) q[1];
sx q[1];
rz(-0.47398403) q[1];
rz(0.75700883) q[3];
sx q[3];
rz(-1.2881713) q[3];
sx q[3];
rz(2.4429537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47161237) q[2];
sx q[2];
rz(-0.70576224) q[2];
sx q[2];
rz(0.82028779) q[2];
rz(-2.7465076) q[3];
sx q[3];
rz(-2.3494651) q[3];
sx q[3];
rz(2.5274966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87913269) q[0];
sx q[0];
rz(-1.4985871) q[0];
sx q[0];
rz(-0.60687989) q[0];
rz(0.34582368) q[1];
sx q[1];
rz(-1.9551829) q[1];
sx q[1];
rz(0.80723673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8223752) q[0];
sx q[0];
rz(-1.2273047) q[0];
sx q[0];
rz(-0.10256711) q[0];
rz(-pi) q[1];
x q[1];
rz(1.973602) q[2];
sx q[2];
rz(-1.8181846) q[2];
sx q[2];
rz(1.1198662) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9620558) q[1];
sx q[1];
rz(-2.7946094) q[1];
sx q[1];
rz(-2.4605303) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5446078) q[3];
sx q[3];
rz(-1.0506949) q[3];
sx q[3];
rz(2.6900684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.81749934) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(0.49638805) q[2];
rz(3.056622) q[3];
sx q[3];
rz(-2.7111263) q[3];
sx q[3];
rz(2.5607204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896511) q[0];
sx q[0];
rz(-1.713151) q[0];
sx q[0];
rz(0.38452837) q[0];
rz(-2.8893068) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(-1.5834454) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2551832) q[0];
sx q[0];
rz(-0.46570581) q[0];
sx q[0];
rz(-0.021115594) q[0];
rz(2.8674312) q[2];
sx q[2];
rz(-1.6841075) q[2];
sx q[2];
rz(-0.99565143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6114735) q[1];
sx q[1];
rz(-0.58936469) q[1];
sx q[1];
rz(1.5477647) q[1];
rz(-pi) q[2];
rz(-2.4676314) q[3];
sx q[3];
rz(-2.1836851) q[3];
sx q[3];
rz(-2.6496668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(1.1400878) q[2];
rz(1.358486) q[3];
sx q[3];
rz(-1.8890817) q[3];
sx q[3];
rz(-2.8235249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7110905) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(-2.7464113) q[0];
rz(-1.2218062) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(2.2538259) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632653) q[0];
sx q[0];
rz(-0.76416956) q[0];
sx q[0];
rz(-0.40212888) q[0];
rz(1.5776921) q[2];
sx q[2];
rz(-2.6665932) q[2];
sx q[2];
rz(-1.4003225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5612167) q[1];
sx q[1];
rz(-1.9264915) q[1];
sx q[1];
rz(3.1144193) q[1];
x q[2];
rz(1.2492845) q[3];
sx q[3];
rz(-0.32673937) q[3];
sx q[3];
rz(-0.48130166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.44783529) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(-2.5698404) q[2];
rz(1.5104431) q[3];
sx q[3];
rz(-1.7084728) q[3];
sx q[3];
rz(1.3408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0157264) q[0];
sx q[0];
rz(-1.8342352) q[0];
sx q[0];
rz(-2.9974708) q[0];
rz(-1.9563328) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(-2.588829) q[2];
sx q[2];
rz(-2.6949099) q[2];
sx q[2];
rz(0.75134122) q[2];
rz(1.8923459) q[3];
sx q[3];
rz(-0.95148181) q[3];
sx q[3];
rz(-0.8368809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
