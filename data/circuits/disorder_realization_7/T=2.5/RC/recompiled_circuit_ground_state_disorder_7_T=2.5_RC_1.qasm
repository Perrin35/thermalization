OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11243842) q[0];
sx q[0];
rz(4.772679) q[0];
sx q[0];
rz(9.7909238) q[0];
rz(2.0868299) q[1];
sx q[1];
rz(-1.8572448) q[1];
sx q[1];
rz(-0.73959124) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407362) q[0];
sx q[0];
rz(-1.9240409) q[0];
sx q[0];
rz(2.059444) q[0];
rz(-pi) q[1];
rz(2.5049091) q[2];
sx q[2];
rz(-1.0708059) q[2];
sx q[2];
rz(-1.9477109) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89902564) q[1];
sx q[1];
rz(-2.0787313) q[1];
sx q[1];
rz(-2.6655156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69460709) q[3];
sx q[3];
rz(-0.41037729) q[3];
sx q[3];
rz(0.17363901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4313878) q[2];
sx q[2];
rz(-1.4375765) q[2];
sx q[2];
rz(2.0598742) q[2];
rz(1.9234575) q[3];
sx q[3];
rz(-1.5798605) q[3];
sx q[3];
rz(1.6703828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0064938) q[0];
sx q[0];
rz(-0.83469892) q[0];
sx q[0];
rz(0.83719069) q[0];
rz(-0.93765014) q[1];
sx q[1];
rz(-1.7559914) q[1];
sx q[1];
rz(0.12486501) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9287024) q[0];
sx q[0];
rz(-0.68395185) q[0];
sx q[0];
rz(0.70080832) q[0];
x q[1];
rz(1.2697095) q[2];
sx q[2];
rz(-1.0137179) q[2];
sx q[2];
rz(0.89424342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0016922) q[1];
sx q[1];
rz(-0.6065953) q[1];
sx q[1];
rz(2.9942375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6065793) q[3];
sx q[3];
rz(-1.8763308) q[3];
sx q[3];
rz(1.7747674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79895926) q[2];
sx q[2];
rz(-1.066076) q[2];
sx q[2];
rz(0.40668818) q[2];
rz(2.0090571) q[3];
sx q[3];
rz(-1.5187902) q[3];
sx q[3];
rz(-0.86743814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7750074) q[0];
sx q[0];
rz(-0.58554119) q[0];
sx q[0];
rz(-1.0803692) q[0];
rz(0.22251546) q[1];
sx q[1];
rz(-2.3718926) q[1];
sx q[1];
rz(-0.57058191) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208698) q[0];
sx q[0];
rz(-1.1900702) q[0];
sx q[0];
rz(-2.8080567) q[0];
x q[1];
rz(-2.5466444) q[2];
sx q[2];
rz(-2.1052268) q[2];
sx q[2];
rz(-2.7723412) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1751082) q[1];
sx q[1];
rz(-2.3382332) q[1];
sx q[1];
rz(3.1109555) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5665652) q[3];
sx q[3];
rz(-1.417727) q[3];
sx q[3];
rz(2.0401772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83501619) q[2];
sx q[2];
rz(-1.4811652) q[2];
sx q[2];
rz(1.906685) q[2];
rz(-2.0645449) q[3];
sx q[3];
rz(-1.5826694) q[3];
sx q[3];
rz(2.7243015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1398337) q[0];
sx q[0];
rz(-1.3707021) q[0];
sx q[0];
rz(2.4613001) q[0];
rz(1.744005) q[1];
sx q[1];
rz(-2.0060507) q[1];
sx q[1];
rz(-0.23522338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5176217) q[0];
sx q[0];
rz(-0.67254449) q[0];
sx q[0];
rz(0.60136232) q[0];
x q[1];
rz(-0.3346457) q[2];
sx q[2];
rz(-2.2586933) q[2];
sx q[2];
rz(0.84977023) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4724303) q[1];
sx q[1];
rz(-2.1390657) q[1];
sx q[1];
rz(-1.8316339) q[1];
rz(-2.8743478) q[3];
sx q[3];
rz(-1.4660204) q[3];
sx q[3];
rz(-1.9418093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.816421) q[2];
sx q[2];
rz(-1.07594) q[2];
sx q[2];
rz(1.761033) q[2];
rz(2.374968) q[3];
sx q[3];
rz(-2.4216757) q[3];
sx q[3];
rz(-0.64002526) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0945053) q[0];
sx q[0];
rz(-2.6386059) q[0];
sx q[0];
rz(-1.267953) q[0];
rz(-2.2044334) q[1];
sx q[1];
rz(-2.2068534) q[1];
sx q[1];
rz(-2.8840205) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1510096) q[0];
sx q[0];
rz(-2.5611612) q[0];
sx q[0];
rz(1.849017) q[0];
rz(2.1721971) q[2];
sx q[2];
rz(-1.3139408) q[2];
sx q[2];
rz(2.1208626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3993381) q[1];
sx q[1];
rz(-0.4961001) q[1];
sx q[1];
rz(-2.6861486) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0004992) q[3];
sx q[3];
rz(-0.52489108) q[3];
sx q[3];
rz(-2.0980245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.81689721) q[2];
sx q[2];
rz(-0.92528737) q[2];
sx q[2];
rz(-0.84358215) q[2];
rz(-2.2364565) q[3];
sx q[3];
rz(-0.89591187) q[3];
sx q[3];
rz(-0.26346537) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0819241) q[0];
sx q[0];
rz(-2.7984239) q[0];
sx q[0];
rz(-1.7053509) q[0];
rz(1.2417271) q[1];
sx q[1];
rz(-1.714434) q[1];
sx q[1];
rz(2.5232975) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.623659) q[0];
sx q[0];
rz(-2.111064) q[0];
sx q[0];
rz(-2.1884657) q[0];
rz(-pi) q[1];
rz(-1.773387) q[2];
sx q[2];
rz(-0.69249047) q[2];
sx q[2];
rz(2.4634374) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3028812) q[1];
sx q[1];
rz(-0.5190604) q[1];
sx q[1];
rz(-2.5005591) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6446436) q[3];
sx q[3];
rz(-0.54936159) q[3];
sx q[3];
rz(3.0086111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5228086) q[2];
sx q[2];
rz(-1.2455384) q[2];
sx q[2];
rz(1.1946542) q[2];
rz(1.4905802) q[3];
sx q[3];
rz(-1.4158764) q[3];
sx q[3];
rz(-1.2913316) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8959344) q[0];
sx q[0];
rz(-0.36933649) q[0];
sx q[0];
rz(-2.9781407) q[0];
rz(-2.5864736) q[1];
sx q[1];
rz(-1.6682383) q[1];
sx q[1];
rz(2.3177573) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4651983) q[0];
sx q[0];
rz(-1.4270743) q[0];
sx q[0];
rz(-2.8222047) q[0];
rz(-pi) q[1];
rz(-1.3029609) q[2];
sx q[2];
rz(-2.3402956) q[2];
sx q[2];
rz(-1.4943701) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8870315) q[1];
sx q[1];
rz(-1.7486835) q[1];
sx q[1];
rz(2.6807941) q[1];
rz(-pi) q[2];
rz(-1.9013405) q[3];
sx q[3];
rz(-1.179266) q[3];
sx q[3];
rz(-0.82465224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0865563) q[2];
sx q[2];
rz(-1.8359416) q[2];
sx q[2];
rz(-1.4646863) q[2];
rz(-3.0875409) q[3];
sx q[3];
rz(-0.87658221) q[3];
sx q[3];
rz(2.7189253) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5257877) q[0];
sx q[0];
rz(-2.7532888) q[0];
sx q[0];
rz(1.4317321) q[0];
rz(-1.8852385) q[1];
sx q[1];
rz(-1.4579371) q[1];
sx q[1];
rz(-1.2725405) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9894858) q[0];
sx q[0];
rz(-0.66755921) q[0];
sx q[0];
rz(-2.9015673) q[0];
rz(-pi) q[1];
rz(1.854004) q[2];
sx q[2];
rz(-2.0618366) q[2];
sx q[2];
rz(-1.1591737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9965197) q[1];
sx q[1];
rz(-0.47059083) q[1];
sx q[1];
rz(0.7820635) q[1];
x q[2];
rz(2.7719487) q[3];
sx q[3];
rz(-2.3497006) q[3];
sx q[3];
rz(-1.3835075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5566179) q[2];
sx q[2];
rz(-1.40404) q[2];
sx q[2];
rz(-0.72648826) q[2];
rz(0.11897421) q[3];
sx q[3];
rz(-1.3930895) q[3];
sx q[3];
rz(0.52072853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.93100905) q[0];
sx q[0];
rz(-1.9887661) q[0];
sx q[0];
rz(2.7880805) q[0];
rz(1.1734236) q[1];
sx q[1];
rz(-1.8559772) q[1];
sx q[1];
rz(-1.2849503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1160082) q[0];
sx q[0];
rz(-1.1366397) q[0];
sx q[0];
rz(0.15196073) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4305306) q[2];
sx q[2];
rz(-2.572142) q[2];
sx q[2];
rz(-2.0650149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28390005) q[1];
sx q[1];
rz(-1.5491747) q[1];
sx q[1];
rz(-1.8649072) q[1];
rz(-pi) q[2];
rz(-1.4015163) q[3];
sx q[3];
rz(-0.96099412) q[3];
sx q[3];
rz(-2.9823401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53359199) q[2];
sx q[2];
rz(-2.3054391) q[2];
sx q[2];
rz(0.81544915) q[2];
rz(2.4943374) q[3];
sx q[3];
rz(-1.3784958) q[3];
sx q[3];
rz(-0.7816202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4746998) q[0];
sx q[0];
rz(-0.26645461) q[0];
sx q[0];
rz(1.4315963) q[0];
rz(2.73009) q[1];
sx q[1];
rz(-1.9930379) q[1];
sx q[1];
rz(-3.0130951) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9978265) q[0];
sx q[0];
rz(-2.3989005) q[0];
sx q[0];
rz(-2.6608442) q[0];
x q[1];
rz(-0.13984404) q[2];
sx q[2];
rz(-2.2763414) q[2];
sx q[2];
rz(-1.2024513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.576927) q[1];
sx q[1];
rz(-0.91178545) q[1];
sx q[1];
rz(-0.66522775) q[1];
rz(-pi) q[2];
rz(2.6144876) q[3];
sx q[3];
rz(-1.2087421) q[3];
sx q[3];
rz(-1.3339588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0633885) q[2];
sx q[2];
rz(-1.3271164) q[2];
sx q[2];
rz(-0.77888954) q[2];
rz(1.286346) q[3];
sx q[3];
rz(-1.0948007) q[3];
sx q[3];
rz(-1.26545) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.417199) q[0];
sx q[0];
rz(-1.2204285) q[0];
sx q[0];
rz(1.3577419) q[0];
rz(-2.4667274) q[1];
sx q[1];
rz(-1.4874896) q[1];
sx q[1];
rz(1.0543324) q[1];
rz(2.8985698) q[2];
sx q[2];
rz(-2.7293548) q[2];
sx q[2];
rz(-1.7133452) q[2];
rz(1.5065083) q[3];
sx q[3];
rz(-0.74623204) q[3];
sx q[3];
rz(-0.85109477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
