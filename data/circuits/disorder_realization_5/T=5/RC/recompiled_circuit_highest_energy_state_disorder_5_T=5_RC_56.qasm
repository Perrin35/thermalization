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
rz(-2.8992803) q[0];
sx q[0];
rz(-0.87376142) q[0];
sx q[0];
rz(-1.4886966) q[0];
rz(-2.0714662) q[1];
sx q[1];
rz(-0.5572328) q[1];
sx q[1];
rz(-2.6723523) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0419641) q[0];
sx q[0];
rz(-1.6828186) q[0];
sx q[0];
rz(2.4449181) q[0];
x q[1];
rz(-0.69122423) q[2];
sx q[2];
rz(-0.35688734) q[2];
sx q[2];
rz(-1.7920902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3011613) q[1];
sx q[1];
rz(-0.60493785) q[1];
sx q[1];
rz(2.147232) q[1];
rz(-pi) q[2];
x q[2];
rz(1.038299) q[3];
sx q[3];
rz(-0.91404741) q[3];
sx q[3];
rz(0.19124243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0327518) q[2];
sx q[2];
rz(-0.27507541) q[2];
sx q[2];
rz(-1.8185505) q[2];
rz(2.1950586) q[3];
sx q[3];
rz(-0.60775477) q[3];
sx q[3];
rz(-0.058145903) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35689795) q[0];
sx q[0];
rz(-2.8234443) q[0];
sx q[0];
rz(1.1784877) q[0];
rz(-2.5937041) q[1];
sx q[1];
rz(-1.9375487) q[1];
sx q[1];
rz(-1.2335802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6797076) q[0];
sx q[0];
rz(-1.2872445) q[0];
sx q[0];
rz(-0.028282982) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.080670653) q[2];
sx q[2];
rz(-1.7364795) q[2];
sx q[2];
rz(2.5625429) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11352053) q[1];
sx q[1];
rz(-1.517761) q[1];
sx q[1];
rz(3.0854532) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5286851) q[3];
sx q[3];
rz(-0.41348413) q[3];
sx q[3];
rz(0.72847937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2445688) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(-1.3429674) q[2];
rz(0.025394214) q[3];
sx q[3];
rz(-1.7829021) q[3];
sx q[3];
rz(3.0620745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19567604) q[0];
sx q[0];
rz(-1.8657277) q[0];
sx q[0];
rz(2.6329686) q[0];
rz(1.4397844) q[1];
sx q[1];
rz(-1.6248117) q[1];
sx q[1];
rz(-1.6280828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6354409) q[0];
sx q[0];
rz(-2.1685528) q[0];
sx q[0];
rz(1.3805375) q[0];
x q[1];
rz(-1.7449101) q[2];
sx q[2];
rz(-1.422494) q[2];
sx q[2];
rz(-2.1789511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95437223) q[1];
sx q[1];
rz(-0.037211671) q[1];
sx q[1];
rz(0.12204917) q[1];
rz(-pi) q[2];
rz(-0.83714385) q[3];
sx q[3];
rz(-0.52609936) q[3];
sx q[3];
rz(-0.21080454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14579183) q[2];
sx q[2];
rz(-2.3172816) q[2];
sx q[2];
rz(-0.15131797) q[2];
rz(2.1794836) q[3];
sx q[3];
rz(-1.6302707) q[3];
sx q[3];
rz(0.37402737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91604084) q[0];
sx q[0];
rz(-2.2489838) q[0];
sx q[0];
rz(1.4501866) q[0];
rz(-2.080503) q[1];
sx q[1];
rz(-2.7174945) q[1];
sx q[1];
rz(2.012097) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0267222) q[0];
sx q[0];
rz(-0.78512275) q[0];
sx q[0];
rz(-2.8644187) q[0];
rz(0.53217066) q[2];
sx q[2];
rz(-1.2304365) q[2];
sx q[2];
rz(2.3305691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0982588) q[1];
sx q[1];
rz(-1.4318649) q[1];
sx q[1];
rz(-0.85161441) q[1];
x q[2];
rz(-3.1067209) q[3];
sx q[3];
rz(-1.4599055) q[3];
sx q[3];
rz(0.66005303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.16047655) q[2];
sx q[2];
rz(-0.99157292) q[2];
sx q[2];
rz(1.1264832) q[2];
rz(2.9142761) q[3];
sx q[3];
rz(-1.7132297) q[3];
sx q[3];
rz(0.050392438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.8166167) q[0];
sx q[0];
rz(-0.80642527) q[0];
sx q[0];
rz(2.4438013) q[0];
rz(-1.6119488) q[1];
sx q[1];
rz(-0.35747129) q[1];
sx q[1];
rz(2.7059817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39640204) q[0];
sx q[0];
rz(-1.309309) q[0];
sx q[0];
rz(3.0130062) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7694791) q[2];
sx q[2];
rz(-2.1796694) q[2];
sx q[2];
rz(-0.13932589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.322239) q[1];
sx q[1];
rz(-1.2585009) q[1];
sx q[1];
rz(0.38375591) q[1];
rz(-0.27542476) q[3];
sx q[3];
rz(-1.6425624) q[3];
sx q[3];
rz(-0.83930086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39516732) q[2];
sx q[2];
rz(-0.94234157) q[2];
sx q[2];
rz(-0.64685267) q[2];
rz(-2.4344889) q[3];
sx q[3];
rz(-3.0349858) q[3];
sx q[3];
rz(-0.99335563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3498722) q[0];
sx q[0];
rz(-2.7339022) q[0];
sx q[0];
rz(0.83734018) q[0];
rz(2.1912241) q[1];
sx q[1];
rz(-1.8888387) q[1];
sx q[1];
rz(-2.9577435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1342574) q[0];
sx q[0];
rz(-2.7401972) q[0];
sx q[0];
rz(-2.4069277) q[0];
rz(0.99509652) q[2];
sx q[2];
rz(-0.3873567) q[2];
sx q[2];
rz(-0.97706276) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5279416) q[1];
sx q[1];
rz(-1.1897478) q[1];
sx q[1];
rz(-0.77186976) q[1];
rz(-pi) q[2];
rz(0.23294874) q[3];
sx q[3];
rz(-1.6701352) q[3];
sx q[3];
rz(0.24391045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4652555) q[2];
sx q[2];
rz(-2.4691171) q[2];
sx q[2];
rz(0.23784168) q[2];
rz(0.075720876) q[3];
sx q[3];
rz(-0.12756158) q[3];
sx q[3];
rz(0.41847509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521562) q[0];
sx q[0];
rz(-2.6223923) q[0];
sx q[0];
rz(0.12553781) q[0];
rz(-2.708066) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(-1.4406406) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214366) q[0];
sx q[0];
rz(-1.6188863) q[0];
sx q[0];
rz(-1.5096971) q[0];
rz(0.48831482) q[2];
sx q[2];
rz(-1.4291233) q[2];
sx q[2];
rz(-0.087208793) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10769612) q[1];
sx q[1];
rz(-1.0748861) q[1];
sx q[1];
rz(2.0350083) q[1];
x q[2];
rz(0.85236211) q[3];
sx q[3];
rz(-0.96489513) q[3];
sx q[3];
rz(-2.5123799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0354249) q[2];
sx q[2];
rz(-1.9751534) q[2];
sx q[2];
rz(-2.9336145) q[2];
rz(0.68317991) q[3];
sx q[3];
rz(-2.4222789) q[3];
sx q[3];
rz(-1.6600367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55027562) q[0];
sx q[0];
rz(-2.2324201) q[0];
sx q[0];
rz(-0.43347484) q[0];
rz(0.47099653) q[1];
sx q[1];
rz(-0.29312557) q[1];
sx q[1];
rz(2.8004144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8467186) q[0];
sx q[0];
rz(-3.0187384) q[0];
sx q[0];
rz(-2.3289178) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8128692) q[2];
sx q[2];
rz(-1.717766) q[2];
sx q[2];
rz(0.16260553) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9097482) q[1];
sx q[1];
rz(-0.97601402) q[1];
sx q[1];
rz(0.21071319) q[1];
rz(-pi) q[2];
rz(-2.9753486) q[3];
sx q[3];
rz(-0.59985929) q[3];
sx q[3];
rz(-1.8458008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53182536) q[2];
sx q[2];
rz(-0.34588459) q[2];
sx q[2];
rz(0.38066614) q[2];
rz(-0.62764132) q[3];
sx q[3];
rz(-2.6142879) q[3];
sx q[3];
rz(0.92831534) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8220383) q[0];
sx q[0];
rz(-1.1975937) q[0];
sx q[0];
rz(2.8354216) q[0];
rz(0.51433688) q[1];
sx q[1];
rz(-2.3898333) q[1];
sx q[1];
rz(2.9591282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5222591) q[0];
sx q[0];
rz(-2.3292208) q[0];
sx q[0];
rz(1.7343397) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9391871) q[2];
sx q[2];
rz(-2.6037964) q[2];
sx q[2];
rz(2.0356195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.811552) q[1];
sx q[1];
rz(-2.7024058) q[1];
sx q[1];
rz(-1.9387092) q[1];
rz(-pi) q[2];
rz(0.45643198) q[3];
sx q[3];
rz(-1.8805537) q[3];
sx q[3];
rz(-2.8995958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3719486) q[2];
sx q[2];
rz(-0.19522788) q[2];
sx q[2];
rz(2.7455184) q[2];
rz(1.1560446) q[3];
sx q[3];
rz(-2.5542185) q[3];
sx q[3];
rz(0.42266947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3892155) q[0];
sx q[0];
rz(-0.14425819) q[0];
sx q[0];
rz(-0.17917646) q[0];
rz(1.3009118) q[1];
sx q[1];
rz(-2.6866899) q[1];
sx q[1];
rz(0.5295583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0146739) q[0];
sx q[0];
rz(-1.3171985) q[0];
sx q[0];
rz(0.34832592) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1346369) q[2];
sx q[2];
rz(-1.8514625) q[2];
sx q[2];
rz(-0.2427189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5770108) q[1];
sx q[1];
rz(-1.8576) q[1];
sx q[1];
rz(-0.69820349) q[1];
rz(2.8488068) q[3];
sx q[3];
rz(-0.55582992) q[3];
sx q[3];
rz(1.0852927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35310465) q[2];
sx q[2];
rz(-0.64211923) q[2];
sx q[2];
rz(-2.8118964) q[2];
rz(1.6126136) q[3];
sx q[3];
rz(-1.2639045) q[3];
sx q[3];
rz(0.32040709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223087) q[0];
sx q[0];
rz(-1.6884463) q[0];
sx q[0];
rz(1.6519188) q[0];
rz(-2.655401) q[1];
sx q[1];
rz(-1.9959027) q[1];
sx q[1];
rz(-2.1015658) q[1];
rz(1.4522105) q[2];
sx q[2];
rz(-2.2439742) q[2];
sx q[2];
rz(-3.0912666) q[2];
rz(-1.7320932) q[3];
sx q[3];
rz(-2.0031569) q[3];
sx q[3];
rz(-1.835837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
