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
rz(-0.22696683) q[0];
sx q[0];
rz(3.7530898) q[0];
sx q[0];
rz(8.8663085) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(-1.7668626) q[1];
sx q[1];
rz(0.078628063) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1100548) q[0];
sx q[0];
rz(-2.1221836) q[0];
sx q[0];
rz(-2.4999775) q[0];
x q[1];
rz(-0.67875454) q[2];
sx q[2];
rz(-1.4732483) q[2];
sx q[2];
rz(-2.7979948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6650162) q[1];
sx q[1];
rz(-1.4890097) q[1];
sx q[1];
rz(-0.25340124) q[1];
rz(-pi) q[2];
rz(-0.73689245) q[3];
sx q[3];
rz(-2.6813563) q[3];
sx q[3];
rz(0.62084197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2416396) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(2.107035) q[2];
rz(-0.11040802) q[3];
sx q[3];
rz(-1.6610828) q[3];
sx q[3];
rz(-2.835527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.13650525) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.9553631) q[0];
rz(1.8738481) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(1.3892106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0076375) q[0];
sx q[0];
rz(-1.2651772) q[0];
sx q[0];
rz(2.1307178) q[0];
rz(1.4591501) q[2];
sx q[2];
rz(-2.0497339) q[2];
sx q[2];
rz(-3.0618596) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8858365) q[1];
sx q[1];
rz(-2.6955018) q[1];
sx q[1];
rz(2.2913833) q[1];
x q[2];
rz(-0.095757397) q[3];
sx q[3];
rz(-2.1928722) q[3];
sx q[3];
rz(0.84703895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17239751) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(1.4364852) q[2];
rz(-1.2596333) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(-3.1148541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7773975) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(1.9245603) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.2478849) q[1];
sx q[1];
rz(1.7020114) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3453044) q[0];
sx q[0];
rz(-0.6718967) q[0];
sx q[0];
rz(-0.20033269) q[0];
rz(0.093482253) q[2];
sx q[2];
rz(-1.0614402) q[2];
sx q[2];
rz(-0.87967726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7172723) q[1];
sx q[1];
rz(-1.3417146) q[1];
sx q[1];
rz(3.0113264) q[1];
rz(-pi) q[2];
rz(-2.2398364) q[3];
sx q[3];
rz(-2.7983411) q[3];
sx q[3];
rz(-0.95860976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8000468) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(-0.81614196) q[2];
rz(3.1331565) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(1.5661904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208045) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(2.8724443) q[0];
rz(1.3828269) q[1];
sx q[1];
rz(-0.56126422) q[1];
sx q[1];
rz(0.47163481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021243367) q[0];
sx q[0];
rz(-2.2902852) q[0];
sx q[0];
rz(-0.18622744) q[0];
rz(3.0582171) q[2];
sx q[2];
rz(-2.7365587) q[2];
sx q[2];
rz(2.0138559) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9593352) q[1];
sx q[1];
rz(-1.8885837) q[1];
sx q[1];
rz(-2.5885568) q[1];
x q[2];
rz(0.90507953) q[3];
sx q[3];
rz(-1.8536012) q[3];
sx q[3];
rz(1.0156516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4857594) q[2];
sx q[2];
rz(-1.4071608) q[2];
sx q[2];
rz(-1.887623) q[2];
rz(0.29821011) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(1.6037174) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1042079) q[0];
sx q[0];
rz(-0.91008121) q[0];
sx q[0];
rz(0.75697672) q[0];
rz(1.0589927) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(0.34245488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.379102) q[0];
sx q[0];
rz(-2.6978793) q[0];
sx q[0];
rz(1.1136454) q[0];
x q[1];
rz(-1.9497237) q[2];
sx q[2];
rz(-2.4034326) q[2];
sx q[2];
rz(1.9110695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96794242) q[1];
sx q[1];
rz(-1.0999803) q[1];
sx q[1];
rz(-1.1556666) q[1];
rz(0.60635546) q[3];
sx q[3];
rz(-2.6487051) q[3];
sx q[3];
rz(2.3541401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4464438) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(-2.5214419) q[2];
rz(0.54689637) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(-2.0641522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81724375) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(1.2212344) q[0];
rz(2.0381894) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(0.81072909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87830905) q[0];
sx q[0];
rz(-1.4118402) q[0];
sx q[0];
rz(0.25860383) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9378433) q[2];
sx q[2];
rz(-1.019291) q[2];
sx q[2];
rz(-1.0288951) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6891046) q[1];
sx q[1];
rz(-2.9924722) q[1];
sx q[1];
rz(1.1959056) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8437988) q[3];
sx q[3];
rz(-1.9341365) q[3];
sx q[3];
rz(-2.9693574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3409884) q[2];
sx q[2];
rz(-1.9438513) q[2];
sx q[2];
rz(1.1654589) q[2];
rz(-3.0321339) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(-3.0808595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8282181) q[0];
sx q[0];
rz(-3.1309541) q[0];
sx q[0];
rz(0.64205545) q[0];
rz(-2.8984046) q[1];
sx q[1];
rz(-2.7311192) q[1];
sx q[1];
rz(0.64116716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579906) q[0];
sx q[0];
rz(-0.80785368) q[0];
sx q[0];
rz(-0.74961787) q[0];
x q[1];
rz(-1.9320609) q[2];
sx q[2];
rz(-0.85339499) q[2];
sx q[2];
rz(0.016005767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6536246) q[1];
sx q[1];
rz(-1.89978) q[1];
sx q[1];
rz(2.3423561) q[1];
x q[2];
rz(-1.9690597) q[3];
sx q[3];
rz(-1.7855111) q[3];
sx q[3];
rz(2.1977148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0082561) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(-0.94375098) q[2];
rz(2.4933695) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(2.4751723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7159395) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-2.7124523) q[0];
rz(0.40402135) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(2.2228352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0149834) q[0];
sx q[0];
rz(-2.1169239) q[0];
sx q[0];
rz(1.2864497) q[0];
rz(0.12523052) q[2];
sx q[2];
rz(-1.3117783) q[2];
sx q[2];
rz(0.7746402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1707216) q[1];
sx q[1];
rz(-1.7890463) q[1];
sx q[1];
rz(-1.8884601) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.346896) q[3];
sx q[3];
rz(-0.32156518) q[3];
sx q[3];
rz(2.0617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1227485) q[2];
sx q[2];
rz(-2.4825373) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(-0.33251479) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(1.7727218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.74828446) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(-1.4211753) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(-0.054904003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6643256) q[0];
sx q[0];
rz(-0.83819637) q[0];
sx q[0];
rz(-3.102836) q[0];
rz(-pi) q[1];
rz(-2.3308067) q[2];
sx q[2];
rz(-0.56466757) q[2];
sx q[2];
rz(-1.1155903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12704472) q[1];
sx q[1];
rz(-0.53847067) q[1];
sx q[1];
rz(1.1380859) q[1];
rz(-pi) q[2];
rz(0.53284455) q[3];
sx q[3];
rz(-2.5379764) q[3];
sx q[3];
rz(-0.6288022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2533337) q[2];
sx q[2];
rz(-0.71037018) q[2];
sx q[2];
rz(-0.19676512) q[2];
rz(-1.0737859) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(-2.4055068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1573023) q[0];
sx q[0];
rz(-2.2881916) q[0];
sx q[0];
rz(0.89944696) q[0];
rz(-2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(1.4003632) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.78681) q[0];
sx q[0];
rz(-1.3616832) q[0];
sx q[0];
rz(2.1378921) q[0];
rz(-1.0850995) q[2];
sx q[2];
rz(-1.7862831) q[2];
sx q[2];
rz(2.4485916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0728112) q[1];
sx q[1];
rz(-0.97807074) q[1];
sx q[1];
rz(-1.0981506) q[1];
rz(-pi) q[2];
rz(-0.66026779) q[3];
sx q[3];
rz(-0.83383027) q[3];
sx q[3];
rz(1.6136606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0805936) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(2.0628085) q[3];
sx q[3];
rz(-2.1507542) q[3];
sx q[3];
rz(2.259528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494736) q[0];
sx q[0];
rz(-1.619512) q[0];
sx q[0];
rz(-1.6721538) q[0];
rz(0.1334162) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(-3.1405814) q[2];
sx q[2];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5569729) q[2];
rz(-2.6399163) q[3];
sx q[3];
rz(-1.5080843) q[3];
sx q[3];
rz(2.9894473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
