OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(0.75403655) q[0];
rz(0.79020483) q[1];
sx q[1];
rz(-1.2269998) q[1];
sx q[1];
rz(-1.1608646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0275118) q[0];
sx q[0];
rz(-2.049276) q[0];
sx q[0];
rz(-1.5374684) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16205807) q[2];
sx q[2];
rz(-1.098212) q[2];
sx q[2];
rz(-0.87640793) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7293284) q[1];
sx q[1];
rz(-2.4630479) q[1];
sx q[1];
rz(-2.2581984) q[1];
rz(-2.394862) q[3];
sx q[3];
rz(-0.24198469) q[3];
sx q[3];
rz(-2.5065968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(0.50981057) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(-1.1126888) q[0];
rz(-0.15377046) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(1.3382834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7897658) q[0];
sx q[0];
rz(-1.7625945) q[0];
sx q[0];
rz(2.4353794) q[0];
rz(-pi) q[1];
rz(-2.3035405) q[2];
sx q[2];
rz(-1.2417972) q[2];
sx q[2];
rz(-0.38603644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2087304) q[1];
sx q[1];
rz(-1.4534229) q[1];
sx q[1];
rz(2.7670303) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58740939) q[3];
sx q[3];
rz(-1.7055159) q[3];
sx q[3];
rz(-0.87232529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(-2.1792049) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96175471) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(-2.8082074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51334914) q[0];
sx q[0];
rz(-1.7772563) q[0];
sx q[0];
rz(2.2702361) q[0];
x q[1];
rz(-0.1336735) q[2];
sx q[2];
rz(-1.954827) q[2];
sx q[2];
rz(3.0093699) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15409878) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(-0.091847329) q[1];
rz(-1.7942579) q[3];
sx q[3];
rz(-0.94933214) q[3];
sx q[3];
rz(-2.7623451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12144111) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(-0.80231673) q[2];
rz(-2.9004167) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(-0.50813466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6619686) q[0];
sx q[0];
rz(-1.1495674) q[0];
sx q[0];
rz(2.7022916) q[0];
rz(-2.5035985) q[2];
sx q[2];
rz(-1.5533414) q[2];
sx q[2];
rz(0.88776112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29624629) q[1];
sx q[1];
rz(-2.2987662) q[1];
sx q[1];
rz(-1.9248357) q[1];
x q[2];
rz(1.8482314) q[3];
sx q[3];
rz(-1.395503) q[3];
sx q[3];
rz(0.41269916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(1.6413123) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-0.19609837) q[0];
rz(-1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-2.0702147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51273726) q[0];
sx q[0];
rz(-0.78539408) q[0];
sx q[0];
rz(-0.94731713) q[0];
x q[1];
rz(-0.41197889) q[2];
sx q[2];
rz(-0.88468555) q[2];
sx q[2];
rz(2.0737322) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6043678) q[1];
sx q[1];
rz(-0.51465263) q[1];
sx q[1];
rz(-2.1746219) q[1];
rz(1.1133744) q[3];
sx q[3];
rz(-1.955982) q[3];
sx q[3];
rz(-2.6089463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(0.98199797) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.8416587) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(-2.9249654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3743065) q[0];
sx q[0];
rz(-1.8348798) q[0];
sx q[0];
rz(-0.81861511) q[0];
x q[1];
rz(-0.40041133) q[2];
sx q[2];
rz(-2.4715804) q[2];
sx q[2];
rz(2.3695721) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10303282) q[1];
sx q[1];
rz(-1.235504) q[1];
sx q[1];
rz(-0.0086987728) q[1];
rz(-pi) q[2];
rz(0.90029193) q[3];
sx q[3];
rz(-1.4196463) q[3];
sx q[3];
rz(0.33526648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(1.023863) q[2];
rz(2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(1.2715626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6525456) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(0.52893692) q[0];
rz(1.5286998) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-2.0524009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80597164) q[0];
sx q[0];
rz(-1.5651337) q[0];
sx q[0];
rz(-1.2732768) q[0];
x q[1];
rz(-0.5577001) q[2];
sx q[2];
rz(-0.93765646) q[2];
sx q[2];
rz(-1.5932839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.293922) q[1];
sx q[1];
rz(-1.5643745) q[1];
sx q[1];
rz(-0.033172219) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96157162) q[3];
sx q[3];
rz(-2.9839091) q[3];
sx q[3];
rz(-0.79573436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65687031) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-2.2739676) q[0];
rz(-3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(0.19518383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0941489) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(0.6028428) q[0];
rz(0.40114258) q[2];
sx q[2];
rz(-1.5578798) q[2];
sx q[2];
rz(-2.0502979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36181297) q[1];
sx q[1];
rz(-2.7098512) q[1];
sx q[1];
rz(-1.4713431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35555367) q[3];
sx q[3];
rz(-2.0204633) q[3];
sx q[3];
rz(2.928283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(2.2360738) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(-1.0409521) q[0];
rz(-0.078661593) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(2.7862766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32896341) q[0];
sx q[0];
rz(-1.9374498) q[0];
sx q[0];
rz(-3.0366412) q[0];
x q[1];
rz(-2.9012868) q[2];
sx q[2];
rz(-1.086364) q[2];
sx q[2];
rz(2.3920609) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91868672) q[1];
sx q[1];
rz(-0.32864535) q[1];
sx q[1];
rz(1.7463643) q[1];
x q[2];
rz(2.8730632) q[3];
sx q[3];
rz(-0.61754698) q[3];
sx q[3];
rz(-2.7919046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67655247) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(-1.7379649) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86826098) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(-1.2982752) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8200127) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(2.8768455) q[0];
rz(1.1624596) q[2];
sx q[2];
rz(-2.0795341) q[2];
sx q[2];
rz(-3.0131154) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4423351) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(-1.5878116) q[1];
rz(-pi) q[2];
rz(-0.91046393) q[3];
sx q[3];
rz(-1.904389) q[3];
sx q[3];
rz(-3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(1.1395617) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.448485) q[0];
sx q[0];
rz(-1.2379452) q[0];
sx q[0];
rz(-2.2647279) q[0];
rz(1.7383472) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(-0.030056603) q[2];
sx q[2];
rz(-3.0855784) q[2];
sx q[2];
rz(-2.2655178) q[2];
rz(1.9308405) q[3];
sx q[3];
rz(-2.4891709) q[3];
sx q[3];
rz(-1.8579033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];