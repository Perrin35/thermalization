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
rz(-2.3875561) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(1.1608646) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55863419) q[0];
sx q[0];
rz(-1.6003803) q[0];
sx q[0];
rz(-2.6628859) q[0];
rz(-pi) q[1];
rz(1.8765175) q[2];
sx q[2];
rz(-2.6439878) q[2];
sx q[2];
rz(-1.2212317) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7293284) q[1];
sx q[1];
rz(-0.67854478) q[1];
sx q[1];
rz(-0.88339427) q[1];
x q[2];
rz(2.394862) q[3];
sx q[3];
rz(-2.899608) q[3];
sx q[3];
rz(0.63499588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51241088) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(1.1126888) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(-1.8033093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35182686) q[0];
sx q[0];
rz(-1.3789982) q[0];
sx q[0];
rz(-0.70621323) q[0];
rz(-pi) q[1];
rz(-2.3035405) q[2];
sx q[2];
rz(-1.2417972) q[2];
sx q[2];
rz(2.7555562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4900134) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(2.8298122) q[1];
rz(-pi) q[2];
rz(-0.23985858) q[3];
sx q[3];
rz(-0.60088241) q[3];
sx q[3];
rz(-0.89752737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(1.9632957) q[2];
rz(2.1792049) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(-1.5270365) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(0.33338526) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51334914) q[0];
sx q[0];
rz(-1.7772563) q[0];
sx q[0];
rz(-0.87135656) q[0];
x q[1];
rz(0.1336735) q[2];
sx q[2];
rz(-1.954827) q[2];
sx q[2];
rz(-3.0093699) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9874939) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(-3.0497453) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63342996) q[3];
sx q[3];
rz(-1.389635) q[3];
sx q[3];
rz(1.0599979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(-2.3392759) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(0.10087092) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(-2.5634649) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(2.633458) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.5882512) q[2];
sx q[2];
rz(2.2538315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9307738) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(-0.37097431) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99677892) q[3];
sx q[3];
rz(-0.32696163) q[3];
sx q[3];
rz(0.60861482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4529139) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-2.508146) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863581) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-2.9454943) q[0];
rz(-1.8798937) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-1.0713779) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8349985) q[0];
sx q[0];
rz(-0.95933611) q[0];
sx q[0];
rz(2.6131265) q[0];
rz(-0.41197889) q[2];
sx q[2];
rz(-2.2569071) q[2];
sx q[2];
rz(-2.0737322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(0.96697076) q[1];
rz(1.1133744) q[3];
sx q[3];
rz(-1.1856106) q[3];
sx q[3];
rz(-0.53264632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(-0.98199797) q[2];
rz(2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380602) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-0.21662724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075561698) q[0];
sx q[0];
rz(-0.78853411) q[0];
sx q[0];
rz(-1.1939474) q[0];
x q[1];
rz(-1.8703307) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.2671721) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0121289) q[1];
sx q[1];
rz(-2.8061917) q[1];
sx q[1];
rz(-1.5458376) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19200237) q[3];
sx q[3];
rz(-0.90930206) q[3];
sx q[3];
rz(-1.1166752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6525456) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(-2.6126557) q[0];
rz(-1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-2.0524009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.335621) q[0];
sx q[0];
rz(-1.5651337) q[0];
sx q[0];
rz(1.8683158) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1955397) q[2];
sx q[2];
rz(-2.3240528) q[2];
sx q[2];
rz(2.3600876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.91428077) q[1];
sx q[1];
rz(-0.033787878) q[1];
sx q[1];
rz(2.950331) q[1];
rz(-pi) q[2];
rz(1.7004622) q[3];
sx q[3];
rz(-1.4808169) q[3];
sx q[3];
rz(0.17168301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(-1.4922173) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65687031) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(-2.2739676) q[0];
rz(-3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(0.19518383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282043) q[0];
sx q[0];
rz(-1.2767316) q[0];
sx q[0];
rz(1.3608027) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40114258) q[2];
sx q[2];
rz(-1.5837129) q[2];
sx q[2];
rz(1.0912947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6703549) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(3.0958789) q[1];
rz(1.0953426) q[3];
sx q[3];
rz(-1.8896777) q[3];
sx q[3];
rz(1.9441324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8937257) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72682196) q[0];
sx q[0];
rz(-1.0960217) q[0];
sx q[0];
rz(2.1006405) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-2.9610596) q[1];
sx q[1];
rz(-2.7862766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61475635) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(1.3044796) q[0];
rz(2.9012868) q[2];
sx q[2];
rz(-2.0552286) q[2];
sx q[2];
rz(2.3920609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81845835) q[1];
sx q[1];
rz(-1.6272021) q[1];
sx q[1];
rz(-1.8947381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5411685) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(3.1353531) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86826098) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(1.2982752) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(2.840852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8200127) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(2.8768455) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9791331) q[2];
sx q[2];
rz(-1.0620585) q[2];
sx q[2];
rz(-0.12847729) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86665323) q[1];
sx q[1];
rz(-1.5870952) q[1];
sx q[1];
rz(-0.29124041) q[1];
rz(-pi) q[2];
rz(0.91046393) q[3];
sx q[3];
rz(-1.2372036) q[3];
sx q[3];
rz(-3.0433082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5131502) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(1.7906174) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6931077) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(1.4032455) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(-0.055988978) q[2];
sx q[2];
rz(-1.5724788) q[2];
sx q[2];
rz(2.4168617) q[2];
rz(-0.26294796) q[3];
sx q[3];
rz(-2.1750952) q[3];
sx q[3];
rz(1.7261214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
