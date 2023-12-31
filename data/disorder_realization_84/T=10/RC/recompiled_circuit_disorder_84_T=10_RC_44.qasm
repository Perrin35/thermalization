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
rz(2.6063483) q[0];
sx q[0];
rz(8.6707414) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(-1.9807281) q[1];
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
rz(1.6041243) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9795346) q[2];
sx q[2];
rz(-2.0433807) q[2];
sx q[2];
rz(-2.2651847) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-pi) q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(-2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(-2.0289039) q[0];
rz(-2.9878222) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(-1.3382834) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3801549) q[0];
sx q[0];
rz(-0.88012154) q[0];
sx q[0];
rz(1.3209016) q[0];
rz(-0.83805214) q[2];
sx q[2];
rz(-1.8997955) q[2];
sx q[2];
rz(2.7555562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4900134) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(-2.8298122) q[1];
rz(-1.409378) q[3];
sx q[3];
rz(-0.98940778) q[3];
sx q[3];
rz(2.5323109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.9632957) q[2];
rz(-0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798379) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-2.8082074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3231455) q[0];
sx q[0];
rz(-2.4172757) q[0];
sx q[0];
rz(-1.2562654) q[0];
x q[1];
rz(1.2522167) q[2];
sx q[2];
rz(-0.40553667) q[2];
sx q[2];
rz(2.6647654) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15409878) q[1];
sx q[1];
rz(-2.0425659) q[1];
sx q[1];
rz(-3.0497453) q[1];
rz(-1.7942579) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(-0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(2.3392759) q[2];
rz(-0.24117593) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872831) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(2.5774082) q[0];
rz(-0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(2.633458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6619686) q[0];
sx q[0];
rz(-1.9920252) q[0];
sx q[0];
rz(-2.7022916) q[0];
x q[1];
rz(1.5490683) q[2];
sx q[2];
rz(-0.93291514) q[2];
sx q[2];
rz(2.4456172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9307738) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(-0.37097431) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9594801) q[3];
sx q[3];
rz(-1.84387) q[3];
sx q[3];
rz(-1.9338716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(2.508146) q[2];
rz(2.5417035) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.0713779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288554) q[0];
sx q[0];
rz(-2.3561986) q[0];
sx q[0];
rz(-2.1942755) q[0];
rz(-0.84153701) q[2];
sx q[2];
rz(-1.8857937) q[2];
sx q[2];
rz(2.3685761) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(-0.96697076) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7171633) q[3];
sx q[3];
rz(-1.1491346) q[3];
sx q[3];
rz(-1.2210341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(-0.98199797) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-0.53034267) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(0.21662724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3743065) q[0];
sx q[0];
rz(-1.3067129) q[0];
sx q[0];
rz(-0.81861511) q[0];
rz(-pi) q[1];
rz(-0.63032052) q[2];
sx q[2];
rz(-1.3263055) q[2];
sx q[2];
rz(2.663161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10303282) q[1];
sx q[1];
rz(-1.235504) q[1];
sx q[1];
rz(-0.0086987728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8111749) q[3];
sx q[3];
rz(-2.4568395) q[3];
sx q[3];
rz(1.7184337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(2.1177297) q[2];
rz(-0.8762382) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4890471) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(-0.52893692) q[0];
rz(-1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-2.0524009) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76656094) q[0];
sx q[0];
rz(-1.2732817) q[0];
sx q[0];
rz(0.0059228063) q[0];
rz(-pi) q[1];
rz(-0.85765526) q[2];
sx q[2];
rz(-2.0115888) q[2];
sx q[2];
rz(2.8105274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2273119) q[1];
sx q[1];
rz(-3.1078048) q[1];
sx q[1];
rz(0.19126161) q[1];
rz(-0.090737061) q[3];
sx q[3];
rz(-1.6999348) q[3];
sx q[3];
rz(-1.7307626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.9160697) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65687031) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-2.2739676) q[0];
rz(-0.067226974) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(2.9464088) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.09571) q[0];
sx q[0];
rz(-1.3699431) q[0];
sx q[0];
rz(2.8413089) q[0];
x q[1];
rz(-0.033069177) q[2];
sx q[2];
rz(-0.40133921) q[2];
sx q[2];
rz(-2.6925342) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36181297) q[1];
sx q[1];
rz(-0.43174141) q[1];
sx q[1];
rz(-1.4713431) q[1];
rz(-0.35555367) q[3];
sx q[3];
rz(-2.0204633) q[3];
sx q[3];
rz(2.928283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.247867) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(2.2360738) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72682196) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(-3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(0.35531607) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126292) q[0];
sx q[0];
rz(-1.2041429) q[0];
sx q[0];
rz(-0.10495149) q[0];
x q[1];
rz(-2.9012868) q[2];
sx q[2];
rz(-1.086364) q[2];
sx q[2];
rz(2.3920609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3231343) q[1];
sx q[1];
rz(-1.5143906) q[1];
sx q[1];
rz(-1.8947381) q[1];
rz(-0.60042419) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67655247) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.5374373) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86826098) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(0.4459933) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83082643) q[0];
sx q[0];
rz(-2.5529773) q[0];
sx q[0];
rz(-1.1525843) q[0];
rz(-pi) q[1];
rz(2.5955574) q[2];
sx q[2];
rz(-1.2166426) q[2];
sx q[2];
rz(1.491577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7584596) q[1];
sx q[1];
rz(-2.8499095) q[1];
sx q[1];
rz(-3.0848857) q[1];
rz(-pi) q[2];
rz(-2.2311287) q[3];
sx q[3];
rz(-1.904389) q[3];
sx q[3];
rz(-0.098284483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(-1.7906174) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(-1.9679507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6931077) q[0];
sx q[0];
rz(-1.2379452) q[0];
sx q[0];
rz(-2.2647279) q[0];
rz(-1.7383472) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(1.5724814) q[2];
sx q[2];
rz(-1.6267852) q[2];
sx q[2];
rz(0.84597107) q[2];
rz(-2.8786447) q[3];
sx q[3];
rz(-0.96649747) q[3];
sx q[3];
rz(-1.4154712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
