OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8755655) q[0];
sx q[0];
rz(-2.6063483) q[0];
sx q[0];
rz(2.3875561) q[0];
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(1.1608646) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1863659) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(0.064155302) q[0];
rz(-2.0487469) q[2];
sx q[2];
rz(-1.7149601) q[2];
sx q[2];
rz(2.521487) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7418993) q[1];
sx q[1];
rz(-2.0772935) q[1];
sx q[1];
rz(2.6687117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4046895) q[3];
sx q[3];
rz(-1.7475834) q[3];
sx q[3];
rz(-0.12648957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(0.66649246) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.8681017) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(1.1126888) q[0];
rz(-2.9878222) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.3382834) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3801549) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(1.3209016) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43054994) q[2];
sx q[2];
rz(-0.88532788) q[2];
sx q[2];
rz(1.4677043) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6515793) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(-2.8298122) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9017341) q[3];
sx q[3];
rz(-0.60088241) q[3];
sx q[3];
rz(0.89752737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.9632957) q[2];
rz(-2.1792049) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(-2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(-1.6145561) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(2.8082074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51334914) q[0];
sx q[0];
rz(-1.3643364) q[0];
sx q[0];
rz(-2.2702361) q[0];
rz(-pi) q[1];
rz(-1.9579499) q[2];
sx q[2];
rz(-1.6946812) q[2];
sx q[2];
rz(-1.3882335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4585321) q[1];
sx q[1];
rz(-1.6525869) q[1];
sx q[1];
rz(2.044278) q[1];
x q[2];
rz(-2.5081627) q[3];
sx q[3];
rz(-1.7519577) q[3];
sx q[3];
rz(-2.0815947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(-0.80231673) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-0.10087092) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0872831) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(-0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-0.50813466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0429338) q[0];
sx q[0];
rz(-1.172195) q[0];
sx q[0];
rz(-1.1111141) q[0];
rz(-pi) q[1];
rz(3.1122909) q[2];
sx q[2];
rz(-2.5033931) q[2];
sx q[2];
rz(2.4820941) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21081884) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(0.37097431) q[1];
rz(-0.99677892) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(-2.5329778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4529139) q[2];
sx q[2];
rz(-2.8082509) q[2];
sx q[2];
rz(2.508146) q[2];
rz(2.5417035) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(2.0702147) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5539861) q[0];
sx q[0];
rz(-1.1452132) q[0];
sx q[0];
rz(2.2527184) q[0];
rz(-pi) q[1];
rz(0.41197889) q[2];
sx q[2];
rz(-2.2569071) q[2];
sx q[2];
rz(2.0737322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86712671) q[1];
sx q[1];
rz(-1.1536088) q[1];
sx q[1];
rz(-2.8309114) q[1];
rz(-pi) q[2];
rz(2.7171633) q[3];
sx q[3];
rz(-1.1491346) q[3];
sx q[3];
rz(-1.2210341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035325) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.2999339) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-2.9249654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43603555) q[0];
sx q[0];
rz(-0.85058054) q[0];
sx q[0];
rz(-2.7869422) q[0];
x q[1];
rz(1.8703307) q[2];
sx q[2];
rz(-2.1795863) q[2];
sx q[2];
rz(-1.2671721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10303282) q[1];
sx q[1];
rz(-1.235504) q[1];
sx q[1];
rz(-0.0086987728) q[1];
rz(-pi) q[2];
rz(-0.19200237) q[3];
sx q[3];
rz(-0.90930206) q[3];
sx q[3];
rz(-2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(0.8762382) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(-1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(0.52893692) q[0];
rz(-1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-2.0524009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74635909) q[0];
sx q[0];
rz(-0.2975718) q[0];
sx q[0];
rz(1.5901106) q[0];
rz(-pi) q[1];
rz(-2.5838926) q[2];
sx q[2];
rz(-2.2039362) q[2];
sx q[2];
rz(1.5483088) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.84767064) q[1];
sx q[1];
rz(-1.5643745) q[1];
sx q[1];
rz(-3.1084204) q[1];
rz(-pi) q[2];
rz(0.090737061) q[3];
sx q[3];
rz(-1.4416579) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(-0.15587458) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(0.86762506) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(-2.9464088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.09571) q[0];
sx q[0];
rz(-1.7716496) q[0];
sx q[0];
rz(2.8413089) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40114258) q[2];
sx q[2];
rz(-1.5837129) q[2];
sx q[2];
rz(1.0912947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36181297) q[1];
sx q[1];
rz(-2.7098512) q[1];
sx q[1];
rz(-1.6702495) q[1];
x q[2];
rz(1.0953426) q[3];
sx q[3];
rz(-1.2519149) q[3];
sx q[3];
rz(-1.9441324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.247867) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(0.90551886) q[2];
rz(1.1577822) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-2.9610596) q[1];
sx q[1];
rz(0.35531607) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5268363) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(-1.8371131) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24030587) q[2];
sx q[2];
rz(-2.0552286) q[2];
sx q[2];
rz(2.3920609) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81845835) q[1];
sx q[1];
rz(-1.6272021) q[1];
sx q[1];
rz(-1.2468546) q[1];
rz(-0.26852946) q[3];
sx q[3];
rz(-0.61754698) q[3];
sx q[3];
rz(0.34968801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67655247) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(1.7379649) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86826098) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(-1.8433174) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-2.840852) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3107662) q[0];
sx q[0];
rz(-2.5529773) q[0];
sx q[0];
rz(-1.1525843) q[0];
x q[1];
rz(0.54603521) q[2];
sx q[2];
rz(-1.9249501) q[2];
sx q[2];
rz(-1.6500157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3831331) q[1];
sx q[1];
rz(-0.29168318) q[1];
sx q[1];
rz(-3.0848857) q[1];
rz(-pi) q[2];
rz(-1.0565287) q[3];
sx q[3];
rz(-0.72838718) q[3];
sx q[3];
rz(-1.2700833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62844244) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-1.1395617) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448485) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(-1.7383472) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(-1.5691112) q[2];
sx q[2];
rz(-1.6267852) q[2];
sx q[2];
rz(0.84597107) q[2];
rz(-2.1915477) q[3];
sx q[3];
rz(-1.3552356) q[3];
sx q[3];
rz(-3.1380359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
