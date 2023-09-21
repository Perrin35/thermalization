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
rz(-1.9807281) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55863419) q[0];
sx q[0];
rz(-1.5412124) q[0];
sx q[0];
rz(0.47870676) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2650752) q[2];
sx q[2];
rz(-0.49760488) q[2];
sx q[2];
rz(-1.9203609) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7272211) q[1];
sx q[1];
rz(-1.1611658) q[1];
sx q[1];
rz(-1.0135256) q[1];
rz(-2.394862) q[3];
sx q[3];
rz(-0.24198469) q[3];
sx q[3];
rz(0.63499588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(0.66649246) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(-1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3614685) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(-1.8033093) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7614377) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(-1.3209016) q[0];
rz(-2.3035405) q[2];
sx q[2];
rz(-1.8997955) q[2];
sx q[2];
rz(0.38603644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6515793) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(-2.8298122) q[1];
x q[2];
rz(-2.9017341) q[3];
sx q[3];
rz(-2.5407102) q[3];
sx q[3];
rz(2.2440653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(-0.96238771) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(-0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96175471) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(0.33338526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282368) q[0];
sx q[0];
rz(-2.2524999) q[0];
sx q[0];
rz(-0.26716726) q[0];
rz(-0.1336735) q[2];
sx q[2];
rz(-1.954827) q[2];
sx q[2];
rz(3.0093699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0957291) q[1];
sx q[1];
rz(-2.6616272) q[1];
sx q[1];
rz(-1.392925) q[1];
rz(-pi) q[2];
rz(2.8414855) q[3];
sx q[3];
rz(-2.4862053) q[3];
sx q[3];
rz(0.75138498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(-0.80231673) q[2];
rz(-2.9004167) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(-0.10087092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(0.56418443) q[0];
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
rz(3.0429338) q[0];
sx q[0];
rz(-1.9693976) q[0];
sx q[0];
rz(-1.1111141) q[0];
rz(-pi) q[1];
rz(-1.5925243) q[2];
sx q[2];
rz(-0.93291514) q[2];
sx q[2];
rz(2.4456172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29624629) q[1];
sx q[1];
rz(-0.84282649) q[1];
sx q[1];
rz(-1.216757) q[1];
rz(-2.1448137) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(2.5329778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4529139) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(-1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(2.9454943) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-1.0713779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51273726) q[0];
sx q[0];
rz(-0.78539408) q[0];
sx q[0];
rz(0.94731713) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1159665) q[2];
sx q[2];
rz(-0.78274667) q[2];
sx q[2];
rz(-2.6775529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5672946) q[1];
sx q[1];
rz(-1.287536) q[1];
sx q[1];
rz(-2.0064555) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3134872) q[3];
sx q[3];
rz(-2.5525186) q[3];
sx q[3];
rz(2.7554054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(2.9249654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43603555) q[0];
sx q[0];
rz(-0.85058054) q[0];
sx q[0];
rz(2.7869422) q[0];
rz(2.7411813) q[2];
sx q[2];
rz(-0.67001221) q[2];
sx q[2];
rz(-2.3695721) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0385598) q[1];
sx q[1];
rz(-1.235504) q[1];
sx q[1];
rz(-3.1328939) q[1];
rz(1.3304177) q[3];
sx q[3];
rz(-0.68475311) q[3];
sx q[3];
rz(1.7184337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(1.2715626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6525456) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(-2.6126557) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-1.0891917) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76656094) q[0];
sx q[0];
rz(-1.2732817) q[0];
sx q[0];
rz(0.0059228063) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1955397) q[2];
sx q[2];
rz(-2.3240528) q[2];
sx q[2];
rz(0.78150502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91428077) q[1];
sx q[1];
rz(-3.1078048) q[1];
sx q[1];
rz(0.19126161) q[1];
rz(3.0508556) q[3];
sx q[3];
rz(-1.6999348) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12525325) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.9160697) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(-2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1104689) q[1];
sx q[1];
rz(-0.19518383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0458826) q[0];
sx q[0];
rz(-1.3699431) q[0];
sx q[0];
rz(2.8413089) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5567661) q[2];
sx q[2];
rz(-1.1696891) q[2];
sx q[2];
rz(0.48497981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6703549) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(-0.045713748) q[1];
x q[2];
rz(2.1956452) q[3];
sx q[3];
rz(-2.5759856) q[3];
sx q[3];
rz(2.2212976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.247867) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(-0.90551886) q[2];
rz(-1.9838105) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72682196) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(2.1006405) q[0];
rz(-0.078661593) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(-0.35531607) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32896341) q[0];
sx q[0];
rz(-1.9374498) q[0];
sx q[0];
rz(-3.0366412) q[0];
rz(2.0673429) q[2];
sx q[2];
rz(-1.3585919) q[2];
sx q[2];
rz(2.4339536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4081817) q[1];
sx q[1];
rz(-1.8942041) q[1];
sx q[1];
rz(0.059493382) q[1];
rz(0.60042419) q[3];
sx q[3];
rz(-1.416559) q[3];
sx q[3];
rz(1.0004117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67655247) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(1.7379649) q[2];
rz(-3.1353531) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86826098) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(1.8433174) q[0];
rz(-2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
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
rz(-2.1029148) q[0];
sx q[0];
rz(0.26474712) q[0];
rz(-pi) q[1];
rz(-0.61873318) q[2];
sx q[2];
rz(-0.64090568) q[2];
sx q[2];
rz(-0.59782019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86665323) q[1];
sx q[1];
rz(-1.5544974) q[1];
sx q[1];
rz(2.8503522) q[1];
rz(-pi) q[2];
rz(2.7281076) q[3];
sx q[3];
rz(-2.1889912) q[3];
sx q[3];
rz(1.9181044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62844244) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-1.1395617) q[2];
rz(1.7906174) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.9679507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6931077) q[0];
sx q[0];
rz(-1.9036475) q[0];
sx q[0];
rz(0.87686476) q[0];
rz(-1.4032455) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(-0.055988978) q[2];
sx q[2];
rz(-1.5724788) q[2];
sx q[2];
rz(2.4168617) q[2];
rz(0.95004497) q[3];
sx q[3];
rz(-1.3552356) q[3];
sx q[3];
rz(-3.1380359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];