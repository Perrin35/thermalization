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
rz(-5.4929805) q[1];
sx q[1];
rz(5.0561855) q[1];
sx q[1];
rz(8.2639134) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863659) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(3.0774374) q[0];
rz(-pi) q[1];
rz(1.8765175) q[2];
sx q[2];
rz(-2.6439878) q[2];
sx q[2];
rz(-1.2212317) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4143715) q[1];
sx q[1];
rz(-1.9804269) q[1];
sx q[1];
rz(-1.0135256) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17920223) q[3];
sx q[3];
rz(-1.407302) q[3];
sx q[3];
rz(1.6678099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(2.4751002) q[2];
rz(-0.50981057) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(-1.8681017) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(-2.9878222) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(-1.3382834) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35182686) q[0];
sx q[0];
rz(-1.3789982) q[0];
sx q[0];
rz(2.4353794) q[0];
rz(-2.0427225) q[2];
sx q[2];
rz(-0.79052351) q[2];
sx q[2];
rz(-2.3015442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6515793) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(-2.8298122) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23985858) q[3];
sx q[3];
rz(-2.5407102) q[3];
sx q[3];
rz(-2.2440653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(1.178297) q[2];
rz(0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.1798379) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(1.5270365) q[0];
rz(-2.4987192) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(0.33338526) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3231455) q[0];
sx q[0];
rz(-2.4172757) q[0];
sx q[0];
rz(1.2562654) q[0];
x q[1];
rz(1.9579499) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(-1.3882335) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.15409878) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(-3.0497453) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8414855) q[3];
sx q[3];
rz(-2.4862053) q[3];
sx q[3];
rz(2.3902077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.12144111) q[2];
sx q[2];
rz(-1.3139498) q[2];
sx q[2];
rz(-0.80231673) q[2];
rz(-0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-0.10087092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0543095) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(0.50813466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80721891) q[0];
sx q[0];
rz(-2.5426572) q[0];
sx q[0];
rz(2.3301624) q[0];
x q[1];
rz(1.5925243) q[2];
sx q[2];
rz(-2.2086775) q[2];
sx q[2];
rz(2.4456172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21081884) q[1];
sx q[1];
rz(-2.3465119) q[1];
sx q[1];
rz(0.37097431) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18211256) q[3];
sx q[3];
rz(-1.84387) q[3];
sx q[3];
rz(-1.207721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-2.508146) q[2];
rz(0.59988919) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-0.19609837) q[0];
rz(-1.8798937) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-1.0713779) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51273726) q[0];
sx q[0];
rz(-2.3561986) q[0];
sx q[0];
rz(-0.94731713) q[0];
rz(-pi) q[1];
rz(-0.41197889) q[2];
sx q[2];
rz(-0.88468555) q[2];
sx q[2];
rz(-1.0678604) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86712671) q[1];
sx q[1];
rz(-1.9879838) q[1];
sx q[1];
rz(2.8309114) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0282182) q[3];
sx q[3];
rz(-1.1856106) q[3];
sx q[3];
rz(-0.53264632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(0.98199797) q[2];
rz(-0.18520959) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.265032) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380602) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(0.21662724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075561698) q[0];
sx q[0];
rz(-0.78853411) q[0];
sx q[0];
rz(1.1939474) q[0];
rz(-pi) q[1];
rz(1.2712619) q[2];
sx q[2];
rz(-2.1795863) q[2];
sx q[2];
rz(1.2671721) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0121289) q[1];
sx q[1];
rz(-2.8061917) q[1];
sx q[1];
rz(1.5458376) q[1];
x q[2];
rz(-2.9495903) q[3];
sx q[3];
rz(-2.2322906) q[3];
sx q[3];
rz(-2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(2.1177297) q[2];
rz(-2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(-0.52893692) q[0];
rz(1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-1.0891917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3750317) q[0];
sx q[0];
rz(-1.8683109) q[0];
sx q[0];
rz(-3.1356698) q[0];
rz(0.94605298) q[2];
sx q[2];
rz(-0.81753987) q[2];
sx q[2];
rz(0.78150502) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84767064) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(-0.033172219) q[1];
x q[2];
rz(-0.090737061) q[3];
sx q[3];
rz(-1.6999348) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12525325) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(-1.2255229) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
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
rz(-2.9464088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0474437) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(0.6028428) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5848265) q[2];
sx q[2];
rz(-1.9719035) q[2];
sx q[2];
rz(-0.48497981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47123779) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(-0.045713748) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.786039) q[3];
sx q[3];
rz(-2.0204633) q[3];
sx q[3];
rz(-2.928283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.247867) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(-0.90551886) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(2.982443) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(1.0409521) q[0];
rz(-0.078661593) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(-0.35531607) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5268363) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(-1.8371131) q[0];
rz(-pi) q[1];
rz(-1.9955194) q[2];
sx q[2];
rz(-2.6051084) q[2];
sx q[2];
rz(1.2338961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73341093) q[1];
sx q[1];
rz(-1.8942041) q[1];
sx q[1];
rz(0.059493382) q[1];
rz(-2.5411685) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(-1.0004117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.7379649) q[2];
rz(3.1353531) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733317) q[0];
sx q[0];
rz(-1.9737759) q[0];
sx q[0];
rz(1.2982752) q[0];
rz(-0.4459933) q[1];
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
rz(1.3215799) q[0];
sx q[0];
rz(-2.1029148) q[0];
sx q[0];
rz(-0.26474712) q[0];
rz(1.1624596) q[2];
sx q[2];
rz(-2.0795341) q[2];
sx q[2];
rz(-3.0131154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69925752) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(1.553781) q[1];
rz(-pi) q[2];
rz(2.0850639) q[3];
sx q[3];
rz(-2.4132055) q[3];
sx q[3];
rz(-1.8715093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62844244) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.2107522) q[3];
sx q[3];
rz(-2.4891709) q[3];
sx q[3];
rz(-1.8579033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
