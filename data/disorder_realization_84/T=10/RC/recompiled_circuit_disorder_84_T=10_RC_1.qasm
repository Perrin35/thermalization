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
rz(0.55863419) q[0];
sx q[0];
rz(-1.5412124) q[0];
sx q[0];
rz(0.47870676) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0928458) q[2];
sx q[2];
rz(-1.4266326) q[2];
sx q[2];
rz(0.62010566) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41226422) q[1];
sx q[1];
rz(-0.67854478) q[1];
sx q[1];
rz(0.88339427) q[1];
x q[2];
rz(-0.74673064) q[3];
sx q[3];
rz(-2.899608) q[3];
sx q[3];
rz(0.63499588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(1.1126888) q[0];
rz(-0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.8033093) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99909821) q[0];
sx q[0];
rz(-2.4141443) q[0];
sx q[0];
rz(2.8508458) q[0];
rz(-pi) q[1];
rz(-1.0988702) q[2];
sx q[2];
rz(-0.79052351) q[2];
sx q[2];
rz(-0.84004842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8255303) q[1];
sx q[1];
rz(-1.9426553) q[1];
sx q[1];
rz(-1.6968246) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7322147) q[3];
sx q[3];
rz(-2.1521849) q[3];
sx q[3];
rz(2.5323109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(1.178297) q[2];
rz(-0.96238771) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(-0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(1.6145561) q[0];
rz(-2.4987192) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(2.8082074) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51334914) q[0];
sx q[0];
rz(-1.7772563) q[0];
sx q[0];
rz(0.87135656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9579499) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(1.7533592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15409878) q[1];
sx q[1];
rz(-2.0425659) q[1];
sx q[1];
rz(-3.0497453) q[1];
rz(-pi) q[2];
rz(-2.8414855) q[3];
sx q[3];
rz(-2.4862053) q[3];
sx q[3];
rz(2.3902077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(2.3392759) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-0.10087092) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-2.5774082) q[0];
rz(-2.5634649) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(2.633458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0429338) q[0];
sx q[0];
rz(-1.9693976) q[0];
sx q[0];
rz(-1.1111141) q[0];
rz(-0.63799413) q[2];
sx q[2];
rz(-1.5882512) q[2];
sx q[2];
rz(0.88776112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.21081884) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(-0.37097431) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1448137) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(0.60861482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4529139) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(2.5417035) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.6413123) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7552345) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(-2.9454943) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-1.0713779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51273726) q[0];
sx q[0];
rz(-2.3561986) q[0];
sx q[0];
rz(0.94731713) q[0];
x q[1];
rz(-1.1159665) q[2];
sx q[2];
rz(-0.78274667) q[2];
sx q[2];
rz(2.6775529) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86712671) q[1];
sx q[1];
rz(-1.1536088) q[1];
sx q[1];
rz(-0.31068128) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7171633) q[3];
sx q[3];
rz(-1.9924581) q[3];
sx q[3];
rz(1.9205586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(0.18520959) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(1.265032) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(0.53034267) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(-0.21662724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075561698) q[0];
sx q[0];
rz(-0.78853411) q[0];
sx q[0];
rz(1.9476452) q[0];
rz(-1.8703307) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.2671721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10303282) q[1];
sx q[1];
rz(-1.9060887) q[1];
sx q[1];
rz(0.0086987728) q[1];
x q[2];
rz(1.3304177) q[3];
sx q[3];
rz(-0.68475311) q[3];
sx q[3];
rz(1.7184337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(-0.8762382) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(-0.52893692) q[0];
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
rz(2.3952336) q[0];
sx q[0];
rz(-0.2975718) q[0];
sx q[0];
rz(1.5901106) q[0];
rz(0.5577001) q[2];
sx q[2];
rz(-2.2039362) q[2];
sx q[2];
rz(1.5483088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72291259) q[1];
sx q[1];
rz(-1.5376248) q[1];
sx q[1];
rz(-1.564371) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0508556) q[3];
sx q[3];
rz(-1.6999348) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65687031) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(0.067226974) q[1];
sx q[1];
rz(-2.1104689) q[1];
sx q[1];
rz(2.9464088) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0474437) q[0];
sx q[0];
rz(-2.7820245) q[0];
sx q[0];
rz(-0.6028428) q[0];
x q[1];
rz(-1.5567661) q[2];
sx q[2];
rz(-1.9719035) q[2];
sx q[2];
rz(-0.48497981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6703549) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(-3.0958789) q[1];
x q[2];
rz(-0.94594749) q[3];
sx q[3];
rz(-0.56560707) q[3];
sx q[3];
rz(0.92029508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8937257) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(2.2360738) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72682196) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(-2.1006405) q[0];
rz(0.078661593) q[1];
sx q[1];
rz(-2.9610596) q[1];
sx q[1];
rz(-0.35531607) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5268363) q[0];
sx q[0];
rz(-2.7608681) q[0];
sx q[0];
rz(-1.3044796) q[0];
x q[1];
rz(1.0742498) q[2];
sx q[2];
rz(-1.3585919) q[2];
sx q[2];
rz(-2.4339536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73341093) q[1];
sx q[1];
rz(-1.2473885) q[1];
sx q[1];
rz(-0.059493382) q[1];
rz(-pi) q[2];
rz(2.5411685) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(-1.7379649) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(-1.6041554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733317) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(-1.2982752) q[0];
rz(-0.4459933) q[1];
sx q[1];
rz(-1.1152209) q[1];
sx q[1];
rz(-2.840852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83082643) q[0];
sx q[0];
rz(-2.5529773) q[0];
sx q[0];
rz(1.1525843) q[0];
rz(-pi) q[1];
rz(-2.5955574) q[2];
sx q[2];
rz(-1.2166426) q[2];
sx q[2];
rz(1.6500157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3831331) q[1];
sx q[1];
rz(-2.8499095) q[1];
sx q[1];
rz(3.0848857) q[1];
rz(0.41348503) q[3];
sx q[3];
rz(-2.1889912) q[3];
sx q[3];
rz(1.2234883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(1.1395617) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6931077) q[0];
sx q[0];
rz(-1.2379452) q[0];
sx q[0];
rz(-2.2647279) q[0];
rz(1.4032455) q[1];
sx q[1];
rz(-1.9156024) q[1];
sx q[1];
rz(1.8617873) q[1];
rz(0.030056603) q[2];
sx q[2];
rz(-0.056014225) q[2];
sx q[2];
rz(0.87607486) q[2];
rz(2.1915477) q[3];
sx q[3];
rz(-1.786357) q[3];
sx q[3];
rz(0.0035567457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
