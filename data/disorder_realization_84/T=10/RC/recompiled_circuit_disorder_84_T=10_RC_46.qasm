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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55863419) q[0];
sx q[0];
rz(-1.6003803) q[0];
sx q[0];
rz(-0.47870676) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16205807) q[2];
sx q[2];
rz(-2.0433807) q[2];
sx q[2];
rz(0.87640793) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4143715) q[1];
sx q[1];
rz(-1.1611658) q[1];
sx q[1];
rz(2.1280671) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4046895) q[3];
sx q[3];
rz(-1.7475834) q[3];
sx q[3];
rz(0.12648957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6291818) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(2.4751002) q[2];
rz(2.6317821) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(2.1424944) q[0];
sx q[0];
rz(-0.72744838) q[0];
sx q[0];
rz(-2.8508458) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43054994) q[2];
sx q[2];
rz(-2.2562648) q[2];
sx q[2];
rz(-1.4677043) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8255303) q[1];
sx q[1];
rz(-1.1989374) q[1];
sx q[1];
rz(1.4447681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.409378) q[3];
sx q[3];
rz(-2.1521849) q[3];
sx q[3];
rz(-0.60928173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(1.178297) q[2];
rz(2.1792049) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(-1.6145561) q[0];
rz(0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(0.33338526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51334914) q[0];
sx q[0];
rz(-1.7772563) q[0];
sx q[0];
rz(0.87135656) q[0];
x q[1];
rz(-0.1336735) q[2];
sx q[2];
rz(-1.954827) q[2];
sx q[2];
rz(-0.13222279) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4585321) q[1];
sx q[1];
rz(-1.6525869) q[1];
sx q[1];
rz(1.0973147) q[1];
rz(-1.7942579) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(2.7623451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.0872831) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(-2.5774082) q[0];
rz(2.5634649) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(2.633458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098658853) q[0];
sx q[0];
rz(-1.9693976) q[0];
sx q[0];
rz(-1.1111141) q[0];
rz(-pi) q[1];
rz(3.1122909) q[2];
sx q[2];
rz(-2.5033931) q[2];
sx q[2];
rz(-0.65949856) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9307738) q[1];
sx q[1];
rz(-0.79508077) q[1];
sx q[1];
rz(2.7706183) q[1];
x q[2];
rz(-1.2933613) q[3];
sx q[3];
rz(-1.7460896) q[3];
sx q[3];
rz(2.7288935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4529139) q[2];
sx q[2];
rz(-2.8082509) q[2];
sx q[2];
rz(-2.508146) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.3863581) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(-0.19609837) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(1.0713779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3065942) q[0];
sx q[0];
rz(-0.95933611) q[0];
sx q[0];
rz(2.6131265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7296138) q[2];
sx q[2];
rz(-0.88468555) q[2];
sx q[2];
rz(-2.0737322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(-2.1746219) q[1];
rz(-pi) q[2];
rz(2.0282182) q[3];
sx q[3];
rz(-1.1856106) q[3];
sx q[3];
rz(-2.6089463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.034996899) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(0.18520959) q[3];
sx q[3];
rz(-0.84398142) q[3];
sx q[3];
rz(-1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075561698) q[0];
sx q[0];
rz(-0.78853411) q[0];
sx q[0];
rz(1.9476452) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40041133) q[2];
sx q[2];
rz(-0.67001221) q[2];
sx q[2];
rz(-2.3695721) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0385598) q[1];
sx q[1];
rz(-1.9060887) q[1];
sx q[1];
rz(-3.1328939) q[1];
rz(2.9495903) q[3];
sx q[3];
rz(-2.2322906) q[3];
sx q[3];
rz(2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65016046) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(1.023863) q[2];
rz(-0.8762382) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(1.8700301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6525456) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(-0.52893692) q[0];
rz(1.6128929) q[1];
sx q[1];
rz(-1.1922319) q[1];
sx q[1];
rz(-1.0891917) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3750317) q[0];
sx q[0];
rz(-1.2732817) q[0];
sx q[0];
rz(-0.0059228063) q[0];
rz(0.5577001) q[2];
sx q[2];
rz(-2.2039362) q[2];
sx q[2];
rz(1.5483088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91428077) q[1];
sx q[1];
rz(-0.033787878) q[1];
sx q[1];
rz(-0.19126161) q[1];
x q[2];
rz(-2.180021) q[3];
sx q[3];
rz(-0.15768356) q[3];
sx q[3];
rz(0.79573436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(-1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(-2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41338838) q[0];
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
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6703549) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(-3.0958789) q[1];
rz(-pi) q[2];
rz(2.04625) q[3];
sx q[3];
rz(-1.8896777) q[3];
sx q[3];
rz(-1.9441324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8937257) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(0.90551886) q[2];
rz(-1.9838105) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(0.15914966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.72682196) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(-1.0409521) q[0];
rz(3.0629311) q[1];
sx q[1];
rz(-2.9610596) q[1];
sx q[1];
rz(0.35531607) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5268363) q[0];
sx q[0];
rz(-2.7608681) q[0];
sx q[0];
rz(1.8371131) q[0];
x q[1];
rz(2.0673429) q[2];
sx q[2];
rz(-1.7830007) q[2];
sx q[2];
rz(-2.4339536) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3231343) q[1];
sx q[1];
rz(-1.5143906) q[1];
sx q[1];
rz(1.2468546) q[1];
rz(-0.60042419) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.4036277) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86826098) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(1.8433174) q[0];
rz(-0.4459933) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(2.840852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38590955) q[0];
sx q[0];
rz(-1.3433546) q[0];
sx q[0];
rz(1.0230416) q[0];
rz(-pi) q[1];
rz(2.5955574) q[2];
sx q[2];
rz(-1.9249501) q[2];
sx q[2];
rz(-1.491577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4423351) q[1];
sx q[1];
rz(-1.2795957) q[1];
sx q[1];
rz(1.5878116) q[1];
x q[2];
rz(1.0565287) q[3];
sx q[3];
rz(-2.4132055) q[3];
sx q[3];
rz(-1.2700833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(2.002031) q[2];
rz(-1.7906174) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.5724814) q[2];
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
