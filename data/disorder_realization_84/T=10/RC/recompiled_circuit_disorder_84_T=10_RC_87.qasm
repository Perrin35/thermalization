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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5829585) q[0];
sx q[0];
rz(-1.5412124) q[0];
sx q[0];
rz(-2.6628859) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9795346) q[2];
sx q[2];
rz(-1.098212) q[2];
sx q[2];
rz(-2.2651847) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3996934) q[1];
sx q[1];
rz(-1.0642991) q[1];
sx q[1];
rz(0.472881) q[1];
x q[2];
rz(1.4046895) q[3];
sx q[3];
rz(-1.3940092) q[3];
sx q[3];
rz(3.0151031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.2966195) q[2];
sx q[2];
rz(-2.4751002) q[2];
rz(-2.6317821) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(-1.8681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(-1.1126888) q[0];
rz(2.9878222) q[1];
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
rz(-0.35182686) q[0];
sx q[0];
rz(-1.3789982) q[0];
sx q[0];
rz(2.4353794) q[0];
x q[1];
rz(0.83805214) q[2];
sx q[2];
rz(-1.2417972) q[2];
sx q[2];
rz(2.7555562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4900134) q[1];
sx q[1];
rz(-0.39169185) q[1];
sx q[1];
rz(-0.3117805) q[1];
rz(-pi) q[2];
rz(-1.7322147) q[3];
sx q[3];
rz(-0.98940778) q[3];
sx q[3];
rz(-2.5323109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.082836941) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(1.178297) q[2];
rz(-2.1792049) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(2.4760831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(1.5270365) q[0];
rz(-2.4987192) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-2.8082074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51334914) q[0];
sx q[0];
rz(-1.3643364) q[0];
sx q[0];
rz(2.2702361) q[0];
rz(-1.1836428) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(-1.3882335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9874939) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(-0.091847329) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3473347) q[3];
sx q[3];
rz(-2.1922605) q[3];
sx q[3];
rz(-0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(0.80231673) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-2.4427588) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(2.633458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80721891) q[0];
sx q[0];
rz(-0.59893543) q[0];
sx q[0];
rz(2.3301624) q[0];
rz(-0.63799413) q[2];
sx q[2];
rz(-1.5882512) q[2];
sx q[2];
rz(-2.2538315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5156887) q[1];
sx q[1];
rz(-1.3090033) q[1];
sx q[1];
rz(2.3817251) q[1];
rz(-pi) q[2];
rz(2.1448137) q[3];
sx q[3];
rz(-2.814631) q[3];
sx q[3];
rz(-2.5329778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4529139) q[2];
sx q[2];
rz(-2.8082509) q[2];
sx q[2];
rz(-0.63344669) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(0.19609837) q[0];
rz(-1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(1.0713779) q[1];
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
x q[1];
rz(0.41197889) q[2];
sx q[2];
rz(-2.2569071) q[2];
sx q[2];
rz(2.0737322) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2744659) q[1];
sx q[1];
rz(-1.1536088) q[1];
sx q[1];
rz(-2.8309114) q[1];
rz(-2.0282182) q[3];
sx q[3];
rz(-1.955982) q[3];
sx q[3];
rz(0.53264632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(-0.98199797) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-2.2976112) q[3];
sx q[3];
rz(1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035325) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(-0.53034267) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(2.9249654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075561698) q[0];
sx q[0];
rz(-2.3530585) q[0];
sx q[0];
rz(-1.1939474) q[0];
rz(-1.8703307) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.2671721) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12946373) q[1];
sx q[1];
rz(-2.8061917) q[1];
sx q[1];
rz(-1.5957551) q[1];
x q[2];
rz(2.9495903) q[3];
sx q[3];
rz(-0.90930206) q[3];
sx q[3];
rz(-2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(-2.1177297) q[2];
rz(-2.2653545) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(-1.2715626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4890471) q[0];
sx q[0];
rz(-1.1627731) q[0];
sx q[0];
rz(-0.52893692) q[0];
rz(1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(-1.0891917) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3750317) q[0];
sx q[0];
rz(-1.2732817) q[0];
sx q[0];
rz(-0.0059228063) q[0];
rz(-0.85765526) q[2];
sx q[2];
rz(-1.1300039) q[2];
sx q[2];
rz(0.33106523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91428077) q[1];
sx q[1];
rz(-3.1078048) q[1];
sx q[1];
rz(0.19126161) q[1];
rz(-0.96157162) q[3];
sx q[3];
rz(-0.15768356) q[3];
sx q[3];
rz(-0.79573436) q[3];
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
rz(-1.2255229) q[2];
rz(-1.4922173) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(2.9857181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65687031) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-0.86762506) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(2.9464088) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41338838) q[0];
sx q[0];
rz(-1.2767316) q[0];
sx q[0];
rz(1.3608027) q[0];
rz(-3.1085235) q[2];
sx q[2];
rz(-0.40133921) q[2];
sx q[2];
rz(-0.44905845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47123779) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(0.045713748) q[1];
rz(-pi) q[2];
x q[2];
rz(2.04625) q[3];
sx q[3];
rz(-1.2519149) q[3];
sx q[3];
rz(1.9441324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8937257) q[2];
sx q[2];
rz(-2.9276431) q[2];
sx q[2];
rz(2.2360738) q[2];
rz(1.9838105) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.18053308) q[1];
sx q[1];
rz(2.7862766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61475635) q[0];
sx q[0];
rz(-2.7608681) q[0];
sx q[0];
rz(-1.8371131) q[0];
rz(2.9012868) q[2];
sx q[2];
rz(-2.0552286) q[2];
sx q[2];
rz(2.3920609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4081817) q[1];
sx q[1];
rz(-1.2473885) q[1];
sx q[1];
rz(0.059493382) q[1];
x q[2];
rz(-2.8730632) q[3];
sx q[3];
rz(-2.5240457) q[3];
sx q[3];
rz(-2.7919046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.7379649) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
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
rz(-1.1152209) q[1];
sx q[1];
rz(0.30074063) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3215799) q[0];
sx q[0];
rz(-2.1029148) q[0];
sx q[0];
rz(0.26474712) q[0];
rz(-pi) q[1];
rz(-2.5228595) q[2];
sx q[2];
rz(-2.500687) q[2];
sx q[2];
rz(-0.59782019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2749394) q[1];
sx q[1];
rz(-1.5544974) q[1];
sx q[1];
rz(-0.29124041) q[1];
rz(2.7281076) q[3];
sx q[3];
rz(-2.1889912) q[3];
sx q[3];
rz(-1.2234883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62844244) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(1.9679507) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-0.030056603) q[2];
sx q[2];
rz(-3.0855784) q[2];
sx q[2];
rz(-2.2655178) q[2];
rz(2.8786447) q[3];
sx q[3];
rz(-2.1750952) q[3];
sx q[3];
rz(1.7261214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];