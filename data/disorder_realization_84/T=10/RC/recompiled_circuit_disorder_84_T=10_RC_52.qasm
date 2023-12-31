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
rz(-2.3513878) q[1];
sx q[1];
rz(-1.9145929) q[1];
sx q[1];
rz(-1.9807281) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863659) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(-0.064155302) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2650752) q[2];
sx q[2];
rz(-2.6439878) q[2];
sx q[2];
rz(1.9203609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4143715) q[1];
sx q[1];
rz(-1.9804269) q[1];
sx q[1];
rz(-1.0135256) q[1];
x q[2];
rz(-2.9623904) q[3];
sx q[3];
rz(-1.7342907) q[3];
sx q[3];
rz(-1.4737827) q[3];
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
rz(1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(-1.1126888) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(-1.8033093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1424944) q[0];
sx q[0];
rz(-2.4141443) q[0];
sx q[0];
rz(-0.29074685) q[0];
rz(-0.43054994) q[2];
sx q[2];
rz(-0.88532788) q[2];
sx q[2];
rz(1.6738883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2087304) q[1];
sx q[1];
rz(-1.6881697) q[1];
sx q[1];
rz(2.7670303) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.409378) q[3];
sx q[3];
rz(-2.1521849) q[3];
sx q[3];
rz(-2.5323109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0587557) q[2];
sx q[2];
rz(-0.78410167) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(-2.1792049) q[3];
sx q[3];
rz(-1.0932086) q[3];
sx q[3];
rz(0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96175471) q[0];
sx q[0];
rz(-0.50757718) q[0];
sx q[0];
rz(-1.6145561) q[0];
rz(-2.4987192) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(2.8082074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3231455) q[0];
sx q[0];
rz(-2.4172757) q[0];
sx q[0];
rz(-1.2562654) q[0];
rz(-pi) q[1];
rz(1.9579499) q[2];
sx q[2];
rz(-1.4469115) q[2];
sx q[2];
rz(-1.3882335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4585321) q[1];
sx q[1];
rz(-1.6525869) q[1];
sx q[1];
rz(2.044278) q[1];
x q[2];
rz(1.3473347) q[3];
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
rz(0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872831) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-2.5774082) q[0];
rz(0.57812771) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(2.633458) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0429338) q[0];
sx q[0];
rz(-1.9693976) q[0];
sx q[0];
rz(2.0304785) q[0];
x q[1];
rz(-0.63799413) q[2];
sx q[2];
rz(-1.5882512) q[2];
sx q[2];
rz(-2.2538315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8453464) q[1];
sx q[1];
rz(-2.2987662) q[1];
sx q[1];
rz(-1.216757) q[1];
x q[2];
rz(-0.99677892) q[3];
sx q[3];
rz(-0.32696163) q[3];
sx q[3];
rz(2.5329778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(0.63344669) q[2];
rz(-2.5417035) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(-0.19609837) q[0];
rz(1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(2.0702147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58760658) q[0];
sx q[0];
rz(-1.1452132) q[0];
sx q[0];
rz(-2.2527184) q[0];
rz(2.0256261) q[2];
sx q[2];
rz(-0.78274667) q[2];
sx q[2];
rz(-0.46403971) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57429806) q[1];
sx q[1];
rz(-1.287536) q[1];
sx q[1];
rz(1.1351372) q[1];
rz(-pi) q[2];
rz(-2.7171633) q[3];
sx q[3];
rz(-1.9924581) q[3];
sx q[3];
rz(-1.2210341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035325) q[0];
sx q[0];
rz(-0.91941994) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(-1.2999339) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(-2.9249654) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43603555) q[0];
sx q[0];
rz(-2.2910121) q[0];
sx q[0];
rz(-0.35465045) q[0];
x q[1];
rz(1.2712619) q[2];
sx q[2];
rz(-2.1795863) q[2];
sx q[2];
rz(-1.8744206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0385598) q[1];
sx q[1];
rz(-1.9060887) q[1];
sx q[1];
rz(3.1328939) q[1];
rz(-pi) q[2];
rz(-2.9495903) q[3];
sx q[3];
rz(-2.2322906) q[3];
sx q[3];
rz(-2.0249174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.65016046) q[2];
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
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-1.9493608) q[1];
sx q[1];
rz(-1.0891917) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80597164) q[0];
sx q[0];
rz(-1.5651337) q[0];
sx q[0];
rz(1.2732768) q[0];
rz(-pi) q[1];
rz(-2.2839374) q[2];
sx q[2];
rz(-2.0115888) q[2];
sx q[2];
rz(-2.8105274) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72291259) q[1];
sx q[1];
rz(-1.6039679) q[1];
sx q[1];
rz(1.564371) q[1];
x q[2];
rz(1.4411304) q[3];
sx q[3];
rz(-1.6607758) q[3];
sx q[3];
rz(0.17168301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0163394) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(-1.9160697) q[2];
rz(1.6493753) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(2.9857181) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-0.06198922) q[0];
sx q[0];
rz(2.2739676) q[0];
rz(-3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(-2.9464088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.09571) q[0];
sx q[0];
rz(-1.3699431) q[0];
sx q[0];
rz(-0.30028371) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1085235) q[2];
sx q[2];
rz(-2.7402534) q[2];
sx q[2];
rz(0.44905845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7797797) q[1];
sx q[1];
rz(-0.43174141) q[1];
sx q[1];
rz(1.4713431) q[1];
rz(-pi) q[2];
rz(-0.94594749) q[3];
sx q[3];
rz(-0.56560707) q[3];
sx q[3];
rz(0.92029508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.247867) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(-1.1577822) q[3];
sx q[3];
rz(-1.9730622) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(-2.1006405) q[0];
rz(0.078661593) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(0.35531607) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5268363) q[0];
sx q[0];
rz(-0.38072452) q[0];
sx q[0];
rz(-1.3044796) q[0];
x q[1];
rz(1.1460733) q[2];
sx q[2];
rz(-2.6051084) q[2];
sx q[2];
rz(1.2338961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73341093) q[1];
sx q[1];
rz(-1.2473885) q[1];
sx q[1];
rz(-0.059493382) q[1];
rz(1.7570417) q[3];
sx q[3];
rz(-0.97847853) q[3];
sx q[3];
rz(-2.4663962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-2.5341454) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(-0.0062395652) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-1.2982752) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(-0.30074063) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38590955) q[0];
sx q[0];
rz(-1.798238) q[0];
sx q[0];
rz(-1.0230416) q[0];
rz(-pi) q[1];
rz(-0.54603521) q[2];
sx q[2];
rz(-1.2166426) q[2];
sx q[2];
rz(-1.6500157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69925752) q[1];
sx q[1];
rz(-1.2795957) q[1];
sx q[1];
rz(1.553781) q[1];
x q[2];
rz(2.0850639) q[3];
sx q[3];
rz(-2.4132055) q[3];
sx q[3];
rz(-1.8715093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5131502) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(1.1395617) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-2.5438178) q[3];
sx q[3];
rz(1.9679507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-1.4032455) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(-1.5724814) q[2];
sx q[2];
rz(-1.5148074) q[2];
sx q[2];
rz(-2.2956216) q[2];
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
