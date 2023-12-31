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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863659) q[0];
sx q[0];
rz(-0.47954924) q[0];
sx q[0];
rz(-0.064155302) q[0];
rz(1.2650752) q[2];
sx q[2];
rz(-2.6439878) q[2];
sx q[2];
rz(-1.9203609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3996934) q[1];
sx q[1];
rz(-2.0772935) q[1];
sx q[1];
rz(-2.6687117) q[1];
rz(0.17920223) q[3];
sx q[3];
rz(-1.7342907) q[3];
sx q[3];
rz(1.6678099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6291818) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(-0.50981057) q[3];
sx q[3];
rz(-1.9911659) q[3];
sx q[3];
rz(-1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.3614685) q[0];
sx q[0];
rz(-1.7391917) q[0];
sx q[0];
rz(2.0289039) q[0];
rz(0.15377046) q[1];
sx q[1];
rz(-0.91111168) q[1];
sx q[1];
rz(1.3382834) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3801549) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(1.820691) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0427225) q[2];
sx q[2];
rz(-2.3510691) q[2];
sx q[2];
rz(2.3015442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3160623) q[1];
sx q[1];
rz(-1.9426553) q[1];
sx q[1];
rz(1.4447681) q[1];
x q[2];
rz(-1.7322147) q[3];
sx q[3];
rz(-2.1521849) q[3];
sx q[3];
rz(-0.60928173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.082836941) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(1.9632957) q[2];
rz(-0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96175471) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(1.5270365) q[0];
rz(2.4987192) q[1];
sx q[1];
rz(-2.0650654) q[1];
sx q[1];
rz(0.33338526) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81844717) q[0];
sx q[0];
rz(-2.4172757) q[0];
sx q[0];
rz(-1.8853272) q[0];
x q[1];
rz(-1.2522167) q[2];
sx q[2];
rz(-2.736056) q[2];
sx q[2];
rz(-0.47682724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15409878) q[1];
sx q[1];
rz(-2.0425659) q[1];
sx q[1];
rz(-3.0497453) q[1];
rz(1.3473347) q[3];
sx q[3];
rz(-0.94933214) q[3];
sx q[3];
rz(0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0201515) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(2.3392759) q[2];
rz(0.24117593) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-3.0407217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0872831) q[0];
sx q[0];
rz(-1.6182951) q[0];
sx q[0];
rz(-0.56418443) q[0];
rz(-0.57812771) q[1];
sx q[1];
rz(-1.6491978) q[1];
sx q[1];
rz(-2.633458) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.479624) q[0];
sx q[0];
rz(-1.9920252) q[0];
sx q[0];
rz(-0.43930102) q[0];
rz(2.5035985) q[2];
sx q[2];
rz(-1.5533414) q[2];
sx q[2];
rz(2.2538315) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6259039) q[1];
sx q[1];
rz(-1.8325894) q[1];
sx q[1];
rz(0.7598676) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9594801) q[3];
sx q[3];
rz(-1.2977227) q[3];
sx q[3];
rz(1.207721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-0.33334175) q[2];
sx q[2];
rz(-2.508146) q[2];
rz(2.5417035) q[3];
sx q[3];
rz(-1.1497295) q[3];
sx q[3];
rz(-1.6413123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7552345) q[0];
sx q[0];
rz(-0.30325493) q[0];
sx q[0];
rz(2.9454943) q[0];
rz(-1.261699) q[1];
sx q[1];
rz(-2.321107) q[1];
sx q[1];
rz(2.0702147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3065942) q[0];
sx q[0];
rz(-2.1822565) q[0];
sx q[0];
rz(-0.5284662) q[0];
rz(-pi) q[1];
rz(-2.3000556) q[2];
sx q[2];
rz(-1.255799) q[2];
sx q[2];
rz(2.3685761) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6043678) q[1];
sx q[1];
rz(-2.62694) q[1];
sx q[1];
rz(0.96697076) q[1];
x q[2];
rz(-0.82810546) q[3];
sx q[3];
rz(-2.5525186) q[3];
sx q[3];
rz(-2.7554054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1065958) q[2];
sx q[2];
rz(-1.3369766) q[2];
sx q[2];
rz(-2.1595947) q[2];
rz(-0.18520959) q[3];
sx q[3];
rz(-2.2976112) q[3];
sx q[3];
rz(-1.8765607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035325) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(0.53034267) q[0];
rz(1.8416587) q[1];
sx q[1];
rz(-1.8118186) q[1];
sx q[1];
rz(2.9249654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075561698) q[0];
sx q[0];
rz(-2.3530585) q[0];
sx q[0];
rz(-1.9476452) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2712619) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.2671721) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10303282) q[1];
sx q[1];
rz(-1.9060887) q[1];
sx q[1];
rz(3.1328939) q[1];
rz(-pi) q[2];
rz(1.8111749) q[3];
sx q[3];
rz(-0.68475311) q[3];
sx q[3];
rz(1.4231589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4914322) q[2];
sx q[2];
rz(-1.5101134) q[2];
sx q[2];
rz(2.1177297) q[2];
rz(0.8762382) q[3];
sx q[3];
rz(-2.439308) q[3];
sx q[3];
rz(-1.2715626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-1.9493608) q[1];
sx q[1];
rz(1.0891917) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74635909) q[0];
sx q[0];
rz(-0.2975718) q[0];
sx q[0];
rz(-1.5901106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5838926) q[2];
sx q[2];
rz(-0.93765646) q[2];
sx q[2];
rz(-1.5932839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91428077) q[1];
sx q[1];
rz(-0.033787878) q[1];
sx q[1];
rz(2.950331) q[1];
rz(-0.96157162) q[3];
sx q[3];
rz(-0.15768356) q[3];
sx q[3];
rz(-0.79573436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0163394) q[2];
sx q[2];
rz(-2.3374127) q[2];
sx q[2];
rz(1.2255229) q[2];
rz(-1.4922173) q[3];
sx q[3];
rz(-1.4368613) q[3];
sx q[3];
rz(-2.9857181) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(2.9464088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0941489) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(2.5387499) q[0];
rz(-pi) q[1];
x q[1];
rz(0.033069177) q[2];
sx q[2];
rz(-0.40133921) q[2];
sx q[2];
rz(2.6925342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47123779) q[1];
sx q[1];
rz(-2.0002623) q[1];
sx q[1];
rz(3.0958789) q[1];
rz(-pi) q[2];
rz(-0.35555367) q[3];
sx q[3];
rz(-2.0204633) q[3];
sx q[3];
rz(-0.21330968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.247867) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-2.2360738) q[2];
rz(-1.9838105) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(2.982443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
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
rz(-2.9610596) q[1];
sx q[1];
rz(-0.35531607) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2040877) q[0];
sx q[0];
rz(-1.6687487) q[0];
sx q[0];
rz(1.9393001) q[0];
x q[1];
rz(2.9012868) q[2];
sx q[2];
rz(-2.0552286) q[2];
sx q[2];
rz(-0.74953178) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.81845835) q[1];
sx q[1];
rz(-1.5143906) q[1];
sx q[1];
rz(-1.2468546) q[1];
rz(-0.60042419) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(1.0004117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4650402) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(-1.4036277) q[2];
rz(0.0062395652) q[3];
sx q[3];
rz(-1.3429567) q[3];
sx q[3];
rz(1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.1152209) q[1];
sx q[1];
rz(2.840852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7556831) q[0];
sx q[0];
rz(-1.3433546) q[0];
sx q[0];
rz(2.118551) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9791331) q[2];
sx q[2];
rz(-2.0795341) q[2];
sx q[2];
rz(-3.0131154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4423351) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(-1.5878116) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0565287) q[3];
sx q[3];
rz(-0.72838718) q[3];
sx q[3];
rz(-1.8715093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5131502) q[2];
sx q[2];
rz(-1.0978038) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(-1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.9679507) q[3];
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
rz(pi/2) q[2];
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
rz(0.055988978) q[2];
sx q[2];
rz(-1.5691139) q[2];
sx q[2];
rz(-0.72473095) q[2];
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
