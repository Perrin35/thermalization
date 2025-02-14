OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.46718207) q[0];
sx q[0];
rz(-1.2527569) q[0];
sx q[0];
rz(7.0897515) q[0];
rz(-0.84713495) q[1];
sx q[1];
rz(-0.48589125) q[1];
sx q[1];
rz(-2.2810305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63634625) q[0];
sx q[0];
rz(-2.5274694) q[0];
sx q[0];
rz(-1.3259802) q[0];
x q[1];
rz(-1.2803308) q[2];
sx q[2];
rz(-0.72735255) q[2];
sx q[2];
rz(-2.1216105) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5845421) q[1];
sx q[1];
rz(-1.8452541) q[1];
sx q[1];
rz(1.518979) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99765649) q[3];
sx q[3];
rz(-1.9966148) q[3];
sx q[3];
rz(-0.65090685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2412771) q[2];
sx q[2];
rz(-0.52854717) q[2];
sx q[2];
rz(2.6701374) q[2];
rz(2.6491162) q[3];
sx q[3];
rz(-1.2329279) q[3];
sx q[3];
rz(-1.2537778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-0.27552283) q[0];
sx q[0];
rz(-2.7560784) q[0];
sx q[0];
rz(-2.8634014) q[0];
rz(-1.1909852) q[1];
sx q[1];
rz(-1.8089801) q[1];
sx q[1];
rz(1.8461548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14034056) q[0];
sx q[0];
rz(-2.3352726) q[0];
sx q[0];
rz(-2.4693523) q[0];
rz(-pi) q[1];
rz(-1.5191742) q[2];
sx q[2];
rz(-1.6502909) q[2];
sx q[2];
rz(2.5148066) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5767314) q[1];
sx q[1];
rz(-2.7612491) q[1];
sx q[1];
rz(1.5768135) q[1];
x q[2];
rz(0.94280394) q[3];
sx q[3];
rz(-2.0217388) q[3];
sx q[3];
rz(1.0058927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.94464716) q[2];
sx q[2];
rz(-1.6702007) q[2];
sx q[2];
rz(2.6317281) q[2];
rz(-1.4465796) q[3];
sx q[3];
rz(-2.6810985) q[3];
sx q[3];
rz(-0.50486008) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6079717) q[0];
sx q[0];
rz(-1.6225659) q[0];
sx q[0];
rz(-0.98440379) q[0];
rz(-3.1254752) q[1];
sx q[1];
rz(-2.0033483) q[1];
sx q[1];
rz(1.3498397) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4181217) q[0];
sx q[0];
rz(-1.7364565) q[0];
sx q[0];
rz(0.70933527) q[0];
rz(-pi) q[1];
rz(1.0060293) q[2];
sx q[2];
rz(-1.2899219) q[2];
sx q[2];
rz(1.8233521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78245263) q[1];
sx q[1];
rz(-1.5305007) q[1];
sx q[1];
rz(-2.4423679) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4593342) q[3];
sx q[3];
rz(-2.3334399) q[3];
sx q[3];
rz(0.34592849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71780378) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(-0.88199893) q[2];
rz(-1.1841904) q[3];
sx q[3];
rz(-2.2004674) q[3];
sx q[3];
rz(-2.4052896) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475567) q[0];
sx q[0];
rz(-1.9147669) q[0];
sx q[0];
rz(3.0645698) q[0];
rz(2.9432964) q[1];
sx q[1];
rz(-1.501333) q[1];
sx q[1];
rz(2.95453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27320751) q[0];
sx q[0];
rz(-1.7284231) q[0];
sx q[0];
rz(-2.8013743) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0346626) q[2];
sx q[2];
rz(-0.38883801) q[2];
sx q[2];
rz(-1.0146245) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9976207) q[1];
sx q[1];
rz(-1.5117497) q[1];
sx q[1];
rz(-0.70934341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8238769) q[3];
sx q[3];
rz(-2.6926929) q[3];
sx q[3];
rz(1.0300058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5932811) q[2];
sx q[2];
rz(-0.45479861) q[2];
sx q[2];
rz(-1.456267) q[2];
rz(2.4233387) q[3];
sx q[3];
rz(-1.7536438) q[3];
sx q[3];
rz(2.3072402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2497571) q[0];
sx q[0];
rz(-1.8064073) q[0];
sx q[0];
rz(0.43858132) q[0];
rz(-2.7843685) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(-1.2507778) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27173247) q[0];
sx q[0];
rz(-1.9969059) q[0];
sx q[0];
rz(-0.10978384) q[0];
x q[1];
rz(1.3322796) q[2];
sx q[2];
rz(-2.8217689) q[2];
sx q[2];
rz(0.25898283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73356556) q[1];
sx q[1];
rz(-1.1981315) q[1];
sx q[1];
rz(1.1421475) q[1];
rz(-pi) q[2];
rz(-2.6762371) q[3];
sx q[3];
rz(-2.8520003) q[3];
sx q[3];
rz(-2.6526439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9095416) q[2];
sx q[2];
rz(-0.48018685) q[2];
sx q[2];
rz(1.5117744) q[2];
rz(0.016228598) q[3];
sx q[3];
rz(-1.0954233) q[3];
sx q[3];
rz(0.79286638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69415724) q[0];
sx q[0];
rz(-1.2801535) q[0];
sx q[0];
rz(1.9919027) q[0];
rz(-0.42090526) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(3.0973869) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6290508) q[0];
sx q[0];
rz(-0.78475941) q[0];
sx q[0];
rz(2.6757338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.030613) q[2];
sx q[2];
rz(-1.1655131) q[2];
sx q[2];
rz(2.0342397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.845712) q[1];
sx q[1];
rz(-2.7428877) q[1];
sx q[1];
rz(-1.7804342) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8474691) q[3];
sx q[3];
rz(-1.0545316) q[3];
sx q[3];
rz(0.54126213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9629024) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(2.1964591) q[2];
rz(-1.4611698) q[3];
sx q[3];
rz(-1.1472568) q[3];
sx q[3];
rz(0.51819658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1951676) q[0];
sx q[0];
rz(-0.31038809) q[0];
sx q[0];
rz(-0.84939605) q[0];
rz(-1.7646344) q[1];
sx q[1];
rz(-2.5701249) q[1];
sx q[1];
rz(1.2845385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234069) q[0];
sx q[0];
rz(-2.0571097) q[0];
sx q[0];
rz(1.125) q[0];
rz(-pi) q[1];
rz(-0.30751245) q[2];
sx q[2];
rz(-1.4685681) q[2];
sx q[2];
rz(-2.6999465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.095605222) q[1];
sx q[1];
rz(-0.52667499) q[1];
sx q[1];
rz(0.070835872) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0128117) q[3];
sx q[3];
rz(-1.9922755) q[3];
sx q[3];
rz(-0.95343219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4756502) q[2];
sx q[2];
rz(-1.3694171) q[2];
sx q[2];
rz(-1.5671889) q[2];
rz(2.808908) q[3];
sx q[3];
rz(-1.0654457) q[3];
sx q[3];
rz(-2.1596597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629267) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(-1.0885619) q[0];
rz(2.2125878) q[1];
sx q[1];
rz(-1.4938415) q[1];
sx q[1];
rz(0.30379024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740363) q[0];
sx q[0];
rz(-2.1308194) q[0];
sx q[0];
rz(-2.014888) q[0];
x q[1];
rz(2.1985198) q[2];
sx q[2];
rz(-2.94063) q[2];
sx q[2];
rz(-2.497884) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15081295) q[1];
sx q[1];
rz(-2.7013198) q[1];
sx q[1];
rz(-1.3516462) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61875069) q[3];
sx q[3];
rz(-1.9721974) q[3];
sx q[3];
rz(1.5681745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8888237) q[2];
sx q[2];
rz(-1.4759651) q[2];
sx q[2];
rz(-0.26407537) q[2];
rz(1.1325599) q[3];
sx q[3];
rz(-0.3370291) q[3];
sx q[3];
rz(-1.0549217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3573414) q[0];
sx q[0];
rz(-0.34264523) q[0];
sx q[0];
rz(1.0145048) q[0];
rz(-2.2134589) q[1];
sx q[1];
rz(-1.721563) q[1];
sx q[1];
rz(2.6511505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3532012) q[0];
sx q[0];
rz(-1.4239131) q[0];
sx q[0];
rz(0.066195458) q[0];
x q[1];
rz(1.956067) q[2];
sx q[2];
rz(-1.2495572) q[2];
sx q[2];
rz(2.3209077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52940581) q[1];
sx q[1];
rz(-1.8711963) q[1];
sx q[1];
rz(3.0985188) q[1];
rz(2.9005592) q[3];
sx q[3];
rz(-2.0969982) q[3];
sx q[3];
rz(1.0722425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.32533112) q[2];
sx q[2];
rz(-2.0897431) q[2];
sx q[2];
rz(-0.10247792) q[2];
rz(-1.2517733) q[3];
sx q[3];
rz(-1.3636369) q[3];
sx q[3];
rz(-1.4832835) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74806279) q[0];
sx q[0];
rz(-3.0952125) q[0];
sx q[0];
rz(1.8512132) q[0];
rz(0.45404592) q[1];
sx q[1];
rz(-1.3244649) q[1];
sx q[1];
rz(2.9581199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35724872) q[0];
sx q[0];
rz(-1.5518477) q[0];
sx q[0];
rz(1.3878893) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.543599) q[2];
sx q[2];
rz(-1.5584297) q[2];
sx q[2];
rz(1.7747161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.390229) q[1];
sx q[1];
rz(-1.5469264) q[1];
sx q[1];
rz(0.49225251) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7137035) q[3];
sx q[3];
rz(-0.57051728) q[3];
sx q[3];
rz(2.4874668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7196677) q[2];
sx q[2];
rz(-2.0823961) q[2];
sx q[2];
rz(-0.025040778) q[2];
rz(-2.5981564) q[3];
sx q[3];
rz(-2.4380324) q[3];
sx q[3];
rz(0.27633015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80171361) q[0];
sx q[0];
rz(-1.0152974) q[0];
sx q[0];
rz(0.80831084) q[0];
rz(0.33357757) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(-1.2407718) q[2];
sx q[2];
rz(-1.6920964) q[2];
sx q[2];
rz(1.7597711) q[2];
rz(2.2441545) q[3];
sx q[3];
rz(-1.031395) q[3];
sx q[3];
rz(-0.54678834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
