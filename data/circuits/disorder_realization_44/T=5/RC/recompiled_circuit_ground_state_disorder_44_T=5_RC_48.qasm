OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(1.692481) q[0];
sx q[0];
rz(11.115885) q[0];
rz(-2.1679572) q[1];
sx q[1];
rz(-1.4373625) q[1];
sx q[1];
rz(0.91926423) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1184753) q[0];
sx q[0];
rz(-0.035041172) q[0];
sx q[0];
rz(-2.6316597) q[0];
rz(-pi) q[1];
rz(1.6617352) q[2];
sx q[2];
rz(-1.6906066) q[2];
sx q[2];
rz(-0.38209846) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.011292) q[1];
sx q[1];
rz(-0.71349547) q[1];
sx q[1];
rz(-1.7111227) q[1];
rz(-pi) q[2];
rz(-0.84897016) q[3];
sx q[3];
rz(-2.5458286) q[3];
sx q[3];
rz(-1.1536382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8093402) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(1.7285041) q[2];
rz(0.20279065) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(-0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8973812) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(2.3731025) q[1];
sx q[1];
rz(-2.0748506) q[1];
sx q[1];
rz(-2.1220727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8038254) q[0];
sx q[0];
rz(-0.44887421) q[0];
sx q[0];
rz(-0.091218936) q[0];
rz(-pi) q[1];
rz(-1.4065845) q[2];
sx q[2];
rz(-0.43638849) q[2];
sx q[2];
rz(1.1498888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72538917) q[1];
sx q[1];
rz(-0.37230154) q[1];
sx q[1];
rz(-1.4854234) q[1];
x q[2];
rz(-0.083262308) q[3];
sx q[3];
rz(-2.1234649) q[3];
sx q[3];
rz(-0.083831122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76356137) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(-1.9192609) q[2];
rz(-1.9289121) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(1.9236176) q[0];
rz(2.6257264) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2005916) q[0];
sx q[0];
rz(-0.070644826) q[0];
sx q[0];
rz(-2.5289422) q[0];
x q[1];
rz(1.0664682) q[2];
sx q[2];
rz(-2.2309003) q[2];
sx q[2];
rz(0.65048993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.565971) q[1];
sx q[1];
rz(-1.0069443) q[1];
sx q[1];
rz(-2.313068) q[1];
x q[2];
rz(-2.6785128) q[3];
sx q[3];
rz(-2.7768917) q[3];
sx q[3];
rz(-2.0228342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.018365232) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(-1.3398735) q[2];
rz(0.54723048) q[3];
sx q[3];
rz(-0.42357835) q[3];
sx q[3];
rz(2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.055534) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(-0.15922971) q[0];
rz(-0.010146443) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(2.8841282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5437357) q[0];
sx q[0];
rz(-2.275995) q[0];
sx q[0];
rz(2.0446834) q[0];
x q[1];
rz(-3.0765216) q[2];
sx q[2];
rz(-2.3461968) q[2];
sx q[2];
rz(0.8000904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60002335) q[1];
sx q[1];
rz(-2.5065055) q[1];
sx q[1];
rz(2.9734008) q[1];
x q[2];
rz(-0.93521714) q[3];
sx q[3];
rz(-0.87163371) q[3];
sx q[3];
rz(0.88672598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7633535) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(2.6687458) q[2];
rz(0.7044479) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(-0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064297) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(-0.065486431) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(-1.4261036) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1714892) q[0];
sx q[0];
rz(-1.7894151) q[0];
sx q[0];
rz(-1.6660369) q[0];
rz(-pi) q[1];
rz(2.9272396) q[2];
sx q[2];
rz(-2.7446236) q[2];
sx q[2];
rz(3.1052542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9551505) q[1];
sx q[1];
rz(-1.3480061) q[1];
sx q[1];
rz(-0.28526116) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0992005) q[3];
sx q[3];
rz(-2.8082153) q[3];
sx q[3];
rz(2.7957145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9466729) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(-1.3067513) q[2];
rz(-1.9715747) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(-3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829247) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(0.40147993) q[0];
rz(0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(2.2875517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9977048) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(-1.3530988) q[0];
rz(-0.67777216) q[2];
sx q[2];
rz(-2.9052264) q[2];
sx q[2];
rz(1.4217699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8484162) q[1];
sx q[1];
rz(-2.5418315) q[1];
sx q[1];
rz(-0.17677115) q[1];
x q[2];
rz(2.0186508) q[3];
sx q[3];
rz(-2.1338226) q[3];
sx q[3];
rz(-0.26859586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4868769) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(-2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(0.83561713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198782) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(-1.1623435) q[0];
rz(-1.989919) q[1];
sx q[1];
rz(-1.78777) q[1];
sx q[1];
rz(2.2231359) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8480523) q[0];
sx q[0];
rz(-1.2402724) q[0];
sx q[0];
rz(-2.8577515) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26084857) q[2];
sx q[2];
rz(-2.2493304) q[2];
sx q[2];
rz(1.4229753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20275252) q[1];
sx q[1];
rz(-1.5781286) q[1];
sx q[1];
rz(-1.5640902) q[1];
rz(-pi) q[2];
rz(-1.1226467) q[3];
sx q[3];
rz(-2.7635241) q[3];
sx q[3];
rz(0.44071769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-0.50298634) q[2];
rz(2.3668187) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(-1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69304943) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(1.5013129) q[0];
rz(0.34128183) q[1];
sx q[1];
rz(-1.7106067) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687896) q[0];
sx q[0];
rz(-2.4478292) q[0];
sx q[0];
rz(0.60141464) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0985435) q[2];
sx q[2];
rz(-2.611428) q[2];
sx q[2];
rz(-1.4250172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4367366) q[1];
sx q[1];
rz(-0.67187998) q[1];
sx q[1];
rz(0.19499548) q[1];
x q[2];
rz(1.7057034) q[3];
sx q[3];
rz(-1.1793609) q[3];
sx q[3];
rz(1.5086482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1559653) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(-0.36925527) q[2];
rz(2.8333832) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(1.5306028) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592598) q[0];
sx q[0];
rz(-1.3066602) q[0];
sx q[0];
rz(-0.34307137) q[0];
rz(-1.169091) q[1];
sx q[1];
rz(-1.3645376) q[1];
sx q[1];
rz(-1.7128568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37821445) q[0];
sx q[0];
rz(-1.9665446) q[0];
sx q[0];
rz(1.2984492) q[0];
x q[1];
rz(0.83447225) q[2];
sx q[2];
rz(-2.1205466) q[2];
sx q[2];
rz(2.2817734) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6862168) q[1];
sx q[1];
rz(-1.1955402) q[1];
sx q[1];
rz(-3.0757559) q[1];
rz(0.87782209) q[3];
sx q[3];
rz(-1.1980499) q[3];
sx q[3];
rz(-0.28723785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2009361) q[2];
sx q[2];
rz(-0.20874615) q[2];
sx q[2];
rz(1.3915871) q[2];
rz(-0.40677795) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(2.2627635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(1.9039924) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(-0.79992574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7163776) q[0];
sx q[0];
rz(-2.3122462) q[0];
sx q[0];
rz(2.8275376) q[0];
rz(2.5913057) q[2];
sx q[2];
rz(-1.3105416) q[2];
sx q[2];
rz(1.257892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25617304) q[1];
sx q[1];
rz(-0.90682632) q[1];
sx q[1];
rz(-1.7682942) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25000817) q[3];
sx q[3];
rz(-2.0431314) q[3];
sx q[3];
rz(1.8369305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37821975) q[2];
sx q[2];
rz(-1.4241445) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(0.25589219) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(-1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.273461) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.3399667) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(3.1392787) q[2];
sx q[2];
rz(-0.70941596) q[2];
sx q[2];
rz(1.8298168) q[2];
rz(-1.003391) q[3];
sx q[3];
rz(-1.0091253) q[3];
sx q[3];
rz(-0.71970018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
