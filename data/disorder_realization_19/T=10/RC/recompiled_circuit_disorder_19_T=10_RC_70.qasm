OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1146381) q[0];
sx q[0];
rz(-1.4517598) q[0];
sx q[0];
rz(2.4960158) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(-1.3728377) q[1];
sx q[1];
rz(1.6436613) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1323224) q[0];
sx q[0];
rz(-1.9267123) q[0];
sx q[0];
rz(0.15412553) q[0];
rz(1.475392) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(0.22104095) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9420535) q[1];
sx q[1];
rz(-2.374211) q[1];
sx q[1];
rz(-2.8295423) q[1];
x q[2];
rz(2.6333991) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(-2.3085964) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-2.5966068) q[0];
rz(-0.90822059) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-0.82495904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23404113) q[0];
sx q[0];
rz(-1.7579161) q[0];
sx q[0];
rz(1.9812816) q[0];
rz(-1.5603793) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(-1.8889697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35926871) q[1];
sx q[1];
rz(-1.4546848) q[1];
sx q[1];
rz(0.5639204) q[1];
rz(-pi) q[2];
rz(1.2689781) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(0.9580108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(-0.13088626) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-0.46974716) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(0.038539561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24667106) q[0];
sx q[0];
rz(-1.5785494) q[0];
sx q[0];
rz(0.00064108032) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4171717) q[2];
sx q[2];
rz(-1.3426174) q[2];
sx q[2];
rz(1.4128026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.027027834) q[1];
sx q[1];
rz(-1.5092883) q[1];
sx q[1];
rz(1.0807178) q[1];
x q[2];
rz(-2.1707118) q[3];
sx q[3];
rz(-1.0661134) q[3];
sx q[3];
rz(-1.8463299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(2.9464338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168266) q[0];
sx q[0];
rz(-2.989528) q[0];
sx q[0];
rz(-0.91465871) q[0];
x q[1];
rz(-1.0267369) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.7550215) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86707592) q[1];
sx q[1];
rz(-1.80596) q[1];
sx q[1];
rz(-1.1067252) q[1];
rz(-pi) q[2];
rz(0.6391265) q[3];
sx q[3];
rz(-1.5577321) q[3];
sx q[3];
rz(0.55707219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(-1.0686482) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-0.1396133) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(-0.37809125) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.970801) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(2.8118954) q[0];
x q[1];
rz(2.4079608) q[2];
sx q[2];
rz(-2.1918104) q[2];
sx q[2];
rz(2.6526407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2851583) q[1];
sx q[1];
rz(-2.563623) q[1];
sx q[1];
rz(-2.2846089) q[1];
rz(2.7630745) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(0.62687031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(2.2536229) q[2];
rz(-2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.75509214) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(-0.91517085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5650425) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(0.25175005) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2890131) q[2];
sx q[2];
rz(-2.4396509) q[2];
sx q[2];
rz(-2.728087) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23657654) q[1];
sx q[1];
rz(-1.502115) q[1];
sx q[1];
rz(-2.828039) q[1];
rz(0.69453199) q[3];
sx q[3];
rz(-1.4834705) q[3];
sx q[3];
rz(-0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6293684) q[0];
sx q[0];
rz(-2.6656796) q[0];
sx q[0];
rz(-0.41643629) q[0];
x q[1];
rz(0.44600365) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(0.29689483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9219626) q[1];
sx q[1];
rz(-1.9486517) q[1];
sx q[1];
rz(0.65803836) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5951049) q[3];
sx q[3];
rz(-0.50384249) q[3];
sx q[3];
rz(0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(1.210775) q[2];
rz(-0.0018421729) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-0.75235596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90826666) q[0];
sx q[0];
rz(-1.255863) q[0];
sx q[0];
rz(2.6106735) q[0];
x q[1];
rz(-0.92213995) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(2.5285072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1250549) q[1];
sx q[1];
rz(-2.9978831) q[1];
sx q[1];
rz(2.2347666) q[1];
rz(-pi) q[2];
rz(-0.58952491) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(2.7713283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2032808) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(1.0060271) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.8267652) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6515503) q[0];
sx q[0];
rz(-2.393912) q[0];
sx q[0];
rz(0.3726532) q[0];
rz(-0.33816955) q[2];
sx q[2];
rz(-0.38376946) q[2];
sx q[2];
rz(1.2107809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41373738) q[1];
sx q[1];
rz(-1.9836042) q[1];
sx q[1];
rz(-0.84819838) q[1];
rz(-pi) q[2];
rz(-2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(2.1237802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7018147) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(2.7108575) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(-2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4196639) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(0.28668177) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(2.749696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6673198) q[0];
sx q[0];
rz(-2.5840238) q[0];
sx q[0];
rz(-1.5865109) q[0];
x q[1];
rz(-3.0060843) q[2];
sx q[2];
rz(-0.91943179) q[2];
sx q[2];
rz(-1.1609921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1230159) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(0.8155483) q[1];
rz(-2.5819671) q[3];
sx q[3];
rz(-1.2979753) q[3];
sx q[3];
rz(0.19459693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75523238) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(-0.77970589) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(-0.34655456) q[3];
sx q[3];
rz(-2.0934436) q[3];
sx q[3];
rz(2.2128076) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
