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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0092702) q[0];
sx q[0];
rz(-1.2148804) q[0];
sx q[0];
rz(-0.15412553) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6662007) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(-0.22104095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77817569) q[1];
sx q[1];
rz(-0.84901224) q[1];
sx q[1];
rz(-1.282882) q[1];
rz(-2.6333991) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(3.1014266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(2.3085964) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-2.5966068) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-0.82495904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23404113) q[0];
sx q[0];
rz(-1.7579161) q[0];
sx q[0];
rz(-1.1603111) q[0];
rz(-pi) q[1];
rz(1.5812133) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(-1.8889697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35926871) q[1];
sx q[1];
rz(-1.6869079) q[1];
sx q[1];
rz(-0.5639204) q[1];
rz(-pi) q[2];
rz(1.8726146) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(-0.17641243) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-0.038539561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24667106) q[0];
sx q[0];
rz(-1.5785494) q[0];
sx q[0];
rz(3.1409516) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23080319) q[2];
sx q[2];
rz(-1.7204086) q[2];
sx q[2];
rz(-0.19300592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7125394) q[1];
sx q[1];
rz(-0.4936115) q[1];
sx q[1];
rz(1.7008971) q[1];
rz(-pi) q[2];
rz(-0.58979804) q[3];
sx q[3];
rz(-2.0876948) q[3];
sx q[3];
rz(0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(-1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53699025) q[0];
sx q[0];
rz(-1.450481) q[0];
sx q[0];
rz(3.0483732) q[0];
rz(-pi) q[1];
rz(-2.1148557) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.3865711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2745167) q[1];
sx q[1];
rz(-1.80596) q[1];
sx q[1];
rz(-2.0348674) q[1];
rz(1.5870729) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(-1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-2.7635014) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3867214) q[0];
sx q[0];
rz(-2.5140962) q[0];
sx q[0];
rz(-1.0794712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81850448) q[2];
sx q[2];
rz(-2.2193925) q[2];
sx q[2];
rz(1.6550145) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34199076) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(1.1127383) q[1];
x q[2];
rz(-1.3729172) q[3];
sx q[3];
rz(-2.0298925) q[3];
sx q[3];
rz(-1.0517694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.4168926) q[3];
sx q[3];
rz(1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3865005) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(-0.37762541) q[0];
rz(-2.8185484) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(2.2264218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.299303) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(-1.4847941) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0038063) q[2];
sx q[2];
rz(-1.1319455) q[2];
sx q[2];
rz(-1.3958508) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.542791) q[1];
sx q[1];
rz(-0.32074499) q[1];
sx q[1];
rz(2.9221605) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1359453) q[3];
sx q[3];
rz(-0.69909401) q[3];
sx q[3];
rz(1.618315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6293684) q[0];
sx q[0];
rz(-2.6656796) q[0];
sx q[0];
rz(-0.41643629) q[0];
rz(-2.695589) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(0.29689483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62896171) q[1];
sx q[1];
rz(-2.1753864) q[1];
sx q[1];
rz(-1.1058034) q[1];
rz(-pi) q[2];
rz(-3.1281934) q[3];
sx q[3];
rz(-1.0671167) q[3];
sx q[3];
rz(-2.3434533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(-1.210775) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(-0.67529768) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(2.3892367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1479748) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(2.5698635) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3812177) q[2];
sx q[2];
rz(-2.484212) q[2];
sx q[2];
rz(2.0331403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68554316) q[1];
sx q[1];
rz(-1.4577663) q[1];
sx q[1];
rz(3.0526524) q[1];
rz(-1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(1.0866144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2032808) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(2.032062) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(-1.8267652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1413404) q[0];
sx q[0];
rz(-2.25657) q[0];
sx q[0];
rz(1.89639) q[0];
rz(-2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(-0.044791128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6445551) q[1];
sx q[1];
rz(-2.2215448) q[1];
sx q[1];
rz(0.52849309) q[1];
x q[2];
rz(-0.0074578961) q[3];
sx q[3];
rz(-2.3141626) q[3];
sx q[3];
rz(-0.5479365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(-2.7108575) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-2.2626256) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4742728) q[0];
sx q[0];
rz(-0.55756888) q[0];
sx q[0];
rz(-1.5865109) q[0];
rz(-1.395412) q[2];
sx q[2];
rz(-2.4782964) q[2];
sx q[2];
rz(-1.3822008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1230159) q[1];
sx q[1];
rz(-1.7565109) q[1];
sx q[1];
rz(-2.3260444) q[1];
rz(-pi) q[2];
rz(1.2519022) q[3];
sx q[3];
rz(-1.0341757) q[3];
sx q[3];
rz(1.5981789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5320756) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-0.77970589) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(-2.1203534) q[3];
sx q[3];
rz(-1.2720576) q[3];
sx q[3];
rz(-2.6779327) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
