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
rz(-0.64557689) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(1.4979314) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(1.9306246) q[0];
rz(-pi) q[1];
x q[1];
rz(1.475392) q[2];
sx q[2];
rz(-2.5089052) q[2];
sx q[2];
rz(-0.22104095) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.363417) q[1];
sx q[1];
rz(-0.84901224) q[1];
sx q[1];
rz(-1.8587106) q[1];
x q[2];
rz(1.1708158) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(0.67584544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92007414) q[2];
sx q[2];
rz(-1.5204844) q[2];
sx q[2];
rz(-2.8011838) q[2];
rz(0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-0.54498589) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-2.3166336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2087814) q[0];
sx q[0];
rz(-2.6926846) q[0];
sx q[0];
rz(-1.1277898) q[0];
rz(-pi) q[1];
rz(-1.5812133) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(-1.252623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35926871) q[1];
sx q[1];
rz(-1.6869079) q[1];
sx q[1];
rz(-0.5639204) q[1];
rz(-pi) q[2];
rz(2.9460658) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(0.60002458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(-0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(0.46974716) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-3.1030531) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8174724) q[0];
sx q[0];
rz(-1.5701553) q[0];
sx q[0];
rz(-1.5785494) q[0];
x q[1];
rz(1.724421) q[2];
sx q[2];
rz(-1.7989752) q[2];
sx q[2];
rz(-1.4128026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.027027834) q[1];
sx q[1];
rz(-1.5092883) q[1];
sx q[1];
rz(-1.0807178) q[1];
rz(-pi) q[2];
rz(-2.5517946) q[3];
sx q[3];
rz(-1.0538978) q[3];
sx q[3];
rz(0.043881744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(0.91252404) q[2];
rz(1.8330666) q[3];
sx q[3];
rz(-1.0037183) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(-1.2448467) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0965658) q[0];
sx q[0];
rz(-1.66334) q[0];
sx q[0];
rz(-1.6916313) q[0];
rz(-pi) q[1];
rz(2.6606576) q[2];
sx q[2];
rz(-2.0630884) q[2];
sx q[2];
rz(0.42602691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8198204) q[1];
sx q[1];
rz(-2.0211377) q[1];
sx q[1];
rz(-2.8798073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5545198) q[3];
sx q[3];
rz(-2.2098594) q[3];
sx q[3];
rz(1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.267743) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-2.7635014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75487126) q[0];
sx q[0];
rz(-2.5140962) q[0];
sx q[0];
rz(-2.0621215) q[0];
x q[1];
rz(-2.3374704) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(-0.59876195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0873449) q[1];
sx q[1];
rz(-1.1451045) q[1];
sx q[1];
rz(-0.40360968) q[1];
rz(-pi) q[2];
rz(-0.46697163) q[3];
sx q[3];
rz(-1.7479556) q[3];
sx q[3];
rz(0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(2.7639672) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(2.2264218) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422896) q[0];
sx q[0];
rz(-1.248484) q[0];
sx q[0];
rz(1.6567985) q[0];
rz(-pi) q[1];
rz(2.1377863) q[2];
sx q[2];
rz(-1.1319455) q[2];
sx q[2];
rz(1.3958508) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23657654) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(0.31355365) q[1];
rz(-1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.3464348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(0.79052314) q[2];
rz(-0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(-0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(2.7745568) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.7162011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7084301) q[0];
sx q[0];
rz(-1.7571974) q[0];
sx q[0];
rz(0.44048803) q[0];
rz(-pi) q[1];
rz(2.695589) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(-0.29689483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.047445) q[1];
sx q[1];
rz(-2.3970251) q[1];
sx q[1];
rz(-0.57569699) q[1];
rz(1.5464877) q[3];
sx q[3];
rz(-2.6377502) q[3];
sx q[3];
rz(-0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(-1.210775) q[2];
rz(0.0018421729) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(0.75235596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589688) q[0];
sx q[0];
rz(-1.0685182) q[0];
sx q[0];
rz(1.9320022) q[0];
rz(-pi) q[1];
rz(-0.92213995) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(-0.61308544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68554316) q[1];
sx q[1];
rz(-1.4577663) q[1];
sx q[1];
rz(0.088940253) q[1];
rz(-1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(1.0866144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(3.1414462) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(2.1355656) q[0];
rz(0.26793119) q[1];
sx q[1];
rz(-1.8012828) q[1];
sx q[1];
rz(1.3148274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49004236) q[0];
sx q[0];
rz(-0.74768066) q[0];
sx q[0];
rz(2.7689395) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8034231) q[2];
sx q[2];
rz(-2.7578232) q[2];
sx q[2];
rz(-1.9308117) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6445551) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(0.52849309) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5626838) q[3];
sx q[3];
rz(-2.3981961) q[3];
sx q[3];
rz(-0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43977794) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4742728) q[0];
sx q[0];
rz(-0.55756888) q[0];
sx q[0];
rz(-1.5550818) q[0];
x q[1];
rz(3.0060843) q[2];
sx q[2];
rz(-2.2221609) q[2];
sx q[2];
rz(-1.1609921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74589409) q[1];
sx q[1];
rz(-0.77334009) q[1];
sx q[1];
rz(1.3032773) q[1];
x q[2];
rz(2.6565353) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(-2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(0.55084294) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(0.77970589) q[2];
sx q[2];
rz(-1.4479326) q[2];
sx q[2];
rz(-0.89276199) q[2];
rz(-2.7950381) q[3];
sx q[3];
rz(-1.0481491) q[3];
sx q[3];
rz(-0.92878503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
