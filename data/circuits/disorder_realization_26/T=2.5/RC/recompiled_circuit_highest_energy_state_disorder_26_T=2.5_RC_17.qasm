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
rz(0.8910203) q[0];
sx q[0];
rz(-1.2863337) q[0];
sx q[0];
rz(-2.2556055) q[0];
rz(0.52357829) q[1];
sx q[1];
rz(3.7012586) q[1];
sx q[1];
rz(8.1044365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.490075) q[0];
sx q[0];
rz(-1.0371672) q[0];
sx q[0];
rz(0.20488157) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34446005) q[2];
sx q[2];
rz(-2.3340324) q[2];
sx q[2];
rz(1.5429614) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7976171) q[1];
sx q[1];
rz(-1.6317751) q[1];
sx q[1];
rz(-1.0073376) q[1];
x q[2];
rz(-2.376972) q[3];
sx q[3];
rz(-2.21661) q[3];
sx q[3];
rz(-0.82987204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3628799) q[2];
sx q[2];
rz(-1.3211297) q[2];
sx q[2];
rz(-0.92817456) q[2];
rz(-1.5859531) q[3];
sx q[3];
rz(-2.0944244) q[3];
sx q[3];
rz(-1.1699404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8856119) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(2.0804491) q[0];
rz(-1.9288918) q[1];
sx q[1];
rz(-2.3586528) q[1];
sx q[1];
rz(1.3139668) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43976682) q[0];
sx q[0];
rz(-1.9226388) q[0];
sx q[0];
rz(-2.3934796) q[0];
rz(-pi) q[1];
rz(2.2167614) q[2];
sx q[2];
rz(-1.2154757) q[2];
sx q[2];
rz(2.3424847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1873389) q[1];
sx q[1];
rz(-2.0994774) q[1];
sx q[1];
rz(-2.9462255) q[1];
rz(-2.2960194) q[3];
sx q[3];
rz(-1.9872287) q[3];
sx q[3];
rz(1.6628979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5709915) q[2];
sx q[2];
rz(-0.97390318) q[2];
sx q[2];
rz(0.904733) q[2];
rz(1.5915126) q[3];
sx q[3];
rz(-1.855987) q[3];
sx q[3];
rz(-1.2929644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631113) q[0];
sx q[0];
rz(-3.0516629) q[0];
sx q[0];
rz(0.91823804) q[0];
rz(0.69084424) q[1];
sx q[1];
rz(-1.7450688) q[1];
sx q[1];
rz(0.7965368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7762866) q[0];
sx q[0];
rz(-1.7788634) q[0];
sx q[0];
rz(3.0479504) q[0];
rz(-pi) q[1];
rz(3.0987691) q[2];
sx q[2];
rz(-1.6302675) q[2];
sx q[2];
rz(-0.24569337) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59350508) q[1];
sx q[1];
rz(-1.6789762) q[1];
sx q[1];
rz(-0.87203474) q[1];
rz(-pi) q[2];
rz(-0.91751601) q[3];
sx q[3];
rz(-1.2923354) q[3];
sx q[3];
rz(2.3115932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.010926509) q[2];
sx q[2];
rz(-0.9504168) q[2];
sx q[2];
rz(-1.8939023) q[2];
rz(2.583336) q[3];
sx q[3];
rz(-1.6853251) q[3];
sx q[3];
rz(-3.1202417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.438544) q[0];
sx q[0];
rz(-0.049228638) q[0];
sx q[0];
rz(-0.7274279) q[0];
rz(-2.9797331) q[1];
sx q[1];
rz(-1.1715803) q[1];
sx q[1];
rz(1.9836099) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25140554) q[0];
sx q[0];
rz(-1.5047538) q[0];
sx q[0];
rz(-1.4736045) q[0];
rz(1.4314992) q[2];
sx q[2];
rz(-1.5698395) q[2];
sx q[2];
rz(0.19569163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0643081) q[1];
sx q[1];
rz(-1.3798215) q[1];
sx q[1];
rz(0.125566) q[1];
rz(-pi) q[2];
rz(2.6154989) q[3];
sx q[3];
rz(-1.4988314) q[3];
sx q[3];
rz(-1.9561033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27266476) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(1.391601) q[2];
rz(-2.9113655) q[3];
sx q[3];
rz(-1.1743436) q[3];
sx q[3];
rz(-0.19190425) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72628438) q[0];
sx q[0];
rz(-0.10157651) q[0];
sx q[0];
rz(-2.0368077) q[0];
rz(3.0305908) q[1];
sx q[1];
rz(-1.1155201) q[1];
sx q[1];
rz(1.312779) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26129042) q[0];
sx q[0];
rz(-1.748607) q[0];
sx q[0];
rz(-1.7755204) q[0];
x q[1];
rz(-3.0022125) q[2];
sx q[2];
rz(-2.5733893) q[2];
sx q[2];
rz(0.85124337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1089563) q[1];
sx q[1];
rz(-1.7708888) q[1];
sx q[1];
rz(1.97831) q[1];
rz(1.1024226) q[3];
sx q[3];
rz(-2.3658731) q[3];
sx q[3];
rz(-0.44743809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.50292) q[2];
sx q[2];
rz(-2.0827677) q[2];
sx q[2];
rz(-1.0271094) q[2];
rz(2.1549759) q[3];
sx q[3];
rz(-0.28237453) q[3];
sx q[3];
rz(-1.7178887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8091549) q[0];
sx q[0];
rz(-1.2323392) q[0];
sx q[0];
rz(2.3409081) q[0];
rz(0.77146161) q[1];
sx q[1];
rz(-2.0636676) q[1];
sx q[1];
rz(1.7652184) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9405917) q[0];
sx q[0];
rz(-1.9671114) q[0];
sx q[0];
rz(0.70756377) q[0];
rz(-1.6402836) q[2];
sx q[2];
rz(-2.67454) q[2];
sx q[2];
rz(-0.75331538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47471646) q[1];
sx q[1];
rz(-2.2762269) q[1];
sx q[1];
rz(1.5882467) q[1];
rz(-pi) q[2];
rz(-0.82628754) q[3];
sx q[3];
rz(-2.1231696) q[3];
sx q[3];
rz(-0.9780405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2996404) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(0.017814962) q[2];
rz(2.8282015) q[3];
sx q[3];
rz(-1.0217977) q[3];
sx q[3];
rz(2.4221086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3571091) q[0];
sx q[0];
rz(-1.5672368) q[0];
sx q[0];
rz(0.16726476) q[0];
rz(2.3978865) q[1];
sx q[1];
rz(-2.2552762) q[1];
sx q[1];
rz(-3.0053265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7881476) q[0];
sx q[0];
rz(-2.2290725) q[0];
sx q[0];
rz(0.35231715) q[0];
x q[1];
rz(1.0616779) q[2];
sx q[2];
rz(-0.3442685) q[2];
sx q[2];
rz(-1.6580251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4140461) q[1];
sx q[1];
rz(-0.36484499) q[1];
sx q[1];
rz(-2.9529497) q[1];
rz(1.6266009) q[3];
sx q[3];
rz(-1.9161092) q[3];
sx q[3];
rz(1.2962411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3027489) q[2];
sx q[2];
rz(-2.9911797) q[2];
sx q[2];
rz(-1.0813084) q[2];
rz(2.9980764) q[3];
sx q[3];
rz(-1.84294) q[3];
sx q[3];
rz(-1.4474086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3926587) q[0];
sx q[0];
rz(-1.8086139) q[0];
sx q[0];
rz(0.64312154) q[0];
rz(-0.92203036) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(1.6304852) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47925835) q[0];
sx q[0];
rz(-1.7500203) q[0];
sx q[0];
rz(-1.2249806) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11428516) q[2];
sx q[2];
rz(-1.8488036) q[2];
sx q[2];
rz(-0.68113995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0550168) q[1];
sx q[1];
rz(-2.9447703) q[1];
sx q[1];
rz(1.5301401) q[1];
rz(-1.2503406) q[3];
sx q[3];
rz(-2.0600187) q[3];
sx q[3];
rz(-1.3716404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.94613218) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(1.3004318) q[2];
rz(2.2293633) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(-0.31360489) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6533971) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(0.61757863) q[0];
rz(0.081534475) q[1];
sx q[1];
rz(-2.640994) q[1];
sx q[1];
rz(2.4542782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9195557) q[0];
sx q[0];
rz(-1.3341781) q[0];
sx q[0];
rz(3.0927004) q[0];
rz(0.388889) q[2];
sx q[2];
rz(-0.28770471) q[2];
sx q[2];
rz(-3.1294745) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8246714) q[1];
sx q[1];
rz(-1.8536398) q[1];
sx q[1];
rz(-2.7088195) q[1];
rz(-pi) q[2];
rz(1.4908904) q[3];
sx q[3];
rz(-0.54466313) q[3];
sx q[3];
rz(-0.93526697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11757892) q[2];
sx q[2];
rz(-1.1774096) q[2];
sx q[2];
rz(0.071361072) q[2];
rz(1.9060382) q[3];
sx q[3];
rz(-2.8046298) q[3];
sx q[3];
rz(0.56628847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.5736893) q[0];
sx q[0];
rz(-0.25147831) q[0];
sx q[0];
rz(0.13791826) q[0];
rz(2.5714286) q[1];
sx q[1];
rz(-2.0591996) q[1];
sx q[1];
rz(-3.1261442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.100525) q[0];
sx q[0];
rz(-2.3599632) q[0];
sx q[0];
rz(-2.0360721) q[0];
rz(-1.1006946) q[2];
sx q[2];
rz(-0.86720556) q[2];
sx q[2];
rz(0.64794618) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5970832) q[1];
sx q[1];
rz(-0.5209777) q[1];
sx q[1];
rz(-0.3358261) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7543162) q[3];
sx q[3];
rz(-1.6257111) q[3];
sx q[3];
rz(-0.6764937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7323759) q[2];
sx q[2];
rz(-0.78745431) q[2];
sx q[2];
rz(0.43238762) q[2];
rz(-0.89808291) q[3];
sx q[3];
rz(-2.2663074) q[3];
sx q[3];
rz(1.557365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6569923) q[0];
sx q[0];
rz(-1.1285755) q[0];
sx q[0];
rz(-2.3791671) q[0];
rz(-2.5905329) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(-1.5253992) q[2];
sx q[2];
rz(-1.9337966) q[2];
sx q[2];
rz(-0.36317229) q[2];
rz(-2.2021709) q[3];
sx q[3];
rz(-2.4236854) q[3];
sx q[3];
rz(-3.051306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
