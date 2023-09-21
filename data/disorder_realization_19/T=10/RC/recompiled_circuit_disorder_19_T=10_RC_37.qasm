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
rz(0.37880701) q[1];
sx q[1];
rz(-1.3728377) q[1];
sx q[1];
rz(1.6436613) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4286715) q[0];
sx q[0];
rz(-0.38654583) q[0];
sx q[0];
rz(-1.1791694) q[0];
x q[1];
rz(-0.069734863) q[2];
sx q[2];
rz(-0.94143922) q[2];
sx q[2];
rz(2.8024408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5421996) q[1];
sx q[1];
rz(-1.3560055) q[1];
sx q[1];
rz(-0.7426803) q[1];
rz(-2.6333991) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(-0.040166044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-0.32749367) q[1];
sx q[1];
rz(-0.82495904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23404113) q[0];
sx q[0];
rz(-1.3836765) q[0];
sx q[0];
rz(1.9812816) q[0];
rz(-1.5603793) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(1.252623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7488885) q[1];
sx q[1];
rz(-0.57447937) q[1];
sx q[1];
rz(-2.9267465) q[1];
x q[2];
rz(0.19552688) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(-0.60002458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5543582) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(1.6167971) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-3.1030531) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32917133) q[0];
sx q[0];
rz(-3.1338131) q[0];
sx q[0];
rz(1.4882985) q[0];
x q[1];
rz(-0.58263393) q[2];
sx q[2];
rz(-0.27432549) q[2];
sx q[2];
rz(0.81253101) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7125394) q[1];
sx q[1];
rz(-0.4936115) q[1];
sx q[1];
rz(1.4406956) q[1];
x q[2];
rz(2.1707118) q[3];
sx q[3];
rz(-1.0661134) q[3];
sx q[3];
rz(-1.2952627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-0.91252404) q[2];
rz(1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(-2.7918949) q[0];
rz(1.2448467) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(-0.19515881) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6046024) q[0];
sx q[0];
rz(-1.6911117) q[0];
sx q[0];
rz(3.0483732) q[0];
x q[1];
rz(-0.85907016) q[2];
sx q[2];
rz(-0.67407437) q[2];
sx q[2];
rz(-2.7328343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3217722) q[1];
sx q[1];
rz(-2.0211377) q[1];
sx q[1];
rz(2.8798073) q[1];
rz(-0.6391265) q[3];
sx q[3];
rz(-1.5577321) q[3];
sx q[3];
rz(2.5845205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.267743) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-3.0019794) q[0];
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
rz(2.3867214) q[0];
sx q[0];
rz(-0.62749642) q[0];
sx q[0];
rz(1.0794712) q[0];
rz(-pi) q[1];
rz(2.3230882) q[2];
sx q[2];
rz(-0.92220014) q[2];
sx q[2];
rz(-1.6550145) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2851583) q[1];
sx q[1];
rz(-0.57796961) q[1];
sx q[1];
rz(-2.2846089) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37851815) q[3];
sx q[3];
rz(-0.49711984) q[3];
sx q[3];
rz(-2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-2.2536229) q[2];
rz(0.97638431) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(-2.2264218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5765502) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(2.8898426) q[0];
x q[1];
rz(2.2890131) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(0.41350565) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8296216) q[1];
sx q[1];
rz(-1.8835856) q[1];
sx q[1];
rz(-1.642986) q[1];
x q[2];
rz(3.0056474) q[3];
sx q[3];
rz(-2.4424986) q[3];
sx q[3];
rz(-1.618315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(0.61736068) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73615605) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(-1.2332747) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(1.4253915) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7084301) q[0];
sx q[0];
rz(-1.7571974) q[0];
sx q[0];
rz(0.44048803) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0076214) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(-1.1531032) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9219626) q[1];
sx q[1];
rz(-1.1929409) q[1];
sx q[1];
rz(-0.65803836) q[1];
x q[2];
rz(-1.5464877) q[3];
sx q[3];
rz(-0.50384249) q[3];
sx q[3];
rz(-0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.210775) q[2];
rz(-3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-0.5843038) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9936179) q[0];
sx q[0];
rz(-2.5320842) q[0];
sx q[0];
rz(2.5698635) q[0];
rz(1.760375) q[2];
sx q[2];
rz(-0.65738064) q[2];
sx q[2];
rz(-2.0331403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4560495) q[1];
sx q[1];
rz(-1.6838264) q[1];
sx q[1];
rz(3.0526524) q[1];
x q[2];
rz(-1.3648871) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(-1.0866144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(2.1355656) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.8267652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00025230322) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(1.2452026) q[0];
rz(-pi) q[1];
rz(2.8034231) q[2];
sx q[2];
rz(-0.38376946) q[2];
sx q[2];
rz(-1.9308117) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73003522) q[1];
sx q[1];
rz(-2.3282603) q[1];
sx q[1];
rz(2.1557393) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0074578961) q[3];
sx q[3];
rz(-2.3141626) q[3];
sx q[3];
rz(-0.5479365) q[3];
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
rz(-1.0916748) q[3];
sx q[3];
rz(-2.664264) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-2.749696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4927917) q[0];
sx q[0];
rz(-2.1282882) q[0];
sx q[0];
rz(-3.1317943) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2266065) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(0.3273302) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7614903) q[1];
sx q[1];
rz(-0.83161608) q[1];
sx q[1];
rz(2.889061) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6565353) q[3];
sx q[3];
rz(-2.5254446) q[3];
sx q[3];
rz(-0.96998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75523238) q[0];
sx q[0];
rz(-1.3437143) q[0];
sx q[0];
rz(1.2585826) q[0];
rz(-1.8023087) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-2.9677283) q[2];
sx q[2];
rz(-0.78730135) q[2];
sx q[2];
rz(0.80136328) q[2];
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
