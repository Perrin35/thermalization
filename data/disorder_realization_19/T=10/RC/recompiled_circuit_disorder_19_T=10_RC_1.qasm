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
rz(-1.6436613) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1323224) q[0];
sx q[0];
rz(-1.9267123) q[0];
sx q[0];
rz(2.9874671) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0718578) q[2];
sx q[2];
rz(-0.94143922) q[2];
sx q[2];
rz(-0.33915181) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5421996) q[1];
sx q[1];
rz(-1.3560055) q[1];
sx q[1];
rz(-2.3989124) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64925361) q[3];
sx q[3];
rz(-1.2459727) q[3];
sx q[3];
rz(-2.0089825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(-0.82495904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.724052) q[0];
sx q[0];
rz(-1.1678956) q[0];
sx q[0];
rz(-0.20362644) q[0];
rz(3.1171563) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(-1.2791866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1383914) q[1];
sx q[1];
rz(-2.1304641) q[1];
sx q[1];
rz(-1.707934) q[1];
rz(1.8726146) q[3];
sx q[3];
rz(-2.5538553) q[3];
sx q[3];
rz(2.1835818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58723441) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-0.17641243) q[2];
rz(0.13088626) q[3];
sx q[3];
rz(-1.9661048) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(-2.6718455) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(-3.1030531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8949216) q[0];
sx q[0];
rz(-1.5785494) q[0];
sx q[0];
rz(-0.00064108032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23080319) q[2];
sx q[2];
rz(-1.421184) q[2];
sx q[2];
rz(2.9485867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4290532) q[1];
sx q[1];
rz(-2.6479811) q[1];
sx q[1];
rz(-1.4406956) q[1];
x q[2];
rz(2.3452957) q[3];
sx q[3];
rz(-2.3781804) q[3];
sx q[3];
rz(2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5588351) q[2];
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
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(1.2448467) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(0.19515881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6046024) q[0];
sx q[0];
rz(-1.6911117) q[0];
sx q[0];
rz(3.0483732) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0267369) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.7550215) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.86707592) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(2.0348674) q[1];
x q[2];
rz(2.5024662) q[3];
sx q[3];
rz(-1.5577321) q[3];
sx q[3];
rz(-0.55707219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23094709) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(-3.0019794) q[0];
rz(-2.0647678) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3867214) q[0];
sx q[0];
rz(-2.5140962) q[0];
sx q[0];
rz(-1.0794712) q[0];
rz(-0.73363186) q[2];
sx q[2];
rz(-0.94978226) q[2];
sx q[2];
rz(0.48895198) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7996019) q[1];
sx q[1];
rz(-1.9365891) q[1];
sx q[1];
rz(2.0288543) q[1];
rz(-2.7630745) q[3];
sx q[3];
rz(-2.6444728) q[3];
sx q[3];
rz(0.62687031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(2.2536229) q[2];
rz(0.97638431) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(1.005727) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75509214) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(2.7639672) q[0];
rz(2.8185484) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(0.91517085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.897402) q[0];
sx q[0];
rz(-1.4892329) q[0];
sx q[0];
rz(0.32342644) q[0];
x q[1];
rz(2.2890131) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(-2.728087) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.542791) q[1];
sx q[1];
rz(-0.32074499) q[1];
sx q[1];
rz(2.9221605) q[1];
x q[2];
rz(1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53720981) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(2.3510695) q[2];
rz(-2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(-1.908318) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.4253915) q[1];
rz(-pi/2) q[2];
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
rz(-1.5064429) q[2];
sx q[2];
rz(-1.1255985) q[2];
sx q[2];
rz(1.3016456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5126309) q[1];
sx q[1];
rz(-0.96620622) q[1];
sx q[1];
rz(-1.1058034) q[1];
rz(-0.013399259) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(-2.3434533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(-0.0018421729) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90826666) q[0];
sx q[0];
rz(-1.8857297) q[0];
sx q[0];
rz(2.6106735) q[0];
x q[1];
rz(-0.1444507) q[2];
sx q[2];
rz(-2.2144197) q[2];
sx q[2];
rz(-0.87063906) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1250549) q[1];
sx q[1];
rz(-0.14370951) q[1];
sx q[1];
rz(-2.2347666) q[1];
rz(-pi) q[2];
rz(0.58952491) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(-2.7713283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-0.00014649815) q[2];
rz(1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-0.29763597) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14934854) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.3148274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7816759) q[0];
sx q[0];
rz(-1.8209551) q[0];
sx q[0];
rz(-0.71235384) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8034231) q[2];
sx q[2];
rz(-2.7578232) q[2];
sx q[2];
rz(1.2107809) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.73003522) q[1];
sx q[1];
rz(-0.81333232) q[1];
sx q[1];
rz(0.98585339) q[1];
x q[2];
rz(1.5626838) q[3];
sx q[3];
rz(-2.3981961) q[3];
sx q[3];
rz(0.5589561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(2.7108575) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(-0.47732863) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(0.39189664) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083188699) q[0];
sx q[0];
rz(-1.579111) q[0];
sx q[0];
rz(1.0132829) q[0];
rz(-pi) q[1];
rz(2.2266065) q[2];
sx q[2];
rz(-1.4631541) q[2];
sx q[2];
rz(2.8142625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1230159) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(-0.8155483) q[1];
rz(-pi) q[2];
rz(0.48505731) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5320756) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(-1.1575451) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.8023087) q[1];
sx q[1];
rz(-0.61359275) q[1];
sx q[1];
rz(0.35992122) q[1];
rz(-1.3988613) q[2];
sx q[2];
rz(-2.3430765) q[2];
sx q[2];
rz(0.55745468) q[2];
rz(1.0380575) q[3];
sx q[3];
rz(-0.61809117) q[3];
sx q[3];
rz(1.5865159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
