OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(-1.0993212) q[1];
sx q[1];
rz(-1.6834747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320839) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-0.30755933) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61383944) q[2];
sx q[2];
rz(-1.5547353) q[2];
sx q[2];
rz(-1.2889372) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2394489) q[1];
sx q[1];
rz(-1.4414756) q[1];
sx q[1];
rz(-0.37954482) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59431608) q[3];
sx q[3];
rz(-1.672097) q[3];
sx q[3];
rz(3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-0.17949417) q[2];
rz(1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.1546086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5886473) q[0];
sx q[0];
rz(-1.5536904) q[0];
sx q[0];
rz(0.018789142) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0915394) q[2];
sx q[2];
rz(-2.4467231) q[2];
sx q[2];
rz(0.75887242) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7696015) q[1];
sx q[1];
rz(-0.76906119) q[1];
sx q[1];
rz(0.12188697) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9566544) q[3];
sx q[3];
rz(-1.3122845) q[3];
sx q[3];
rz(-2.1188494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-2.2581805) q[2];
rz(-0.47131053) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(1.6261684) q[0];
rz(-2.5405163) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(-2.0498958) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1611623) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(-0.79361332) q[0];
rz(-pi) q[1];
rz(-1.0513564) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.5067593) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97869067) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(0.49180007) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.6347173) q[3];
sx q[3];
rz(-1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-2.8312347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045517) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(-1.4287352) q[0];
x q[1];
rz(1.4448302) q[2];
sx q[2];
rz(-0.8896041) q[2];
sx q[2];
rz(1.3931526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.51018184) q[1];
sx q[1];
rz(-1.5893755) q[1];
sx q[1];
rz(-1.2127962) q[1];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(-1.8087169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0115396) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-0.87019428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0349883) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(-2.5614221) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-2.4325271) q[2];
sx q[2];
rz(2.2774334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40904564) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(-1.9802666) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3051885) q[3];
sx q[3];
rz(-2.5634273) q[3];
sx q[3];
rz(2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(2.4678521) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314635) q[0];
sx q[0];
rz(-1.1420113) q[0];
sx q[0];
rz(-1.6786806) q[0];
rz(-pi) q[1];
rz(-0.67955534) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7266453) q[1];
sx q[1];
rz(-2.2989095) q[1];
sx q[1];
rz(-2.3689518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9353988) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(0.091726124) q[2];
rz(-2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-3.1076028) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2798529) q[0];
sx q[0];
rz(-1.4089157) q[0];
sx q[0];
rz(-2.3539761) q[0];
rz(0.54079536) q[2];
sx q[2];
rz(-2.135709) q[2];
sx q[2];
rz(-2.5164547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9719203) q[1];
sx q[1];
rz(-1.5202513) q[1];
sx q[1];
rz(-2.2488942) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9339318) q[3];
sx q[3];
rz(-2.5163979) q[3];
sx q[3];
rz(0.085426424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(0.87654385) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1720393) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(-2.0041549) q[0];
rz(-pi) q[1];
rz(2.5567899) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(-2.1540097) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3280914) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(-2.8495795) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94533841) q[3];
sx q[3];
rz(-1.006554) q[3];
sx q[3];
rz(-0.64627796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(-1.5173222) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354934) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(-1.8898213) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1447434) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(-1.4247308) q[0];
x q[1];
rz(-1.5179713) q[2];
sx q[2];
rz(-0.17715684) q[2];
sx q[2];
rz(-1.0747386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.138315) q[1];
sx q[1];
rz(-1.7764846) q[1];
sx q[1];
rz(-0.059271952) q[1];
rz(-pi) q[2];
rz(0.59659776) q[3];
sx q[3];
rz(-1.3571686) q[3];
sx q[3];
rz(1.240977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4831055) q[0];
sx q[0];
rz(-2.8327836) q[0];
sx q[0];
rz(-1.2312066) q[0];
rz(-pi) q[1];
rz(2.3775616) q[2];
sx q[2];
rz(-1.6198297) q[2];
sx q[2];
rz(-1.3438091) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2227576) q[1];
sx q[1];
rz(-1.7822052) q[1];
sx q[1];
rz(-2.0445776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5851192) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(-2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88400921) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(0.94669) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(2.3161841) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
rz(2.0675038) q[3];
sx q[3];
rz(-0.91377331) q[3];
sx q[3];
rz(2.5627315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
