OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(4.5933525) q[0];
sx q[0];
rz(10.070355) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(-1.6436613) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6490373) q[0];
sx q[0];
rz(-1.4264002) q[0];
sx q[0];
rz(-1.210968) q[0];
rz(-pi) q[1];
rz(-2.2013118) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(-1.8688569) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9420535) q[1];
sx q[1];
rz(-0.76738165) q[1];
sx q[1];
rz(2.8295423) q[1];
x q[2];
rz(-1.1708158) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(2.4657472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(2.5966068) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(0.82495904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93281125) q[0];
sx q[0];
rz(-2.6926846) q[0];
sx q[0];
rz(2.0138028) q[0];
rz(3.1171563) q[2];
sx q[2];
rz(-0.40301286) q[2];
sx q[2];
rz(-1.2791866) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7488885) q[1];
sx q[1];
rz(-0.57447937) q[1];
sx q[1];
rz(0.21484612) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0042079) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(-2.2752938) q[3];
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
rz(0.13088626) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(0.032657284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3787057) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(-0.038539561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8174724) q[0];
sx q[0];
rz(-1.5714374) q[0];
sx q[0];
rz(1.5785494) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58263393) q[2];
sx q[2];
rz(-2.8672672) q[2];
sx q[2];
rz(0.81253101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5650428) q[1];
sx q[1];
rz(-1.0817263) q[1];
sx q[1];
rz(3.0719041) q[1];
x q[2];
rz(0.79629691) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(0.91252404) q[2];
rz(-1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.4413888) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7581166) q[0];
sx q[0];
rz(-1.3154727) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(1.2448467) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6046024) q[0];
sx q[0];
rz(-1.6911117) q[0];
sx q[0];
rz(-3.0483732) q[0];
rz(1.0267369) q[2];
sx q[2];
rz(-1.1508905) q[2];
sx q[2];
rz(-1.3865711) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2745167) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(-1.1067252) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1196932) q[3];
sx q[3];
rz(-0.63924131) q[3];
sx q[3];
rz(2.1454449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9106456) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(1.267743) q[2];
rz(-1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(0.1396133) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-0.37809125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17079167) q[0];
sx q[0];
rz(-1.0266725) q[0];
sx q[0];
rz(-0.32969726) q[0];
rz(0.80412229) q[2];
sx q[2];
rz(-2.1466549) q[2];
sx q[2];
rz(2.5428307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0873449) q[1];
sx q[1];
rz(-1.1451045) q[1];
sx q[1];
rz(-0.40360968) q[1];
x q[2];
rz(-2.674621) q[3];
sx q[3];
rz(-1.3936371) q[3];
sx q[3];
rz(0.60764473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.7639672) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(2.2264218) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897402) q[0];
sx q[0];
rz(-1.6523598) q[0];
sx q[0];
rz(0.32342644) q[0];
x q[1];
rz(-2.1377863) q[2];
sx q[2];
rz(-2.0096471) q[2];
sx q[2];
rz(-1.7457419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5988016) q[1];
sx q[1];
rz(-2.8208477) q[1];
sx q[1];
rz(-0.21943211) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69453199) q[3];
sx q[3];
rz(-1.6581222) q[3];
sx q[3];
rz(-2.9897523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53720981) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(0.61736068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.7162011) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51222425) q[0];
sx q[0];
rz(-0.47591305) q[0];
sx q[0];
rz(2.7251564) q[0];
rz(1.5064429) q[2];
sx q[2];
rz(-1.1255985) q[2];
sx q[2];
rz(1.839947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9219626) q[1];
sx q[1];
rz(-1.1929409) q[1];
sx q[1];
rz(-2.4835543) q[1];
x q[2];
rz(2.0745139) q[3];
sx q[3];
rz(-1.5590612) q[3];
sx q[3];
rz(0.77912441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(-1.9308176) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(3.0650744) q[0];
rz(-0.67529768) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(2.3892367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479748) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(0.57172914) q[0];
x q[1];
rz(-1.760375) q[2];
sx q[2];
rz(-0.65738064) q[2];
sx q[2];
rz(2.0331403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4560495) q[1];
sx q[1];
rz(-1.6838264) q[1];
sx q[1];
rz(-3.0526524) q[1];
rz(-2.5520677) q[3];
sx q[3];
rz(-1.7426963) q[3];
sx q[3];
rz(0.37026438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.3148274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1413404) q[0];
sx q[0];
rz(-0.88502266) q[0];
sx q[0];
rz(1.89639) q[0];
rz(-pi) q[1];
rz(-2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(-0.044791128) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4970376) q[1];
sx q[1];
rz(-2.2215448) q[1];
sx q[1];
rz(-0.52849309) q[1];
rz(2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(1.0178125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(-2.8549109) q[0];
rz(0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(0.39189664) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4927917) q[0];
sx q[0];
rz(-1.0133044) q[0];
sx q[0];
rz(-3.1317943) q[0];
rz(-pi) q[1];
rz(-0.13550831) q[2];
sx q[2];
rz(-0.91943179) q[2];
sx q[2];
rz(1.1609921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38010234) q[1];
sx q[1];
rz(-2.3099766) q[1];
sx q[1];
rz(0.2525316) q[1];
rz(-0.48505731) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(-2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(1.1575451) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(-1.9633861) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(-1.8023087) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-0.17386439) q[2];
sx q[2];
rz(-2.3542913) q[2];
sx q[2];
rz(-2.3402294) q[2];
rz(-2.1035351) q[3];
sx q[3];
rz(-0.61809117) q[3];
sx q[3];
rz(1.5865159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
