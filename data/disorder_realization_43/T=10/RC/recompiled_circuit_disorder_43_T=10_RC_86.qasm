OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2655576) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(2.5493456) q[0];
x q[1];
rz(1.0693597) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(-1.8471579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0426892) q[1];
sx q[1];
rz(-0.96682036) q[1];
sx q[1];
rz(2.9788115) q[1];
x q[2];
rz(-1.1638155) q[3];
sx q[3];
rz(-0.80899901) q[3];
sx q[3];
rz(2.80796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(2.3577918) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3621688) q[0];
sx q[0];
rz(-2.9075025) q[0];
sx q[0];
rz(-0.97658821) q[0];
x q[1];
rz(2.2750521) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(1.3726335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9176863) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(0.19885893) q[1];
rz(-2.5793521) q[3];
sx q[3];
rz(-1.7842818) q[3];
sx q[3];
rz(1.7052887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(0.0058962065) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(2.7745461) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(-0.095741622) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.30004) q[0];
sx q[0];
rz(-1.7471572) q[0];
sx q[0];
rz(-0.15533133) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.1635457) q[2];
sx q[2];
rz(0.72088748) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9853471) q[1];
sx q[1];
rz(-0.74838446) q[1];
sx q[1];
rz(-2.2845539) q[1];
rz(-pi) q[2];
rz(1.6060353) q[3];
sx q[3];
rz(-1.0198776) q[3];
sx q[3];
rz(-2.2743724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(-0.52571458) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(1.0338763) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047065145) q[0];
sx q[0];
rz(-2.3405502) q[0];
sx q[0];
rz(-2.4497776) q[0];
rz(-2.6150377) q[2];
sx q[2];
rz(-2.0806081) q[2];
sx q[2];
rz(1.1916135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8139207) q[1];
sx q[1];
rz(-1.76183) q[1];
sx q[1];
rz(-3.0753067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86446188) q[3];
sx q[3];
rz(-2.7221788) q[3];
sx q[3];
rz(0.97233397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(-2.7205617) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-2.72686) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(0.57410747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78050437) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(0.9719406) q[0];
x q[1];
rz(-2.4262869) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(-0.41358435) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.638006) q[1];
sx q[1];
rz(-1.2212911) q[1];
sx q[1];
rz(-1.867678) q[1];
rz(-2.2361034) q[3];
sx q[3];
rz(-1.4169766) q[3];
sx q[3];
rz(-0.32115768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85270143) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-2.4093157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9376611) q[0];
sx q[0];
rz(-1.4121778) q[0];
sx q[0];
rz(-0.21477867) q[0];
x q[1];
rz(0.90768355) q[2];
sx q[2];
rz(-2.216202) q[2];
sx q[2];
rz(-1.6931319) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69387586) q[1];
sx q[1];
rz(-0.69732053) q[1];
sx q[1];
rz(0.35481528) q[1];
rz(-pi) q[2];
rz(2.5712588) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(-1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.1614655) q[2];
rz(2.8325864) q[3];
sx q[3];
rz(-1.8892663) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(-2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-0.77004534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55243385) q[0];
sx q[0];
rz(-1.3981515) q[0];
sx q[0];
rz(-2.5347559) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1427878) q[2];
sx q[2];
rz(-2.4704128) q[2];
sx q[2];
rz(-0.83516781) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8416482) q[1];
sx q[1];
rz(-2.3023459) q[1];
sx q[1];
rz(0.70920918) q[1];
rz(-pi) q[2];
rz(-1.9023864) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-0.36639211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893433) q[0];
sx q[0];
rz(-2.0848668) q[0];
sx q[0];
rz(2.423717) q[0];
rz(-pi) q[1];
rz(1.1364469) q[2];
sx q[2];
rz(-0.98698101) q[2];
sx q[2];
rz(2.7455612) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7494292) q[1];
sx q[1];
rz(-0.87595075) q[1];
sx q[1];
rz(-1.7872582) q[1];
rz(-pi) q[2];
rz(3.0931926) q[3];
sx q[3];
rz(-1.0194922) q[3];
sx q[3];
rz(-1.6325523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1239132) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095187) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(2.2576387) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.489483) q[0];
sx q[0];
rz(-1.2533422) q[0];
sx q[0];
rz(0.86274685) q[0];
rz(-2.6762166) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(-2.344775) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8225122) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(-2.6507069) q[1];
x q[2];
rz(-1.0320372) q[3];
sx q[3];
rz(-1.7146535) q[3];
sx q[3];
rz(-1.1080351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(-2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.7609319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7847583) q[0];
sx q[0];
rz(-2.1958302) q[0];
sx q[0];
rz(0.16108315) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5153377) q[2];
sx q[2];
rz(-2.7055253) q[2];
sx q[2];
rz(-0.01576327) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41439357) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(0.90331932) q[1];
x q[2];
rz(0.72762604) q[3];
sx q[3];
rz(-2.8907223) q[3];
sx q[3];
rz(2.6550456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3537139) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(2.805368) q[2];
rz(-2.0119038) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(0.54429383) q[1];
sx q[1];
rz(-1.9017362) q[1];
sx q[1];
rz(-1.5128296) q[1];
rz(2.3168646) q[2];
sx q[2];
rz(-1.8246973) q[2];
sx q[2];
rz(1.4551103) q[2];
rz(0.23562283) q[3];
sx q[3];
rz(-1.0141254) q[3];
sx q[3];
rz(0.49401382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];