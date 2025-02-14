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
rz(0.81646252) q[0];
sx q[0];
rz(-3.0397968) q[0];
sx q[0];
rz(0.53959227) q[0];
rz(3.6379023) q[1];
sx q[1];
rz(3.451347) q[1];
sx q[1];
rz(5.7440905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236299) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(2.7001803) q[0];
rz(1.7221872) q[2];
sx q[2];
rz(-1.2149723) q[2];
sx q[2];
rz(2.9265938) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0338933) q[1];
sx q[1];
rz(-1.4310208) q[1];
sx q[1];
rz(1.8962533) q[1];
rz(-pi) q[2];
x q[2];
rz(0.081955628) q[3];
sx q[3];
rz(-2.2181619) q[3];
sx q[3];
rz(-2.820197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3794136) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(1.1890821) q[2];
rz(-1.1842229) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14520833) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(-1.3231963) q[0];
rz(0.48201758) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(0.97420305) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486748) q[0];
sx q[0];
rz(-2.4705624) q[0];
sx q[0];
rz(-1.2811518) q[0];
rz(-1.6840132) q[2];
sx q[2];
rz(-1.6694476) q[2];
sx q[2];
rz(-1.1455918) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1434553) q[1];
sx q[1];
rz(-2.1051198) q[1];
sx q[1];
rz(1.2972531) q[1];
x q[2];
rz(0.5663381) q[3];
sx q[3];
rz(-1.3674595) q[3];
sx q[3];
rz(1.7254064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(-2.3919487) q[2];
rz(-0.43198112) q[3];
sx q[3];
rz(-2.0260729) q[3];
sx q[3];
rz(1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5169446) q[0];
sx q[0];
rz(-2.8693146) q[0];
sx q[0];
rz(2.5352449) q[0];
rz(1.3336522) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(-0.25951728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0344275) q[0];
sx q[0];
rz(-1.8176862) q[0];
sx q[0];
rz(-0.16040032) q[0];
x q[1];
rz(0.17510842) q[2];
sx q[2];
rz(-1.3710183) q[2];
sx q[2];
rz(-1.7410884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3644581) q[1];
sx q[1];
rz(-1.196047) q[1];
sx q[1];
rz(-1.115429) q[1];
rz(2.1331779) q[3];
sx q[3];
rz(-1.1964238) q[3];
sx q[3];
rz(-1.238747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8848662) q[2];
sx q[2];
rz(-1.7776411) q[2];
sx q[2];
rz(1.5474896) q[2];
rz(-2.3571842) q[3];
sx q[3];
rz(-1.514785) q[3];
sx q[3];
rz(1.5756395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49482685) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(-0.25949091) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-0.44993284) q[1];
sx q[1];
rz(-1.4422013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7176162) q[0];
sx q[0];
rz(-2.0956796) q[0];
sx q[0];
rz(-0.54846008) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97564189) q[2];
sx q[2];
rz(-1.6791428) q[2];
sx q[2];
rz(-0.25749046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3842786) q[1];
sx q[1];
rz(-2.7065249) q[1];
sx q[1];
rz(-2.5237094) q[1];
rz(-pi) q[2];
rz(1.6569225) q[3];
sx q[3];
rz(-0.68400506) q[3];
sx q[3];
rz(2.4620584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0356902) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-2.7316459) q[2];
rz(-1.2116872) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(1.3338026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7164417) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(-0.43824276) q[1];
sx q[1];
rz(-1.848105) q[1];
sx q[1];
rz(-0.73572198) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05984027) q[0];
sx q[0];
rz(-0.81938374) q[0];
sx q[0];
rz(-1.6900586) q[0];
rz(-0.14938905) q[2];
sx q[2];
rz(-1.427703) q[2];
sx q[2];
rz(0.18994513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.84668881) q[1];
sx q[1];
rz(-1.3219993) q[1];
sx q[1];
rz(-1.873068) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92145958) q[3];
sx q[3];
rz(-1.9917352) q[3];
sx q[3];
rz(-2.4025847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2030486) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(-2.5277444) q[2];
rz(2.7326873) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(-2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530389) q[0];
sx q[0];
rz(-0.63816324) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(1.6963814) q[1];
sx q[1];
rz(-0.45672363) q[1];
sx q[1];
rz(-0.17527418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5492903) q[0];
sx q[0];
rz(-1.3978551) q[0];
sx q[0];
rz(-0.025385234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1138328) q[2];
sx q[2];
rz(-0.98688417) q[2];
sx q[2];
rz(-1.2500545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3289017) q[1];
sx q[1];
rz(-1.1393424) q[1];
sx q[1];
rz(-2.575909) q[1];
x q[2];
rz(2.3168269) q[3];
sx q[3];
rz(-2.3769393) q[3];
sx q[3];
rz(0.95461707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1767629) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(2.1691587) q[2];
rz(1.6644299) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.2870188) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(2.7884685) q[0];
rz(2.1137721) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(2.1305398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42732692) q[0];
sx q[0];
rz(-1.2135226) q[0];
sx q[0];
rz(-1.5554122) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89247668) q[2];
sx q[2];
rz(-0.8304285) q[2];
sx q[2];
rz(2.3811946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2680929) q[1];
sx q[1];
rz(-0.38485369) q[1];
sx q[1];
rz(1.1420239) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8500438) q[3];
sx q[3];
rz(-1.8729775) q[3];
sx q[3];
rz(0.19552375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(1.774452) q[2];
rz(-1.8732871) q[3];
sx q[3];
rz(-1.4788078) q[3];
sx q[3];
rz(-2.2936599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276412) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(2.0116346) q[0];
rz(0.21895151) q[1];
sx q[1];
rz(-1.4066701) q[1];
sx q[1];
rz(2.2311282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6235891) q[0];
sx q[0];
rz(-1.9027038) q[0];
sx q[0];
rz(-2.5123572) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47825079) q[2];
sx q[2];
rz(-1.646981) q[2];
sx q[2];
rz(-2.8752799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76669756) q[1];
sx q[1];
rz(-1.8378705) q[1];
sx q[1];
rz(-1.0407991) q[1];
rz(-0.72515709) q[3];
sx q[3];
rz(-1.8817543) q[3];
sx q[3];
rz(-0.5288611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8255446) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.8749219) q[3];
sx q[3];
rz(-2.6958444) q[3];
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
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6398741) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(0.83183944) q[0];
rz(-2.4141451) q[1];
sx q[1];
rz(-1.7214382) q[1];
sx q[1];
rz(3.0453392) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1548658) q[0];
sx q[0];
rz(-0.71457246) q[0];
sx q[0];
rz(-2.3113219) q[0];
rz(-pi) q[1];
rz(0.25419323) q[2];
sx q[2];
rz(-1.5113748) q[2];
sx q[2];
rz(2.3583902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61273328) q[1];
sx q[1];
rz(-1.7816356) q[1];
sx q[1];
rz(3.0490521) q[1];
x q[2];
rz(-1.7670125) q[3];
sx q[3];
rz(-1.4150146) q[3];
sx q[3];
rz(0.84583144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7353797) q[2];
sx q[2];
rz(-1.4344183) q[2];
sx q[2];
rz(1.5926788) q[2];
rz(0.63273543) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(-0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525472) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(-2.7885875) q[0];
rz(2.8189335) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(-0.37193146) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3561365) q[0];
sx q[0];
rz(-1.8710941) q[0];
sx q[0];
rz(0.2072643) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3916799) q[2];
sx q[2];
rz(-1.8955911) q[2];
sx q[2];
rz(1.3472404) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10807762) q[1];
sx q[1];
rz(-1.4107499) q[1];
sx q[1];
rz(2.8041502) q[1];
x q[2];
rz(-1.5266034) q[3];
sx q[3];
rz(-2.4165476) q[3];
sx q[3];
rz(2.8738662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17452621) q[2];
sx q[2];
rz(-1.9517978) q[2];
sx q[2];
rz(1.2971499) q[2];
rz(-0.8477115) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2321155) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(1.7535946) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(0.83609381) q[2];
sx q[2];
rz(-2.2420364) q[2];
sx q[2];
rz(-1.4902761) q[2];
rz(-0.41107117) q[3];
sx q[3];
rz(-0.84079327) q[3];
sx q[3];
rz(-0.83415915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
