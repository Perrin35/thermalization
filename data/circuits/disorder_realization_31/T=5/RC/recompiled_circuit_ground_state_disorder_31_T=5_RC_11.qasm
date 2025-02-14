OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6346729) q[0];
sx q[0];
rz(1.4022175) q[0];
sx q[0];
rz(13.316857) q[0];
rz(-3.0783202) q[1];
sx q[1];
rz(-0.81875357) q[1];
sx q[1];
rz(-1.0190581) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4472103) q[0];
sx q[0];
rz(-1.3827494) q[0];
sx q[0];
rz(0.5699372) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86111492) q[2];
sx q[2];
rz(-1.7758596) q[2];
sx q[2];
rz(2.4697863) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99996725) q[1];
sx q[1];
rz(-2.5427175) q[1];
sx q[1];
rz(-0.24354045) q[1];
rz(-pi) q[2];
rz(-0.5818855) q[3];
sx q[3];
rz(-0.62633461) q[3];
sx q[3];
rz(-0.7382381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84746209) q[2];
sx q[2];
rz(-2.1475466) q[2];
sx q[2];
rz(1.4685941) q[2];
rz(0.20644203) q[3];
sx q[3];
rz(-1.3279746) q[3];
sx q[3];
rz(0.65202057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9701397) q[0];
sx q[0];
rz(-1.411922) q[0];
sx q[0];
rz(-0.10301244) q[0];
rz(0.39762321) q[1];
sx q[1];
rz(-0.42639521) q[1];
sx q[1];
rz(2.2893589) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6838339) q[0];
sx q[0];
rz(-1.5105828) q[0];
sx q[0];
rz(-2.1046066) q[0];
x q[1];
rz(-1.2322517) q[2];
sx q[2];
rz(-1.8391092) q[2];
sx q[2];
rz(1.2635096) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5928015) q[1];
sx q[1];
rz(-2.617814) q[1];
sx q[1];
rz(1.5984867) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9485365) q[3];
sx q[3];
rz(-1.5279433) q[3];
sx q[3];
rz(-0.022006527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6323866) q[2];
sx q[2];
rz(-2.1337324) q[2];
sx q[2];
rz(-2.3243135) q[2];
rz(-2.9362074) q[3];
sx q[3];
rz(-0.21586625) q[3];
sx q[3];
rz(-0.5425905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83194724) q[0];
sx q[0];
rz(-1.0878071) q[0];
sx q[0];
rz(-0.52016869) q[0];
rz(-0.7236411) q[1];
sx q[1];
rz(-1.2836722) q[1];
sx q[1];
rz(2.9626194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72523141) q[0];
sx q[0];
rz(-2.0746914) q[0];
sx q[0];
rz(2.079062) q[0];
x q[1];
rz(-0.51836022) q[2];
sx q[2];
rz(-0.39769618) q[2];
sx q[2];
rz(2.8373371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6295669) q[1];
sx q[1];
rz(-2.8993052) q[1];
sx q[1];
rz(2.971388) q[1];
x q[2];
rz(-2.7951205) q[3];
sx q[3];
rz(-0.68744171) q[3];
sx q[3];
rz(0.32780743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2554864) q[2];
sx q[2];
rz(-0.28033689) q[2];
sx q[2];
rz(-0.96957668) q[2];
rz(0.3420091) q[3];
sx q[3];
rz(-1.0229144) q[3];
sx q[3];
rz(-1.3268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372357) q[0];
sx q[0];
rz(-2.3396753) q[0];
sx q[0];
rz(-0.075415762) q[0];
rz(0.30741179) q[1];
sx q[1];
rz(-1.1630029) q[1];
sx q[1];
rz(-0.42617646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722011) q[0];
sx q[0];
rz(-1.0474023) q[0];
sx q[0];
rz(-0.70573576) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46838872) q[2];
sx q[2];
rz(-1.8231744) q[2];
sx q[2];
rz(0.017367432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0016617) q[1];
sx q[1];
rz(-2.5298018) q[1];
sx q[1];
rz(-1.2242713) q[1];
rz(-pi) q[2];
rz(-0.40554096) q[3];
sx q[3];
rz(-1.5933005) q[3];
sx q[3];
rz(0.5675494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9093466) q[2];
sx q[2];
rz(-1.1541977) q[2];
sx q[2];
rz(-2.2171891) q[2];
rz(2.038548) q[3];
sx q[3];
rz(-2.0274935) q[3];
sx q[3];
rz(-1.2801142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1396609) q[0];
sx q[0];
rz(-2.6502471) q[0];
sx q[0];
rz(0.77144462) q[0];
rz(-1.3123243) q[1];
sx q[1];
rz(-1.7701365) q[1];
sx q[1];
rz(-1.9171453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3712922) q[0];
sx q[0];
rz(-1.5427569) q[0];
sx q[0];
rz(-1.4337149) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.13891) q[2];
sx q[2];
rz(-0.61842881) q[2];
sx q[2];
rz(-1.9069049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6836732) q[1];
sx q[1];
rz(-2.7453642) q[1];
sx q[1];
rz(1.8093682) q[1];
x q[2];
rz(1.1125074) q[3];
sx q[3];
rz(-0.84064129) q[3];
sx q[3];
rz(2.2793913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9839342) q[2];
sx q[2];
rz(-0.57571405) q[2];
sx q[2];
rz(-1.0945339) q[2];
rz(0.076722773) q[3];
sx q[3];
rz(-0.50944296) q[3];
sx q[3];
rz(0.69242394) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43129608) q[0];
sx q[0];
rz(-2.4831979) q[0];
sx q[0];
rz(-0.21507138) q[0];
rz(-0.98945016) q[1];
sx q[1];
rz(-0.96885252) q[1];
sx q[1];
rz(-1.6542124) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3475018) q[0];
sx q[0];
rz(-2.5985744) q[0];
sx q[0];
rz(2.8229982) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8645668) q[2];
sx q[2];
rz(-1.9508298) q[2];
sx q[2];
rz(-0.88949163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3645059) q[1];
sx q[1];
rz(-1.987127) q[1];
sx q[1];
rz(-2.7440666) q[1];
rz(-pi) q[2];
rz(0.66296799) q[3];
sx q[3];
rz(-1.9687395) q[3];
sx q[3];
rz(2.8090973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6626849) q[2];
sx q[2];
rz(-1.4143133) q[2];
sx q[2];
rz(1.5642081) q[2];
rz(2.1936737) q[3];
sx q[3];
rz(-1.0736059) q[3];
sx q[3];
rz(1.4028367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3941512) q[0];
sx q[0];
rz(-1.489137) q[0];
sx q[0];
rz(-0.58962756) q[0];
rz(1.6472752) q[1];
sx q[1];
rz(-1.3811771) q[1];
sx q[1];
rz(-3.0628915) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.323228) q[0];
sx q[0];
rz(-1.8591037) q[0];
sx q[0];
rz(-2.7298729) q[0];
rz(-pi) q[1];
rz(-2.6888618) q[2];
sx q[2];
rz(-2.6331278) q[2];
sx q[2];
rz(1.1217211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5509531) q[1];
sx q[1];
rz(-1.7788243) q[1];
sx q[1];
rz(-0.25352884) q[1];
rz(-0.77459333) q[3];
sx q[3];
rz(-2.4586185) q[3];
sx q[3];
rz(-2.3247064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9239203) q[2];
sx q[2];
rz(-2.2675026) q[2];
sx q[2];
rz(3.0295642) q[2];
rz(2.7364386) q[3];
sx q[3];
rz(-1.743789) q[3];
sx q[3];
rz(-0.20598327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9236295) q[0];
sx q[0];
rz(-2.1681652) q[0];
sx q[0];
rz(1.4189036) q[0];
rz(1.0796684) q[1];
sx q[1];
rz(-2.7579312) q[1];
sx q[1];
rz(1.4071646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1315707) q[0];
sx q[0];
rz(-0.52531201) q[0];
sx q[0];
rz(2.4729687) q[0];
rz(-pi) q[1];
rz(-2.7101369) q[2];
sx q[2];
rz(-1.4979216) q[2];
sx q[2];
rz(-0.88957722) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0041173) q[1];
sx q[1];
rz(-1.1761929) q[1];
sx q[1];
rz(2.7156107) q[1];
rz(-0.24875764) q[3];
sx q[3];
rz(-1.837656) q[3];
sx q[3];
rz(-2.998179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.63274038) q[2];
sx q[2];
rz(-0.036157046) q[2];
sx q[2];
rz(2.4031694) q[2];
rz(2.9449055) q[3];
sx q[3];
rz(-2.5236712) q[3];
sx q[3];
rz(2.5014307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38615534) q[0];
sx q[0];
rz(-2.8897987) q[0];
sx q[0];
rz(0.039903076) q[0];
rz(-0.2208605) q[1];
sx q[1];
rz(-1.618229) q[1];
sx q[1];
rz(-1.1562645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98469767) q[0];
sx q[0];
rz(-1.8348719) q[0];
sx q[0];
rz(-0.15021245) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9258181) q[2];
sx q[2];
rz(-1.7606944) q[2];
sx q[2];
rz(0.18958651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58312341) q[1];
sx q[1];
rz(-1.9052747) q[1];
sx q[1];
rz(2.1534647) q[1];
x q[2];
rz(-1.6764033) q[3];
sx q[3];
rz(-1.6855414) q[3];
sx q[3];
rz(-2.4345943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7463189) q[2];
sx q[2];
rz(-2.7776182) q[2];
sx q[2];
rz(-3.086997) q[2];
rz(0.78628457) q[3];
sx q[3];
rz(-1.9449284) q[3];
sx q[3];
rz(-1.1597077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7603124) q[0];
sx q[0];
rz(-0.78802839) q[0];
sx q[0];
rz(-0.81083167) q[0];
rz(0.52931085) q[1];
sx q[1];
rz(-1.7050754) q[1];
sx q[1];
rz(2.0307821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7082695) q[0];
sx q[0];
rz(-2.4305211) q[0];
sx q[0];
rz(2.5266886) q[0];
x q[1];
rz(2.2718515) q[2];
sx q[2];
rz(-0.68863622) q[2];
sx q[2];
rz(2.5272126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6880092) q[1];
sx q[1];
rz(-0.78026672) q[1];
sx q[1];
rz(2.9095786) q[1];
x q[2];
rz(0.30382321) q[3];
sx q[3];
rz(-2.15832) q[3];
sx q[3];
rz(-1.2161906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79798737) q[2];
sx q[2];
rz(-0.56861773) q[2];
sx q[2];
rz(-0.24151754) q[2];
rz(0.24164116) q[3];
sx q[3];
rz(-1.6499949) q[3];
sx q[3];
rz(-2.9901166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.3342313) q[0];
sx q[0];
rz(-1.7932899) q[0];
sx q[0];
rz(0.67247969) q[0];
rz(-0.45202759) q[1];
sx q[1];
rz(-0.94257911) q[1];
sx q[1];
rz(-0.51530757) q[1];
rz(-1.4448901) q[2];
sx q[2];
rz(-0.44582146) q[2];
sx q[2];
rz(1.3498342) q[2];
rz(-0.47361278) q[3];
sx q[3];
rz(-1.0992186) q[3];
sx q[3];
rz(0.038224955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
