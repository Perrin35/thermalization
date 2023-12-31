OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(0.89755091) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7449887) q[0];
sx q[0];
rz(-0.94192266) q[0];
sx q[0];
rz(-3.0552342) q[0];
rz(-pi) q[1];
rz(-0.93809442) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(-3.0245568) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(-0.27889241) q[1];
x q[2];
rz(0.1944794) q[3];
sx q[3];
rz(-0.76768657) q[3];
sx q[3];
rz(-3.1377813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(0.73195362) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(2.2251341) q[0];
rz(-2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40936138) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(2.2344927) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5231045) q[2];
sx q[2];
rz(-1.7881219) q[2];
sx q[2];
rz(2.6478812) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0658873) q[1];
sx q[1];
rz(-1.9888708) q[1];
sx q[1];
rz(-1.651152) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2611748) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(-2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(2.0203967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(2.0545161) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7240702) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(-3.0561471) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9417291) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(-0.73168879) q[1];
rz(2.6477473) q[3];
sx q[3];
rz(-1.6742799) q[3];
sx q[3];
rz(2.8007357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1040092) q[0];
sx q[0];
rz(-1.596367) q[0];
sx q[0];
rz(2.5081162) q[0];
x q[1];
rz(2.7599081) q[2];
sx q[2];
rz(-1.0978062) q[2];
sx q[2];
rz(-2.4406529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1981922) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(1.1073768) q[1];
x q[2];
rz(1.6223525) q[3];
sx q[3];
rz(-1.0815074) q[3];
sx q[3];
rz(-1.4435022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(-1.1486357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38265739) q[0];
sx q[0];
rz(-1.4969345) q[0];
sx q[0];
rz(-2.1942684) q[0];
rz(2.1158475) q[2];
sx q[2];
rz(-0.75136853) q[2];
sx q[2];
rz(0.92912208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3682813) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.3694622) q[1];
rz(-2.4162021) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1308412) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(1.032069) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5500568) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(0.39598879) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-2.9352303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9128742) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(-2.3939783) q[0];
rz(0.66531078) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(-2.2523508) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98112647) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(-1.4658982) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2431074) q[3];
sx q[3];
rz(-1.7565691) q[3];
sx q[3];
rz(0.9542619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-2.5653429) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.9304088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241087) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(-0.80362513) q[0];
rz(0.25522916) q[2];
sx q[2];
rz(-0.50349871) q[2];
sx q[2];
rz(2.2358759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2562099) q[1];
sx q[1];
rz(-1.4146283) q[1];
sx q[1];
rz(0.068710879) q[1];
rz(2.1738449) q[3];
sx q[3];
rz(-1.839404) q[3];
sx q[3];
rz(-0.3097765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3367735) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.221009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2894665) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(0.37129398) q[0];
rz(0.17692716) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(0.71080506) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11215969) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(0.68316858) q[1];
x q[2];
rz(-0.82151316) q[3];
sx q[3];
rz(-1.8213846) q[3];
sx q[3];
rz(-0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-2.4364046) q[2];
rz(2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(3.0158214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14995689) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(1.3791023) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.618082) q[2];
sx q[2];
rz(-2.7611809) q[2];
sx q[2];
rz(-1.7720122) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8217433) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-2.7601932) q[1];
rz(-1.3689234) q[3];
sx q[3];
rz(-1.2807506) q[3];
sx q[3];
rz(-1.2253075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.003309) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(-2.8740846) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(2.0013924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32894293) q[0];
sx q[0];
rz(-0.2812627) q[0];
sx q[0];
rz(-1.9360696) q[0];
x q[1];
rz(-2.0613725) q[2];
sx q[2];
rz(-0.51418257) q[2];
sx q[2];
rz(1.0487923) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(1.6608095) q[1];
x q[2];
rz(0.79348989) q[3];
sx q[3];
rz(-1.7975382) q[3];
sx q[3];
rz(1.2236809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(-1.3557419) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-0.5207516) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-2.5614212) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
