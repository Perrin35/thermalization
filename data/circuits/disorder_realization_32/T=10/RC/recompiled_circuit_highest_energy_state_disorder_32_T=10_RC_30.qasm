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
rz(-1.929317) q[0];
sx q[0];
rz(-1.1317929) q[0];
sx q[0];
rz(-1.3728859) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(4.5166587) q[1];
sx q[1];
rz(5.4969129) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4117811) q[0];
sx q[0];
rz(-1.775911) q[0];
sx q[0];
rz(1.1954444) q[0];
x q[1];
rz(0.9144056) q[2];
sx q[2];
rz(-1.1160276) q[2];
sx q[2];
rz(0.17780534) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3481625) q[1];
sx q[1];
rz(-1.600666) q[1];
sx q[1];
rz(1.8601599) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0065828) q[3];
sx q[3];
rz(-0.67970961) q[3];
sx q[3];
rz(-0.93384472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11975153) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(1.8915668) q[2];
rz(0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64330548) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(0.071320891) q[0];
rz(-1.9885063) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(2.6128795) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3246177) q[0];
sx q[0];
rz(-1.8864838) q[0];
sx q[0];
rz(-1.0782918) q[0];
rz(0.18888338) q[2];
sx q[2];
rz(-1.802465) q[2];
sx q[2];
rz(-2.846039) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9897044) q[1];
sx q[1];
rz(-0.80422771) q[1];
sx q[1];
rz(-2.3249435) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1558481) q[3];
sx q[3];
rz(-1.364721) q[3];
sx q[3];
rz(1.8725439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(-3.1129692) q[2];
rz(-2.7875767) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(-2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(-2.2502374) q[0];
rz(0.50093961) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0863773) q[0];
sx q[0];
rz(-0.098345938) q[0];
sx q[0];
rz(2.6340061) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0979314) q[2];
sx q[2];
rz(-1.322515) q[2];
sx q[2];
rz(0.52094007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8938287) q[1];
sx q[1];
rz(-2.2732452) q[1];
sx q[1];
rz(0.6599627) q[1];
x q[2];
rz(-2.9611582) q[3];
sx q[3];
rz(-1.5984319) q[3];
sx q[3];
rz(2.7829426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3567051) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(0.2365665) q[2];
rz(-2.2208354) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9193566) q[0];
sx q[0];
rz(-1.8574497) q[0];
sx q[0];
rz(2.2663569) q[0];
rz(-1.5454166) q[1];
sx q[1];
rz(-0.86219209) q[1];
sx q[1];
rz(-0.045305591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52462477) q[0];
sx q[0];
rz(-3.0201206) q[0];
sx q[0];
rz(-2.139702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1854109) q[2];
sx q[2];
rz(-0.76757694) q[2];
sx q[2];
rz(1.1078579) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2341174) q[1];
sx q[1];
rz(-1.6576515) q[1];
sx q[1];
rz(-2.240827) q[1];
x q[2];
rz(-1.8604467) q[3];
sx q[3];
rz(-2.5389606) q[3];
sx q[3];
rz(-1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8594325) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(1.947594) q[2];
rz(-0.89030877) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10876656) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(2.4591675) q[0];
rz(-1.8531063) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.3312181) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82560655) q[0];
sx q[0];
rz(-2.2974112) q[0];
sx q[0];
rz(0.021922317) q[0];
x q[1];
rz(-1.17679) q[2];
sx q[2];
rz(-1.8117684) q[2];
sx q[2];
rz(-2.7888128) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9115548) q[1];
sx q[1];
rz(-1.7111519) q[1];
sx q[1];
rz(2.7799003) q[1];
rz(2.6909037) q[3];
sx q[3];
rz(-2.0540049) q[3];
sx q[3];
rz(0.86408981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1149301) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(2.4746573) q[2];
rz(-2.5866348) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(-2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.32894593) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(0.70372787) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(-0.15883787) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1554398) q[0];
sx q[0];
rz(-1.6912529) q[0];
sx q[0];
rz(-1.7140165) q[0];
rz(-pi) q[1];
rz(-2.9877404) q[2];
sx q[2];
rz(-1.6922608) q[2];
sx q[2];
rz(0.50869321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68137121) q[1];
sx q[1];
rz(-1.4386144) q[1];
sx q[1];
rz(-0.48355196) q[1];
rz(-0.30384003) q[3];
sx q[3];
rz(-0.66038495) q[3];
sx q[3];
rz(-2.8840898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3682897) q[2];
sx q[2];
rz(-2.1998019) q[2];
sx q[2];
rz(0.85757315) q[2];
rz(1.9722021) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9661949) q[0];
sx q[0];
rz(-2.3747787) q[0];
sx q[0];
rz(2.4229557) q[0];
rz(2.8596558) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(-0.66863543) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40067264) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(2.7223865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.093542) q[2];
sx q[2];
rz(-2.4562533) q[2];
sx q[2];
rz(2.6743741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9397647) q[1];
sx q[1];
rz(-1.9376905) q[1];
sx q[1];
rz(0.35122996) q[1];
rz(2.5103522) q[3];
sx q[3];
rz(-1.5641912) q[3];
sx q[3];
rz(-0.94830482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64688993) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(1.0364214) q[2];
rz(-2.5070665) q[3];
sx q[3];
rz(-0.98792881) q[3];
sx q[3];
rz(2.3488267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6805639) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(2.4155937) q[0];
rz(1.1622102) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(1.9451709) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0045753) q[0];
sx q[0];
rz(-1.5992303) q[0];
sx q[0];
rz(1.7445376) q[0];
x q[1];
rz(0.047646626) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(-0.36587151) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6641621) q[1];
sx q[1];
rz(-2.4242085) q[1];
sx q[1];
rz(2.391361) q[1];
rz(0.4287339) q[3];
sx q[3];
rz(-2.4247243) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90765816) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(0.93351239) q[2];
rz(2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1652949) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(-0.053255178) q[0];
rz(-0.28757295) q[1];
sx q[1];
rz(-2.2122999) q[1];
sx q[1];
rz(-2.8823749) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142562) q[0];
sx q[0];
rz(-1.5467318) q[0];
sx q[0];
rz(1.5177478) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0854112) q[2];
sx q[2];
rz(-0.26517235) q[2];
sx q[2];
rz(-2.3725703) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6280243) q[1];
sx q[1];
rz(-1.6322989) q[1];
sx q[1];
rz(2.3223367) q[1];
rz(1.7227371) q[3];
sx q[3];
rz(-0.81063089) q[3];
sx q[3];
rz(-1.5347655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5549434) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-0.1813691) q[2];
rz(0.3950611) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(-0.980353) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3996537) q[0];
sx q[0];
rz(-1.5927277) q[0];
sx q[0];
rz(-2.7594866) q[0];
rz(0.29640472) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(1.872725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.554396) q[0];
sx q[0];
rz(-0.66080873) q[0];
sx q[0];
rz(-2.3095678) q[0];
rz(1.7364794) q[2];
sx q[2];
rz(-2.4615917) q[2];
sx q[2];
rz(-2.9939637) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2603392) q[1];
sx q[1];
rz(-1.9104382) q[1];
sx q[1];
rz(3.1189671) q[1];
rz(-pi) q[2];
rz(-2.6878854) q[3];
sx q[3];
rz(-1.6804916) q[3];
sx q[3];
rz(0.060893313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84539139) q[2];
sx q[2];
rz(-0.56075823) q[2];
sx q[2];
rz(0.38983795) q[2];
rz(1.2975533) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(2.6806954) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(1.233558) q[3];
sx q[3];
rz(-2.3261286) q[3];
sx q[3];
rz(0.89589768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
