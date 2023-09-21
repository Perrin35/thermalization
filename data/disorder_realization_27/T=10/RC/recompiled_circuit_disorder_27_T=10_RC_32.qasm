OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(5.0737557) q[1];
sx q[1];
rz(4.3901246) q[1];
sx q[1];
rz(7.6683383) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614852) q[0];
sx q[0];
rz(-0.84851096) q[0];
sx q[0];
rz(-1.1397584) q[0];
rz(2.7117549) q[2];
sx q[2];
rz(-2.5463383) q[2];
sx q[2];
rz(2.0855479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79917158) q[1];
sx q[1];
rz(-2.85596) q[1];
sx q[1];
rz(-2.530453) q[1];
rz(-pi) q[2];
rz(0.6119421) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(2.2236168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40760621) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(-2.6876887) q[0];
rz(2.1167963) q[1];
sx q[1];
rz(-0.4075993) q[1];
sx q[1];
rz(1.227238) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62118976) q[0];
sx q[0];
rz(-1.495201) q[0];
sx q[0];
rz(1.7002506) q[0];
rz(-pi) q[1];
rz(2.4279459) q[2];
sx q[2];
rz(-0.64672856) q[2];
sx q[2];
rz(1.5681859) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7397346) q[1];
sx q[1];
rz(-1.2793853) q[1];
sx q[1];
rz(-1.8406364) q[1];
x q[2];
rz(-2.9407528) q[3];
sx q[3];
rz(-1.4174995) q[3];
sx q[3];
rz(-2.4917345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(2.5741637) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(-0.89865249) q[0];
rz(-2.1458416) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-2.8083037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22526564) q[0];
sx q[0];
rz(-2.3258381) q[0];
sx q[0];
rz(-2.4832721) q[0];
rz(-pi) q[1];
rz(-1.7477481) q[2];
sx q[2];
rz(-1.4743917) q[2];
sx q[2];
rz(0.63444885) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.052735141) q[1];
sx q[1];
rz(-2.5087068) q[1];
sx q[1];
rz(-1.8295374) q[1];
rz(1.0201449) q[3];
sx q[3];
rz(-1.9198951) q[3];
sx q[3];
rz(2.246644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4553392) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(0.84093705) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(-1.1545198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.44160098) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(-1.0355863) q[1];
sx q[1];
rz(-2.6338449) q[1];
sx q[1];
rz(-1.9365786) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.372615) q[0];
sx q[0];
rz(-0.87991558) q[0];
sx q[0];
rz(3.0983503) q[0];
rz(1.8393458) q[2];
sx q[2];
rz(-2.4797202) q[2];
sx q[2];
rz(2.4115987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0282064) q[1];
sx q[1];
rz(-0.75062597) q[1];
sx q[1];
rz(-1.1651462) q[1];
rz(-pi) q[2];
rz(0.071675008) q[3];
sx q[3];
rz(-2.274401) q[3];
sx q[3];
rz(-0.091761656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-2.7704346) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054759653) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(-0.99331028) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(-0.55508074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13984891) q[0];
sx q[0];
rz(-1.886133) q[0];
sx q[0];
rz(0.01339162) q[0];
rz(-1.416989) q[2];
sx q[2];
rz(-2.7222735) q[2];
sx q[2];
rz(0.16659444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86075393) q[1];
sx q[1];
rz(-2.6959531) q[1];
sx q[1];
rz(-1.526236) q[1];
x q[2];
rz(0.080901905) q[3];
sx q[3];
rz(-1.5187851) q[3];
sx q[3];
rz(-2.4018283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.83539) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(-3.0026657) q[2];
rz(-2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(-2.561835) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.6019843) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67028763) q[0];
sx q[0];
rz(-1.7738713) q[0];
sx q[0];
rz(2.2905486) q[0];
x q[1];
rz(-3.0110353) q[2];
sx q[2];
rz(-1.6160384) q[2];
sx q[2];
rz(-3.0322078) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.08338883) q[1];
sx q[1];
rz(-1.997588) q[1];
sx q[1];
rz(0.31174387) q[1];
rz(-pi) q[2];
rz(-0.88605373) q[3];
sx q[3];
rz(-1.7985117) q[3];
sx q[3];
rz(-0.39623228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(-2.8721151) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(0.68429464) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(-0.51876846) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25012384) q[0];
sx q[0];
rz(-2.2995673) q[0];
sx q[0];
rz(0.17983371) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3709626) q[2];
sx q[2];
rz(-1.4118886) q[2];
sx q[2];
rz(-2.4385914) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3541975) q[1];
sx q[1];
rz(-1.1095611) q[1];
sx q[1];
rz(-2.4128777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026168907) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(-2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8873022) q[2];
sx q[2];
rz(-1.3672978) q[2];
sx q[2];
rz(-2.5781412) q[2];
rz(0.051579483) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(-2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.1812487) q[0];
sx q[0];
rz(-0.5287756) q[0];
sx q[0];
rz(1.7425591) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(0.74434892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5520185) q[0];
sx q[0];
rz(-2.9840143) q[0];
sx q[0];
rz(-0.35430674) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8366351) q[2];
sx q[2];
rz(-2.611534) q[2];
sx q[2];
rz(-2.838138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5163755) q[1];
sx q[1];
rz(-0.55621925) q[1];
sx q[1];
rz(2.7079627) q[1];
x q[2];
rz(-3.0934422) q[3];
sx q[3];
rz(-2.6865494) q[3];
sx q[3];
rz(-2.9330848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7156334) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(2.5320833) q[2];
rz(-2.4842747) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1881926) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(2.3666568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79561728) q[0];
sx q[0];
rz(-0.91751639) q[0];
sx q[0];
rz(-0.4972636) q[0];
rz(2.0140892) q[2];
sx q[2];
rz(-1.7121676) q[2];
sx q[2];
rz(-3.1183426) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1800268) q[1];
sx q[1];
rz(-2.0225836) q[1];
sx q[1];
rz(1.7801442) q[1];
x q[2];
rz(-1.7899412) q[3];
sx q[3];
rz(-1.431576) q[3];
sx q[3];
rz(0.62300357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9541786) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(-0.15787086) q[2];
rz(1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(-2.9272595) q[0];
rz(-0.65746039) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(1.0459895) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5770618) q[0];
sx q[0];
rz(-2.7663167) q[0];
sx q[0];
rz(-0.29348404) q[0];
rz(-0.59098737) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(-1.1514593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54088456) q[1];
sx q[1];
rz(-0.67544671) q[1];
sx q[1];
rz(-0.9687959) q[1];
x q[2];
rz(0.40542116) q[3];
sx q[3];
rz(-1.5463721) q[3];
sx q[3];
rz(-0.072211857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7252698) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-2.4889448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.6090341) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(-1.2621461) q[2];
sx q[2];
rz(-1.9625004) q[2];
sx q[2];
rz(0.99448268) q[2];
rz(3.0388721) q[3];
sx q[3];
rz(-2.5892047) q[3];
sx q[3];
rz(1.8451167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
