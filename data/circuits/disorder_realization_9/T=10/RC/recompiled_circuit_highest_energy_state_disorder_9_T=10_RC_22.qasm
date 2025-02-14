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
rz(0.23586805) q[0];
sx q[0];
rz(-1.1455102) q[0];
sx q[0];
rz(0.048576485) q[0];
rz(-1.6254758) q[1];
sx q[1];
rz(-1.0928417) q[1];
sx q[1];
rz(0.71576524) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019139394) q[0];
sx q[0];
rz(-1.4000443) q[0];
sx q[0];
rz(-3.1117361) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2907998) q[2];
sx q[2];
rz(-1.7699307) q[2];
sx q[2];
rz(3.0906349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29808334) q[1];
sx q[1];
rz(-2.2560826) q[1];
sx q[1];
rz(2.09649) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51718398) q[3];
sx q[3];
rz(-1.9265012) q[3];
sx q[3];
rz(-2.9773447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6526661) q[2];
sx q[2];
rz(-2.581614) q[2];
sx q[2];
rz(1.0555335) q[2];
rz(0.40600285) q[3];
sx q[3];
rz(-1.8279653) q[3];
sx q[3];
rz(-0.82160151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15776289) q[0];
sx q[0];
rz(-1.0394179) q[0];
sx q[0];
rz(-2.5854172) q[0];
rz(-1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(0.083273085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10262415) q[0];
sx q[0];
rz(-2.4100465) q[0];
sx q[0];
rz(-0.49517531) q[0];
rz(-pi) q[1];
rz(2.3131392) q[2];
sx q[2];
rz(-0.76702416) q[2];
sx q[2];
rz(1.7491064) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1494807) q[1];
sx q[1];
rz(-0.50051033) q[1];
sx q[1];
rz(0.47187279) q[1];
rz(-0.49442716) q[3];
sx q[3];
rz(-1.0599679) q[3];
sx q[3];
rz(-1.1100681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1978153) q[2];
sx q[2];
rz(-1.752172) q[2];
sx q[2];
rz(-0.49276349) q[2];
rz(2.1150151) q[3];
sx q[3];
rz(-0.30288282) q[3];
sx q[3];
rz(-1.7290285) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76348412) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(0.43877959) q[0];
rz(-1.6890866) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(-2.7486393) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753112) q[0];
sx q[0];
rz(-1.5428758) q[0];
sx q[0];
rz(-2.0620538) q[0];
rz(-pi) q[1];
rz(-0.032756373) q[2];
sx q[2];
rz(-2.0989145) q[2];
sx q[2];
rz(-0.32450766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71686983) q[1];
sx q[1];
rz(-2.7174414) q[1];
sx q[1];
rz(1.39759) q[1];
x q[2];
rz(1.4912823) q[3];
sx q[3];
rz(-0.92303716) q[3];
sx q[3];
rz(-0.11707725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4517333) q[2];
sx q[2];
rz(-2.2971575) q[2];
sx q[2];
rz(-2.8680657) q[2];
rz(-0.61255974) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(0.85649049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3278811) q[0];
sx q[0];
rz(-0.59309816) q[0];
sx q[0];
rz(-2.3408422) q[0];
rz(-0.71209359) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(-0.48316479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8034604) q[0];
sx q[0];
rz(-2.3134941) q[0];
sx q[0];
rz(-1.5716193) q[0];
rz(-pi) q[1];
rz(-2.4078385) q[2];
sx q[2];
rz(-0.68063191) q[2];
sx q[2];
rz(1.5891927) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8193244) q[1];
sx q[1];
rz(-1.5419648) q[1];
sx q[1];
rz(-2.4954456) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9075059) q[3];
sx q[3];
rz(-0.90970618) q[3];
sx q[3];
rz(-2.8291836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98372769) q[2];
sx q[2];
rz(-2.4309776) q[2];
sx q[2];
rz(1.857081) q[2];
rz(2.6537248) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(1.454486) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81216413) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(-0.46464768) q[0];
rz(2.1931785) q[1];
sx q[1];
rz(-0.83810884) q[1];
sx q[1];
rz(1.4533739) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68935716) q[0];
sx q[0];
rz(-0.043137155) q[0];
sx q[0];
rz(-2.1998911) q[0];
rz(-pi) q[1];
rz(2.4011432) q[2];
sx q[2];
rz(-1.3794823) q[2];
sx q[2];
rz(-2.180661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5504811) q[1];
sx q[1];
rz(-1.5036991) q[1];
sx q[1];
rz(-0.61136742) q[1];
rz(-pi) q[2];
rz(2.0876168) q[3];
sx q[3];
rz(-2.2272416) q[3];
sx q[3];
rz(2.2051728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98443085) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(-2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30094576) q[0];
sx q[0];
rz(-0.41579682) q[0];
sx q[0];
rz(2.0796602) q[0];
rz(-0.95147079) q[1];
sx q[1];
rz(-0.87363344) q[1];
sx q[1];
rz(-1.5207312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5439107) q[0];
sx q[0];
rz(-1.8476227) q[0];
sx q[0];
rz(3.0863161) q[0];
rz(-pi) q[1];
rz(2.4302519) q[2];
sx q[2];
rz(-0.78186505) q[2];
sx q[2];
rz(0.021325354) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6013612) q[1];
sx q[1];
rz(-1.4077606) q[1];
sx q[1];
rz(-2.3705179) q[1];
rz(0.49455182) q[3];
sx q[3];
rz(-2.7003717) q[3];
sx q[3];
rz(-1.0426723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1157397) q[2];
sx q[2];
rz(-1.226475) q[2];
sx q[2];
rz(-2.8958877) q[2];
rz(-1.4541516) q[3];
sx q[3];
rz(-1.5116296) q[3];
sx q[3];
rz(-2.3112442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65569735) q[0];
sx q[0];
rz(-2.4035154) q[0];
sx q[0];
rz(0.6947211) q[0];
rz(-3.0409536) q[1];
sx q[1];
rz(-1.0357608) q[1];
sx q[1];
rz(2.2763841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14510205) q[0];
sx q[0];
rz(-1.1528413) q[0];
sx q[0];
rz(1.6161902) q[0];
rz(-0.91036441) q[2];
sx q[2];
rz(-1.5857134) q[2];
sx q[2];
rz(2.2098324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6859436) q[1];
sx q[1];
rz(-2.5091617) q[1];
sx q[1];
rz(-1.3157822) q[1];
rz(-pi) q[2];
rz(2.4860704) q[3];
sx q[3];
rz(-0.78981863) q[3];
sx q[3];
rz(1.569455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.18409099) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(-0.91020477) q[2];
rz(1.6893859) q[3];
sx q[3];
rz(-1.2319177) q[3];
sx q[3];
rz(-0.46428251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24377395) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(2.5183103) q[0];
rz(0.15577623) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(0.26330858) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0000738) q[0];
sx q[0];
rz(-2.9263068) q[0];
sx q[0];
rz(-1.5207284) q[0];
x q[1];
rz(0.66923334) q[2];
sx q[2];
rz(-0.097429052) q[2];
sx q[2];
rz(-0.92684666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.840082) q[1];
sx q[1];
rz(-0.58159822) q[1];
sx q[1];
rz(0.14029293) q[1];
rz(0.54375387) q[3];
sx q[3];
rz(-2.2787146) q[3];
sx q[3];
rz(1.5073564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89173633) q[2];
sx q[2];
rz(-0.24954924) q[2];
sx q[2];
rz(0.88805324) q[2];
rz(0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(2.1685062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71858281) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(2.6408559) q[0];
rz(-2.0210733) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(2.3401071) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716176) q[0];
sx q[0];
rz(-1.5688251) q[0];
sx q[0];
rz(-0.070974102) q[0];
rz(0.5653462) q[2];
sx q[2];
rz(-1.1156811) q[2];
sx q[2];
rz(-2.3836977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7940143) q[1];
sx q[1];
rz(-1.9982575) q[1];
sx q[1];
rz(-2.4540837) q[1];
x q[2];
rz(-2.9467907) q[3];
sx q[3];
rz(-1.8240989) q[3];
sx q[3];
rz(1.7326488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72006172) q[2];
sx q[2];
rz(-1.820182) q[2];
sx q[2];
rz(0.23019543) q[2];
rz(-3.1318393) q[3];
sx q[3];
rz(-1.516927) q[3];
sx q[3];
rz(-0.15596381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941876) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(0.76960027) q[0];
rz(0.96254483) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(-2.1544971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1302764) q[0];
sx q[0];
rz(-0.53142953) q[0];
sx q[0];
rz(1.7087144) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9883698) q[2];
sx q[2];
rz(-1.330687) q[2];
sx q[2];
rz(-0.78165141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0278575) q[1];
sx q[1];
rz(-3.0887071) q[1];
sx q[1];
rz(-2.2478893) q[1];
rz(2.3651667) q[3];
sx q[3];
rz(-1.2033495) q[3];
sx q[3];
rz(1.1212136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7139682) q[2];
sx q[2];
rz(-0.93728137) q[2];
sx q[2];
rz(3.0734708) q[2];
rz(-2.6918329) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(-1.4382039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.149067) q[0];
sx q[0];
rz(-1.5328007) q[0];
sx q[0];
rz(-1.3971064) q[0];
rz(-2.8515011) q[1];
sx q[1];
rz(-0.42872226) q[1];
sx q[1];
rz(0.52679481) q[1];
rz(0.59103237) q[2];
sx q[2];
rz(-2.5606511) q[2];
sx q[2];
rz(-2.64369) q[2];
rz(-1.5519334) q[3];
sx q[3];
rz(-1.1723779) q[3];
sx q[3];
rz(1.0642024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
