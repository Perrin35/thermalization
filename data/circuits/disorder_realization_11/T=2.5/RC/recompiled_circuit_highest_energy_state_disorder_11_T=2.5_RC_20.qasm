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
rz(0.053057916) q[0];
sx q[0];
rz(-2.7795656) q[0];
sx q[0];
rz(1.9757353) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(1.4282164) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7390763) q[0];
sx q[0];
rz(-1.5105643) q[0];
sx q[0];
rz(-2.6291558) q[0];
rz(-2.8901454) q[2];
sx q[2];
rz(-1.5000952) q[2];
sx q[2];
rz(1.155702) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62318351) q[1];
sx q[1];
rz(-1.3819547) q[1];
sx q[1];
rz(-0.13407003) q[1];
rz(-pi) q[2];
rz(0.26784874) q[3];
sx q[3];
rz(-1.0976657) q[3];
sx q[3];
rz(0.79896636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4521788) q[2];
sx q[2];
rz(-0.016760085) q[2];
sx q[2];
rz(0.20081271) q[2];
rz(2.9928442) q[3];
sx q[3];
rz(-3.1368308) q[3];
sx q[3];
rz(-2.8301921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1397322) q[0];
sx q[0];
rz(-0.59356028) q[0];
sx q[0];
rz(-2.1019905) q[0];
rz(3.1270341) q[1];
sx q[1];
rz(-1.2332375) q[1];
sx q[1];
rz(-1.5537517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56180219) q[0];
sx q[0];
rz(-2.1952711) q[0];
sx q[0];
rz(-0.66656163) q[0];
rz(-pi) q[1];
rz(-1.3755685) q[2];
sx q[2];
rz(-0.075610925) q[2];
sx q[2];
rz(1.6603254) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5841516) q[1];
sx q[1];
rz(-1.5366035) q[1];
sx q[1];
rz(-1.8299915) q[1];
rz(1.3782937) q[3];
sx q[3];
rz(-1.8979591) q[3];
sx q[3];
rz(0.006803676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62850922) q[2];
sx q[2];
rz(-1.6078948) q[2];
sx q[2];
rz(1.3879363) q[2];
rz(1.7657109) q[3];
sx q[3];
rz(-1.0466156) q[3];
sx q[3];
rz(-2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7498216) q[0];
sx q[0];
rz(-2.9048558) q[0];
sx q[0];
rz(-2.5340875) q[0];
rz(-1.5433743) q[1];
sx q[1];
rz(-0.18074712) q[1];
sx q[1];
rz(0.9651331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6032219) q[0];
sx q[0];
rz(-0.020909034) q[0];
sx q[0];
rz(1.7666817) q[0];
x q[1];
rz(2.7982462) q[2];
sx q[2];
rz(-1.1330686) q[2];
sx q[2];
rz(0.9870607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0697433) q[1];
sx q[1];
rz(-2.9925821) q[1];
sx q[1];
rz(-1.8144375) q[1];
rz(-pi) q[2];
rz(-1.4749281) q[3];
sx q[3];
rz(-1.6536668) q[3];
sx q[3];
rz(1.7516983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8049916) q[2];
sx q[2];
rz(-2.470863) q[2];
sx q[2];
rz(0.86656183) q[2];
rz(1.1066412) q[3];
sx q[3];
rz(-1.5525147) q[3];
sx q[3];
rz(1.6710056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5690145) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(-1.6160075) q[0];
rz(3.1309879) q[1];
sx q[1];
rz(-3.1378101) q[1];
sx q[1];
rz(-2.3912281) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.01975) q[0];
sx q[0];
rz(-1.6656309) q[0];
sx q[0];
rz(-0.66480831) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32230538) q[2];
sx q[2];
rz(-2.5111329) q[2];
sx q[2];
rz(0.96994627) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4669686) q[1];
sx q[1];
rz(-1.4556938) q[1];
sx q[1];
rz(-0.50927598) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31171215) q[3];
sx q[3];
rz(-1.4501713) q[3];
sx q[3];
rz(-1.3159304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7287207) q[2];
sx q[2];
rz(-2.0527288) q[2];
sx q[2];
rz(1.9000165) q[2];
rz(0.0056754644) q[3];
sx q[3];
rz(-2.3311876) q[3];
sx q[3];
rz(-2.2976105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4960957) q[0];
sx q[0];
rz(-3.0912919) q[0];
sx q[0];
rz(0.92329931) q[0];
rz(-0.80054379) q[1];
sx q[1];
rz(-0.0034595483) q[1];
sx q[1];
rz(-2.961535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6430017) q[0];
sx q[0];
rz(-0.24234903) q[0];
sx q[0];
rz(1.3968299) q[0];
rz(-pi) q[1];
rz(-3.0434612) q[2];
sx q[2];
rz(-1.4906825) q[2];
sx q[2];
rz(2.0126041) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1504768) q[1];
sx q[1];
rz(-1.2186945) q[1];
sx q[1];
rz(1.8262509) q[1];
x q[2];
rz(1.8534142) q[3];
sx q[3];
rz(-0.37943951) q[3];
sx q[3];
rz(-1.5490378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8397612) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(1.4361471) q[2];
rz(-1.885421) q[3];
sx q[3];
rz(-1.4830282) q[3];
sx q[3];
rz(-3.0459611) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65211463) q[0];
sx q[0];
rz(-0.57179946) q[0];
sx q[0];
rz(2.6982464) q[0];
rz(1.1706785) q[1];
sx q[1];
rz(-3.1405293) q[1];
sx q[1];
rz(-2.7098157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74878377) q[0];
sx q[0];
rz(-1.8816299) q[0];
sx q[0];
rz(0.75324599) q[0];
rz(-pi) q[1];
rz(2.0197996) q[2];
sx q[2];
rz(-2.9509955) q[2];
sx q[2];
rz(-2.7191741) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1237632) q[1];
sx q[1];
rz(-2.4430954) q[1];
sx q[1];
rz(-2.5269954) q[1];
x q[2];
rz(2.5434142) q[3];
sx q[3];
rz(-2.7477187) q[3];
sx q[3];
rz(-1.0306213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7410437) q[2];
sx q[2];
rz(-2.2811175) q[2];
sx q[2];
rz(-1.7575556) q[2];
rz(2.4436229) q[3];
sx q[3];
rz(-0.75708404) q[3];
sx q[3];
rz(0.32148263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31239241) q[0];
sx q[0];
rz(-1.7990524) q[0];
sx q[0];
rz(-0.37305748) q[0];
rz(-2.864605) q[1];
sx q[1];
rz(-0.00033683446) q[1];
sx q[1];
rz(-2.3854947) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38967237) q[0];
sx q[0];
rz(-0.7451267) q[0];
sx q[0];
rz(0.57655277) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.937708) q[2];
sx q[2];
rz(-1.9607301) q[2];
sx q[2];
rz(-2.6904047) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4226729) q[1];
sx q[1];
rz(-0.70168272) q[1];
sx q[1];
rz(0.91928457) q[1];
x q[2];
rz(-0.30310615) q[3];
sx q[3];
rz(-1.167093) q[3];
sx q[3];
rz(-1.2991326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21195352) q[2];
sx q[2];
rz(-0.55277199) q[2];
sx q[2];
rz(2.2201404) q[2];
rz(-0.067342162) q[3];
sx q[3];
rz(-1.9332956) q[3];
sx q[3];
rz(1.6578081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9901554) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(0.13023278) q[0];
rz(0.79365927) q[1];
sx q[1];
rz(-0.001611324) q[1];
sx q[1];
rz(-0.28009716) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1156688) q[0];
sx q[0];
rz(-1.6243132) q[0];
sx q[0];
rz(1.6497067) q[0];
x q[1];
rz(2.75861) q[2];
sx q[2];
rz(-2.9071701) q[2];
sx q[2];
rz(-1.926926) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62781266) q[1];
sx q[1];
rz(-2.4901721) q[1];
sx q[1];
rz(-0.86077229) q[1];
x q[2];
rz(1.0956826) q[3];
sx q[3];
rz(-0.91714782) q[3];
sx q[3];
rz(-0.16330367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.085999504) q[2];
sx q[2];
rz(-1.4693825) q[2];
sx q[2];
rz(-2.3386193) q[2];
rz(1.5351852) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(0.83139658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0292173) q[0];
sx q[0];
rz(-3.1387098) q[0];
sx q[0];
rz(3.032384) q[0];
rz(-2.7681007) q[1];
sx q[1];
rz(-1.9544574) q[1];
sx q[1];
rz(-2.5901897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5057988) q[0];
sx q[0];
rz(-2.3058545) q[0];
sx q[0];
rz(-2.4018025) q[0];
rz(-0.18420561) q[2];
sx q[2];
rz(-1.6937243) q[2];
sx q[2];
rz(2.9254928) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.046519) q[1];
sx q[1];
rz(-2.1035026) q[1];
sx q[1];
rz(-2.2979876) q[1];
x q[2];
rz(1.3250454) q[3];
sx q[3];
rz(-1.1781404) q[3];
sx q[3];
rz(2.806854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6138844) q[2];
sx q[2];
rz(-1.8324499) q[2];
sx q[2];
rz(-1.323553) q[2];
rz(1.2644794) q[3];
sx q[3];
rz(-1.8529961) q[3];
sx q[3];
rz(0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5962113) q[0];
sx q[0];
rz(-0.63544202) q[0];
sx q[0];
rz(2.4052461) q[0];
rz(-2.9296854) q[1];
sx q[1];
rz(-1.0286464) q[1];
sx q[1];
rz(-1.5440936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6074267) q[0];
sx q[0];
rz(-1.5527927) q[0];
sx q[0];
rz(0.27915633) q[0];
x q[1];
rz(-3.1112291) q[2];
sx q[2];
rz(-1.5734451) q[2];
sx q[2];
rz(-2.9752258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7791833) q[1];
sx q[1];
rz(-1.1115555) q[1];
sx q[1];
rz(-0.68024858) q[1];
rz(-2.6968217) q[3];
sx q[3];
rz(-1.8666779) q[3];
sx q[3];
rz(0.55287655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3046917) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(1.2738073) q[2];
rz(1.4443719) q[3];
sx q[3];
rz(-3.0675409) q[3];
sx q[3];
rz(1.4950289) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16784167) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(1.5372859) q[1];
sx q[1];
rz(-0.91265408) q[1];
sx q[1];
rz(0.18462054) q[1];
rz(0.040097728) q[2];
sx q[2];
rz(-1.5798777) q[2];
sx q[2];
rz(-0.29691534) q[2];
rz(1.7459695) q[3];
sx q[3];
rz(-0.61476954) q[3];
sx q[3];
rz(0.7634306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
