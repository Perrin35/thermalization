OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2663015) q[0];
sx q[0];
rz(-0.40691352) q[0];
sx q[0];
rz(0.20242515) q[0];
rz(-1.3099194) q[1];
sx q[1];
rz(-1.1453495) q[1];
sx q[1];
rz(-2.1529012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3920445) q[0];
sx q[0];
rz(-0.77717669) q[0];
sx q[0];
rz(1.2991608) q[0];
x q[1];
rz(-1.7086171) q[2];
sx q[2];
rz(-1.4933961) q[2];
sx q[2];
rz(3.1397506) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9548077) q[1];
sx q[1];
rz(-2.4802172) q[1];
sx q[1];
rz(-1.9438882) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8396138) q[3];
sx q[3];
rz(-0.43987396) q[3];
sx q[3];
rz(-2.0206106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53628355) q[2];
sx q[2];
rz(-0.80159694) q[2];
sx q[2];
rz(-2.0901704) q[2];
rz(1.5031523) q[3];
sx q[3];
rz(-1.4132615) q[3];
sx q[3];
rz(-2.9063291) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1992092) q[0];
sx q[0];
rz(-1.0301882) q[0];
sx q[0];
rz(0.27401608) q[0];
rz(-0.38385299) q[1];
sx q[1];
rz(-0.63261384) q[1];
sx q[1];
rz(1.3207818) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21652554) q[0];
sx q[0];
rz(-1.4122047) q[0];
sx q[0];
rz(0.065404371) q[0];
x q[1];
rz(1.3396367) q[2];
sx q[2];
rz(-1.0826687) q[2];
sx q[2];
rz(1.7423992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40055028) q[1];
sx q[1];
rz(-2.2930751) q[1];
sx q[1];
rz(0.063401179) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76791538) q[3];
sx q[3];
rz(-2.1817099) q[3];
sx q[3];
rz(-3.0141941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5968093) q[2];
sx q[2];
rz(-1.4027255) q[2];
sx q[2];
rz(-1.7421494) q[2];
rz(2.5341189) q[3];
sx q[3];
rz(-1.4676658) q[3];
sx q[3];
rz(-2.7510711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52943277) q[0];
sx q[0];
rz(-2.2859892) q[0];
sx q[0];
rz(-0.31923077) q[0];
rz(0.051479738) q[1];
sx q[1];
rz(-1.1795283) q[1];
sx q[1];
rz(1.0565588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.257911) q[0];
sx q[0];
rz(-1.5192832) q[0];
sx q[0];
rz(1.7345588) q[0];
rz(-pi) q[1];
x q[1];
rz(2.276307) q[2];
sx q[2];
rz(-1.6076536) q[2];
sx q[2];
rz(1.8344022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1358976) q[1];
sx q[1];
rz(-2.4881004) q[1];
sx q[1];
rz(1.280927) q[1];
rz(-pi) q[2];
rz(3.0491054) q[3];
sx q[3];
rz(-1.1690946) q[3];
sx q[3];
rz(0.71850151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1928548) q[2];
sx q[2];
rz(-1.4317908) q[2];
sx q[2];
rz(1.7970435) q[2];
rz(-2.2236842) q[3];
sx q[3];
rz(-2.5036006) q[3];
sx q[3];
rz(1.8857672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767839) q[0];
sx q[0];
rz(-0.5905686) q[0];
sx q[0];
rz(-1.508924) q[0];
rz(0.50679874) q[1];
sx q[1];
rz(-1.9281887) q[1];
sx q[1];
rz(2.7510344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67462639) q[0];
sx q[0];
rz(-1.5635365) q[0];
sx q[0];
rz(-1.566559) q[0];
x q[1];
rz(1.3803595) q[2];
sx q[2];
rz(-1.9445297) q[2];
sx q[2];
rz(0.98183364) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0259121) q[1];
sx q[1];
rz(-1.5424219) q[1];
sx q[1];
rz(2.1898153) q[1];
x q[2];
rz(-2.9831577) q[3];
sx q[3];
rz(-1.6601175) q[3];
sx q[3];
rz(1.4591372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.249923) q[2];
sx q[2];
rz(-0.93081433) q[2];
sx q[2];
rz(-2.2451952) q[2];
rz(2.2253288) q[3];
sx q[3];
rz(-0.93583411) q[3];
sx q[3];
rz(2.1808482) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498098) q[0];
sx q[0];
rz(-0.35583219) q[0];
sx q[0];
rz(0.47019666) q[0];
rz(1.4037941) q[1];
sx q[1];
rz(-2.4411968) q[1];
sx q[1];
rz(-1.2276924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.550023) q[0];
sx q[0];
rz(-1.8409022) q[0];
sx q[0];
rz(0.0065992532) q[0];
rz(-0.63666261) q[2];
sx q[2];
rz(-0.65386727) q[2];
sx q[2];
rz(0.9034397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0684546) q[1];
sx q[1];
rz(-2.2139991) q[1];
sx q[1];
rz(0.094631852) q[1];
rz(2.6484102) q[3];
sx q[3];
rz(-2.5458284) q[3];
sx q[3];
rz(-2.7920134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1710743) q[2];
sx q[2];
rz(-2.4288364) q[2];
sx q[2];
rz(-1.3842545) q[2];
rz(-0.61128831) q[3];
sx q[3];
rz(-1.7620112) q[3];
sx q[3];
rz(0.39613625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1432081) q[0];
sx q[0];
rz(-1.0118326) q[0];
sx q[0];
rz(-0.55214733) q[0];
rz(-2.4119008) q[1];
sx q[1];
rz(-0.98760664) q[1];
sx q[1];
rz(2.139835) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14267966) q[0];
sx q[0];
rz(-2.985654) q[0];
sx q[0];
rz(-1.0264615) q[0];
rz(0.92879734) q[2];
sx q[2];
rz(-1.9812407) q[2];
sx q[2];
rz(-0.14308077) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.525108) q[1];
sx q[1];
rz(-1.387038) q[1];
sx q[1];
rz(0.77145578) q[1];
rz(2.0902315) q[3];
sx q[3];
rz(-0.95963235) q[3];
sx q[3];
rz(0.93944235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6843162) q[2];
sx q[2];
rz(-1.3490973) q[2];
sx q[2];
rz(-1.3433733) q[2];
rz(0.7243048) q[3];
sx q[3];
rz(-1.9380006) q[3];
sx q[3];
rz(-0.30266416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90172076) q[0];
sx q[0];
rz(-1.8093103) q[0];
sx q[0];
rz(-2.4161762) q[0];
rz(-1.2984917) q[1];
sx q[1];
rz(-2.6893171) q[1];
sx q[1];
rz(-2.8628912) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1567247) q[0];
sx q[0];
rz(-1.4188123) q[0];
sx q[0];
rz(-2.5990309) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33228525) q[2];
sx q[2];
rz(-1.9481619) q[2];
sx q[2];
rz(-1.2720296) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6141521) q[1];
sx q[1];
rz(-2.2598055) q[1];
sx q[1];
rz(-0.39871881) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6288644) q[3];
sx q[3];
rz(-0.89539274) q[3];
sx q[3];
rz(1.9844296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7574888) q[2];
sx q[2];
rz(-2.8748547) q[2];
sx q[2];
rz(-0.98814803) q[2];
rz(-0.10495505) q[3];
sx q[3];
rz(-1.2982488) q[3];
sx q[3];
rz(2.9102563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4526378) q[0];
sx q[0];
rz(-0.80113688) q[0];
sx q[0];
rz(-0.6193921) q[0];
rz(0.014892189) q[1];
sx q[1];
rz(-2.2512524) q[1];
sx q[1];
rz(-2.4598918) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084472585) q[0];
sx q[0];
rz(-2.1051877) q[0];
sx q[0];
rz(0.33584331) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4686844) q[2];
sx q[2];
rz(-2.4739304) q[2];
sx q[2];
rz(0.84102902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3577963) q[1];
sx q[1];
rz(-1.2192084) q[1];
sx q[1];
rz(1.5635673) q[1];
rz(-0.19187627) q[3];
sx q[3];
rz(-0.8423281) q[3];
sx q[3];
rz(1.6863562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90998489) q[2];
sx q[2];
rz(-2.3384428) q[2];
sx q[2];
rz(-2.8019359) q[2];
rz(-0.61399442) q[3];
sx q[3];
rz(-1.3795779) q[3];
sx q[3];
rz(1.5617255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851819) q[0];
sx q[0];
rz(-2.7934533) q[0];
sx q[0];
rz(-2.2871995) q[0];
rz(-2.5482381) q[1];
sx q[1];
rz(-2.1095468) q[1];
sx q[1];
rz(-2.8462483) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2345384) q[0];
sx q[0];
rz(-3.0772644) q[0];
sx q[0];
rz(-0.5366908) q[0];
rz(-1.1530442) q[2];
sx q[2];
rz(-1.5502073) q[2];
sx q[2];
rz(-2.0143687) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1251245) q[1];
sx q[1];
rz(-1.09519) q[1];
sx q[1];
rz(0.7866025) q[1];
rz(1.8804824) q[3];
sx q[3];
rz(-1.4638374) q[3];
sx q[3];
rz(2.2992319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9796543) q[2];
sx q[2];
rz(-1.6831968) q[2];
sx q[2];
rz(0.68183199) q[2];
rz(-1.9409059) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(0.43584263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1107263) q[0];
sx q[0];
rz(-1.8616572) q[0];
sx q[0];
rz(-2.6357292) q[0];
rz(-2.1234296) q[1];
sx q[1];
rz(-1.70645) q[1];
sx q[1];
rz(0.86046576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8448526) q[0];
sx q[0];
rz(-0.89549795) q[0];
sx q[0];
rz(1.6313511) q[0];
rz(-0.14218743) q[2];
sx q[2];
rz(-1.280613) q[2];
sx q[2];
rz(-2.1423401) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9221895) q[1];
sx q[1];
rz(-1.8199395) q[1];
sx q[1];
rz(-1.5369028) q[1];
x q[2];
rz(2.1846049) q[3];
sx q[3];
rz(-1.0928983) q[3];
sx q[3];
rz(-2.1769644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0597509) q[2];
sx q[2];
rz(-1.1885208) q[2];
sx q[2];
rz(2.5583978) q[2];
rz(-0.0010631337) q[3];
sx q[3];
rz(-0.93102264) q[3];
sx q[3];
rz(-2.6482436) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4790333) q[0];
sx q[0];
rz(-1.722535) q[0];
sx q[0];
rz(-2.2395635) q[0];
rz(0.61530151) q[1];
sx q[1];
rz(-1.6533783) q[1];
sx q[1];
rz(1.3396214) q[1];
rz(1.3134708) q[2];
sx q[2];
rz(-2.0791292) q[2];
sx q[2];
rz(2.8544478) q[2];
rz(1.2801805) q[3];
sx q[3];
rz(-2.5463085) q[3];
sx q[3];
rz(0.13704722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
