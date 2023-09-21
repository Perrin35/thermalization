OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9665943) q[0];
sx q[0];
rz(-2.7881665) q[0];
sx q[0];
rz(2.0768291) q[0];
rz(3.937768) q[1];
sx q[1];
rz(1.9328971) q[1];
sx q[1];
rz(9.9608496) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85815) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(1.0648849) q[0];
x q[1];
rz(0.86059086) q[2];
sx q[2];
rz(-1.5194367) q[2];
sx q[2];
rz(-2.0757338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8592035) q[1];
sx q[1];
rz(-2.5232362) q[1];
sx q[1];
rz(0.32787852) q[1];
rz(-pi) q[2];
rz(-1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(-2.1171452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(-0.25201592) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3268711) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(2.4844556) q[0];
rz(-pi) q[1];
rz(-0.96946851) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(-3.0960992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(2.3398188) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80531081) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(-3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(2.5699239) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088061995) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(0.17266973) q[0];
rz(-0.89102913) q[2];
sx q[2];
rz(-0.53855145) q[2];
sx q[2];
rz(-1.0559168) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.09517041) q[1];
sx q[1];
rz(-1.7032402) q[1];
sx q[1];
rz(2.4080647) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27922697) q[3];
sx q[3];
rz(-2.673827) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-0.70880115) q[2];
rz(0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(3.0917621) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(0.18049151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1946963) q[0];
sx q[0];
rz(-1.6616271) q[0];
sx q[0];
rz(0.38740654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7454342) q[2];
sx q[2];
rz(-1.2067814) q[2];
sx q[2];
rz(-0.82496914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6199477) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(-2.820693) q[1];
rz(-pi) q[2];
rz(1.5482076) q[3];
sx q[3];
rz(-2.6456607) q[3];
sx q[3];
rz(-3.0043234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-0.93262514) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809526) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(0.25156897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895296) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(-2.2625838) q[0];
rz(-pi) q[1];
rz(2.6816363) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(2.3098582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.014616414) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(1.5429392) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1523347) q[3];
sx q[3];
rz(-2.3495418) q[3];
sx q[3];
rz(-2.0056412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(2.1110558) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-2.5103536) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0133936) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(1.2584932) q[0];
x q[1];
rz(0.82627798) q[2];
sx q[2];
rz(-1.3958566) q[2];
sx q[2];
rz(0.28184055) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6720851) q[1];
sx q[1];
rz(-0.91963327) q[1];
sx q[1];
rz(0.0068411946) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29720184) q[3];
sx q[3];
rz(-1.8330049) q[3];
sx q[3];
rz(-2.8924243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(-0.97314107) q[2];
rz(2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(0.0059676776) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(1.1447385) q[0];
rz(-pi) q[1];
rz(-1.1184095) q[2];
sx q[2];
rz(-0.87086073) q[2];
sx q[2];
rz(-1.6391022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2318864) q[1];
sx q[1];
rz(-1.5071553) q[1];
sx q[1];
rz(-1.015889) q[1];
x q[2];
rz(1.5238477) q[3];
sx q[3];
rz(-1.1991683) q[3];
sx q[3];
rz(-0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(0.19006426) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-2.1933864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873916) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(-0.76822922) q[0];
rz(0.87585978) q[2];
sx q[2];
rz(-2.5033853) q[2];
sx q[2];
rz(1.1832331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5612948) q[1];
sx q[1];
rz(-1.6333688) q[1];
sx q[1];
rz(-2.9194174) q[1];
rz(2.7847399) q[3];
sx q[3];
rz(-1.6423422) q[3];
sx q[3];
rz(-0.94716351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-2.7622973) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(2.7889263) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1044554) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(-3.0260968) q[0];
rz(-pi) q[1];
rz(-1.7151095) q[2];
sx q[2];
rz(-1.8679973) q[2];
sx q[2];
rz(2.4852963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7949617) q[1];
sx q[1];
rz(-1.7966086) q[1];
sx q[1];
rz(0.56983106) q[1];
rz(-pi) q[2];
rz(2.6412233) q[3];
sx q[3];
rz(-2.2420068) q[3];
sx q[3];
rz(-2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(-2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47637981) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(-0.54540821) q[0];
rz(1.0087183) q[2];
sx q[2];
rz(-2.808411) q[2];
sx q[2];
rz(2.9269232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5156538) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(-2.291912) q[1];
rz(2.0274721) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(0.014952095) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99075714) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(-1.5785718) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(1.8216495) q[2];
sx q[2];
rz(-1.7095507) q[2];
sx q[2];
rz(1.6301353) q[2];
rz(0.43185497) q[3];
sx q[3];
rz(-0.57341822) q[3];
sx q[3];
rz(1.740406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
