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
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28344261) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(2.0767077) q[0];
rz(-pi) q[1];
rz(-0.067692368) q[2];
sx q[2];
rz(-0.86172416) q[2];
sx q[2];
rz(-2.6807705) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2823892) q[1];
sx q[1];
rz(-0.61835641) q[1];
sx q[1];
rz(2.8137141) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53341289) q[3];
sx q[3];
rz(-2.3331113) q[3];
sx q[3];
rz(-1.4097139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.082829647) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-0.70835152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3268711) q[0];
sx q[0];
rz(-0.46778361) q[0];
sx q[0];
rz(-2.4844556) q[0];
rz(1.7827665) q[2];
sx q[2];
rz(-0.6119298) q[2];
sx q[2];
rz(-1.3509392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-2.0611144) q[1];
sx q[1];
rz(-2.3398188) q[1];
rz(-2.6459341) q[3];
sx q[3];
rz(-2.0303876) q[3];
sx q[3];
rz(1.3575777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11671242) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(0.29176816) q[2];
rz(3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3290688) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-2.5699239) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5105195) q[0];
sx q[0];
rz(-1.6594995) q[0];
sx q[0];
rz(-2.6064794) q[0];
rz(-pi) q[1];
x q[1];
rz(1.13582) q[2];
sx q[2];
rz(-1.8990714) q[2];
sx q[2];
rz(-0.091718397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.09517041) q[1];
sx q[1];
rz(-1.7032402) q[1];
sx q[1];
rz(2.4080647) q[1];
x q[2];
rz(-0.45205558) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(3.0878386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-2.4327915) q[2];
rz(0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(0.79950142) q[0];
rz(-3.0917621) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(2.9611011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364396) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(-2.9050164) q[0];
rz(-1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(0.82496914) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52164493) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(-0.32089969) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1293731) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(1.3752939) q[3];
sx q[3];
rz(-1.8245274) q[3];
sx q[3];
rz(-0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3809526) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-0.25156897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0753206) q[0];
sx q[0];
rz(-2.2499654) q[0];
sx q[0];
rz(-2.9156296) q[0];
rz(2.9676874) q[2];
sx q[2];
rz(-0.46602962) q[2];
sx q[2];
rz(0.58338651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1269762) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(-1.5429392) q[1];
rz(-pi) q[2];
rz(-0.86815636) q[3];
sx q[3];
rz(-1.9725102) q[3];
sx q[3];
rz(-2.274184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(2.1110558) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(-1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(-0.39302557) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5315006) q[0];
sx q[0];
rz(-1.2628265) q[0];
sx q[0];
rz(-0.17225762) q[0];
rz(1.3156462) q[2];
sx q[2];
rz(-0.7609376) q[2];
sx q[2];
rz(-1.4756502) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6833718) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(1.5797735) q[1];
x q[2];
rz(-1.2971446) q[3];
sx q[3];
rz(-1.2840464) q[3];
sx q[3];
rz(1.8991889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3605911) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(2.1684516) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-0.0059676776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3403444) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(-1.9968541) q[0];
rz(-pi) q[1];
rz(2.3890424) q[2];
sx q[2];
rz(-1.9117022) q[2];
sx q[2];
rz(-2.7698851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90970627) q[1];
sx q[1];
rz(-1.5071553) q[1];
sx q[1];
rz(-1.015889) q[1];
x q[2];
rz(-3.0217516) q[3];
sx q[3];
rz(-0.37444515) q[3];
sx q[3];
rz(-0.10711174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(2.4068508) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(0.94820625) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78806879) q[0];
sx q[0];
rz(-0.81809645) q[0];
sx q[0];
rz(-1.3226932) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44342946) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(0.37920096) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1369689) q[1];
sx q[1];
rz(-1.7925295) q[1];
sx q[1];
rz(1.6349413) q[1];
rz(0.35685278) q[3];
sx q[3];
rz(-1.4992504) q[3];
sx q[3];
rz(-0.94716351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.3482288) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(0.25094029) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-2.7889263) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34516016) q[0];
sx q[0];
rz(-0.16819084) q[0];
sx q[0];
rz(-2.3229984) q[0];
rz(-1.4264832) q[2];
sx q[2];
rz(-1.8679973) q[2];
sx q[2];
rz(0.65629634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36665146) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(-1.8371546) q[1];
rz(2.6412233) q[3];
sx q[3];
rz(-2.2420068) q[3];
sx q[3];
rz(-2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(0.68332589) q[2];
rz(1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(-2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-0.93651071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47637981) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(2.5961844) q[0];
rz(2.1328743) q[2];
sx q[2];
rz(-2.808411) q[2];
sx q[2];
rz(-2.9269232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9156956) q[1];
sx q[1];
rz(-0.86663336) q[1];
sx q[1];
rz(-2.8812863) q[1];
rz(-pi) q[2];
rz(1.1141206) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(-0.041681899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(0.014952095) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1508355) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(-1.5785718) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(-1.8216495) q[2];
sx q[2];
rz(-1.432042) q[2];
sx q[2];
rz(-1.5114573) q[2];
rz(-2.7097377) q[3];
sx q[3];
rz(-0.57341822) q[3];
sx q[3];
rz(1.740406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
