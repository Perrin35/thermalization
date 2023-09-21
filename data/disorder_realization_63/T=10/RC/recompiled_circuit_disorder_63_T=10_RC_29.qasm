OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(-0.35342616) q[0];
sx q[0];
rz(1.0647635) q[0];
rz(0.79617533) q[1];
sx q[1];
rz(-1.9328971) q[1];
sx q[1];
rz(-2.605521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85815) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(2.0767077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4921161) q[2];
sx q[2];
rz(-0.71173758) q[2];
sx q[2];
rz(-0.56456883) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2823892) q[1];
sx q[1];
rz(-2.5232362) q[1];
sx q[1];
rz(0.32787852) q[1];
rz(1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(-1.0244474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7321695) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(0.25201592) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(0.70835152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3591374) q[0];
sx q[0];
rz(-1.8498427) q[0];
sx q[0];
rz(-0.38048394) q[0];
x q[1];
rz(-0.14658908) q[2];
sx q[2];
rz(-2.1671038) q[2];
sx q[2];
rz(1.6080315) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(0.80177387) q[1];
rz(-pi) q[2];
rz(-0.80531081) q[3];
sx q[3];
rz(-0.66262965) q[3];
sx q[3];
rz(-0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(0.29176816) q[2];
rz(-0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(0.30971757) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(-0.57166878) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11274352) q[0];
sx q[0];
rz(-2.1035806) q[0];
sx q[0];
rz(1.6738196) q[0];
rz(2.0057726) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(-0.091718397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5204822) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(2.945167) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8623657) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(-1.8750909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16033515) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4335094) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-0.18049151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4806992) q[0];
sx q[0];
rz(-1.9565214) q[0];
sx q[0];
rz(-1.4727403) q[0];
rz(-2.7137202) q[2];
sx q[2];
rz(-2.7395436) q[2];
sx q[2];
rz(-2.7767162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6199477) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(2.820693) q[1];
rz(-pi) q[2];
rz(3.1293731) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(-0.11158768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.583741) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-2.8900237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35206301) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(-0.87900889) q[0];
rz(1.4839843) q[2];
sx q[2];
rz(-2.0292536) q[2];
sx q[2];
rz(0.77755962) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5800036) q[1];
sx q[1];
rz(-1.5434693) q[1];
sx q[1];
rz(2.9461529) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5079481) q[3];
sx q[3];
rz(-0.93379279) q[3];
sx q[3];
rz(0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128199) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(1.2584932) q[0];
x q[1];
rz(-2.9057301) q[2];
sx q[2];
rz(-2.3013407) q[2];
sx q[2];
rz(-2.0116531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46950754) q[1];
sx q[1];
rz(-0.91963327) q[1];
sx q[1];
rz(-0.0068411946) q[1];
x q[2];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.8330049) q[3];
sx q[3];
rz(-0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(3.135625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8864266) q[0];
sx q[0];
rz(-0.92029858) q[0];
sx q[0];
rz(-0.35264539) q[0];
rz(2.0231831) q[2];
sx q[2];
rz(-0.87086073) q[2];
sx q[2];
rz(1.5024904) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5827427) q[1];
sx q[1];
rz(-0.55816459) q[1];
sx q[1];
rz(1.4504257) q[1];
rz(-1.617745) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(-2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(2.1933864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43319234) q[0];
sx q[0];
rz(-0.78484479) q[0];
sx q[0];
rz(2.8851896) q[0];
rz(-2.2657329) q[2];
sx q[2];
rz(-2.5033853) q[2];
sx q[2];
rz(-1.9583595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0046237) q[1];
sx q[1];
rz(-1.7925295) q[1];
sx q[1];
rz(-1.5066513) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4944581) q[3];
sx q[3];
rz(-1.2148972) q[3];
sx q[3];
rz(2.5446041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(0.24027696) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.3482288) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76294476) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6610722) q[0];
sx q[0];
rz(-1.6854223) q[0];
sx q[0];
rz(1.6941487) q[0];
x q[1];
rz(-2.8414531) q[2];
sx q[2];
rz(-1.7087414) q[2];
sx q[2];
rz(-2.2696242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7749412) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(-1.304438) q[1];
rz(-1.0274067) q[3];
sx q[3];
rz(-0.81334844) q[3];
sx q[3];
rz(-0.020997626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(2.4582668) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(-1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.8413683) q[0];
rz(2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(0.93651071) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721601) q[0];
sx q[0];
rz(-2.5494808) q[0];
sx q[0];
rz(-0.44606146) q[0];
rz(1.8556701) q[2];
sx q[2];
rz(-1.7459918) q[2];
sx q[2];
rz(-0.81923649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.3733555) q[1];
sx q[1];
rz(2.291912) q[1];
rz(-pi) q[2];
rz(-2.1533222) q[3];
sx q[3];
rz(-1.8347782) q[3];
sx q[3];
rz(-1.9895944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-0.014952095) q[2];
rz(0.22710083) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.2623513) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1508355) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(-1.5630209) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(-2.0832534) q[2];
sx q[2];
rz(-0.28596157) q[2];
sx q[2];
rz(2.7059976) q[2];
rz(2.7097377) q[3];
sx q[3];
rz(-2.5681744) q[3];
sx q[3];
rz(-1.4011866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];