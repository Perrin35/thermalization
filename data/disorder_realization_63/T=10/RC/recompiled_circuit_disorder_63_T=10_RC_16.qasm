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
rz(3.937768) q[1];
sx q[1];
rz(1.9328971) q[1];
sx q[1];
rz(9.9608496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28344261) q[0];
sx q[0];
rz(-1.3387696) q[0];
sx q[0];
rz(-2.0767077) q[0];
rz(-pi) q[1];
rz(1.4921161) q[2];
sx q[2];
rz(-2.4298551) q[2];
sx q[2];
rz(0.56456883) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2546665) q[1];
sx q[1];
rz(-0.98985043) q[1];
sx q[1];
rz(1.7960153) q[1];
x q[2];
rz(1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(-1.0244474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(1.0268964) q[2];
rz(0.25201592) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10193292) q[0];
sx q[0];
rz(-1.9358557) q[0];
sx q[0];
rz(-1.2714766) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9950036) q[2];
sx q[2];
rz(-2.1671038) q[2];
sx q[2];
rz(-1.6080315) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0185512) q[1];
sx q[1];
rz(-2.0611144) q[1];
sx q[1];
rz(-0.80177387) q[1];
x q[2];
rz(-2.0833011) q[3];
sx q[3];
rz(-2.0111492) q[3];
sx q[3];
rz(0.022170598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(3.0388888) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(0.57166878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088061995) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(-0.17266973) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2505635) q[2];
sx q[2];
rz(-2.6030412) q[2];
sx q[2];
rz(-2.0856759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.09517041) q[1];
sx q[1];
rz(-1.7032402) q[1];
sx q[1];
rz(-2.4080647) q[1];
rz(1.709135) q[3];
sx q[3];
rz(-1.1225015) q[3];
sx q[3];
rz(1.5773147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(2.4327915) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.4806992) q[0];
sx q[0];
rz(-1.9565214) q[0];
sx q[0];
rz(-1.6688523) q[0];
rz(-pi) q[1];
rz(-0.42787243) q[2];
sx q[2];
rz(-2.7395436) q[2];
sx q[2];
rz(-0.36487647) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0536249) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(1.4474523) q[1];
rz(-pi) q[2];
rz(-1.5933851) q[3];
sx q[3];
rz(-0.49593192) q[3];
sx q[3];
rz(-0.13726928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.583741) q[2];
rz(1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-2.8900237) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895296) q[0];
sx q[0];
rz(-1.7460199) q[0];
sx q[0];
rz(2.2625838) q[0];
x q[1];
rz(1.6576084) q[2];
sx q[2];
rz(-1.112339) q[2];
sx q[2];
rz(-2.364033) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.014616414) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(1.5986534) q[1];
rz(-0.5079481) q[3];
sx q[3];
rz(-0.93379279) q[3];
sx q[3];
rz(0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.183737) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(1.0305369) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(-0.39302557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0133936) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(-1.8830995) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9057301) q[2];
sx q[2];
rz(-2.3013407) q[2];
sx q[2];
rz(-2.0116531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1054354) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(-0.91962199) q[1];
rz(-pi) q[2];
rz(-0.7418467) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(1.1175931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3605911) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(-2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(3.135625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(1.9968541) q[0];
rz(-0.75255021) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(-0.37170751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55885) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(-1.4504257) q[1];
x q[2];
rz(-1.5238477) q[3];
sx q[3];
rz(-1.1991683) q[3];
sx q[3];
rz(0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-2.1933864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873916) q[0];
sx q[0];
rz(-1.750995) q[0];
sx q[0];
rz(-2.3733634) q[0];
x q[1];
rz(0.44342946) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(-2.7623917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4211385) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(0.2770263) q[1];
rz(-pi) q[2];
rz(2.9392397) q[3];
sx q[3];
rz(-0.36365299) q[3];
sx q[3];
rz(2.3285151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(0.24027696) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(0.25094029) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-2.7889263) q[1];
rz(-pi/2) q[2];
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
rz(-pi) q[1];
rz(-0.43899957) q[2];
sx q[2];
rz(-0.3294496) q[2];
sx q[2];
rz(-1.1169369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0295804) q[1];
sx q[1];
rz(-0.60830804) q[1];
sx q[1];
rz(0.40257247) q[1];
rz(-pi) q[2];
rz(-2.3064763) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(-1.9851828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76715604) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1922558) q[0];
sx q[0];
rz(-1.0431457) q[0];
sx q[0];
rz(-1.2883745) q[0];
rz(-pi) q[1];
rz(2.1328743) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(2.9269232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5156538) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(-2.291912) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0274721) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(0.041681899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(0.014952095) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(1.3199432) q[2];
sx q[2];
rz(-1.432042) q[2];
sx q[2];
rz(-1.5114573) q[2];
rz(-1.3068009) q[3];
sx q[3];
rz(-2.0859857) q[3];
sx q[3];
rz(1.2386238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
