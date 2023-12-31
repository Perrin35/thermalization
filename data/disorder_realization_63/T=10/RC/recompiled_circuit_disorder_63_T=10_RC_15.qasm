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
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4608085) q[0];
sx q[0];
rz(-0.55235282) q[0];
sx q[0];
rz(2.0244563) q[0];
rz(0.86059086) q[2];
sx q[2];
rz(-1.6221559) q[2];
sx q[2];
rz(2.0757338) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8592035) q[1];
sx q[1];
rz(-0.61835641) q[1];
sx q[1];
rz(0.32787852) q[1];
rz(-pi) q[2];
rz(1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(2.1171452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-0.70835152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0396597) q[0];
sx q[0];
rz(-1.205737) q[0];
sx q[0];
rz(-1.870116) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3588261) q[2];
sx q[2];
rz(-0.6119298) q[2];
sx q[2];
rz(1.7906534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9006151) q[1];
sx q[1];
rz(-2.2574189) q[1];
sx q[1];
rz(2.2254506) q[1];
rz(-pi) q[2];
rz(-0.49565855) q[3];
sx q[3];
rz(-1.111205) q[3];
sx q[3];
rz(1.3575777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(3.0388888) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(-1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(-1.6614871) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(-2.5699239) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11274352) q[0];
sx q[0];
rz(-2.1035806) q[0];
sx q[0];
rz(-1.6738196) q[0];
rz(0.35927202) q[2];
sx q[2];
rz(-1.9810988) q[2];
sx q[2];
rz(1.3303732) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.09517041) q[1];
sx q[1];
rz(-1.4383525) q[1];
sx q[1];
rz(-0.73352791) q[1];
rz(0.45205558) q[3];
sx q[3];
rz(-1.6953903) q[3];
sx q[3];
rz(-0.053754036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4335094) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-0.79950142) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(-2.9611011) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1946963) q[0];
sx q[0];
rz(-1.6616271) q[0];
sx q[0];
rz(-0.38740654) q[0];
rz(-pi) q[1];
rz(1.7454342) q[2];
sx q[2];
rz(-1.2067814) q[2];
sx q[2];
rz(-2.3166235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4304639) q[1];
sx q[1];
rz(-2.8006878) q[1];
sx q[1];
rz(2.7845963) q[1];
rz(2.0666215) q[3];
sx q[3];
rz(-1.5600481) q[3];
sx q[3];
rz(1.6881975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(1.7662988) q[3];
sx q[3];
rz(-1.8245274) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-0.25156897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7154327) q[0];
sx q[0];
rz(-0.71007219) q[0];
sx q[0];
rz(1.8415113) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45995633) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(2.3098582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0136452) q[1];
sx q[1];
rz(-0.19731678) q[1];
sx q[1];
rz(3.0017588) q[1];
rz(0.5079481) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(-2.7579443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(2.1110558) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-2.5103536) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(-0.39302557) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61009208) q[0];
sx q[0];
rz(-1.2628265) q[0];
sx q[0];
rz(-0.17225762) q[0];
x q[1];
rz(-0.82627798) q[2];
sx q[2];
rz(-1.3958566) q[2];
sx q[2];
rz(2.8597521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1054354) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(-0.91962199) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.844448) q[3];
sx q[3];
rz(-1.2840464) q[3];
sx q[3];
rz(1.2424038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3605911) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(-2.1684516) q[2];
rz(0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.4457552) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8864266) q[0];
sx q[0];
rz(-2.2212941) q[0];
sx q[0];
rz(-0.35264539) q[0];
x q[1];
rz(-2.3890424) q[2];
sx q[2];
rz(-1.9117022) q[2];
sx q[2];
rz(2.7698851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.55885) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(-1.4504257) q[1];
rz(-pi) q[2];
rz(1.617745) q[3];
sx q[3];
rz(-1.1991683) q[3];
sx q[3];
rz(0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(-2.7251785) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(2.4068508) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(-0.94820625) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873916) q[0];
sx q[0];
rz(-1.750995) q[0];
sx q[0];
rz(0.76822922) q[0];
rz(-pi) q[1];
rz(2.0886705) q[2];
sx q[2];
rz(-1.1793943) q[2];
sx q[2];
rz(-0.97757593) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1369689) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(1.5066513) q[1];
rz(-0.20235297) q[3];
sx q[3];
rz(-0.36365299) q[3];
sx q[3];
rz(2.3285151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(-2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-2.7889263) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371373) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(3.0260968) q[0];
x q[1];
rz(2.7025931) q[2];
sx q[2];
rz(-0.3294496) q[2];
sx q[2];
rz(2.0246558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7749412) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(1.304438) q[1];
x q[2];
rz(-0.8351164) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(-1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.3002243) q[0];
rz(0.70872778) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6694326) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(-2.6955312) q[0];
rz(-1.2859225) q[2];
sx q[2];
rz(-1.3956009) q[2];
sx q[2];
rz(0.81923649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.3733555) q[1];
sx q[1];
rz(2.291912) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1141206) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(-3.0999108) q[3];
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
rz(2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(-1.5630209) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(-1.8216495) q[2];
sx q[2];
rz(-1.432042) q[2];
sx q[2];
rz(-1.5114573) q[2];
rz(-2.6111504) q[3];
sx q[3];
rz(-1.3417288) q[3];
sx q[3];
rz(-0.19977278) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
