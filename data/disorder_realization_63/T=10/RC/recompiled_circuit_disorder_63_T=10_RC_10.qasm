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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85815) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(-1.0648849) q[0];
rz(-1.4921161) q[2];
sx q[2];
rz(-0.71173758) q[2];
sx q[2];
rz(0.56456883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55878996) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(2.5488528) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4077971) q[3];
sx q[3];
rz(-1.194209) q[3];
sx q[3];
rz(0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-2.4332411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10193292) q[0];
sx q[0];
rz(-1.205737) q[0];
sx q[0];
rz(1.2714766) q[0];
rz(-2.1721241) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(3.0960992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0201455) q[1];
sx q[1];
rz(-2.2312806) q[1];
sx q[1];
rz(-2.5026908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6459341) q[3];
sx q[3];
rz(-1.111205) q[3];
sx q[3];
rz(-1.7840149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(-0.10270384) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(2.5699239) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088061995) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(-2.9689229) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2505635) q[2];
sx q[2];
rz(-2.6030412) q[2];
sx q[2];
rz(-1.0559168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5204822) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(-2.945167) q[1];
rz(-pi) q[2];
rz(1.709135) q[3];
sx q[3];
rz(-2.0190911) q[3];
sx q[3];
rz(-1.5773147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(2.4327915) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70808327) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(2.9611011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9468964) q[0];
sx q[0];
rz(-1.6616271) q[0];
sx q[0];
rz(0.38740654) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7137202) q[2];
sx q[2];
rz(-0.40204907) q[2];
sx q[2];
rz(0.36487647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0879678) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(1.6941403) q[1];
rz(1.0749712) q[3];
sx q[3];
rz(-1.5600481) q[3];
sx q[3];
rz(1.4533952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0522456) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(2.8900237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7895296) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(2.2625838) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4839843) q[2];
sx q[2];
rz(-2.0292536) q[2];
sx q[2];
rz(-2.364033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5800036) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(0.1954397) q[1];
rz(-pi) q[2];
rz(-0.5079481) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(2.7579443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(-2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(2.7485671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61009208) q[0];
sx q[0];
rz(-1.2628265) q[0];
sx q[0];
rz(0.17225762) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9057301) q[2];
sx q[2];
rz(-2.3013407) q[2];
sx q[2];
rz(1.1299396) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46950754) q[1];
sx q[1];
rz(-0.91963327) q[1];
sx q[1];
rz(0.0068411946) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7418467) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(2.0239995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(-0.97314107) q[2];
rz(0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(3.135625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0452506) q[0];
sx q[0];
rz(-1.8492286) q[0];
sx q[0];
rz(2.2521426) q[0];
rz(-pi) q[1];
rz(2.0231831) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(1.6391022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2318864) q[1];
sx q[1];
rz(-1.6344374) q[1];
sx q[1];
rz(1.015889) q[1];
x q[2];
rz(3.0217516) q[3];
sx q[3];
rz(-0.37444515) q[3];
sx q[3];
rz(-3.0344809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(0.19006426) q[2];
rz(0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3535239) q[0];
sx q[0];
rz(-2.3234962) q[0];
sx q[0];
rz(1.3226932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0886705) q[2];
sx q[2];
rz(-1.1793943) q[2];
sx q[2];
rz(-0.97757593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72045418) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(2.8645664) q[1];
rz(-1.4944581) q[3];
sx q[3];
rz(-1.9266955) q[3];
sx q[3];
rz(-2.5446041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(-2.9013157) q[2];
rz(-2.7622973) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3786479) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-2.7889263) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6610722) q[0];
sx q[0];
rz(-1.4561704) q[0];
sx q[0];
rz(1.6941487) q[0];
rz(-pi) q[1];
rz(0.43899957) q[2];
sx q[2];
rz(-2.8121431) q[2];
sx q[2];
rz(2.0246558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36665146) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(1.304438) q[1];
rz(-pi) q[2];
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
rz(2.3744366) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.8413683) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-0.93651071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47637981) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(0.54540821) q[0];
x q[1];
rz(-0.18239393) q[2];
sx q[2];
rz(-1.8511904) q[2];
sx q[2];
rz(0.80255752) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(0.84968062) q[1];
rz(2.1533222) q[3];
sx q[3];
rz(-1.3068145) q[3];
sx q[3];
rz(-1.9895944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-0.014952095) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(-1.5785718) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(-1.3199432) q[2];
sx q[2];
rz(-1.7095507) q[2];
sx q[2];
rz(1.6301353) q[2];
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
