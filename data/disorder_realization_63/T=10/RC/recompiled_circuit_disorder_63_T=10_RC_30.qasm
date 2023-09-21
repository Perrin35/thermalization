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
rz(-2.0768291) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(-0.53607166) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85815) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(-2.0767077) q[0];
rz(-pi) q[1];
rz(-0.86059086) q[2];
sx q[2];
rz(-1.5194367) q[2];
sx q[2];
rz(-1.0658588) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2546665) q[1];
sx q[1];
rz(-0.98985043) q[1];
sx q[1];
rz(-1.7960153) q[1];
rz(-pi) q[2];
rz(2.4077971) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(-0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(2.1146963) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0396597) q[0];
sx q[0];
rz(-1.205737) q[0];
sx q[0];
rz(-1.870116) q[0];
x q[1];
rz(1.3588261) q[2];
sx q[2];
rz(-2.5296629) q[2];
sx q[2];
rz(1.7906534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.12304141) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(-2.3398188) q[1];
x q[2];
rz(-2.6459341) q[3];
sx q[3];
rz(-2.0303876) q[3];
sx q[3];
rz(1.3575777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(3.0388888) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(0.57166878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11274352) q[0];
sx q[0];
rz(-1.0380121) q[0];
sx q[0];
rz(-1.4677731) q[0];
rz(-pi) q[1];
rz(-2.7823206) q[2];
sx q[2];
rz(-1.9810988) q[2];
sx q[2];
rz(-1.8112195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5204822) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(0.1964257) q[1];
rz(-pi) q[2];
rz(-1.4324576) q[3];
sx q[3];
rz(-2.0190911) q[3];
sx q[3];
rz(1.5642779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(2.4327915) q[2];
rz(-0.30250868) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(0.18049151) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4806992) q[0];
sx q[0];
rz(-1.1850712) q[0];
sx q[0];
rz(1.6688523) q[0];
rz(-2.7724491) q[2];
sx q[2];
rz(-1.4077079) q[2];
sx q[2];
rz(0.8085608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0536249) q[1];
sx q[1];
rz(-1.2521724) q[1];
sx q[1];
rz(1.6941403) q[1];
x q[2];
rz(-0.012219592) q[3];
sx q[3];
rz(-2.0665902) q[3];
sx q[3];
rz(3.030005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(1.5578516) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(0.93262514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(2.8900237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7895296) q[0];
sx q[0];
rz(-1.7460199) q[0];
sx q[0];
rz(-2.2625838) q[0];
rz(-pi) q[1];
rz(-0.45995633) q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1269762) q[1];
sx q[1];
rz(-1.3754305) q[1];
sx q[1];
rz(-1.5986534) q[1];
rz(0.5079481) q[3];
sx q[3];
rz(-2.2077999) q[3];
sx q[3];
rz(0.38364832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(0.63123909) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(2.2914698) q[0];
rz(-1.1823581) q[1];
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
rz(-2.5315006) q[0];
sx q[0];
rz(-1.2628265) q[0];
sx q[0];
rz(-2.969335) q[0];
rz(-pi) q[1];
rz(2.3153147) q[2];
sx q[2];
rz(-1.7457361) q[2];
sx q[2];
rz(-2.8597521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0361573) q[1];
sx q[1];
rz(-1.5762377) q[1];
sx q[1];
rz(-0.91962199) q[1];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(0.24916838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-1.1258619) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.5234891) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8864266) q[0];
sx q[0];
rz(-0.92029858) q[0];
sx q[0];
rz(2.7889473) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47875328) q[2];
sx q[2];
rz(-0.81214777) q[2];
sx q[2];
rz(2.2854545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5827427) q[1];
sx q[1];
rz(-0.55816459) q[1];
sx q[1];
rz(1.6911669) q[1];
rz(-pi) q[2];
rz(2.7695914) q[3];
sx q[3];
rz(-1.614538) q[3];
sx q[3];
rz(1.352076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(-2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-0.53623143) q[0];
rz(-0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-2.1933864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1873916) q[0];
sx q[0];
rz(-1.750995) q[0];
sx q[0];
rz(2.3733634) q[0];
rz(-pi) q[1];
rz(-1.0529222) q[2];
sx q[2];
rz(-1.9621984) q[2];
sx q[2];
rz(0.97757593) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4211385) q[1];
sx q[1];
rz(-2.9109143) q[1];
sx q[1];
rz(2.8645664) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20235297) q[3];
sx q[3];
rz(-2.7779397) q[3];
sx q[3];
rz(-0.81307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(0.24027696) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(2.8885686) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1044554) q[0];
sx q[0];
rz(-1.6933352) q[0];
sx q[0];
rz(-0.11549581) q[0];
rz(-0.30013957) q[2];
sx q[2];
rz(-1.7087414) q[2];
sx q[2];
rz(-0.87196841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11201227) q[1];
sx q[1];
rz(-0.60830804) q[1];
sx q[1];
rz(2.7390202) q[1];
rz(-2.3064763) q[3];
sx q[3];
rz(-1.1856688) q[3];
sx q[3];
rz(-1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(-1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(-1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3764573) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.3002243) q[0];
rz(0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-2.2050819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6694326) q[0];
sx q[0];
rz(-2.5494808) q[0];
sx q[0];
rz(0.44606146) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8556701) q[2];
sx q[2];
rz(-1.7459918) q[2];
sx q[2];
rz(0.81923649) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.7682372) q[1];
sx q[1];
rz(0.84968062) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98827045) q[3];
sx q[3];
rz(-1.8347782) q[3];
sx q[3];
rz(1.9895944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(2.0832534) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(-1.8347918) q[3];
sx q[3];
rz(-1.055607) q[3];
sx q[3];
rz(-1.9029688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
