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
rz(0.063989446) q[0];
sx q[0];
rz(4.0539157) q[0];
sx q[0];
rz(11.231448) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27150422) q[0];
sx q[0];
rz(-3.0850924) q[0];
sx q[0];
rz(1.9046049) q[0];
rz(0.71662997) q[2];
sx q[2];
rz(-0.83558768) q[2];
sx q[2];
rz(2.731833) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0483889) q[1];
sx q[1];
rz(-1.2902765) q[1];
sx q[1];
rz(2.2040221) q[1];
x q[2];
rz(-1.3282449) q[3];
sx q[3];
rz(-1.537916) q[3];
sx q[3];
rz(-1.3803079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7923183) q[2];
sx q[2];
rz(-0.26733843) q[2];
sx q[2];
rz(-0.054923687) q[2];
rz(-1.6867636) q[3];
sx q[3];
rz(-2.7456386) q[3];
sx q[3];
rz(-1.5168064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216264) q[0];
sx q[0];
rz(-1.8333789) q[0];
sx q[0];
rz(0.24800214) q[0];
rz(1.8964881) q[1];
sx q[1];
rz(-0.73687941) q[1];
sx q[1];
rz(-1.2287593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7354483) q[0];
sx q[0];
rz(-1.3856674) q[0];
sx q[0];
rz(-1.6044751) q[0];
x q[1];
rz(2.3360148) q[2];
sx q[2];
rz(-0.68289872) q[2];
sx q[2];
rz(0.77308347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1540888) q[1];
sx q[1];
rz(-2.262907) q[1];
sx q[1];
rz(1.8471884) q[1];
rz(-pi) q[2];
rz(1.6307488) q[3];
sx q[3];
rz(-1.9904117) q[3];
sx q[3];
rz(-2.7762716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9241141) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(0.58383101) q[2];
rz(2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(1.0113641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15447021) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(-1.467147) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.4091622) q[1];
sx q[1];
rz(1.0999701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407134) q[0];
sx q[0];
rz(-1.3726808) q[0];
sx q[0];
rz(1.4639356) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8520865) q[2];
sx q[2];
rz(-1.2443674) q[2];
sx q[2];
rz(-3.141186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7138805) q[1];
sx q[1];
rz(-1.1888479) q[1];
sx q[1];
rz(-1.9314585) q[1];
rz(-pi) q[2];
rz(-0.93020029) q[3];
sx q[3];
rz(-0.67852321) q[3];
sx q[3];
rz(0.92032209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.786342) q[2];
sx q[2];
rz(-0.92266551) q[2];
sx q[2];
rz(1.5477017) q[2];
rz(-1.330438) q[3];
sx q[3];
rz(-1.7682313) q[3];
sx q[3];
rz(1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.380577) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(-0.15783489) q[0];
rz(-2.8997968) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(0.48113021) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7173955) q[0];
sx q[0];
rz(-1.6592053) q[0];
sx q[0];
rz(-2.428695) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35901661) q[2];
sx q[2];
rz(-1.8098157) q[2];
sx q[2];
rz(1.1306819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2896636) q[1];
sx q[1];
rz(-0.73411513) q[1];
sx q[1];
rz(1.0490225) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1948482) q[3];
sx q[3];
rz(-2.7928154) q[3];
sx q[3];
rz(-2.6117976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(2.8672186) q[2];
rz(1.4659878) q[3];
sx q[3];
rz(-1.8612739) q[3];
sx q[3];
rz(0.93332851) q[3];
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
rz(pi/2) q[0];
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
rz(0.45193732) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(-2.080132) q[0];
rz(-2.6929216) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(1.7015069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0827918) q[0];
sx q[0];
rz(-0.77547204) q[0];
sx q[0];
rz(-2.916921) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8605804) q[2];
sx q[2];
rz(-2.7216879) q[2];
sx q[2];
rz(-0.29422255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0651814) q[1];
sx q[1];
rz(-1.5440327) q[1];
sx q[1];
rz(-1.3127432) q[1];
rz(-pi) q[2];
rz(1.7251882) q[3];
sx q[3];
rz(-2.5254211) q[3];
sx q[3];
rz(2.7026724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(-1.9484005) q[2];
rz(2.7423972) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(-3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4611918) q[0];
sx q[0];
rz(-1.8541279) q[0];
sx q[0];
rz(0.68049085) q[0];
rz(0.83241278) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(0.15636538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0427317) q[0];
sx q[0];
rz(-1.2953238) q[0];
sx q[0];
rz(-1.5163438) q[0];
x q[1];
rz(0.83145468) q[2];
sx q[2];
rz(-2.105684) q[2];
sx q[2];
rz(1.2818467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.029441) q[1];
sx q[1];
rz(-1.8990771) q[1];
sx q[1];
rz(-2.8480457) q[1];
rz(-pi) q[2];
rz(1.0632319) q[3];
sx q[3];
rz(-0.2476736) q[3];
sx q[3];
rz(1.4192691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6303595) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(2.2516001) q[2];
rz(-0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156859) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(0.61221468) q[0];
rz(-1.7828434) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(-2.8403958) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4801203) q[0];
sx q[0];
rz(-1.8422707) q[0];
sx q[0];
rz(-3.0508811) q[0];
rz(0.37294228) q[2];
sx q[2];
rz(-0.94764793) q[2];
sx q[2];
rz(-0.22280338) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14872486) q[1];
sx q[1];
rz(-2.5588256) q[1];
sx q[1];
rz(0.73281835) q[1];
rz(-pi) q[2];
rz(-2.4047818) q[3];
sx q[3];
rz(-2.7932146) q[3];
sx q[3];
rz(0.37955561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23413868) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(-2.786484) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(-2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2316786) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(0.32383305) q[0];
rz(-2.553885) q[1];
sx q[1];
rz(-2.4617742) q[1];
sx q[1];
rz(-0.30276611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2229109) q[0];
sx q[0];
rz(-2.7245173) q[0];
sx q[0];
rz(0.66584754) q[0];
rz(2.9423334) q[2];
sx q[2];
rz(-2.8997921) q[2];
sx q[2];
rz(1.3346145) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0787133) q[1];
sx q[1];
rz(-2.4643366) q[1];
sx q[1];
rz(-0.96911624) q[1];
x q[2];
rz(-1.15074) q[3];
sx q[3];
rz(-1.3358634) q[3];
sx q[3];
rz(-2.3714424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6451463) q[2];
sx q[2];
rz(-1.030693) q[2];
sx q[2];
rz(2.5808064) q[2];
rz(2.136611) q[3];
sx q[3];
rz(-1.8533665) q[3];
sx q[3];
rz(-2.0463478) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305369) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(0.74991599) q[0];
rz(-0.38052446) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(0.97533018) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0935287) q[0];
sx q[0];
rz(-0.87665999) q[0];
sx q[0];
rz(2.3175236) q[0];
rz(-pi) q[1];
rz(-1.777095) q[2];
sx q[2];
rz(-2.8562284) q[2];
sx q[2];
rz(2.9138164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6187894) q[1];
sx q[1];
rz(-2.0674043) q[1];
sx q[1];
rz(1.9267531) q[1];
rz(-pi) q[2];
rz(-1.2694025) q[3];
sx q[3];
rz(-1.3580672) q[3];
sx q[3];
rz(-1.0476867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19899496) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(-1.0814166) q[2];
rz(-0.069247581) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13114318) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(-2.8358054) q[0];
rz(2.8336613) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(1.1526795) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7869513) q[0];
sx q[0];
rz(-0.85022038) q[0];
sx q[0];
rz(1.2194751) q[0];
rz(-pi) q[1];
rz(-2.9342317) q[2];
sx q[2];
rz(-0.93610686) q[2];
sx q[2];
rz(0.35516741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88886966) q[1];
sx q[1];
rz(-2.6111341) q[1];
sx q[1];
rz(2.2525851) q[1];
rz(1.6341428) q[3];
sx q[3];
rz(-1.5858486) q[3];
sx q[3];
rz(-1.9346332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(-1.6884241) q[2];
rz(0.48062634) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65344812) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(1.0427955) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(1.8389134) q[2];
sx q[2];
rz(-2.0901879) q[2];
sx q[2];
rz(-1.2218634) q[2];
rz(-3.0328209) q[3];
sx q[3];
rz(-0.36009195) q[3];
sx q[3];
rz(2.910955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
