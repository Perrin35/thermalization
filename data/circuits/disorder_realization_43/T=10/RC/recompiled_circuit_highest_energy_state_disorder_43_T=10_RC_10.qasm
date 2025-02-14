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
rz(0.69005203) q[0];
sx q[0];
rz(8.5741841) q[0];
sx q[0];
rz(9.2246715) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(5.7601647) q[1];
sx q[1];
rz(10.756607) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1855186) q[0];
sx q[0];
rz(-1.3393434) q[0];
sx q[0];
rz(-1.1197907) q[0];
rz(-2.5426504) q[2];
sx q[2];
rz(-1.9859196) q[2];
sx q[2];
rz(2.6225363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15439776) q[1];
sx q[1];
rz(-1.618865) q[1];
sx q[1];
rz(0.67832077) q[1];
rz(-pi) q[2];
rz(-3.0588714) q[3];
sx q[3];
rz(-1.0961856) q[3];
sx q[3];
rz(0.079291346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.3145831) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(-2.7956853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697934) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1163569) q[0];
sx q[0];
rz(-1.5486002) q[0];
sx q[0];
rz(1.8216013) q[0];
rz(0.34423265) q[2];
sx q[2];
rz(-1.5308799) q[2];
sx q[2];
rz(1.0913864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.98709244) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(-2.2461056) q[1];
rz(-1.3608906) q[3];
sx q[3];
rz(-0.91448254) q[3];
sx q[3];
rz(2.8950952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0756125) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(-1.570545) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(-1.7861231) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9539255) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-0.060610108) q[0];
rz(-0.68471471) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(-0.33904591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6936172) q[0];
sx q[0];
rz(-2.1537499) q[0];
sx q[0];
rz(-1.5956123) q[0];
rz(-pi) q[1];
rz(2.9525063) q[2];
sx q[2];
rz(-1.6918618) q[2];
sx q[2];
rz(2.6775286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10968929) q[1];
sx q[1];
rz(-1.8432143) q[1];
sx q[1];
rz(1.8930356) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5590844) q[3];
sx q[3];
rz(-2.3475966) q[3];
sx q[3];
rz(0.3670961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(-0.47046146) q[2];
rz(-0.86236924) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(1.580015) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3622793) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-2.733574) q[0];
rz(-1.3511924) q[1];
sx q[1];
rz(-1.8441169) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6141339) q[0];
sx q[0];
rz(-1.8728932) q[0];
sx q[0];
rz(0.52665773) q[0];
x q[1];
rz(-0.91561134) q[2];
sx q[2];
rz(-1.5024937) q[2];
sx q[2];
rz(-3.0223522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3258354) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.6440637) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77797555) q[3];
sx q[3];
rz(-0.56853349) q[3];
sx q[3];
rz(-1.3028627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-0.72285405) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(0.55996672) q[0];
rz(-2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-0.29108873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4093035) q[0];
sx q[0];
rz(-1.8390053) q[0];
sx q[0];
rz(-1.3951673) q[0];
rz(0.79925691) q[2];
sx q[2];
rz(-1.6148668) q[2];
sx q[2];
rz(1.26233) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4704125) q[1];
sx q[1];
rz(-0.7268097) q[1];
sx q[1];
rz(1.2695168) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.970333) q[3];
sx q[3];
rz(-0.75920877) q[3];
sx q[3];
rz(1.9943135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9186972) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(1.1693303) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(-1.0317624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87001291) q[0];
sx q[0];
rz(-1.7450512) q[0];
sx q[0];
rz(-0.016079024) q[0];
rz(-pi) q[1];
rz(-2.6844031) q[2];
sx q[2];
rz(-0.44538272) q[2];
sx q[2];
rz(-3.074844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9578625) q[1];
sx q[1];
rz(-2.0477844) q[1];
sx q[1];
rz(-1.3579923) q[1];
rz(-2.2797645) q[3];
sx q[3];
rz(-2.3579881) q[3];
sx q[3];
rz(-2.6605822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80456698) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(2.2931781) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(3.1019822) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1161716) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(2.1436932) q[0];
rz(2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.3708699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8664198) q[0];
sx q[0];
rz(-0.93979561) q[0];
sx q[0];
rz(-1.5874528) q[0];
x q[1];
rz(-0.15372686) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(1.597126) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5120041) q[1];
sx q[1];
rz(-1.8114212) q[1];
sx q[1];
rz(-2.5179067) q[1];
rz(2.3418535) q[3];
sx q[3];
rz(-0.99501401) q[3];
sx q[3];
rz(-2.6588001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-0.029646309) q[2];
rz(-0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(2.7992547) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(1.7436183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2919162) q[0];
sx q[0];
rz(-1.9789985) q[0];
sx q[0];
rz(1.085039) q[0];
rz(-0.2335642) q[2];
sx q[2];
rz(-1.8478881) q[2];
sx q[2];
rz(-1.8943) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82247558) q[1];
sx q[1];
rz(-2.2386595) q[1];
sx q[1];
rz(-0.0112337) q[1];
rz(-pi) q[2];
rz(-0.21381883) q[3];
sx q[3];
rz(-2.345746) q[3];
sx q[3];
rz(-2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(1.3207377) q[2];
rz(1.0197506) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19972292) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(-0.95440188) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53756489) q[0];
sx q[0];
rz(-1.5286235) q[0];
sx q[0];
rz(2.2794363) q[0];
rz(2.886359) q[2];
sx q[2];
rz(-1.8430147) q[2];
sx q[2];
rz(1.24025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22218765) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(-0.081187831) q[1];
rz(-1.3865269) q[3];
sx q[3];
rz(-0.83067719) q[3];
sx q[3];
rz(2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(1.3433749) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0992391) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(1.4105256) q[0];
rz(0.22077416) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-0.9332307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8659023) q[0];
sx q[0];
rz(-2.3161016) q[0];
sx q[0];
rz(-1.5522214) q[0];
rz(-pi) q[1];
rz(-1.9311229) q[2];
sx q[2];
rz(-2.1585352) q[2];
sx q[2];
rz(2.7992942) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2360742) q[1];
sx q[1];
rz(-0.80414861) q[1];
sx q[1];
rz(2.9574636) q[1];
rz(-pi) q[2];
rz(1.7228863) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(1.9877951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(-0.7190052) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.63856335) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-2.6999264) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(0.031671192) q[2];
sx q[2];
rz(-1.9513672) q[2];
sx q[2];
rz(0.58917108) q[2];
rz(-1.6735703) q[3];
sx q[3];
rz(-1.0036052) q[3];
sx q[3];
rz(-2.6617692) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
