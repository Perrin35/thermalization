OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(2.6240786) q[0];
rz(1.327688) q[1];
sx q[1];
rz(7.1751243) q[1];
sx q[1];
rz(8.8827477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9292277) q[0];
sx q[0];
rz(-0.80885086) q[0];
sx q[0];
rz(1.3332643) q[0];
rz(1.3165242) q[2];
sx q[2];
rz(-0.82077998) q[2];
sx q[2];
rz(-2.9645142) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.65329018) q[1];
sx q[1];
rz(-0.79217109) q[1];
sx q[1];
rz(-1.2912871) q[1];
rz(-pi) q[2];
rz(1.7446506) q[3];
sx q[3];
rz(-1.3411634) q[3];
sx q[3];
rz(-1.053912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1338256) q[2];
sx q[2];
rz(-1.0453036) q[2];
sx q[2];
rz(0.93057752) q[2];
rz(3.0104356) q[3];
sx q[3];
rz(-2.7635062) q[3];
sx q[3];
rz(-1.4906496) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5432878) q[0];
sx q[0];
rz(-2.2853993) q[0];
sx q[0];
rz(1.244586) q[0];
rz(0.92567956) q[1];
sx q[1];
rz(-1.8857748) q[1];
sx q[1];
rz(-1.6474887) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5410575) q[0];
sx q[0];
rz(-0.78406683) q[0];
sx q[0];
rz(-2.0406141) q[0];
rz(0.67333198) q[2];
sx q[2];
rz(-1.8367366) q[2];
sx q[2];
rz(0.8348726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0665613) q[1];
sx q[1];
rz(-1.9157456) q[1];
sx q[1];
rz(-2.6607772) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5794287) q[3];
sx q[3];
rz(-1.8749471) q[3];
sx q[3];
rz(2.9184273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1458448) q[2];
sx q[2];
rz(-2.6971942) q[2];
sx q[2];
rz(0.90163976) q[2];
rz(-1.7835279) q[3];
sx q[3];
rz(-1.3172904) q[3];
sx q[3];
rz(-2.5942514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77483344) q[0];
sx q[0];
rz(-1.9623373) q[0];
sx q[0];
rz(1.1460079) q[0];
rz(2.215812) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(2.944223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1903052) q[0];
sx q[0];
rz(-1.7078551) q[0];
sx q[0];
rz(-0.97676116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8497265) q[2];
sx q[2];
rz(-1.2780683) q[2];
sx q[2];
rz(-1.9809993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2132385) q[1];
sx q[1];
rz(-2.9103807) q[1];
sx q[1];
rz(-1.3998881) q[1];
x q[2];
rz(1.7719384) q[3];
sx q[3];
rz(-2.1052868) q[3];
sx q[3];
rz(0.25132685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7388514) q[2];
sx q[2];
rz(-2.1612942) q[2];
sx q[2];
rz(2.951238) q[2];
rz(-0.061577408) q[3];
sx q[3];
rz(-2.172251) q[3];
sx q[3];
rz(-2.0891345) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15678081) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(-0.62830997) q[0];
rz(1.620694) q[1];
sx q[1];
rz(-1.2077121) q[1];
sx q[1];
rz(-0.77888387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464061) q[0];
sx q[0];
rz(-2.4281326) q[0];
sx q[0];
rz(0.97795217) q[0];
rz(2.4905048) q[2];
sx q[2];
rz(-0.44241239) q[2];
sx q[2];
rz(-0.96786066) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0592093) q[1];
sx q[1];
rz(-1.6184855) q[1];
sx q[1];
rz(-1.5842745) q[1];
rz(-pi) q[2];
rz(2.3301937) q[3];
sx q[3];
rz(-2.2795942) q[3];
sx q[3];
rz(2.3611695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4817619) q[2];
sx q[2];
rz(-2.063282) q[2];
sx q[2];
rz(3.1415494) q[2];
rz(1.0846042) q[3];
sx q[3];
rz(-1.1559887) q[3];
sx q[3];
rz(-2.675975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67149177) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(1.8430365) q[0];
rz(-0.57688722) q[1];
sx q[1];
rz(-1.2737609) q[1];
sx q[1];
rz(-2.6947122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80432207) q[0];
sx q[0];
rz(-1.7511586) q[0];
sx q[0];
rz(0.65482636) q[0];
rz(-pi) q[1];
rz(1.8700897) q[2];
sx q[2];
rz(-1.6422049) q[2];
sx q[2];
rz(1.4324607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3342377) q[1];
sx q[1];
rz(-0.45354834) q[1];
sx q[1];
rz(-1.0940882) q[1];
x q[2];
rz(2.9561012) q[3];
sx q[3];
rz(-1.5122454) q[3];
sx q[3];
rz(1.0492988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7070861) q[2];
sx q[2];
rz(-0.27927566) q[2];
sx q[2];
rz(-2.9217829) q[2];
rz(-1.48742) q[3];
sx q[3];
rz(-1.6916964) q[3];
sx q[3];
rz(-2.4421104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18911067) q[0];
sx q[0];
rz(-2.7087152) q[0];
sx q[0];
rz(0.89749807) q[0];
rz(1.7706361) q[1];
sx q[1];
rz(-1.0400583) q[1];
sx q[1];
rz(-2.2460361) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093256) q[0];
sx q[0];
rz(-2.4054962) q[0];
sx q[0];
rz(-0.21047896) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3562548) q[2];
sx q[2];
rz(-0.78162748) q[2];
sx q[2];
rz(1.8434032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59437737) q[1];
sx q[1];
rz(-1.2378232) q[1];
sx q[1];
rz(0.39526387) q[1];
rz(-pi) q[2];
rz(1.0895133) q[3];
sx q[3];
rz(-1.1466422) q[3];
sx q[3];
rz(-1.4312351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36876496) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(0.15711288) q[2];
rz(-0.67240063) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(-0.11252832) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2812578) q[0];
sx q[0];
rz(-2.4761138) q[0];
sx q[0];
rz(1.459664) q[0];
rz(-2.7809987) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(1.8416539) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324289) q[0];
sx q[0];
rz(-0.43540186) q[0];
sx q[0];
rz(2.8485427) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5333161) q[2];
sx q[2];
rz(-2.424834) q[2];
sx q[2];
rz(-2.3065717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9323401) q[1];
sx q[1];
rz(-2.8289218) q[1];
sx q[1];
rz(-2.8220842) q[1];
rz(-0.61494334) q[3];
sx q[3];
rz(-1.9087026) q[3];
sx q[3];
rz(0.29454052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.931687) q[2];
sx q[2];
rz(-2.1114712) q[2];
sx q[2];
rz(2.9295861) q[2];
rz(2.0906406) q[3];
sx q[3];
rz(-1.3764328) q[3];
sx q[3];
rz(3.0534548) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6050922) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(0.17298803) q[0];
rz(2.0269003) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(-0.10985049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3152471) q[0];
sx q[0];
rz(-1.7638806) q[0];
sx q[0];
rz(2.5772334) q[0];
rz(-pi) q[1];
rz(-2.8809261) q[2];
sx q[2];
rz(-2.3015966) q[2];
sx q[2];
rz(2.276401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.122095) q[1];
sx q[1];
rz(-1.3155342) q[1];
sx q[1];
rz(2.9649516) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55446378) q[3];
sx q[3];
rz(-1.4687456) q[3];
sx q[3];
rz(-0.79342519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7383808) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(-0.66486764) q[2];
rz(0.46857771) q[3];
sx q[3];
rz(-1.5237619) q[3];
sx q[3];
rz(2.328228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68798962) q[0];
sx q[0];
rz(-0.60992321) q[0];
sx q[0];
rz(2.4000121) q[0];
rz(1.954151) q[1];
sx q[1];
rz(-1.6148184) q[1];
sx q[1];
rz(2.5724519) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4837836) q[0];
sx q[0];
rz(-1.8632194) q[0];
sx q[0];
rz(2.800709) q[0];
rz(2.9958565) q[2];
sx q[2];
rz(-0.45782858) q[2];
sx q[2];
rz(-0.65953461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7267159) q[1];
sx q[1];
rz(-2.0316213) q[1];
sx q[1];
rz(-1.0970835) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3027116) q[3];
sx q[3];
rz(-1.0195707) q[3];
sx q[3];
rz(-1.6987926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6471214) q[2];
sx q[2];
rz(-0.98860604) q[2];
sx q[2];
rz(0.76623255) q[2];
rz(-1.9689485) q[3];
sx q[3];
rz(-1.5394883) q[3];
sx q[3];
rz(-0.22181454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2072993) q[0];
sx q[0];
rz(-2.793029) q[0];
sx q[0];
rz(-0.19185129) q[0];
rz(0.81709298) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(2.4339035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548508) q[0];
sx q[0];
rz(-1.8157209) q[0];
sx q[0];
rz(-1.7502511) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0470263) q[2];
sx q[2];
rz(-2.4942538) q[2];
sx q[2];
rz(0.29905427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5194044) q[1];
sx q[1];
rz(-0.71336245) q[1];
sx q[1];
rz(-0.98825561) q[1];
rz(-pi) q[2];
rz(-0.58587425) q[3];
sx q[3];
rz(-2.088024) q[3];
sx q[3];
rz(-0.45211238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2231458) q[2];
sx q[2];
rz(-2.196849) q[2];
sx q[2];
rz(2.4165912) q[2];
rz(-0.21507344) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(0.71896368) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216777) q[0];
sx q[0];
rz(-0.81659962) q[0];
sx q[0];
rz(0.30246977) q[0];
rz(1.1014145) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(2.106059) q[2];
sx q[2];
rz(-1.4695833) q[2];
sx q[2];
rz(1.5493456) q[2];
rz(1.5982895) q[3];
sx q[3];
rz(-2.0381654) q[3];
sx q[3];
rz(0.5034133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
