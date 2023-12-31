OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(2.8013464) q[1];
sx q[1];
rz(10.624788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8631247) q[0];
sx q[0];
rz(-1.5651363) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8067402) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(-1.877117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.66558054) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(-2.0631454) q[1];
rz(0.37813152) q[3];
sx q[3];
rz(-1.5298801) q[3];
sx q[3];
rz(-0.10842987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(-3.0572609) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.3551691) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38793135) q[0];
sx q[0];
rz(-0.61457115) q[0];
sx q[0];
rz(-1.4580926) q[0];
rz(1.6927035) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(-0.09300692) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.720571) q[1];
sx q[1];
rz(-2.6877626) q[1];
sx q[1];
rz(1.717091) q[1];
rz(-1.7543206) q[3];
sx q[3];
rz(-1.6471383) q[3];
sx q[3];
rz(-0.21218382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(2.0764988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051620313) q[0];
sx q[0];
rz(-1.9301156) q[0];
sx q[0];
rz(3.1348455) q[0];
rz(-pi) q[1];
rz(-0.87419072) q[2];
sx q[2];
rz(-1.0993996) q[2];
sx q[2];
rz(-3.0490321) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1467495) q[1];
sx q[1];
rz(-0.61770505) q[1];
sx q[1];
rz(1.5084933) q[1];
rz(2.8912656) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(-1.7771429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(0.24969077) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(3.1304741) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269037) q[0];
sx q[0];
rz(-2.5137797) q[0];
sx q[0];
rz(-2.3054302) q[0];
rz(-pi) q[1];
rz(2.3035994) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(2.879564) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.30918446) q[1];
sx q[1];
rz(-2.5384181) q[1];
sx q[1];
rz(-0.96997728) q[1];
rz(-pi) q[2];
rz(0.53400455) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(-0.66343534) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(1.8146851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793314) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(-0.22236951) q[0];
rz(-pi) q[1];
rz(-2.1660216) q[2];
sx q[2];
rz(-2.274548) q[2];
sx q[2];
rz(0.032035839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9304282) q[1];
sx q[1];
rz(-1.2439562) q[1];
sx q[1];
rz(-1.8153166) q[1];
rz(1.9782449) q[3];
sx q[3];
rz(-1.3031928) q[3];
sx q[3];
rz(2.101055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1422687) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.8959321) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(-2.9340414) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(1.0043762) q[0];
x q[1];
rz(-2.0954779) q[2];
sx q[2];
rz(-2.4002053) q[2];
sx q[2];
rz(1.3337097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11568497) q[1];
sx q[1];
rz(-0.58870643) q[1];
sx q[1];
rz(-3.1097417) q[1];
rz(-pi) q[2];
rz(-2.9168105) q[3];
sx q[3];
rz(-2.4596679) q[3];
sx q[3];
rz(3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(-2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-2.9102944) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092098504) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(-0.90256079) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28568761) q[2];
sx q[2];
rz(-2.7966768) q[2];
sx q[2];
rz(-3.0722741) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.454969) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(1.7186233) q[1];
rz(0.087248487) q[3];
sx q[3];
rz(-2.1656519) q[3];
sx q[3];
rz(-2.0283386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-2.7451519) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.5213535) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0884468) q[0];
sx q[0];
rz(-1.3741115) q[0];
sx q[0];
rz(-1.501207) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0422158) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(-2.9277756) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8926881) q[1];
sx q[1];
rz(-2.3612594) q[1];
sx q[1];
rz(-0.08463879) q[1];
rz(-2.2500854) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(2.7822568) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(0.27063453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04979241) q[0];
sx q[0];
rz(-1.6066178) q[0];
sx q[0];
rz(1.5396176) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80910271) q[2];
sx q[2];
rz(-1.9168233) q[2];
sx q[2];
rz(-1.5835294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7120413) q[1];
sx q[1];
rz(-1.6152641) q[1];
sx q[1];
rz(3.0958423) q[1];
x q[2];
rz(-2.4828033) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-2.4023138) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(0.49490067) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98289821) q[0];
sx q[0];
rz(-1.6086846) q[0];
sx q[0];
rz(1.0844896) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.020936326) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(-0.39839572) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4329263) q[1];
sx q[1];
rz(-2.4927944) q[1];
sx q[1];
rz(-1.7556346) q[1];
rz(-1.6845735) q[3];
sx q[3];
rz(-2.4416231) q[3];
sx q[3];
rz(-1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(-2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-0.83256759) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-0.20559786) q[2];
sx q[2];
rz(-2.1046706) q[2];
sx q[2];
rz(-0.54866366) q[2];
rz(-1.8264063) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
