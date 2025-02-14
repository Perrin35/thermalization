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
rz(0.47705874) q[0];
sx q[0];
rz(2.8395489) q[0];
sx q[0];
rz(10.680847) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(1.8355651) q[1];
sx q[1];
rz(10.553283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1953729) q[0];
sx q[0];
rz(-1.2913307) q[0];
sx q[0];
rz(2.165876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62126764) q[2];
sx q[2];
rz(-1.8227907) q[2];
sx q[2];
rz(-2.9997343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8339631) q[1];
sx q[1];
rz(-1.071573) q[1];
sx q[1];
rz(-2.7021033) q[1];
x q[2];
rz(-2.009274) q[3];
sx q[3];
rz(-0.98377675) q[3];
sx q[3];
rz(-1.349821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6065373) q[2];
sx q[2];
rz(-0.29416072) q[2];
sx q[2];
rz(1.199022) q[2];
rz(-0.27896518) q[3];
sx q[3];
rz(-1.4594892) q[3];
sx q[3];
rz(0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8798384) q[0];
sx q[0];
rz(-2.6848875) q[0];
sx q[0];
rz(2.7754011) q[0];
rz(1.2884033) q[1];
sx q[1];
rz(-2.2346456) q[1];
sx q[1];
rz(0.2624951) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2912882) q[0];
sx q[0];
rz(-1.8666963) q[0];
sx q[0];
rz(1.759191) q[0];
x q[1];
rz(1.8196202) q[2];
sx q[2];
rz(-1.5142528) q[2];
sx q[2];
rz(2.4849934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7028864) q[1];
sx q[1];
rz(-1.2186495) q[1];
sx q[1];
rz(-1.8789324) q[1];
rz(-pi) q[2];
rz(-1.6672584) q[3];
sx q[3];
rz(-0.25611195) q[3];
sx q[3];
rz(0.90463582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75258201) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(-2.9627723) q[2];
rz(-0.016077476) q[3];
sx q[3];
rz(-2.6060846) q[3];
sx q[3];
rz(-2.2107562) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7054684) q[0];
sx q[0];
rz(-3.0146764) q[0];
sx q[0];
rz(0.57417589) q[0];
rz(2.9899959) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(-1.2281598) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41802619) q[0];
sx q[0];
rz(-1.0581338) q[0];
sx q[0];
rz(0.17983564) q[0];
rz(1.6723694) q[2];
sx q[2];
rz(-0.59528643) q[2];
sx q[2];
rz(-0.02073076) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5500726) q[1];
sx q[1];
rz(-1.6603163) q[1];
sx q[1];
rz(0.95687434) q[1];
rz(1.0574041) q[3];
sx q[3];
rz(-2.5988165) q[3];
sx q[3];
rz(-2.0273939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0000618) q[2];
sx q[2];
rz(-0.63566339) q[2];
sx q[2];
rz(2.2504811) q[2];
rz(-2.7530503) q[3];
sx q[3];
rz(-1.6308234) q[3];
sx q[3];
rz(2.8019606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-1.5163088) q[0];
sx q[0];
rz(-0.044965222) q[0];
sx q[0];
rz(-2.6594824) q[0];
rz(0.66857839) q[1];
sx q[1];
rz(-1.6079638) q[1];
sx q[1];
rz(0.63749981) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6808785) q[0];
sx q[0];
rz(-1.3377046) q[0];
sx q[0];
rz(-2.1713269) q[0];
x q[1];
rz(-2.2442071) q[2];
sx q[2];
rz(-1.4801133) q[2];
sx q[2];
rz(2.6081843) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1177534) q[1];
sx q[1];
rz(-1.4938524) q[1];
sx q[1];
rz(-1.6047828) q[1];
rz(2.8907475) q[3];
sx q[3];
rz(-1.224784) q[3];
sx q[3];
rz(3.0437247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1981226) q[2];
sx q[2];
rz(-2.1169457) q[2];
sx q[2];
rz(-0.060297273) q[2];
rz(-0.30334011) q[3];
sx q[3];
rz(-0.077849418) q[3];
sx q[3];
rz(-0.2651324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9918793) q[0];
sx q[0];
rz(-0.35357058) q[0];
sx q[0];
rz(-0.41325945) q[0];
rz(0.39516559) q[1];
sx q[1];
rz(-1.3576077) q[1];
sx q[1];
rz(2.1392335) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43267879) q[0];
sx q[0];
rz(-1.029338) q[0];
sx q[0];
rz(2.0869291) q[0];
x q[1];
rz(3.1055519) q[2];
sx q[2];
rz(-1.7805575) q[2];
sx q[2];
rz(0.4309817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2152948) q[1];
sx q[1];
rz(-1.7855254) q[1];
sx q[1];
rz(0.88820998) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0071383) q[3];
sx q[3];
rz(-2.7652869) q[3];
sx q[3];
rz(0.2175771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5137382) q[2];
sx q[2];
rz(-2.0444874) q[2];
sx q[2];
rz(-1.537079) q[2];
rz(-1.1995992) q[3];
sx q[3];
rz(-1.1211841) q[3];
sx q[3];
rz(-0.77170038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62726778) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(-2.7728873) q[0];
rz(-0.74327028) q[1];
sx q[1];
rz(-0.5629881) q[1];
sx q[1];
rz(2.7875913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10121541) q[0];
sx q[0];
rz(-0.64283744) q[0];
sx q[0];
rz(1.6580216) q[0];
rz(-pi) q[1];
rz(-1.8839624) q[2];
sx q[2];
rz(-2.0764362) q[2];
sx q[2];
rz(-0.87228105) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5716155) q[1];
sx q[1];
rz(-2.1381525) q[1];
sx q[1];
rz(2.9302271) q[1];
x q[2];
rz(-2.5052222) q[3];
sx q[3];
rz(-1.2546347) q[3];
sx q[3];
rz(-3.0077601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.02504286) q[2];
sx q[2];
rz(-1.2094689) q[2];
sx q[2];
rz(-0.0054736007) q[2];
rz(-2.6239851) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(-3.1143809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7356877) q[0];
sx q[0];
rz(-0.55957496) q[0];
sx q[0];
rz(-2.9285808) q[0];
rz(1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(2.2921553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54100688) q[0];
sx q[0];
rz(-1.6141506) q[0];
sx q[0];
rz(2.5227273) q[0];
x q[1];
rz(-0.95291887) q[2];
sx q[2];
rz(-3.0885792) q[2];
sx q[2];
rz(1.2575146) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35104677) q[1];
sx q[1];
rz(-0.67706579) q[1];
sx q[1];
rz(-0.76404567) q[1];
x q[2];
rz(-1.8324424) q[3];
sx q[3];
rz(-0.54352409) q[3];
sx q[3];
rz(-1.8147008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21462943) q[2];
sx q[2];
rz(-1.9081076) q[2];
sx q[2];
rz(-2.5978973) q[2];
rz(2.8804273) q[3];
sx q[3];
rz(-1.5615014) q[3];
sx q[3];
rz(-2.9847667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9549304) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(-2.9539811) q[0];
rz(-2.018351) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(2.979028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.613224) q[0];
sx q[0];
rz(-0.75308164) q[0];
sx q[0];
rz(2.6591419) q[0];
x q[1];
rz(1.2741873) q[2];
sx q[2];
rz(-2.5979418) q[2];
sx q[2];
rz(2.6988876) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2073719) q[1];
sx q[1];
rz(-1.5255063) q[1];
sx q[1];
rz(0.95691935) q[1];
rz(-pi) q[2];
rz(-1.4054589) q[3];
sx q[3];
rz(-2.5389034) q[3];
sx q[3];
rz(-1.7053982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8673914) q[2];
sx q[2];
rz(-2.826639) q[2];
sx q[2];
rz(2.9401722) q[2];
rz(0.47657403) q[3];
sx q[3];
rz(-1.7628935) q[3];
sx q[3];
rz(2.1417638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808712) q[0];
sx q[0];
rz(-0.23715401) q[0];
sx q[0];
rz(-2.5647105) q[0];
rz(-2.8374529) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(1.1041799) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9228684) q[0];
sx q[0];
rz(-1.4822685) q[0];
sx q[0];
rz(1.1909816) q[0];
rz(-0.18851243) q[2];
sx q[2];
rz(-2.5645263) q[2];
sx q[2];
rz(3.0182676) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7930709) q[1];
sx q[1];
rz(-0.99731113) q[1];
sx q[1];
rz(-2.6751233) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9619759) q[3];
sx q[3];
rz(-1.5558387) q[3];
sx q[3];
rz(-2.7103031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1748109) q[2];
sx q[2];
rz(-2.8947783) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(-1.5934058) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(-0.21827503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9504647) q[0];
sx q[0];
rz(-2.2948965) q[0];
sx q[0];
rz(-2.6306613) q[0];
rz(-1.9966985) q[1];
sx q[1];
rz(-1.1868008) q[1];
sx q[1];
rz(1.5085545) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36105267) q[0];
sx q[0];
rz(-1.0298488) q[0];
sx q[0];
rz(-1.279328) q[0];
x q[1];
rz(1.0912446) q[2];
sx q[2];
rz(-2.6065718) q[2];
sx q[2];
rz(-1.0045235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76782334) q[1];
sx q[1];
rz(-1.1255742) q[1];
sx q[1];
rz(1.3740463) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3168961) q[3];
sx q[3];
rz(-1.2798556) q[3];
sx q[3];
rz(-1.8559141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7213584) q[2];
sx q[2];
rz(-0.63174641) q[2];
sx q[2];
rz(-1.4314123) q[2];
rz(3.1259649) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(0.50610745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90841993) q[0];
sx q[0];
rz(-1.40042) q[0];
sx q[0];
rz(-0.9216876) q[0];
rz(-1.5064916) q[1];
sx q[1];
rz(-1.9252621) q[1];
sx q[1];
rz(2.6797934) q[1];
rz(-2.2057381) q[2];
sx q[2];
rz(-0.46347413) q[2];
sx q[2];
rz(-2.2141963) q[2];
rz(1.4931645) q[3];
sx q[3];
rz(-1.5356728) q[3];
sx q[3];
rz(0.079726797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
