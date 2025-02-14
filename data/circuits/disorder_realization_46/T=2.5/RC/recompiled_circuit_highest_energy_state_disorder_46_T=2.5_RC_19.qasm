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
rz(-0.30204371) q[0];
sx q[0];
rz(1.8855236) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(1.8355651) q[1];
sx q[1];
rz(10.553283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23763424) q[0];
sx q[0];
rz(-0.65017831) q[0];
sx q[0];
rz(1.0976492) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41646012) q[2];
sx q[2];
rz(-2.4774733) q[2];
sx q[2];
rz(2.0478319) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0417718) q[1];
sx q[1];
rz(-1.9536293) q[1];
sx q[1];
rz(2.113093) q[1];
rz(2.5079191) q[3];
sx q[3];
rz(-1.2094991) q[3];
sx q[3];
rz(-0.47508815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5350554) q[2];
sx q[2];
rz(-2.8474319) q[2];
sx q[2];
rz(-1.199022) q[2];
rz(-2.8626275) q[3];
sx q[3];
rz(-1.6821034) q[3];
sx q[3];
rz(0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26175427) q[0];
sx q[0];
rz(-2.6848875) q[0];
sx q[0];
rz(0.36619151) q[0];
rz(-1.2884033) q[1];
sx q[1];
rz(-2.2346456) q[1];
sx q[1];
rz(-0.2624951) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703471) q[0];
sx q[0];
rz(-0.34929212) q[0];
sx q[0];
rz(-2.5906934) q[0];
x q[1];
rz(-3.0832564) q[2];
sx q[2];
rz(-1.3223786) q[2];
sx q[2];
rz(2.2417541) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69345638) q[1];
sx q[1];
rz(-2.6779406) q[1];
sx q[1];
rz(-0.69000375) q[1];
rz(-1.4743342) q[3];
sx q[3];
rz(-0.25611195) q[3];
sx q[3];
rz(-0.90463582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75258201) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(-2.9627723) q[2];
rz(-3.1255152) q[3];
sx q[3];
rz(-0.5355081) q[3];
sx q[3];
rz(-2.2107562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43612424) q[0];
sx q[0];
rz(-0.12691623) q[0];
sx q[0];
rz(-0.57417589) q[0];
rz(0.15159675) q[1];
sx q[1];
rz(-0.44436374) q[1];
sx q[1];
rz(1.9134329) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41802619) q[0];
sx q[0];
rz(-2.0834588) q[0];
sx q[0];
rz(-2.961757) q[0];
rz(-pi) q[1];
rz(3.0730293) q[2];
sx q[2];
rz(-0.97899635) q[2];
sx q[2];
rz(-0.14321274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0361402) q[1];
sx q[1];
rz(-0.6195809) q[1];
sx q[1];
rz(-1.416227) q[1];
rz(-2.853573) q[3];
sx q[3];
rz(-1.1041118) q[3];
sx q[3];
rz(2.6096901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0000618) q[2];
sx q[2];
rz(-2.5059293) q[2];
sx q[2];
rz(2.2504811) q[2];
rz(-0.38854232) q[3];
sx q[3];
rz(-1.5107692) q[3];
sx q[3];
rz(2.8019606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.6252839) q[0];
sx q[0];
rz(-3.0966274) q[0];
sx q[0];
rz(-2.6594824) q[0];
rz(-2.4730143) q[1];
sx q[1];
rz(-1.6079638) q[1];
sx q[1];
rz(0.63749981) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26698819) q[0];
sx q[0];
rz(-0.98867304) q[0];
sx q[0];
rz(0.28018392) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89738557) q[2];
sx q[2];
rz(-1.4801133) q[2];
sx q[2];
rz(0.53340837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1177534) q[1];
sx q[1];
rz(-1.6477403) q[1];
sx q[1];
rz(-1.6047828) q[1];
rz(-pi) q[2];
rz(2.8907475) q[3];
sx q[3];
rz(-1.9168087) q[3];
sx q[3];
rz(0.097867997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1981226) q[2];
sx q[2];
rz(-1.024647) q[2];
sx q[2];
rz(-3.0812954) q[2];
rz(0.30334011) q[3];
sx q[3];
rz(-3.0637432) q[3];
sx q[3];
rz(-0.2651324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9918793) q[0];
sx q[0];
rz(-2.7880221) q[0];
sx q[0];
rz(-0.41325945) q[0];
rz(0.39516559) q[1];
sx q[1];
rz(-1.7839849) q[1];
sx q[1];
rz(-2.1392335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7189594) q[0];
sx q[0];
rz(-2.0074685) q[0];
sx q[0];
rz(2.5366001) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4031314) q[2];
sx q[2];
rz(-2.9288026) q[2];
sx q[2];
rz(2.8820702) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9262979) q[1];
sx q[1];
rz(-1.3560672) q[1];
sx q[1];
rz(-0.88820998) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20805863) q[3];
sx q[3];
rz(-1.8866619) q[3];
sx q[3];
rz(-2.7622595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5137382) q[2];
sx q[2];
rz(-1.0971053) q[2];
sx q[2];
rz(-1.6045137) q[2];
rz(1.1995992) q[3];
sx q[3];
rz(-2.0204085) q[3];
sx q[3];
rz(-0.77170038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5143249) q[0];
sx q[0];
rz(-1.5031313) q[0];
sx q[0];
rz(-0.36870536) q[0];
rz(0.74327028) q[1];
sx q[1];
rz(-2.5786046) q[1];
sx q[1];
rz(2.7875913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6021332) q[0];
sx q[0];
rz(-1.51855) q[0];
sx q[0];
rz(-2.2118072) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2576302) q[2];
sx q[2];
rz(-2.0764362) q[2];
sx q[2];
rz(-0.87228105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.94985574) q[1];
sx q[1];
rz(-2.5402148) q[1];
sx q[1];
rz(-1.8888425) q[1];
x q[2];
rz(-2.6383868) q[3];
sx q[3];
rz(-0.70069289) q[3];
sx q[3];
rz(-1.038674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.02504286) q[2];
sx q[2];
rz(-1.9321238) q[2];
sx q[2];
rz(-3.1361191) q[2];
rz(0.51760751) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(0.027211729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590498) q[0];
sx q[0];
rz(-2.5820177) q[0];
sx q[0];
rz(-2.9285808) q[0];
rz(1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(2.2921553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54100688) q[0];
sx q[0];
rz(-1.5274421) q[0];
sx q[0];
rz(-2.5227273) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.527571) q[2];
sx q[2];
rz(-1.6014978) q[2];
sx q[2];
rz(0.30393201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.28038) q[1];
sx q[1];
rz(-2.0191128) q[1];
sx q[1];
rz(-2.6157152) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8324424) q[3];
sx q[3];
rz(-2.5980686) q[3];
sx q[3];
rz(1.8147008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9269632) q[2];
sx q[2];
rz(-1.9081076) q[2];
sx q[2];
rz(-0.54369533) q[2];
rz(0.26116535) q[3];
sx q[3];
rz(-1.5615014) q[3];
sx q[3];
rz(2.9847667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549304) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(-0.18761158) q[0];
rz(-2.018351) q[1];
sx q[1];
rz(-2.3767327) q[1];
sx q[1];
rz(-2.979028) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99062067) q[0];
sx q[0];
rz(-2.2216068) q[0];
sx q[0];
rz(-1.9810173) q[0];
x q[1];
rz(-0.17485042) q[2];
sx q[2];
rz(-1.0533337) q[2];
sx q[2];
rz(2.3558921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29926649) q[1];
sx q[1];
rz(-0.61533058) q[1];
sx q[1];
rz(-1.4922813) q[1];
rz(-pi) q[2];
rz(-3.0288234) q[3];
sx q[3];
rz(-0.9774607) q[3];
sx q[3];
rz(-1.5055553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8673914) q[2];
sx q[2];
rz(-2.826639) q[2];
sx q[2];
rz(2.9401722) q[2];
rz(-2.6650186) q[3];
sx q[3];
rz(-1.3786992) q[3];
sx q[3];
rz(0.99982888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808712) q[0];
sx q[0];
rz(-2.9044386) q[0];
sx q[0];
rz(-0.57688212) q[0];
rz(-2.8374529) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(-2.0374128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21872422) q[0];
sx q[0];
rz(-1.6593242) q[0];
sx q[0];
rz(-1.9506111) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9530802) q[2];
sx q[2];
rz(-2.5645263) q[2];
sx q[2];
rz(0.12332502) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6526281) q[1];
sx q[1];
rz(-1.9581989) q[1];
sx q[1];
rz(-0.94469686) q[1];
rz(-1.6100092) q[3];
sx q[3];
rz(-0.3914507) q[3];
sx q[3];
rz(1.9658364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1748109) q[2];
sx q[2];
rz(-2.8947783) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(-1.5934058) q[3];
sx q[3];
rz(-2.4631409) q[3];
sx q[3];
rz(2.9233176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19112793) q[0];
sx q[0];
rz(-0.8466962) q[0];
sx q[0];
rz(0.51093131) q[0];
rz(-1.9966985) q[1];
sx q[1];
rz(-1.9547918) q[1];
sx q[1];
rz(-1.5085545) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16647969) q[0];
sx q[0];
rz(-0.60752019) q[0];
sx q[0];
rz(0.44618531) q[0];
rz(-1.0867003) q[2];
sx q[2];
rz(-1.8082613) q[2];
sx q[2];
rz(2.9959842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33433769) q[1];
sx q[1];
rz(-0.48408135) q[1];
sx q[1];
rz(-0.38881199) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3168961) q[3];
sx q[3];
rz(-1.861737) q[3];
sx q[3];
rz(-1.2856785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7213584) q[2];
sx q[2];
rz(-2.5098462) q[2];
sx q[2];
rz(1.4314123) q[2];
rz(0.015627705) q[3];
sx q[3];
rz(-1.2651919) q[3];
sx q[3];
rz(2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2331727) q[0];
sx q[0];
rz(-1.40042) q[0];
sx q[0];
rz(-0.9216876) q[0];
rz(1.5064916) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(2.2057381) q[2];
sx q[2];
rz(-2.6781185) q[2];
sx q[2];
rz(0.92739633) q[2];
rz(-3.1063632) q[3];
sx q[3];
rz(-1.4932125) q[3];
sx q[3];
rz(1.6532547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
