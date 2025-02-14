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
rz(-1.0492078) q[0];
sx q[0];
rz(-2.3823491) q[0];
sx q[0];
rz(0.37641755) q[0];
rz(1.5578101) q[1];
sx q[1];
rz(-1.0695142) q[1];
sx q[1];
rz(0.95626107) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.555684) q[0];
sx q[0];
rz(-0.95147419) q[0];
sx q[0];
rz(2.7509718) q[0];
rz(-pi) q[1];
rz(0.42266385) q[2];
sx q[2];
rz(-1.3939121) q[2];
sx q[2];
rz(1.6852578) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1714345) q[1];
sx q[1];
rz(-1.4490607) q[1];
sx q[1];
rz(2.3129169) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3322863) q[3];
sx q[3];
rz(-1.3068891) q[3];
sx q[3];
rz(-1.1562386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3236397) q[2];
sx q[2];
rz(-2.2223667) q[2];
sx q[2];
rz(1.9383355) q[2];
rz(-2.1008927) q[3];
sx q[3];
rz(-0.87472707) q[3];
sx q[3];
rz(0.11999764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1131209) q[0];
sx q[0];
rz(-0.57336837) q[0];
sx q[0];
rz(-2.0290802) q[0];
rz(1.6587229) q[1];
sx q[1];
rz(-0.590938) q[1];
sx q[1];
rz(-2.9842751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33799115) q[0];
sx q[0];
rz(-2.4287619) q[0];
sx q[0];
rz(0.26960752) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6578761) q[2];
sx q[2];
rz(-1.198999) q[2];
sx q[2];
rz(1.7024937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8437905) q[1];
sx q[1];
rz(-0.71759598) q[1];
sx q[1];
rz(1.1949431) q[1];
x q[2];
rz(-2.9187536) q[3];
sx q[3];
rz(-1.6442219) q[3];
sx q[3];
rz(1.8334695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6436254) q[2];
sx q[2];
rz(-2.6622541) q[2];
sx q[2];
rz(0.41110006) q[2];
rz(0.57605612) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(-0.99803734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75329798) q[0];
sx q[0];
rz(-1.0030712) q[0];
sx q[0];
rz(2.4004747) q[0];
rz(1.6916212) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(-2.5695739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2021671) q[0];
sx q[0];
rz(-1.8644445) q[0];
sx q[0];
rz(-2.9333326) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7775615) q[2];
sx q[2];
rz(-1.8305147) q[2];
sx q[2];
rz(2.5842359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.05217543) q[1];
sx q[1];
rz(-1.1686109) q[1];
sx q[1];
rz(1.1074462) q[1];
rz(-2.1628863) q[3];
sx q[3];
rz(-2.7436069) q[3];
sx q[3];
rz(-3.132511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.584562) q[2];
sx q[2];
rz(-1.9494373) q[2];
sx q[2];
rz(-2.3438047) q[2];
rz(-1.3095464) q[3];
sx q[3];
rz(-1.2582425) q[3];
sx q[3];
rz(1.7358739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0107467) q[0];
sx q[0];
rz(-0.89711419) q[0];
sx q[0];
rz(-2.77453) q[0];
rz(-1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(-0.22148111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9629134) q[0];
sx q[0];
rz(-1.7811462) q[0];
sx q[0];
rz(1.5014929) q[0];
rz(0.34299739) q[2];
sx q[2];
rz(-0.38831899) q[2];
sx q[2];
rz(1.266154) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1180956) q[1];
sx q[1];
rz(-2.6129236) q[1];
sx q[1];
rz(-0.65331991) q[1];
x q[2];
rz(-0.52677299) q[3];
sx q[3];
rz(-1.3484517) q[3];
sx q[3];
rz(-2.6759345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14747846) q[2];
sx q[2];
rz(-1.1722112) q[2];
sx q[2];
rz(-1.8032903) q[2];
rz(-0.5021247) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(2.8978735) q[3];
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
rz(-0.62622708) q[0];
sx q[0];
rz(-0.60972649) q[0];
sx q[0];
rz(-1.6060265) q[0];
rz(-3.0961127) q[1];
sx q[1];
rz(-1.3810424) q[1];
sx q[1];
rz(0.46612003) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67959866) q[0];
sx q[0];
rz(-1.3211622) q[0];
sx q[0];
rz(-0.14139463) q[0];
rz(-pi) q[1];
rz(0.32987288) q[2];
sx q[2];
rz(-2.1876642) q[2];
sx q[2];
rz(-2.5170779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22039761) q[1];
sx q[1];
rz(-0.092893727) q[1];
sx q[1];
rz(2.9542406) q[1];
rz(-pi) q[2];
rz(-0.97753559) q[3];
sx q[3];
rz(-1.007742) q[3];
sx q[3];
rz(-0.42773358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2443627) q[2];
sx q[2];
rz(-1.4655317) q[2];
sx q[2];
rz(-0.7828632) q[2];
rz(-0.016294567) q[3];
sx q[3];
rz(-1.3085082) q[3];
sx q[3];
rz(2.079594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29430729) q[0];
sx q[0];
rz(-2.6562302) q[0];
sx q[0];
rz(0.37931994) q[0];
rz(2.8324221) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(2.9177623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4113385) q[0];
sx q[0];
rz(-2.3650432) q[0];
sx q[0];
rz(-2.0221439) q[0];
x q[1];
rz(1.3794704) q[2];
sx q[2];
rz(-1.4606287) q[2];
sx q[2];
rz(0.7502816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8221082) q[1];
sx q[1];
rz(-0.15076877) q[1];
sx q[1];
rz(-0.78298486) q[1];
rz(-2.4357314) q[3];
sx q[3];
rz(-2.0076723) q[3];
sx q[3];
rz(2.049751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0763187) q[2];
sx q[2];
rz(-1.167401) q[2];
sx q[2];
rz(-0.25263146) q[2];
rz(-2.1624055) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(2.7948936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.613649) q[0];
sx q[0];
rz(-0.76512965) q[0];
sx q[0];
rz(-1.8494404) q[0];
rz(2.6076803) q[1];
sx q[1];
rz(-1.6894059) q[1];
sx q[1];
rz(0.36010489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394516) q[0];
sx q[0];
rz(-1.9308254) q[0];
sx q[0];
rz(-0.074010297) q[0];
x q[1];
rz(-0.5582997) q[2];
sx q[2];
rz(-0.55634004) q[2];
sx q[2];
rz(0.46890989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8812027) q[1];
sx q[1];
rz(-0.61169988) q[1];
sx q[1];
rz(-2.7124497) q[1];
x q[2];
rz(-3.1084314) q[3];
sx q[3];
rz(-1.6735184) q[3];
sx q[3];
rz(0.20513137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99122938) q[2];
sx q[2];
rz(-1.6232619) q[2];
sx q[2];
rz(-0.72597996) q[2];
rz(-2.7109072) q[3];
sx q[3];
rz(-2.3613598) q[3];
sx q[3];
rz(-2.6881257) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78366572) q[0];
sx q[0];
rz(-1.6538606) q[0];
sx q[0];
rz(-0.7861535) q[0];
rz(0.17732009) q[1];
sx q[1];
rz(-1.0847849) q[1];
sx q[1];
rz(-1.0160944) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1209512) q[0];
sx q[0];
rz(-2.8403611) q[0];
sx q[0];
rz(-1.2094638) q[0];
x q[1];
rz(2.5943726) q[2];
sx q[2];
rz(-0.24697082) q[2];
sx q[2];
rz(0.17683593) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6888231) q[1];
sx q[1];
rz(-2.4773543) q[1];
sx q[1];
rz(0.027605357) q[1];
x q[2];
rz(1.949259) q[3];
sx q[3];
rz(-2.2871823) q[3];
sx q[3];
rz(-0.56179699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59755406) q[2];
sx q[2];
rz(-0.95557135) q[2];
sx q[2];
rz(-0.7698783) q[2];
rz(0.17635135) q[3];
sx q[3];
rz(-1.6917112) q[3];
sx q[3];
rz(-2.8281853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65799323) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(1.0461079) q[0];
rz(1.5484035) q[1];
sx q[1];
rz(-1.7240588) q[1];
sx q[1];
rz(-1.2215325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0964805) q[0];
sx q[0];
rz(-1.4661705) q[0];
sx q[0];
rz(-1.9278931) q[0];
x q[1];
rz(2.0030973) q[2];
sx q[2];
rz(-2.7932248) q[2];
sx q[2];
rz(-1.5846202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5719749) q[1];
sx q[1];
rz(-1.3571661) q[1];
sx q[1];
rz(-2.9353082) q[1];
x q[2];
rz(0.11193377) q[3];
sx q[3];
rz(-0.043826274) q[3];
sx q[3];
rz(2.6153713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.017642411) q[2];
sx q[2];
rz(-1.4175748) q[2];
sx q[2];
rz(1.9564015) q[2];
rz(-0.3398529) q[3];
sx q[3];
rz(-0.9797107) q[3];
sx q[3];
rz(-3.0237696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7793943) q[0];
sx q[0];
rz(-1.3267936) q[0];
sx q[0];
rz(-0.69778824) q[0];
rz(2.4367874) q[1];
sx q[1];
rz(-1.1230527) q[1];
sx q[1];
rz(2.8885081) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721203) q[0];
sx q[0];
rz(-1.4968431) q[0];
sx q[0];
rz(-1.3373313) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0519876) q[2];
sx q[2];
rz(-0.98271433) q[2];
sx q[2];
rz(1.1747509) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.733787) q[1];
sx q[1];
rz(-1.1152667) q[1];
sx q[1];
rz(-2.2940919) q[1];
x q[2];
rz(-0.28530085) q[3];
sx q[3];
rz(-0.75114512) q[3];
sx q[3];
rz(1.7342664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42085984) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(0.68812686) q[2];
rz(2.1963035) q[3];
sx q[3];
rz(-1.1752081) q[3];
sx q[3];
rz(0.92760408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.59104334) q[0];
sx q[0];
rz(-1.782438) q[0];
sx q[0];
rz(1.8989643) q[0];
rz(0.91724829) q[1];
sx q[1];
rz(-1.4731673) q[1];
sx q[1];
rz(0.57631667) q[1];
rz(1.9333712) q[2];
sx q[2];
rz(-1.0248263) q[2];
sx q[2];
rz(-0.43760763) q[2];
rz(1.3618906) q[3];
sx q[3];
rz(-2.7068797) q[3];
sx q[3];
rz(-2.5754365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
