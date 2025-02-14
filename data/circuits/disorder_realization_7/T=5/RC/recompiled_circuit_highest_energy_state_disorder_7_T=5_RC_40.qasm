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
rz(-3.0149674) q[0];
sx q[0];
rz(-1.5746483) q[0];
sx q[0];
rz(-2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(5.5362267) q[1];
sx q[1];
rz(9.8510392) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3328457) q[0];
sx q[0];
rz(-1.3984183) q[0];
sx q[0];
rz(0.7953978) q[0];
x q[1];
rz(-1.2834163) q[2];
sx q[2];
rz(-1.0848019) q[2];
sx q[2];
rz(-3.0241242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1805686) q[1];
sx q[1];
rz(-2.5371309) q[1];
sx q[1];
rz(-3.040041) q[1];
rz(-pi) q[2];
rz(1.0996885) q[3];
sx q[3];
rz(-0.63835164) q[3];
sx q[3];
rz(-1.9698576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2446186) q[2];
sx q[2];
rz(-1.3099058) q[2];
sx q[2];
rz(0.2429602) q[2];
rz(-2.7729559) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(-1.2322371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.31545562) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(0.36732236) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(-2.499089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272907) q[0];
sx q[0];
rz(-1.1513984) q[0];
sx q[0];
rz(2.6227996) q[0];
x q[1];
rz(1.1682061) q[2];
sx q[2];
rz(-2.085683) q[2];
sx q[2];
rz(1.3001315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52740188) q[1];
sx q[1];
rz(-2.3878015) q[1];
sx q[1];
rz(2.4587653) q[1];
rz(1.8090463) q[3];
sx q[3];
rz(-0.49360156) q[3];
sx q[3];
rz(2.3675338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.502304) q[2];
sx q[2];
rz(-0.93561155) q[2];
sx q[2];
rz(1.7712234) q[2];
rz(0.227452) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(-1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3962536) q[0];
sx q[0];
rz(-2.9662913) q[0];
sx q[0];
rz(1.1981717) q[0];
rz(2.0612969) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(3.0467196) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6300044) q[0];
sx q[0];
rz(-2.5325724) q[0];
sx q[0];
rz(-2.0320351) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40887654) q[2];
sx q[2];
rz(-2.1890702) q[2];
sx q[2];
rz(2.8586819) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.171189) q[1];
sx q[1];
rz(-0.53817777) q[1];
sx q[1];
rz(-1.490905) q[1];
x q[2];
rz(-0.27697825) q[3];
sx q[3];
rz(-1.9414895) q[3];
sx q[3];
rz(2.1408368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0692811) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(0.4064694) q[2];
rz(1.1786002) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.1627722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56226319) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(-0.33600268) q[0];
rz(2.621189) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(3.0373108) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2662243) q[0];
sx q[0];
rz(-1.2422529) q[0];
sx q[0];
rz(1.203241) q[0];
rz(-pi) q[1];
rz(-1.5566795) q[2];
sx q[2];
rz(-1.4143362) q[2];
sx q[2];
rz(0.9597646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49890624) q[1];
sx q[1];
rz(-2.3489174) q[1];
sx q[1];
rz(-1.0426056) q[1];
rz(-0.23747344) q[3];
sx q[3];
rz(-2.6769612) q[3];
sx q[3];
rz(-2.1647418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61863724) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.9970419) q[2];
rz(2.5942904) q[3];
sx q[3];
rz(-1.9132883) q[3];
sx q[3];
rz(-1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(0.8892141) q[0];
sx q[0];
rz(-0.19110282) q[0];
sx q[0];
rz(-0.55437535) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.8781885) q[1];
sx q[1];
rz(-2.8401781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9145935) q[0];
sx q[0];
rz(-2.1560139) q[0];
sx q[0];
rz(-0.27497681) q[0];
rz(-pi) q[1];
rz(-0.38827916) q[2];
sx q[2];
rz(-1.3781606) q[2];
sx q[2];
rz(-1.9658058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6993048) q[1];
sx q[1];
rz(-1.3584474) q[1];
sx q[1];
rz(1.1816756) q[1];
rz(-pi) q[2];
rz(-0.88121342) q[3];
sx q[3];
rz(-0.47166079) q[3];
sx q[3];
rz(-0.72615964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.021412795) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-3.1346698) q[2];
rz(-0.93112469) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-1.1091703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865006) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(1.6078) q[0];
rz(-0.7026698) q[1];
sx q[1];
rz(-2.6103795) q[1];
sx q[1];
rz(-2.839397) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5836499) q[0];
sx q[0];
rz(-1.6725701) q[0];
sx q[0];
rz(-0.90910615) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2185368) q[2];
sx q[2];
rz(-2.2552236) q[2];
sx q[2];
rz(-0.5796488) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.211765) q[1];
sx q[1];
rz(-0.95410871) q[1];
sx q[1];
rz(0.21287983) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4190361) q[3];
sx q[3];
rz(-1.6223063) q[3];
sx q[3];
rz(2.6035978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2493784) q[2];
sx q[2];
rz(-2.0912781) q[2];
sx q[2];
rz(0.92602473) q[2];
rz(-1.9269491) q[3];
sx q[3];
rz(-2.6483783) q[3];
sx q[3];
rz(-0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72325426) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(2.8077937) q[1];
sx q[1];
rz(-1.326694) q[1];
sx q[1];
rz(1.0858067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5949149) q[0];
sx q[0];
rz(-1.4913861) q[0];
sx q[0];
rz(-2.9360442) q[0];
rz(-pi) q[1];
x q[1];
rz(1.169431) q[2];
sx q[2];
rz(-1.4039206) q[2];
sx q[2];
rz(-2.3292975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2183509) q[1];
sx q[1];
rz(-0.78303799) q[1];
sx q[1];
rz(2.3438575) q[1];
rz(-pi) q[2];
rz(2.2074039) q[3];
sx q[3];
rz(-1.459834) q[3];
sx q[3];
rz(2.8257089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51599017) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(1.9773352) q[2];
rz(0.26816756) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5371573) q[0];
sx q[0];
rz(-0.51979655) q[0];
sx q[0];
rz(-1.9926158) q[0];
rz(-0.2941429) q[1];
sx q[1];
rz(-1.5584757) q[1];
sx q[1];
rz(2.099096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54759083) q[0];
sx q[0];
rz(-0.85193513) q[0];
sx q[0];
rz(0.44525624) q[0];
rz(-2.9127321) q[2];
sx q[2];
rz(-2.3386152) q[2];
sx q[2];
rz(2.2265018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5120705) q[1];
sx q[1];
rz(-2.7334556) q[1];
sx q[1];
rz(-0.44993181) q[1];
rz(-pi) q[2];
rz(2.2348079) q[3];
sx q[3];
rz(-1.5935002) q[3];
sx q[3];
rz(-1.5496538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(1.3647122) q[2];
rz(-2.4890066) q[3];
sx q[3];
rz(-2.3502974) q[3];
sx q[3];
rz(-0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.515601) q[0];
sx q[0];
rz(-2.7643272) q[0];
sx q[0];
rz(3.1242477) q[0];
rz(0.01677244) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(-0.15388547) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2543593) q[0];
sx q[0];
rz(-1.6463929) q[0];
sx q[0];
rz(-2.8235954) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4520313) q[2];
sx q[2];
rz(-0.84367263) q[2];
sx q[2];
rz(0.52983701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20275234) q[1];
sx q[1];
rz(-1.2766663) q[1];
sx q[1];
rz(-2.7619656) q[1];
rz(0.023970402) q[3];
sx q[3];
rz(-1.6466738) q[3];
sx q[3];
rz(-0.45437231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1880356) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(0.44450644) q[2];
rz(-0.19017531) q[3];
sx q[3];
rz(-1.5835652) q[3];
sx q[3];
rz(-1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.001215) q[0];
sx q[0];
rz(-1.8920521) q[0];
sx q[0];
rz(2.0709399) q[0];
rz(-3.095678) q[1];
sx q[1];
rz(-1.670198) q[1];
sx q[1];
rz(-2.3505223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2915724) q[0];
sx q[0];
rz(-0.34713161) q[0];
sx q[0];
rz(-1.8308543) q[0];
x q[1];
rz(1.3020817) q[2];
sx q[2];
rz(-1.955532) q[2];
sx q[2];
rz(3.0962944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4262096) q[1];
sx q[1];
rz(-1.9265623) q[1];
sx q[1];
rz(-0.4224311) q[1];
x q[2];
rz(1.3338575) q[3];
sx q[3];
rz(-0.51722368) q[3];
sx q[3];
rz(2.5631529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.644824) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(1.6688639) q[2];
rz(-1.5276927) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(1.9014026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.938217) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(-0.51364246) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(-0.41507872) q[2];
sx q[2];
rz(-0.54176258) q[2];
sx q[2];
rz(-1.1272346) q[2];
rz(-2.0075825) q[3];
sx q[3];
rz(-0.81098771) q[3];
sx q[3];
rz(0.24843957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
