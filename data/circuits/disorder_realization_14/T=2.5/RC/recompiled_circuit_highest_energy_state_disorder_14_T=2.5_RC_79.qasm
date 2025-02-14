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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(-0.60751539) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018463919) q[0];
sx q[0];
rz(-1.0842609) q[0];
sx q[0];
rz(-2.8790265) q[0];
rz(-pi) q[1];
rz(0.11694853) q[2];
sx q[2];
rz(-2.0772572) q[2];
sx q[2];
rz(-3.1053271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7360037) q[1];
sx q[1];
rz(-2.0479255) q[1];
sx q[1];
rz(1.3911584) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.257454) q[3];
sx q[3];
rz(-0.55301412) q[3];
sx q[3];
rz(0.64698863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5241549) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(2.4784135) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8993768) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(-2.2826165) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(-1.4069517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024723) q[0];
sx q[0];
rz(-2.6022096) q[0];
sx q[0];
rz(1.1155737) q[0];
rz(1.2386049) q[2];
sx q[2];
rz(-1.536035) q[2];
sx q[2];
rz(1.3369651) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8415815) q[1];
sx q[1];
rz(-1.661013) q[1];
sx q[1];
rz(-1.0567699) q[1];
rz(1.155726) q[3];
sx q[3];
rz(-1.0336116) q[3];
sx q[3];
rz(1.9946757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0824288) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(1.3272237) q[2];
rz(2.0969157) q[3];
sx q[3];
rz(-0.71377126) q[3];
sx q[3];
rz(2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5663719) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(-2.7161993) q[0];
rz(-1.3771903) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(1.4345217) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4664513) q[0];
sx q[0];
rz(-1.688862) q[0];
sx q[0];
rz(0.64339126) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4332831) q[2];
sx q[2];
rz(-1.5454486) q[2];
sx q[2];
rz(-0.011653221) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7075478) q[1];
sx q[1];
rz(-0.92869379) q[1];
sx q[1];
rz(1.0471859) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9886964) q[3];
sx q[3];
rz(-2.4800081) q[3];
sx q[3];
rz(2.083287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.44125685) q[2];
sx q[2];
rz(-1.2563027) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(-1.357632) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.6443845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212379) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(2.3714016) q[0];
rz(-2.1417446) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(1.3410478) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5504549) q[0];
sx q[0];
rz(-1.3584175) q[0];
sx q[0];
rz(-0.55340931) q[0];
rz(-pi) q[1];
rz(-2.6821939) q[2];
sx q[2];
rz(-0.93743491) q[2];
sx q[2];
rz(2.6448696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25136687) q[1];
sx q[1];
rz(-0.40491762) q[1];
sx q[1];
rz(-2.4144854) q[1];
x q[2];
rz(1.8561268) q[3];
sx q[3];
rz(-1.4754681) q[3];
sx q[3];
rz(2.8872629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8009214) q[2];
sx q[2];
rz(-2.0719216) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(-0.10041222) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(-1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(-1.7452128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76043789) q[0];
sx q[0];
rz(-1.9312857) q[0];
sx q[0];
rz(1.1973778) q[0];
rz(-1.4047844) q[2];
sx q[2];
rz(-1.976227) q[2];
sx q[2];
rz(2.977598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3815617) q[1];
sx q[1];
rz(-2.2027822) q[1];
sx q[1];
rz(1.6595592) q[1];
rz(-3.1110686) q[3];
sx q[3];
rz(-2.3168193) q[3];
sx q[3];
rz(-1.8061226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9153626) q[2];
sx q[2];
rz(-1.1504983) q[2];
sx q[2];
rz(-0.076233141) q[2];
rz(-1.7539615) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(-1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607894) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(1.8060818) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(2.9551771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7714027) q[0];
sx q[0];
rz(-2.0982842) q[0];
sx q[0];
rz(-2.6341088) q[0];
x q[1];
rz(2.8373986) q[2];
sx q[2];
rz(-2.549571) q[2];
sx q[2];
rz(-1.6006058) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.586372) q[1];
sx q[1];
rz(-2.456291) q[1];
sx q[1];
rz(1.1838811) q[1];
rz(-pi) q[2];
rz(2.3535054) q[3];
sx q[3];
rz(-2.4895146) q[3];
sx q[3];
rz(1.2913845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67941252) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(-1.1421674) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(2.5073945) q[0];
rz(-1.0448666) q[1];
sx q[1];
rz(-2.2506782) q[1];
sx q[1];
rz(-2.1814836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7564108) q[0];
sx q[0];
rz(-0.39968458) q[0];
sx q[0];
rz(-1.5002285) q[0];
rz(-1.5788001) q[2];
sx q[2];
rz(-1.8718613) q[2];
sx q[2];
rz(-1.0775527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3893435) q[1];
sx q[1];
rz(-1.5823808) q[1];
sx q[1];
rz(-1.0453792) q[1];
rz(-pi) q[2];
rz(-0.74713196) q[3];
sx q[3];
rz(-0.91643667) q[3];
sx q[3];
rz(0.33734712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94379696) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(1.724285) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5328131) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(0.10803647) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(1.1688165) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625789) q[0];
sx q[0];
rz(-2.4860408) q[0];
sx q[0];
rz(-0.84540875) q[0];
rz(0.76185267) q[2];
sx q[2];
rz(-1.5707865) q[2];
sx q[2];
rz(0.61864432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7119693) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(1.3369249) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5904558) q[3];
sx q[3];
rz(-1.0791313) q[3];
sx q[3];
rz(1.2439072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.03269) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(-0.6558134) q[2];
rz(-0.10410318) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8769237) q[0];
sx q[0];
rz(-1.8380565) q[0];
sx q[0];
rz(1.6498097) q[0];
rz(-2.1022294) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(-2.6731491) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85850785) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(-1.5751189) q[0];
rz(-pi) q[1];
rz(-2.5889187) q[2];
sx q[2];
rz(-0.50853339) q[2];
sx q[2];
rz(0.39297152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8711618) q[1];
sx q[1];
rz(-2.3498145) q[1];
sx q[1];
rz(1.7986078) q[1];
rz(0.78448589) q[3];
sx q[3];
rz(-0.38803369) q[3];
sx q[3];
rz(1.7540384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.09482065) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(0.92791933) q[2];
rz(1.8038484) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(1.6524338) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0874262) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(-1.077865) q[0];
rz(-2.7742591) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.2841388) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7736762) q[0];
sx q[0];
rz(-1.319029) q[0];
sx q[0];
rz(1.0836224) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91725332) q[2];
sx q[2];
rz(-1.6628569) q[2];
sx q[2];
rz(1.1600509) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91032797) q[1];
sx q[1];
rz(-0.18016768) q[1];
sx q[1];
rz(2.641201) q[1];
x q[2];
rz(-0.73478847) q[3];
sx q[3];
rz(-1.7289203) q[3];
sx q[3];
rz(1.2722335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1401356) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(-0.4293116) q[2];
rz(0.15520994) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(1.8163053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406381) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(2.2776729) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(0.39045329) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(-0.94384296) q[3];
sx q[3];
rz(-1.2921492) q[3];
sx q[3];
rz(2.9155801) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
