OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(-1.1251261) q[0];
rz(1.0220802) q[1];
sx q[1];
rz(-0.66749579) q[1];
sx q[1];
rz(0.8134841) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047377) q[0];
sx q[0];
rz(-1.9362402) q[0];
sx q[0];
rz(-0.15113896) q[0];
rz(-pi) q[1];
rz(2.4762597) q[2];
sx q[2];
rz(-0.54845458) q[2];
sx q[2];
rz(-1.0496333) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95619394) q[1];
sx q[1];
rz(-0.46727249) q[1];
sx q[1];
rz(-0.33561349) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4072106) q[3];
sx q[3];
rz(-1.4814875) q[3];
sx q[3];
rz(0.79998868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(-2.2129464) q[2];
rz(1.5422025) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(1.9699875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(-3.0145338) q[0];
rz(0.98310414) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(-0.7712706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3160123) q[0];
sx q[0];
rz(-2.2716652) q[0];
sx q[0];
rz(-0.46937816) q[0];
rz(-pi) q[1];
rz(3.0808582) q[2];
sx q[2];
rz(-1.9078476) q[2];
sx q[2];
rz(-0.13523808) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.88395547) q[1];
sx q[1];
rz(-1.2246545) q[1];
sx q[1];
rz(-0.37948541) q[1];
rz(-1.0326321) q[3];
sx q[3];
rz(-1.0618883) q[3];
sx q[3];
rz(-0.63652388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9942921) q[2];
sx q[2];
rz(-1.9376126) q[2];
sx q[2];
rz(-3.1331983) q[2];
rz(0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10107772) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(-0.44152942) q[0];
rz(0.98835522) q[1];
sx q[1];
rz(-1.1343196) q[1];
sx q[1];
rz(0.13557869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6147185) q[0];
sx q[0];
rz(-1.563094) q[0];
sx q[0];
rz(3.1253424) q[0];
rz(1.1889815) q[2];
sx q[2];
rz(-2.2509607) q[2];
sx q[2];
rz(1.176468) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49040321) q[1];
sx q[1];
rz(-2.5156543) q[1];
sx q[1];
rz(2.2041026) q[1];
rz(-pi) q[2];
rz(0.83643416) q[3];
sx q[3];
rz(-1.7178917) q[3];
sx q[3];
rz(2.7874352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93379891) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(0.03820339) q[2];
rz(-0.52538747) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(-0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4811089) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(0.61087459) q[0];
rz(-1.8065709) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(-2.9023721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4668381) q[0];
sx q[0];
rz(-1.309899) q[0];
sx q[0];
rz(0.99715085) q[0];
x q[1];
rz(1.1785422) q[2];
sx q[2];
rz(-2.0118656) q[2];
sx q[2];
rz(0.46229306) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4799616) q[1];
sx q[1];
rz(-2.0212272) q[1];
sx q[1];
rz(1.2110787) q[1];
rz(1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(-3.0016921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7188344) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(3.139843) q[2];
rz(-2.8478029) q[3];
sx q[3];
rz(-1.9521451) q[3];
sx q[3];
rz(-0.26688117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9002429) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(0.84306651) q[0];
rz(2.8040366) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(-1.6036124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3030515) q[0];
sx q[0];
rz(-1.7619131) q[0];
sx q[0];
rz(2.3421846) q[0];
rz(-pi) q[1];
rz(-0.29422167) q[2];
sx q[2];
rz(-2.9595642) q[2];
sx q[2];
rz(2.1862669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2258721) q[1];
sx q[1];
rz(-2.8166526) q[1];
sx q[1];
rz(-2.48808) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4902275) q[3];
sx q[3];
rz(-1.3617097) q[3];
sx q[3];
rz(0.63864708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7097077) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(-1.0464) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-2.4304978) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043561291) q[0];
sx q[0];
rz(-1.2579608) q[0];
sx q[0];
rz(-0.64796722) q[0];
rz(-1.7421534) q[1];
sx q[1];
rz(-1.6439227) q[1];
sx q[1];
rz(-1.3290149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2752987) q[0];
sx q[0];
rz(-1.4614551) q[0];
sx q[0];
rz(1.7924395) q[0];
x q[1];
rz(-1.7850661) q[2];
sx q[2];
rz(-1.8545215) q[2];
sx q[2];
rz(1.8510712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44101366) q[1];
sx q[1];
rz(-2.3668336) q[1];
sx q[1];
rz(2.0225594) q[1];
rz(-pi) q[2];
rz(-0.51402199) q[3];
sx q[3];
rz(-1.8933305) q[3];
sx q[3];
rz(0.0016101282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7602188) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(-0.94183952) q[0];
rz(-2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9401682) q[0];
sx q[0];
rz(-1.3469101) q[0];
sx q[0];
rz(-3.1046449) q[0];
rz(-pi) q[1];
rz(1.0490369) q[2];
sx q[2];
rz(-1.6055487) q[2];
sx q[2];
rz(-3.1086189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6814179) q[1];
sx q[1];
rz(-0.81528403) q[1];
sx q[1];
rz(0.47487201) q[1];
rz(-0.022054733) q[3];
sx q[3];
rz(-1.0733177) q[3];
sx q[3];
rz(-0.22145222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3840702) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(0.43357098) q[2];
rz(-1.5844257) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(0.43237329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(-1.7973416) q[0];
rz(1.2035707) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(-1.3628091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575131) q[0];
sx q[0];
rz(-1.6005033) q[0];
sx q[0];
rz(-3.1196575) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3154738) q[2];
sx q[2];
rz(-1.6996133) q[2];
sx q[2];
rz(-1.4908189) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40920112) q[1];
sx q[1];
rz(-0.26136569) q[1];
sx q[1];
rz(0.50326668) q[1];
rz(2.8239488) q[3];
sx q[3];
rz(-2.2322725) q[3];
sx q[3];
rz(2.3931488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5726996) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(-0.9838689) q[2];
rz(-2.900506) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(2.4510395) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702061) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(2.8669226) q[0];
rz(1.341691) q[1];
sx q[1];
rz(-2.5119753) q[1];
sx q[1];
rz(-1.0460269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957841) q[0];
sx q[0];
rz(-1.2465451) q[0];
sx q[0];
rz(-0.62600531) q[0];
rz(-pi) q[1];
rz(-1.1860023) q[2];
sx q[2];
rz(-2.9287387) q[2];
sx q[2];
rz(1.720088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0053802) q[1];
sx q[1];
rz(-2.5544205) q[1];
sx q[1];
rz(-0.11759742) q[1];
rz(-pi) q[2];
rz(0.041954354) q[3];
sx q[3];
rz(-1.0681515) q[3];
sx q[3];
rz(-1.6142538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(2.7421303) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(-0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.28867662) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(2.5592819) q[0];
rz(2.9583926) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(-2.3351672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69003478) q[0];
sx q[0];
rz(-1.6584883) q[0];
sx q[0];
rz(-1.466485) q[0];
rz(0.33558947) q[2];
sx q[2];
rz(-1.8238471) q[2];
sx q[2];
rz(-2.0584681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2576221) q[1];
sx q[1];
rz(-2.3720354) q[1];
sx q[1];
rz(0.29410024) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6819477) q[3];
sx q[3];
rz(-1.2254997) q[3];
sx q[3];
rz(2.2636641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23239423) q[2];
sx q[2];
rz(-2.4269203) q[2];
sx q[2];
rz(-0.8052899) q[2];
rz(2.9946839) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(-2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.88937) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(-1.5421142) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-1.7003822) q[2];
sx q[2];
rz(-2.414188) q[2];
sx q[2];
rz(1.3368171) q[2];
rz(-0.44381683) q[3];
sx q[3];
rz(-2.4469821) q[3];
sx q[3];
rz(-1.2367005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
