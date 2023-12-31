OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(-2.864569) q[1];
sx q[1];
rz(-2.6695873) q[1];
sx q[1];
rz(-0.0013874887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479552) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(1.1975343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.087287993) q[2];
sx q[2];
rz(-0.44863551) q[2];
sx q[2];
rz(-2.0729614) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6359771) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(0.14640267) q[1];
x q[2];
rz(0.9790768) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(-1.0867659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(-2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(0.97066561) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(0.81545365) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.095123) q[0];
sx q[0];
rz(-0.9233343) q[0];
sx q[0];
rz(-1.4955273) q[0];
rz(-pi) q[1];
rz(1.843812) q[2];
sx q[2];
rz(-2.2815939) q[2];
sx q[2];
rz(-0.057991512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7533469) q[1];
sx q[1];
rz(-1.2699632) q[1];
sx q[1];
rz(1.5335598) q[1];
x q[2];
rz(-2.8498597) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(-0.093689703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8130428) q[0];
sx q[0];
rz(-1.2721491) q[0];
sx q[0];
rz(2.9850328) q[0];
rz(0.83061647) q[2];
sx q[2];
rz(-0.81095552) q[2];
sx q[2];
rz(-2.9142771) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42400186) q[1];
sx q[1];
rz(-2.058299) q[1];
sx q[1];
rz(1.8600149) q[1];
x q[2];
rz(-0.079654982) q[3];
sx q[3];
rz(-1.0632535) q[3];
sx q[3];
rz(1.5910651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6040566) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(2.5615454) q[2];
rz(-2.3245658) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(-0.91039175) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7461473) q[0];
sx q[0];
rz(-1.4822042) q[0];
sx q[0];
rz(0.079061411) q[0];
rz(2.0312112) q[2];
sx q[2];
rz(-1.5693671) q[2];
sx q[2];
rz(-1.9542076) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6064925) q[1];
sx q[1];
rz(-0.98787687) q[1];
sx q[1];
rz(2.5938354) q[1];
x q[2];
rz(1.7219909) q[3];
sx q[3];
rz(-0.70221838) q[3];
sx q[3];
rz(0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.794902) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.516974) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6277916) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(-0.4517201) q[0];
rz(-pi) q[1];
rz(-0.5337358) q[2];
sx q[2];
rz(-1.9118475) q[2];
sx q[2];
rz(1.5460207) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5114054) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(-1.1955111) q[1];
rz(0.77002854) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8020442) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(-2.8009159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5013803) q[0];
sx q[0];
rz(-3.0684154) q[0];
sx q[0];
rz(-1.120938) q[0];
rz(2.3855626) q[2];
sx q[2];
rz(-1.7871734) q[2];
sx q[2];
rz(-2.0331969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.016547116) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(-0.46334456) q[1];
x q[2];
rz(-2.2716899) q[3];
sx q[3];
rz(-2.5330336) q[3];
sx q[3];
rz(-2.9121947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(2.5781393) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(2.6766434) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15360399) q[0];
sx q[0];
rz(-2.0360332) q[0];
sx q[0];
rz(0.21501712) q[0];
rz(-pi) q[1];
rz(-0.99636997) q[2];
sx q[2];
rz(-1.9949706) q[2];
sx q[2];
rz(-0.043957274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3685776) q[1];
sx q[1];
rz(-1.3524719) q[1];
sx q[1];
rz(-1.6792084) q[1];
x q[2];
rz(0.025860272) q[3];
sx q[3];
rz(-1.5246632) q[3];
sx q[3];
rz(2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-0.31759343) q[2];
rz(0.57146227) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-3.0632339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0067622234) q[0];
sx q[0];
rz(-1.5015331) q[0];
sx q[0];
rz(-0.80471054) q[0];
rz(-1.3252844) q[2];
sx q[2];
rz(-1.3405372) q[2];
sx q[2];
rz(-0.85207176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2625418) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(-2.2294728) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3571635) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(-2.8182639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(2.890214) q[2];
rz(0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79779977) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-0.87402469) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6459991) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(-1.5828703) q[0];
rz(2.2888695) q[2];
sx q[2];
rz(-0.52041473) q[2];
sx q[2];
rz(2.4783217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6898432) q[1];
sx q[1];
rz(-0.50536957) q[1];
sx q[1];
rz(-1.7187353) q[1];
rz(0.47253982) q[3];
sx q[3];
rz(-0.66415411) q[3];
sx q[3];
rz(3.1411375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.7782036) q[3];
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
rz(-1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(-1.9627409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9261242) q[0];
sx q[0];
rz(-1.7755531) q[0];
sx q[0];
rz(-1.8490851) q[0];
rz(-pi) q[1];
rz(-0.12702282) q[2];
sx q[2];
rz(-1.4782895) q[2];
sx q[2];
rz(0.18735838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9465543) q[1];
sx q[1];
rz(-2.185501) q[1];
sx q[1];
rz(2.9123995) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6188649) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(-2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(2.2422092) q[2];
rz(0.87456885) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.62109229) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(-3.1065431) q[2];
sx q[2];
rz(-1.1309584) q[2];
sx q[2];
rz(-2.1396648) q[2];
rz(1.7189797) q[3];
sx q[3];
rz(-1.8437181) q[3];
sx q[3];
rz(-3.0317882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
