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
rz(-2.5042188) q[0];
sx q[0];
rz(-1.0567935) q[0];
sx q[0];
rz(-2.3699397) q[0];
rz(-0.15089384) q[1];
sx q[1];
rz(-0.3124736) q[1];
sx q[1];
rz(-1.6627275) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61687311) q[0];
sx q[0];
rz(-1.2390854) q[0];
sx q[0];
rz(1.2613368) q[0];
rz(-pi) q[1];
rz(-0.25881473) q[2];
sx q[2];
rz(-2.0740899) q[2];
sx q[2];
rz(0.81282114) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87541619) q[1];
sx q[1];
rz(-2.3485942) q[1];
sx q[1];
rz(2.7686053) q[1];
rz(-pi) q[2];
rz(-0.67920864) q[3];
sx q[3];
rz(-1.7082885) q[3];
sx q[3];
rz(-0.42543558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.60407698) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(0.83211952) q[2];
rz(0.65699792) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(2.7895797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81545365) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(0.64999181) q[0];
rz(-3.0187712) q[1];
sx q[1];
rz(-0.42418066) q[1];
sx q[1];
rz(2.4688683) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8433326) q[0];
sx q[0];
rz(-1.6125868) q[0];
sx q[0];
rz(0.01173794) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8459199) q[2];
sx q[2];
rz(-1.8725935) q[2];
sx q[2];
rz(1.9375305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0716463) q[1];
sx q[1];
rz(-0.59161797) q[1];
sx q[1];
rz(-1.2143192) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0510212) q[3];
sx q[3];
rz(-1.5477763) q[3];
sx q[3];
rz(2.1995403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1136721) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(-2.7552354) q[2];
rz(-1.1682642) q[3];
sx q[3];
rz(-1.2315653) q[3];
sx q[3];
rz(-2.0333576) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89762178) q[0];
sx q[0];
rz(-1.4742999) q[0];
sx q[0];
rz(0.058187159) q[0];
rz(2.8367786) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(0.85495368) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8182504) q[0];
sx q[0];
rz(-1.5691681) q[0];
sx q[0];
rz(3.073284) q[0];
rz(0.57010713) q[2];
sx q[2];
rz(-0.7749346) q[2];
sx q[2];
rz(2.3910445) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7185229) q[1];
sx q[1];
rz(-0.81941019) q[1];
sx q[1];
rz(1.2889483) q[1];
rz(-pi) q[2];
rz(-2.1257954) q[3];
sx q[3];
rz(-1.3449752) q[3];
sx q[3];
rz(-0.34103901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.052224) q[2];
sx q[2];
rz(-1.8171909) q[2];
sx q[2];
rz(-1.268187) q[2];
rz(0.44102272) q[3];
sx q[3];
rz(-2.5206168) q[3];
sx q[3];
rz(0.84180251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.12255254) q[0];
sx q[0];
rz(-0.95128107) q[0];
sx q[0];
rz(0.71504354) q[0];
rz(0.41636458) q[1];
sx q[1];
rz(-0.48960296) q[1];
sx q[1];
rz(0.29410902) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.490192) q[0];
sx q[0];
rz(-0.6985025) q[0];
sx q[0];
rz(1.6609037) q[0];
rz(-1.3662228) q[2];
sx q[2];
rz(-1.507221) q[2];
sx q[2];
rz(-0.83067375) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3336368) q[1];
sx q[1];
rz(-0.2273493) q[1];
sx q[1];
rz(0.61509404) q[1];
rz(-pi) q[2];
rz(-0.99331345) q[3];
sx q[3];
rz(-1.0619831) q[3];
sx q[3];
rz(-1.4636544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9010345) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(0.72875363) q[2];
rz(-0.48480222) q[3];
sx q[3];
rz(-1.0032283) q[3];
sx q[3];
rz(-2.9639967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6840376) q[0];
sx q[0];
rz(-1.3223248) q[0];
sx q[0];
rz(-0.27873248) q[0];
rz(0.067525603) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(-0.50594893) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099952817) q[0];
sx q[0];
rz(-2.1258265) q[0];
sx q[0];
rz(1.6567848) q[0];
rz(-pi) q[1];
rz(-1.672411) q[2];
sx q[2];
rz(-1.7676465) q[2];
sx q[2];
rz(-2.5972899) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5488779) q[1];
sx q[1];
rz(-2.1134106) q[1];
sx q[1];
rz(0.9701258) q[1];
rz(1.9002731) q[3];
sx q[3];
rz(-2.4412324) q[3];
sx q[3];
rz(2.2735689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2013596) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(3.039321) q[2];
rz(-0.30397948) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(2.6612018) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893696) q[0];
sx q[0];
rz(-2.3134573) q[0];
sx q[0];
rz(2.4612259) q[0];
rz(2.1807189) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(-0.56112498) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38491098) q[0];
sx q[0];
rz(-1.2088684) q[0];
sx q[0];
rz(0.2189526) q[0];
x q[1];
rz(-2.9343453) q[2];
sx q[2];
rz(-2.2842513) q[2];
sx q[2];
rz(-1.1175454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0099854) q[1];
sx q[1];
rz(-1.6254043) q[1];
sx q[1];
rz(-2.8427678) q[1];
x q[2];
rz(2.0797727) q[3];
sx q[3];
rz(-2.1474995) q[3];
sx q[3];
rz(-2.1901166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9095824) q[2];
sx q[2];
rz(-0.97753111) q[2];
sx q[2];
rz(1.5634465) q[2];
rz(2.1252508) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(-1.1004755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4119754) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(0.49402657) q[0];
rz(1.5644851) q[1];
sx q[1];
rz(-2.1421075) q[1];
sx q[1];
rz(3.0513501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3615389) q[0];
sx q[0];
rz(-2.001273) q[0];
sx q[0];
rz(1.3132877) q[0];
rz(-pi) q[1];
rz(0.97930564) q[2];
sx q[2];
rz(-2.3537209) q[2];
sx q[2];
rz(-2.5798714) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7612865) q[1];
sx q[1];
rz(-0.50380987) q[1];
sx q[1];
rz(-2.1977717) q[1];
rz(-pi) q[2];
rz(-0.97848864) q[3];
sx q[3];
rz(-2.8045125) q[3];
sx q[3];
rz(0.85142985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6287441) q[2];
sx q[2];
rz(-0.57400846) q[2];
sx q[2];
rz(-0.64344704) q[2];
rz(-2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(-3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028320463) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(-1.9182308) q[0];
rz(0.89448294) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(1.9953413) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47726163) q[0];
sx q[0];
rz(-2.1210175) q[0];
sx q[0];
rz(-2.6653637) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38739631) q[2];
sx q[2];
rz(-1.9700389) q[2];
sx q[2];
rz(2.7170167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97630298) q[1];
sx q[1];
rz(-2.7600651) q[1];
sx q[1];
rz(-0.69451992) q[1];
rz(1.9548863) q[3];
sx q[3];
rz(-0.62058228) q[3];
sx q[3];
rz(-0.24031249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.869183) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(1.4746846) q[2];
rz(-2.9278751) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(-0.86764446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2992582) q[0];
sx q[0];
rz(-0.61115757) q[0];
sx q[0];
rz(-2.4364731) q[0];
rz(0.57535386) q[1];
sx q[1];
rz(-1.7333142) q[1];
sx q[1];
rz(-0.56952482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3045159) q[0];
sx q[0];
rz(-1.6479371) q[0];
sx q[0];
rz(1.4894372) q[0];
x q[1];
rz(-0.51409431) q[2];
sx q[2];
rz(-0.55214685) q[2];
sx q[2];
rz(-2.8469126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25334206) q[1];
sx q[1];
rz(-0.98501316) q[1];
sx q[1];
rz(2.3164877) q[1];
rz(1.2678746) q[3];
sx q[3];
rz(-1.1008796) q[3];
sx q[3];
rz(2.5055714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9465955) q[2];
sx q[2];
rz(-0.94299287) q[2];
sx q[2];
rz(-0.97879624) q[2];
rz(-0.43092522) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8579213) q[0];
sx q[0];
rz(-0.019275276) q[0];
sx q[0];
rz(-1.4625782) q[0];
rz(-2.1025533) q[1];
sx q[1];
rz(-1.7580527) q[1];
sx q[1];
rz(-1.9336112) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53762728) q[0];
sx q[0];
rz(-2.6485291) q[0];
sx q[0];
rz(0.8925256) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0459837) q[2];
sx q[2];
rz(-1.5417347) q[2];
sx q[2];
rz(0.20653221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17716781) q[1];
sx q[1];
rz(-2.1806114) q[1];
sx q[1];
rz(-0.17788203) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.966217) q[3];
sx q[3];
rz(-1.5081769) q[3];
sx q[3];
rz(-2.5015802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.155948) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(1.8211625) q[2];
rz(-2.7909347) q[3];
sx q[3];
rz(-2.3242293) q[3];
sx q[3];
rz(-2.5613274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.25528) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(-1.4844004) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(-1.6419353) q[2];
sx q[2];
rz(-1.1924469) q[2];
sx q[2];
rz(3.057657) q[2];
rz(-1.2091533) q[3];
sx q[3];
rz(-2.4222838) q[3];
sx q[3];
rz(-2.1771976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
