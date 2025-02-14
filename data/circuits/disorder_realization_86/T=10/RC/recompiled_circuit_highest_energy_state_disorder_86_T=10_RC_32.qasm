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
rz(1.9011693) q[0];
sx q[0];
rz(-1.9440396) q[0];
sx q[0];
rz(-2.9252606) q[0];
rz(4.0989838) q[1];
sx q[1];
rz(5.6070072) q[1];
sx q[1];
rz(11.054872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0656242) q[0];
sx q[0];
rz(-1.0145717) q[0];
sx q[0];
rz(1.5821304) q[0];
rz(1.897282) q[2];
sx q[2];
rz(-1.1190345) q[2];
sx q[2];
rz(1.2690074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.88684593) q[1];
sx q[1];
rz(-1.0751845) q[1];
sx q[1];
rz(-1.9504471) q[1];
rz(2.7483398) q[3];
sx q[3];
rz(-2.234708) q[3];
sx q[3];
rz(2.8877986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.034885255) q[2];
sx q[2];
rz(-2.7896176) q[2];
sx q[2];
rz(-1.0484288) q[2];
rz(0.18167051) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(0.65339965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.492391) q[0];
sx q[0];
rz(-2.1972456) q[0];
sx q[0];
rz(-2.6948068) q[0];
rz(-0.88042879) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(-2.3588038) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.114417) q[0];
sx q[0];
rz(-1.3265298) q[0];
sx q[0];
rz(-0.063112325) q[0];
rz(-pi) q[1];
rz(-2.5249285) q[2];
sx q[2];
rz(-1.2149802) q[2];
sx q[2];
rz(0.43105506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70322733) q[1];
sx q[1];
rz(-1.2709054) q[1];
sx q[1];
rz(2.6542526) q[1];
rz(-pi) q[2];
rz(-2.2530678) q[3];
sx q[3];
rz(-1.1924713) q[3];
sx q[3];
rz(-0.64680566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.11144) q[2];
sx q[2];
rz(-1.6609265) q[2];
sx q[2];
rz(2.1622369) q[2];
rz(-0.3848981) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(2.9773007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48591831) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(-0.78980494) q[0];
rz(2.9549331) q[1];
sx q[1];
rz(-1.4013545) q[1];
sx q[1];
rz(-0.99367118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8756722) q[0];
sx q[0];
rz(-1.3215995) q[0];
sx q[0];
rz(-0.56206352) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7970071) q[2];
sx q[2];
rz(-2.1717584) q[2];
sx q[2];
rz(0.41659875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7977596) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(0.47618687) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8348654) q[3];
sx q[3];
rz(-2.3969458) q[3];
sx q[3];
rz(-2.4414203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32187244) q[2];
sx q[2];
rz(-2.0237782) q[2];
sx q[2];
rz(2.7703088) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45469859) q[0];
sx q[0];
rz(-1.0201539) q[0];
sx q[0];
rz(1.4042847) q[0];
rz(0.67717254) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(-1.6489395) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7248047) q[0];
sx q[0];
rz(-1.8128403) q[0];
sx q[0];
rz(1.7927732) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8437496) q[2];
sx q[2];
rz(-1.1561511) q[2];
sx q[2];
rz(-2.9306102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4863534) q[1];
sx q[1];
rz(-2.7905373) q[1];
sx q[1];
rz(2.4454861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.221232) q[3];
sx q[3];
rz(-2.3858983) q[3];
sx q[3];
rz(-0.44103482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45381418) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(2.7933534) q[2];
rz(1.6866775) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(1.3454364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763879) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(-2.2494466) q[0];
rz(-0.46420321) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(1.2947327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034445914) q[0];
sx q[0];
rz(-2.4305516) q[0];
sx q[0];
rz(-0.76816316) q[0];
rz(-pi) q[1];
rz(-2.6038997) q[2];
sx q[2];
rz(-1.6811996) q[2];
sx q[2];
rz(-2.3619719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4720104) q[1];
sx q[1];
rz(-1.246844) q[1];
sx q[1];
rz(2.0355909) q[1];
x q[2];
rz(-0.61035778) q[3];
sx q[3];
rz(-0.98516432) q[3];
sx q[3];
rz(0.40025362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8190454) q[2];
sx q[2];
rz(-2.5307541) q[2];
sx q[2];
rz(0.56582212) q[2];
rz(3.0602509) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(-2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.464798) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(-1.5860522) q[0];
rz(-2.0809035) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(-0.63180822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9264797) q[0];
sx q[0];
rz(-0.2479015) q[0];
sx q[0];
rz(2.8126024) q[0];
rz(-2.1366227) q[2];
sx q[2];
rz(-0.7938876) q[2];
sx q[2];
rz(0.53295202) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0785219) q[1];
sx q[1];
rz(-2.1025189) q[1];
sx q[1];
rz(0.39549455) q[1];
rz(0.17978823) q[3];
sx q[3];
rz(-0.69151141) q[3];
sx q[3];
rz(0.19516842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1598728) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(-2.6427606) q[2];
rz(-1.8286797) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(-1.4341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53443921) q[0];
sx q[0];
rz(-2.7767015) q[0];
sx q[0];
rz(-0.69951192) q[0];
rz(-2.6761159) q[1];
sx q[1];
rz(-2.270348) q[1];
sx q[1];
rz(0.93719283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95109601) q[0];
sx q[0];
rz(-1.7825923) q[0];
sx q[0];
rz(1.7429211) q[0];
x q[1];
rz(1.466339) q[2];
sx q[2];
rz(-1.3323931) q[2];
sx q[2];
rz(1.6676211) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7286789) q[1];
sx q[1];
rz(-1.2866486) q[1];
sx q[1];
rz(-2.9212664) q[1];
x q[2];
rz(-2.8392467) q[3];
sx q[3];
rz(-1.8519326) q[3];
sx q[3];
rz(-0.1880364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.297544) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(0.37332264) q[2];
rz(-1.1897872) q[3];
sx q[3];
rz(-0.50896421) q[3];
sx q[3];
rz(2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.3968286) q[0];
sx q[0];
rz(-2.0592392) q[0];
sx q[0];
rz(-0.74791351) q[0];
rz(-0.76639908) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(3.1386197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.317694) q[0];
sx q[0];
rz(-0.71461535) q[0];
sx q[0];
rz(-1.2750285) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1220436) q[2];
sx q[2];
rz(-1.4998933) q[2];
sx q[2];
rz(-0.30724684) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0743979) q[1];
sx q[1];
rz(-1.41681) q[1];
sx q[1];
rz(-0.043937307) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19488867) q[3];
sx q[3];
rz(-1.3759383) q[3];
sx q[3];
rz(0.080435924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11091867) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(1.0670916) q[2];
rz(-0.083960697) q[3];
sx q[3];
rz(-2.6559918) q[3];
sx q[3];
rz(-0.77897227) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79754168) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(3.0352266) q[0];
rz(2.1616409) q[1];
sx q[1];
rz(-1.4915024) q[1];
sx q[1];
rz(-2.3756557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2666616) q[0];
sx q[0];
rz(-2.4342172) q[0];
sx q[0];
rz(2.4445663) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6607051) q[2];
sx q[2];
rz(-1.769763) q[2];
sx q[2];
rz(1.0294017) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3801743) q[1];
sx q[1];
rz(-1.5360218) q[1];
sx q[1];
rz(0.67616391) q[1];
x q[2];
rz(-2.8838653) q[3];
sx q[3];
rz(-3.1270087) q[3];
sx q[3];
rz(0.38250438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.47436675) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-0.70029798) q[2];
rz(-2.4750366) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67548442) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(-1.3960557) q[0];
rz(-1.4736157) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(2.4748763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6121126) q[0];
sx q[0];
rz(-2.1938938) q[0];
sx q[0];
rz(-2.3701131) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6595419) q[2];
sx q[2];
rz(-1.7983984) q[2];
sx q[2];
rz(-1.2510117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.279505) q[1];
sx q[1];
rz(-0.34395978) q[1];
sx q[1];
rz(0.39744795) q[1];
x q[2];
rz(-2.5693043) q[3];
sx q[3];
rz(-1.0977355) q[3];
sx q[3];
rz(-2.8561887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.066976808) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(-2.2508049) q[2];
rz(0.43205076) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(0.74705684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73410949) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(-1.2474077) q[1];
sx q[1];
rz(-2.4220962) q[1];
sx q[1];
rz(-1.9069506) q[1];
rz(2.9802889) q[2];
sx q[2];
rz(-1.2304753) q[2];
sx q[2];
rz(-0.58438042) q[2];
rz(-0.4890781) q[3];
sx q[3];
rz(-1.6001971) q[3];
sx q[3];
rz(-1.8501545) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
