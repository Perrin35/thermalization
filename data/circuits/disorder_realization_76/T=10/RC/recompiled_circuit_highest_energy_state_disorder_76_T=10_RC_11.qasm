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
rz(-0.79987502) q[0];
sx q[0];
rz(2.3246111) q[0];
sx q[0];
rz(9.8819879) q[0];
rz(0.41181052) q[1];
sx q[1];
rz(4.7635912) q[1];
sx q[1];
rz(6.8430321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7871817) q[0];
sx q[0];
rz(-2.5776064) q[0];
sx q[0];
rz(-2.5111879) q[0];
rz(-pi) q[1];
rz(-0.28607762) q[2];
sx q[2];
rz(-2.0564579) q[2];
sx q[2];
rz(2.3229675) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7840883) q[1];
sx q[1];
rz(-0.58114806) q[1];
sx q[1];
rz(-0.52304348) q[1];
rz(1.4234957) q[3];
sx q[3];
rz(-1.788874) q[3];
sx q[3];
rz(-2.9678621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7655699) q[2];
sx q[2];
rz(-2.1799808) q[2];
sx q[2];
rz(1.2961071) q[2];
rz(2.5206595) q[3];
sx q[3];
rz(-0.66578484) q[3];
sx q[3];
rz(-2.2165829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45944443) q[0];
sx q[0];
rz(-2.5542673) q[0];
sx q[0];
rz(-2.4272954) q[0];
rz(-0.722305) q[1];
sx q[1];
rz(-2.0562833) q[1];
sx q[1];
rz(1.3246271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9324397) q[0];
sx q[0];
rz(-0.55190933) q[0];
sx q[0];
rz(-2.1789684) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73280735) q[2];
sx q[2];
rz(-1.6116953) q[2];
sx q[2];
rz(-0.78860229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5058742) q[1];
sx q[1];
rz(-1.7298475) q[1];
sx q[1];
rz(-0.37271865) q[1];
rz(-pi) q[2];
rz(0.88543798) q[3];
sx q[3];
rz(-1.4760222) q[3];
sx q[3];
rz(1.4023905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.033919949) q[2];
sx q[2];
rz(-1.4227957) q[2];
sx q[2];
rz(2.1495492) q[2];
rz(-0.77573675) q[3];
sx q[3];
rz(-0.78191596) q[3];
sx q[3];
rz(1.0824664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-2.7139605) q[0];
sx q[0];
rz(-1.7949224) q[0];
sx q[0];
rz(2.2295075) q[0];
rz(2.7569547) q[1];
sx q[1];
rz(-2.1802528) q[1];
sx q[1];
rz(2.6461163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3289483) q[0];
sx q[0];
rz(-1.5248796) q[0];
sx q[0];
rz(-1.4981235) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3922973) q[2];
sx q[2];
rz(-1.4297419) q[2];
sx q[2];
rz(2.2749449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7621618) q[1];
sx q[1];
rz(-0.84745211) q[1];
sx q[1];
rz(1.5033733) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1285237) q[3];
sx q[3];
rz(-0.94892293) q[3];
sx q[3];
rz(-2.6789846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2428525) q[2];
sx q[2];
rz(-2.0169368) q[2];
sx q[2];
rz(0.43453547) q[2];
rz(0.09856002) q[3];
sx q[3];
rz(-1.4222654) q[3];
sx q[3];
rz(0.77384531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74715215) q[0];
sx q[0];
rz(-0.23500615) q[0];
sx q[0];
rz(2.8727942) q[0];
rz(2.0647743) q[1];
sx q[1];
rz(-1.7348758) q[1];
sx q[1];
rz(-1.9487618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4556325) q[0];
sx q[0];
rz(-1.7436899) q[0];
sx q[0];
rz(-0.25211035) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7014163) q[2];
sx q[2];
rz(-1.4358178) q[2];
sx q[2];
rz(-2.2041952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7320648) q[1];
sx q[1];
rz(-1.3573779) q[1];
sx q[1];
rz(1.7312538) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4402839) q[3];
sx q[3];
rz(-1.4150672) q[3];
sx q[3];
rz(2.6913672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8863135) q[2];
sx q[2];
rz(-2.2463319) q[2];
sx q[2];
rz(-2.8423584) q[2];
rz(-0.29087654) q[3];
sx q[3];
rz(-1.2521005) q[3];
sx q[3];
rz(-2.4860184) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84151477) q[0];
sx q[0];
rz(-0.97989196) q[0];
sx q[0];
rz(-0.078068659) q[0];
rz(0.95371753) q[1];
sx q[1];
rz(-2.6225312) q[1];
sx q[1];
rz(-1.5740707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8918607) q[0];
sx q[0];
rz(-1.390269) q[0];
sx q[0];
rz(-2.2523227) q[0];
rz(-pi) q[1];
rz(-0.40549739) q[2];
sx q[2];
rz(-0.73842305) q[2];
sx q[2];
rz(-2.9905295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5513796) q[1];
sx q[1];
rz(-2.1825983) q[1];
sx q[1];
rz(-3.1019866) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0573274) q[3];
sx q[3];
rz(-2.4467496) q[3];
sx q[3];
rz(-0.39277276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2786402) q[2];
sx q[2];
rz(-2.7095257) q[2];
sx q[2];
rz(0.26818177) q[2];
rz(-2.6523318) q[3];
sx q[3];
rz(-2.0475976) q[3];
sx q[3];
rz(1.4485654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.0972524) q[0];
sx q[0];
rz(-3.1179929) q[0];
sx q[0];
rz(-2.6771255) q[0];
rz(2.6181472) q[1];
sx q[1];
rz(-2.372066) q[1];
sx q[1];
rz(-0.73062599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084186144) q[0];
sx q[0];
rz(-0.75907403) q[0];
sx q[0];
rz(-1.7763863) q[0];
rz(-pi) q[1];
rz(-0.36412698) q[2];
sx q[2];
rz(-2.4362323) q[2];
sx q[2];
rz(-0.10157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82416049) q[1];
sx q[1];
rz(-0.80866122) q[1];
sx q[1];
rz(2.5416994) q[1];
rz(2.6650362) q[3];
sx q[3];
rz(-2.8598197) q[3];
sx q[3];
rz(1.2721541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0940493) q[2];
sx q[2];
rz(-1.0796248) q[2];
sx q[2];
rz(0.13761061) q[2];
rz(2.008647) q[3];
sx q[3];
rz(-1.1762985) q[3];
sx q[3];
rz(1.2631811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9889744) q[0];
sx q[0];
rz(-1.512383) q[0];
sx q[0];
rz(-3.1076987) q[0];
rz(-0.35941091) q[1];
sx q[1];
rz(-1.5243328) q[1];
sx q[1];
rz(2.3072037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2330025) q[0];
sx q[0];
rz(-2.6437573) q[0];
sx q[0];
rz(-0.69186886) q[0];
rz(2.4530386) q[2];
sx q[2];
rz(-1.1821185) q[2];
sx q[2];
rz(-1.2297024) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31454489) q[1];
sx q[1];
rz(-0.27092182) q[1];
sx q[1];
rz(-1.2071868) q[1];
rz(-pi) q[2];
rz(0.078975695) q[3];
sx q[3];
rz(-0.67310909) q[3];
sx q[3];
rz(-2.1977294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72839165) q[2];
sx q[2];
rz(-2.8381556) q[2];
sx q[2];
rz(2.0835853) q[2];
rz(0.47438619) q[3];
sx q[3];
rz(-2.3460903) q[3];
sx q[3];
rz(1.0559731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0742842) q[0];
sx q[0];
rz(-1.5165167) q[0];
sx q[0];
rz(-1.3577331) q[0];
rz(2.7108497) q[1];
sx q[1];
rz(-1.4746702) q[1];
sx q[1];
rz(0.46708435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8170094) q[0];
sx q[0];
rz(-1.6139493) q[0];
sx q[0];
rz(0.0014716455) q[0];
rz(-0.22612382) q[2];
sx q[2];
rz(-2.4259822) q[2];
sx q[2];
rz(-2.0307045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43710762) q[1];
sx q[1];
rz(-0.52981716) q[1];
sx q[1];
rz(0.036661224) q[1];
rz(-2.2086618) q[3];
sx q[3];
rz(-1.6137505) q[3];
sx q[3];
rz(-1.651498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0507386) q[2];
sx q[2];
rz(-0.41768062) q[2];
sx q[2];
rz(-2.0609071) q[2];
rz(-2.4596227) q[3];
sx q[3];
rz(-2.3365648) q[3];
sx q[3];
rz(2.4431156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4526116) q[0];
sx q[0];
rz(-0.51849759) q[0];
sx q[0];
rz(-2.3338351) q[0];
rz(-1.0001596) q[1];
sx q[1];
rz(-0.88614416) q[1];
sx q[1];
rz(2.3449786) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0611872) q[0];
sx q[0];
rz(-0.10594254) q[0];
sx q[0];
rz(-1.1619912) q[0];
x q[1];
rz(-2.8136425) q[2];
sx q[2];
rz(-0.65153507) q[2];
sx q[2];
rz(-0.33199379) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1353242) q[1];
sx q[1];
rz(-0.43244967) q[1];
sx q[1];
rz(-1.2342374) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2563989) q[3];
sx q[3];
rz(-1.8125497) q[3];
sx q[3];
rz(0.35255656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5920068) q[2];
sx q[2];
rz(-1.0909811) q[2];
sx q[2];
rz(0.31527147) q[2];
rz(-1.573805) q[3];
sx q[3];
rz(-1.1476293) q[3];
sx q[3];
rz(-1.1228336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0386117) q[0];
sx q[0];
rz(-2.4760315) q[0];
sx q[0];
rz(-0.82157201) q[0];
rz(0.483825) q[1];
sx q[1];
rz(-1.7211823) q[1];
sx q[1];
rz(2.9169567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0517462) q[0];
sx q[0];
rz(-1.7210074) q[0];
sx q[0];
rz(1.8972904) q[0];
rz(-pi) q[1];
rz(-0.63172747) q[2];
sx q[2];
rz(-1.016482) q[2];
sx q[2];
rz(2.4265223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1146176) q[1];
sx q[1];
rz(-1.3624422) q[1];
sx q[1];
rz(0.88233739) q[1];
x q[2];
rz(3.0167019) q[3];
sx q[3];
rz(-1.7833059) q[3];
sx q[3];
rz(-2.8954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3839174) q[2];
sx q[2];
rz(-3.0249247) q[2];
sx q[2];
rz(-0.095001027) q[2];
rz(-1.4875686) q[3];
sx q[3];
rz(-2.5520958) q[3];
sx q[3];
rz(0.76475638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8944396) q[0];
sx q[0];
rz(-2.7338487) q[0];
sx q[0];
rz(-0.26640531) q[0];
rz(15/(8*pi)) q[1];
sx q[1];
rz(-1.6886371) q[1];
sx q[1];
rz(1.4773038) q[1];
rz(-2.785037) q[2];
sx q[2];
rz(-0.92477634) q[2];
sx q[2];
rz(-1.7348031) q[2];
rz(1.568902) q[3];
sx q[3];
rz(-1.1952778) q[3];
sx q[3];
rz(0.3921685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
