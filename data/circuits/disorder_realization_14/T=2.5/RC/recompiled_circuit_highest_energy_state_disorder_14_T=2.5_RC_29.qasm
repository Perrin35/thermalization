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
rz(-0.7402339) q[0];
sx q[0];
rz(4.8010173) q[0];
sx q[0];
rz(12.231449) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(-1.0022751) q[1];
sx q[1];
rz(0.60751539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6014221) q[0];
sx q[0];
rz(-2.5937732) q[0];
sx q[0];
rz(-1.1146077) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0801699) q[2];
sx q[2];
rz(-1.6730089) q[2];
sx q[2];
rz(1.6639903) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40558896) q[1];
sx q[1];
rz(-2.0479255) q[1];
sx q[1];
rz(-1.7504343) q[1];
rz(-2.7685952) q[3];
sx q[3];
rz(-1.9891051) q[3];
sx q[3];
rz(-1.4137063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6174378) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(2.4784135) q[2];
rz(-0.080862008) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(-1.1801571) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2422159) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(2.2826165) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(1.7346409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55757421) q[0];
sx q[0];
rz(-2.0502591) q[0];
sx q[0];
rz(2.8842501) q[0];
rz(-1.9029877) q[2];
sx q[2];
rz(-1.6055577) q[2];
sx q[2];
rz(1.8046276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3000112) q[1];
sx q[1];
rz(-1.661013) q[1];
sx q[1];
rz(1.0567699) q[1];
rz(-pi) q[2];
rz(-0.57699109) q[3];
sx q[3];
rz(-1.2170047) q[3];
sx q[3];
rz(-0.64567034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0591639) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.814369) q[2];
rz(-2.0969157) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(2.7248342) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(-0.4253934) q[0];
rz(-1.3771903) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(-1.4345217) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4664513) q[0];
sx q[0];
rz(-1.688862) q[0];
sx q[0];
rz(-2.4982014) q[0];
rz(-pi) q[1];
x q[1];
rz(0.038952053) q[2];
sx q[2];
rz(-2.4329081) q[2];
sx q[2];
rz(-1.612029) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9454398) q[1];
sx q[1];
rz(-1.1588133) q[1];
sx q[1];
rz(-0.71228551) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1528963) q[3];
sx q[3];
rz(-0.66158453) q[3];
sx q[3];
rz(-2.083287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7003358) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.4972081) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62035471) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(0.77019101) q[0];
rz(0.99984804) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(-1.3410478) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59113778) q[0];
sx q[0];
rz(-1.3584175) q[0];
sx q[0];
rz(-0.55340931) q[0];
rz(1.0275339) q[2];
sx q[2];
rz(-2.3781666) q[2];
sx q[2];
rz(1.1929407) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.25136687) q[1];
sx q[1];
rz(-0.40491762) q[1];
sx q[1];
rz(-2.4144854) q[1];
x q[2];
rz(0.099319066) q[3];
sx q[3];
rz(-1.8547952) q[3];
sx q[3];
rz(-1.8530396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3406713) q[2];
sx q[2];
rz(-2.0719216) q[2];
sx q[2];
rz(-0.7684024) q[2];
rz(-0.10041222) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860344) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(-2.9146063) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(1.3963799) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5433301) q[0];
sx q[0];
rz(-2.6285183) q[0];
sx q[0];
rz(0.76900478) q[0];
rz(-2.7740741) q[2];
sx q[2];
rz(-0.43635363) q[2];
sx q[2];
rz(-2.9038716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13670838) q[1];
sx q[1];
rz(-1.6423823) q[1];
sx q[1];
rz(2.5077255) q[1];
rz(-pi) q[2];
rz(-1.6038069) q[3];
sx q[3];
rz(-0.74652687) q[3];
sx q[3];
rz(1.7611662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9153626) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(-0.076233141) q[2];
rz(1.7539615) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6607894) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(1.3355108) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(2.9551771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.37019) q[0];
sx q[0];
rz(-2.0982842) q[0];
sx q[0];
rz(-0.50748388) q[0];
rz(-pi) q[1];
rz(2.5711362) q[2];
sx q[2];
rz(-1.738731) q[2];
sx q[2];
rz(2.8569375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32115667) q[1];
sx q[1];
rz(-1.3296491) q[1];
sx q[1];
rz(0.92280573) q[1];
rz(-0.49390467) q[3];
sx q[3];
rz(-1.1260238) q[3];
sx q[3];
rz(-2.1879823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67941252) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(-1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6776176) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(2.5073945) q[0];
rz(-1.0448666) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(2.1814836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2506367) q[0];
sx q[0];
rz(-1.5433558) q[0];
sx q[0];
rz(1.1720042) q[0];
rz(-pi) q[1];
rz(2.8405186) q[2];
sx q[2];
rz(-1.5631526) q[2];
sx q[2];
rz(-2.6507225) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83852488) q[1];
sx q[1];
rz(-0.52553287) q[1];
sx q[1];
rz(-1.5477033) q[1];
rz(-0.74713196) q[3];
sx q[3];
rz(-2.225156) q[3];
sx q[3];
rz(-0.33734712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1977957) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(1.5459527) q[2];
rz(1.4173077) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5328131) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(-0.97484318) q[0];
rz(0.10803647) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(1.9727762) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0371712) q[0];
sx q[0];
rz(-1.9871431) q[0];
sx q[0];
rz(-1.0486616) q[0];
rz(-pi) q[1];
rz(-1.5708099) q[2];
sx q[2];
rz(-0.80894366) q[2];
sx q[2];
rz(2.1894313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9795684) q[1];
sx q[1];
rz(-1.3378394) q[1];
sx q[1];
rz(-3.0515025) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5904558) q[3];
sx q[3];
rz(-2.0624614) q[3];
sx q[3];
rz(1.2439072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.03269) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(-2.4857793) q[2];
rz(3.0374895) q[3];
sx q[3];
rz(-1.7963573) q[3];
sx q[3];
rz(-2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(1.6498097) q[0];
rz(2.1022294) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(-0.46844354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3062564) q[0];
sx q[0];
rz(-0.1876615) q[0];
sx q[0];
rz(-0.022764877) q[0];
rz(-1.8554815) q[2];
sx q[2];
rz(-1.9980901) q[2];
sx q[2];
rz(0.22186771) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13889986) q[1];
sx q[1];
rz(-1.7322092) q[1];
sx q[1];
rz(2.3494865) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28161093) q[3];
sx q[3];
rz(-1.8413895) q[3];
sx q[3];
rz(-2.5789346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.046772) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(2.2136733) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0874262) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(1.077865) q[0];
rz(-0.36733356) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7736762) q[0];
sx q[0];
rz(-1.8225637) q[0];
sx q[0];
rz(2.0579703) q[0];
x q[1];
rz(1.7214891) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(-2.6113752) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4176191) q[1];
sx q[1];
rz(-1.7286757) q[1];
sx q[1];
rz(-1.6579614) q[1];
rz(-pi) q[2];
rz(-0.73478847) q[3];
sx q[3];
rz(-1.4126724) q[3];
sx q[3];
rz(-1.2722335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(-2.7122811) q[2];
rz(-0.15520994) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(1.8163053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406381) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(-2.2776729) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(0.39045329) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(0.33959099) q[3];
sx q[3];
rz(-2.17008) q[3];
sx q[3];
rz(-1.9934987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
