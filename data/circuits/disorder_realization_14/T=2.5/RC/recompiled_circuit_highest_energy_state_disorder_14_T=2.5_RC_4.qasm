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
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(-0.60751539) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54017055) q[0];
sx q[0];
rz(-2.5937732) q[0];
sx q[0];
rz(-2.026985) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11694853) q[2];
sx q[2];
rz(-2.0772572) q[2];
sx q[2];
rz(0.036265515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40558896) q[1];
sx q[1];
rz(-1.0936671) q[1];
sx q[1];
rz(-1.7504343) q[1];
rz(-pi) q[2];
rz(2.257454) q[3];
sx q[3];
rz(-0.55301412) q[3];
sx q[3];
rz(-0.64698863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6174378) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(-0.080862008) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8993768) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(0.85897613) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(1.7346409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039120395) q[0];
sx q[0];
rz(-0.539383) q[0];
sx q[0];
rz(-1.1155737) q[0];
rz(-1.2386049) q[2];
sx q[2];
rz(-1.536035) q[2];
sx q[2];
rz(1.8046276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3000112) q[1];
sx q[1];
rz(-1.4805796) q[1];
sx q[1];
rz(-1.0567699) q[1];
rz(-2.5646016) q[3];
sx q[3];
rz(-1.9245879) q[3];
sx q[3];
rz(-0.64567034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0824288) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.814369) q[2];
rz(-2.0969157) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(-0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5752207) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(0.4253934) q[0];
rz(-1.3771903) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(-1.707071) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67514133) q[0];
sx q[0];
rz(-1.4527307) q[0];
sx q[0];
rz(-0.64339126) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6041669) q[2];
sx q[2];
rz(-2.2788308) q[2];
sx q[2];
rz(-1.58085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19615281) q[1];
sx q[1];
rz(-1.1588133) q[1];
sx q[1];
rz(-0.71228551) q[1];
rz(-pi) q[2];
rz(-1.1528963) q[3];
sx q[3];
rz(-2.4800081) q[3];
sx q[3];
rz(1.0583056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.2563027) q[2];
sx q[2];
rz(1.2909935) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212379) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(0.77019101) q[0];
rz(-2.1417446) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(1.3410478) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59113778) q[0];
sx q[0];
rz(-1.3584175) q[0];
sx q[0];
rz(0.55340931) q[0];
x q[1];
rz(2.6821939) q[2];
sx q[2];
rz(-0.93743491) q[2];
sx q[2];
rz(-2.6448696) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8902258) q[1];
sx q[1];
rz(-0.40491762) q[1];
sx q[1];
rz(0.7271073) q[1];
rz(-pi) q[2];
rz(-3.0422736) q[3];
sx q[3];
rz(-1.2867974) q[3];
sx q[3];
rz(-1.2885531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8009214) q[2];
sx q[2];
rz(-2.0719216) q[2];
sx q[2];
rz(2.3731903) q[2];
rz(-3.0411804) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95555821) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(-0.22698639) q[0];
rz(-1.3817894) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(-1.3963799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5433301) q[0];
sx q[0];
rz(-0.51307438) q[0];
sx q[0];
rz(2.3725879) q[0];
x q[1];
rz(-0.41047217) q[2];
sx q[2];
rz(-1.4183525) q[2];
sx q[2];
rz(1.6688011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5310881) q[1];
sx q[1];
rz(-0.63734326) q[1];
sx q[1];
rz(3.0211041) q[1];
x q[2];
rz(1.5377858) q[3];
sx q[3];
rz(-0.74652687) q[3];
sx q[3];
rz(1.7611662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9153626) q[2];
sx q[2];
rz(-1.1504983) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(-1.3876312) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(-1.7414198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48080322) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(-1.8060818) q[0];
rz(2.238359) q[1];
sx q[1];
rz(-1.3390373) q[1];
sx q[1];
rz(-2.9551771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6064998) q[0];
sx q[0];
rz(-0.71487521) q[0];
sx q[0];
rz(-0.8755279) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.769563) q[2];
sx q[2];
rz(-1.0093401) q[2];
sx q[2];
rz(1.1793062) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.820436) q[1];
sx q[1];
rz(-1.8119436) q[1];
sx q[1];
rz(-0.92280573) q[1];
rz(2.3535054) q[3];
sx q[3];
rz(-0.65207802) q[3];
sx q[3];
rz(-1.2913845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67941252) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(1.3235486) q[2];
rz(1.1421674) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(-2.2511258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46397504) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(-0.63419813) q[0];
rz(-1.0448666) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(2.1814836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38518181) q[0];
sx q[0];
rz(-0.39968458) q[0];
sx q[0];
rz(-1.6413641) q[0];
rz(-0.30107408) q[2];
sx q[2];
rz(-1.5784401) q[2];
sx q[2];
rz(-0.49087015) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83852488) q[1];
sx q[1];
rz(-0.52553287) q[1];
sx q[1];
rz(-1.5938894) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3784901) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(-2.4214937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1977957) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(-1.5459527) q[2];
rz(-1.724285) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(-1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6087795) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(3.0335562) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(-1.1688165) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790137) q[0];
sx q[0];
rz(-0.65555182) q[0];
sx q[0];
rz(2.2961839) q[0];
rz(-pi) q[1];
rz(-3.1415784) q[2];
sx q[2];
rz(-0.76185267) q[2];
sx q[2];
rz(-0.95214168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9795684) q[1];
sx q[1];
rz(-1.8037533) q[1];
sx q[1];
rz(-0.090090171) q[1];
x q[2];
rz(-1.5511369) q[3];
sx q[3];
rz(-1.0791313) q[3];
sx q[3];
rz(-1.2439072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10890266) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(-0.6558134) q[2];
rz(-3.0374895) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(-2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.6498097) q[0];
rz(-2.1022294) q[1];
sx q[1];
rz(-0.84016687) q[1];
sx q[1];
rz(2.6731491) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2830848) q[0];
sx q[0];
rz(-1.7584086) q[0];
sx q[0];
rz(-1.5664738) q[0];
x q[1];
rz(2.5889187) q[2];
sx q[2];
rz(-2.6330593) q[2];
sx q[2];
rz(0.39297152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2704308) q[1];
sx q[1];
rz(-2.3498145) q[1];
sx q[1];
rz(1.7986078) q[1];
rz(-pi) q[2];
rz(-1.8519206) q[3];
sx q[3];
rz(-1.2996965) q[3];
sx q[3];
rz(-2.2106314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.046772) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(-2.2136733) q[2];
rz(-1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(-1.4891589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0874262) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(1.077865) q[0];
rz(0.36733356) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7736762) q[0];
sx q[0];
rz(-1.8225637) q[0];
sx q[0];
rz(-1.0836224) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7214891) q[2];
sx q[2];
rz(-2.4825399) q[2];
sx q[2];
rz(-0.53021741) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.16691567) q[1];
sx q[1];
rz(-1.484718) q[1];
sx q[1];
rz(-2.9831216) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.359109) q[3];
sx q[3];
rz(-0.84722391) q[3];
sx q[3];
rz(0.43989691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(2.7122811) q[2];
rz(2.9863827) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(-1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406381) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(0.86391972) q[1];
sx q[1];
rz(-2.7680631) q[1];
sx q[1];
rz(1.6815129) q[1];
rz(0.39045329) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(2.024509) q[3];
sx q[3];
rz(-0.67839218) q[3];
sx q[3];
rz(-1.4339504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
