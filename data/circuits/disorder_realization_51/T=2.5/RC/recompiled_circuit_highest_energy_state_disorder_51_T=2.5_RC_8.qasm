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
rz(2.3866374) q[0];
sx q[0];
rz(-0.75140262) q[0];
sx q[0];
rz(2.0928535) q[0];
rz(0.41930786) q[1];
sx q[1];
rz(-1.3860621) q[1];
sx q[1];
rz(0.11828932) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6615579) q[0];
sx q[0];
rz(-1.2999417) q[0];
sx q[0];
rz(-2.1476905) q[0];
x q[1];
rz(0.050968214) q[2];
sx q[2];
rz(-1.5545003) q[2];
sx q[2];
rz(-2.7560134) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.58913) q[1];
sx q[1];
rz(-1.5877866) q[1];
sx q[1];
rz(-0.029689992) q[1];
rz(-pi) q[2];
rz(0.30338538) q[3];
sx q[3];
rz(-1.6441364) q[3];
sx q[3];
rz(-2.908542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92363858) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(0.16036073) q[2];
rz(1.1637566) q[3];
sx q[3];
rz(-2.0573503) q[3];
sx q[3];
rz(2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9775951) q[0];
sx q[0];
rz(-1.6544592) q[0];
sx q[0];
rz(0.98597041) q[0];
rz(-1.5522955) q[1];
sx q[1];
rz(-0.28737107) q[1];
sx q[1];
rz(-1.5812965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4673104) q[0];
sx q[0];
rz(-2.5105866) q[0];
sx q[0];
rz(2.254807) q[0];
rz(-pi) q[1];
rz(-0.59448096) q[2];
sx q[2];
rz(-1.5882601) q[2];
sx q[2];
rz(1.524615) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5039798) q[1];
sx q[1];
rz(-0.39800522) q[1];
sx q[1];
rz(-2.7102273) q[1];
x q[2];
rz(1.738882) q[3];
sx q[3];
rz(-2.5017284) q[3];
sx q[3];
rz(0.17639032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61957773) q[2];
sx q[2];
rz(-1.3238944) q[2];
sx q[2];
rz(2.1066693) q[2];
rz(-0.50152913) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(2.1653304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72306776) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(-1.204741) q[0];
rz(-2.095626) q[1];
sx q[1];
rz(-3.0515262) q[1];
sx q[1];
rz(-0.14030309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36049309) q[0];
sx q[0];
rz(-1.3083959) q[0];
sx q[0];
rz(2.187832) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42201747) q[2];
sx q[2];
rz(-2.2520503) q[2];
sx q[2];
rz(1.2764507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.122815) q[1];
sx q[1];
rz(-2.5053456) q[1];
sx q[1];
rz(-0.38738363) q[1];
rz(-pi) q[2];
rz(1.9123593) q[3];
sx q[3];
rz(-1.4384801) q[3];
sx q[3];
rz(1.4190471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77421618) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(0.2224758) q[2];
rz(0.065464822) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-0.84622598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533933) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(2.5945493) q[0];
rz(-2.8833) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(2.8000854) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753174) q[0];
sx q[0];
rz(-1.5721653) q[0];
sx q[0];
rz(-1.5773236) q[0];
rz(-pi) q[1];
rz(-1.1522305) q[2];
sx q[2];
rz(-1.2593049) q[2];
sx q[2];
rz(-0.77889393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66879067) q[1];
sx q[1];
rz(-2.3040132) q[1];
sx q[1];
rz(2.0269143) q[1];
rz(-0.47291748) q[3];
sx q[3];
rz(-2.0295791) q[3];
sx q[3];
rz(-2.6948351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42939886) q[2];
sx q[2];
rz(-1.8455576) q[2];
sx q[2];
rz(-2.1923501) q[2];
rz(-0.74887577) q[3];
sx q[3];
rz(-1.2810992) q[3];
sx q[3];
rz(-0.062189814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267938) q[0];
sx q[0];
rz(-3.1068046) q[0];
sx q[0];
rz(-1.590796) q[0];
rz(-1.3560449) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(-0.063025085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7558799) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(0.066019375) q[0];
x q[1];
rz(2.5777528) q[2];
sx q[2];
rz(-0.83602521) q[2];
sx q[2];
rz(-2.6631402) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.121513) q[1];
sx q[1];
rz(-0.20993349) q[1];
sx q[1];
rz(-1.6990183) q[1];
rz(2.6659417) q[3];
sx q[3];
rz(-2.3010572) q[3];
sx q[3];
rz(0.85640872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4847792) q[2];
sx q[2];
rz(-1.3098837) q[2];
sx q[2];
rz(-2.4978034) q[2];
rz(-2.2667609) q[3];
sx q[3];
rz(-2.8372786) q[3];
sx q[3];
rz(-2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0952045) q[0];
sx q[0];
rz(-3.0859741) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(-0.19861673) q[1];
sx q[1];
rz(-3.1348517) q[1];
sx q[1];
rz(-0.14828646) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4041408) q[0];
sx q[0];
rz(-0.1830398) q[0];
sx q[0];
rz(3.1279081) q[0];
rz(1.4951493) q[2];
sx q[2];
rz(-1.9809857) q[2];
sx q[2];
rz(-0.73034053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10210382) q[1];
sx q[1];
rz(-1.7051656) q[1];
sx q[1];
rz(2.9062382) q[1];
rz(1.8174174) q[3];
sx q[3];
rz(-2.518836) q[3];
sx q[3];
rz(1.0657708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4699012) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(0.12413231) q[2];
rz(-2.5668674) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(-0.1709443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588876) q[0];
sx q[0];
rz(-3.0165065) q[0];
sx q[0];
rz(0.74147725) q[0];
rz(-0.28400907) q[1];
sx q[1];
rz(-0.0037071204) q[1];
sx q[1];
rz(-0.31518087) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3698502) q[0];
sx q[0];
rz(-0.070842243) q[0];
sx q[0];
rz(1.9596582) q[0];
rz(2.9340247) q[2];
sx q[2];
rz(-1.1313442) q[2];
sx q[2];
rz(2.07711) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1629535) q[1];
sx q[1];
rz(-2.1617315) q[1];
sx q[1];
rz(1.7617474) q[1];
x q[2];
rz(-1.7944502) q[3];
sx q[3];
rz(-1.5478494) q[3];
sx q[3];
rz(0.060088559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86889851) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(-0.72186738) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.1451274) q[3];
sx q[3];
rz(2.0731879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.566399) q[0];
sx q[0];
rz(-0.02481758) q[0];
sx q[0];
rz(-1.5665293) q[0];
rz(0.20340915) q[1];
sx q[1];
rz(-1.2982439) q[1];
sx q[1];
rz(-2.4967616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5813371) q[0];
sx q[0];
rz(-1.7783209) q[0];
sx q[0];
rz(-3.1330714) q[0];
x q[1];
rz(-0.036083607) q[2];
sx q[2];
rz(-2.5390194) q[2];
sx q[2];
rz(2.7831603) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3159807) q[1];
sx q[1];
rz(-2.59814) q[1];
sx q[1];
rz(1.7135336) q[1];
rz(-pi) q[2];
rz(0.4911812) q[3];
sx q[3];
rz(-1.1648263) q[3];
sx q[3];
rz(-0.58455672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5440392) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(-0.37975797) q[2];
rz(-1.0455658) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(1.1988962) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3658635) q[0];
sx q[0];
rz(-0.033670306) q[0];
sx q[0];
rz(-1.3566383) q[0];
rz(-0.44048539) q[1];
sx q[1];
rz(-2.0511274) q[1];
sx q[1];
rz(0.7007362) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7788575) q[0];
sx q[0];
rz(-0.88934169) q[0];
sx q[0];
rz(-0.89618857) q[0];
rz(-pi) q[1];
rz(1.3323302) q[2];
sx q[2];
rz(-2.4372413) q[2];
sx q[2];
rz(-0.65504247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4266738) q[1];
sx q[1];
rz(-2.2853984) q[1];
sx q[1];
rz(1.8780519) q[1];
rz(1.5748936) q[3];
sx q[3];
rz(-1.1381686) q[3];
sx q[3];
rz(0.062502472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7880154) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(1.8549982) q[2];
rz(-0.50518099) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(1.9193468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(-1.5420445) q[0];
rz(2.385335) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(-2.8047628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9040457) q[0];
sx q[0];
rz(-1.2531452) q[0];
sx q[0];
rz(-1.5685515) q[0];
x q[1];
rz(-2.9354503) q[2];
sx q[2];
rz(-1.8461602) q[2];
sx q[2];
rz(0.58616591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5440798) q[1];
sx q[1];
rz(-1.4855396) q[1];
sx q[1];
rz(-1.4592749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7261717) q[3];
sx q[3];
rz(-2.0924699) q[3];
sx q[3];
rz(-3.0565302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6280262) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(-2.8619518) q[2];
rz(-0.63129342) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(-2.5173371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30017988) q[0];
sx q[0];
rz(-1.5914088) q[0];
sx q[0];
rz(1.7803022) q[0];
rz(2.367876) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(2.5071267) q[2];
sx q[2];
rz(-0.96426156) q[2];
sx q[2];
rz(-2.5360863) q[2];
rz(1.5480883) q[3];
sx q[3];
rz(-1.9469713) q[3];
sx q[3];
rz(3.136829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
