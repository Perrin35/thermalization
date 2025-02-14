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
rz(-2.2588377) q[0];
sx q[0];
rz(-0.14482276) q[0];
sx q[0];
rz(2.3902399) q[0];
rz(1.524628) q[1];
sx q[1];
rz(-2.1722062) q[1];
sx q[1];
rz(0.17925395) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67231945) q[0];
sx q[0];
rz(-1.0814965) q[0];
sx q[0];
rz(-0.31764046) q[0];
rz(2.4186736) q[2];
sx q[2];
rz(-1.5953879) q[2];
sx q[2];
rz(-1.6484156) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1826659) q[1];
sx q[1];
rz(-1.8799056) q[1];
sx q[1];
rz(1.6462417) q[1];
x q[2];
rz(-1.6514014) q[3];
sx q[3];
rz(-2.5063519) q[3];
sx q[3];
rz(1.0280392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2396607) q[2];
sx q[2];
rz(-1.8472981) q[2];
sx q[2];
rz(2.530976) q[2];
rz(-2.2527952) q[3];
sx q[3];
rz(-2.461268) q[3];
sx q[3];
rz(2.4077267) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69773847) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.3296211) q[0];
rz(-1.656172) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(0.86404538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3812019) q[0];
sx q[0];
rz(-1.4684033) q[0];
sx q[0];
rz(1.824462) q[0];
rz(-1.1077706) q[2];
sx q[2];
rz(-1.290375) q[2];
sx q[2];
rz(-3.1299431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6384008) q[1];
sx q[1];
rz(-1.710689) q[1];
sx q[1];
rz(-1.3963455) q[1];
x q[2];
rz(-2.1883499) q[3];
sx q[3];
rz(-0.88897486) q[3];
sx q[3];
rz(-2.1566331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0835421) q[2];
sx q[2];
rz(-0.6821878) q[2];
sx q[2];
rz(-2.9502499) q[2];
rz(2.8355016) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(0.93650854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8323583) q[0];
sx q[0];
rz(-1.2795871) q[0];
sx q[0];
rz(0.424463) q[0];
rz(1.8008495) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(2.4901966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6474045) q[0];
sx q[0];
rz(-1.5636824) q[0];
sx q[0];
rz(1.4071541) q[0];
rz(-pi) q[1];
rz(1.7537139) q[2];
sx q[2];
rz(-1.9486893) q[2];
sx q[2];
rz(1.479508) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9561852) q[1];
sx q[1];
rz(-0.53403097) q[1];
sx q[1];
rz(2.7772481) q[1];
rz(-pi) q[2];
rz(-1.3440787) q[3];
sx q[3];
rz(-1.4170839) q[3];
sx q[3];
rz(0.24776974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.035630781) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(-2.4448709) q[2];
rz(-0.68909711) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(-2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7032787) q[0];
sx q[0];
rz(-2.8499481) q[0];
sx q[0];
rz(2.1412204) q[0];
rz(1.6814303) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.2132852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54736894) q[0];
sx q[0];
rz(-1.692019) q[0];
sx q[0];
rz(1.2884883) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4126265) q[2];
sx q[2];
rz(-1.2104958) q[2];
sx q[2];
rz(-3.136022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.69566947) q[1];
sx q[1];
rz(-1.4497821) q[1];
sx q[1];
rz(-2.8112429) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0509316) q[3];
sx q[3];
rz(-2.1075776) q[3];
sx q[3];
rz(-1.8859552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0393684) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-0.064112045) q[2];
rz(2.4168849) q[3];
sx q[3];
rz(-1.2084081) q[3];
sx q[3];
rz(-2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.0940014) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-2.3498348) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(1.3349104) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955973) q[0];
sx q[0];
rz(-2.3438518) q[0];
sx q[0];
rz(2.255385) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15891517) q[2];
sx q[2];
rz(-2.5549485) q[2];
sx q[2];
rz(1.212709) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89692749) q[1];
sx q[1];
rz(-1.7411971) q[1];
sx q[1];
rz(-0.076431304) q[1];
rz(-pi) q[2];
rz(-1.8800635) q[3];
sx q[3];
rz(-2.5430508) q[3];
sx q[3];
rz(-2.8054597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2935334) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(2.000957) q[2];
rz(-3.0349777) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(0.30985668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(-0.003224592) q[0];
rz(0.047317304) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.7914194) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4989717) q[0];
sx q[0];
rz(-1.3184217) q[0];
sx q[0];
rz(-3.0809771) q[0];
rz(-pi) q[1];
x q[1];
rz(2.240594) q[2];
sx q[2];
rz(-2.0428847) q[2];
sx q[2];
rz(-2.6028518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1803169) q[1];
sx q[1];
rz(-0.84876981) q[1];
sx q[1];
rz(0.78603334) q[1];
x q[2];
rz(1.7058701) q[3];
sx q[3];
rz(-2.9189442) q[3];
sx q[3];
rz(-0.50375736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99001592) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-3.0604176) q[2];
rz(-0.89933991) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(-1.8544633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.4261632) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(-0.50450605) q[0];
rz(2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(-3.1212433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85914579) q[0];
sx q[0];
rz(-1.9345645) q[0];
sx q[0];
rz(3.0412741) q[0];
rz(1.7791012) q[2];
sx q[2];
rz(-1.4464738) q[2];
sx q[2];
rz(2.1712239) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.523218) q[1];
sx q[1];
rz(-1.3058387) q[1];
sx q[1];
rz(-2.9441903) q[1];
x q[2];
rz(-2.0836104) q[3];
sx q[3];
rz(-0.098976299) q[3];
sx q[3];
rz(1.3266226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4703579) q[2];
sx q[2];
rz(-1.3037668) q[2];
sx q[2];
rz(-1.7669558) q[2];
rz(-1.6124407) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(-3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2328211) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(0.94295162) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97639192) q[0];
sx q[0];
rz(-1.2726674) q[0];
sx q[0];
rz(-2.0701029) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23343691) q[2];
sx q[2];
rz(-0.95165157) q[2];
sx q[2];
rz(0.34282986) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8170191) q[1];
sx q[1];
rz(-0.96265154) q[1];
sx q[1];
rz(-0.84546802) q[1];
rz(-pi) q[2];
rz(-3.037463) q[3];
sx q[3];
rz(-1.2648598) q[3];
sx q[3];
rz(-0.96398523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1711787) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(-1.9937493) q[2];
rz(-2.7796699) q[3];
sx q[3];
rz(-1.3779093) q[3];
sx q[3];
rz(1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4107133) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-0.10243375) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(-2.3238497) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9303904) q[0];
sx q[0];
rz(-1.0387392) q[0];
sx q[0];
rz(-2.3389057) q[0];
rz(-pi) q[1];
rz(-1.7162343) q[2];
sx q[2];
rz(-0.67611968) q[2];
sx q[2];
rz(-2.6349677) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9407258) q[1];
sx q[1];
rz(-1.5897017) q[1];
sx q[1];
rz(2.3000642) q[1];
x q[2];
rz(-0.46157591) q[3];
sx q[3];
rz(-0.99050039) q[3];
sx q[3];
rz(-0.020315276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(-2.8279772) q[2];
rz(0.46401986) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.87026507) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(-0.23649293) q[0];
rz(0.21367167) q[1];
sx q[1];
rz(-2.2387319) q[1];
sx q[1];
rz(0.22458354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86041566) q[0];
sx q[0];
rz(-0.91360859) q[0];
sx q[0];
rz(-2.1967931) q[0];
rz(-pi) q[1];
rz(0.5245536) q[2];
sx q[2];
rz(-0.65215014) q[2];
sx q[2];
rz(2.1422276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27222363) q[1];
sx q[1];
rz(-1.0450796) q[1];
sx q[1];
rz(1.438526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8524695) q[3];
sx q[3];
rz(-1.895088) q[3];
sx q[3];
rz(1.3642023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1647722) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(2.9898047) q[2];
rz(-0.031489059) q[3];
sx q[3];
rz(-1.3969235) q[3];
sx q[3];
rz(0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0781773) q[0];
sx q[0];
rz(-1.6178394) q[0];
sx q[0];
rz(2.7375426) q[0];
rz(1.0833441) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(-1.5189717) q[2];
sx q[2];
rz(-1.3695649) q[2];
sx q[2];
rz(2.8893378) q[2];
rz(-0.32947551) q[3];
sx q[3];
rz(-1.4737045) q[3];
sx q[3];
rz(1.2274689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
