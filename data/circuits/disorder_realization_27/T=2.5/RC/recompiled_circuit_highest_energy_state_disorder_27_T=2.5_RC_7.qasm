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
rz(1.4577515) q[0];
sx q[0];
rz(-3.0384851) q[0];
sx q[0];
rz(-0.74080324) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(3.8771602) q[1];
sx q[1];
rz(9.1280042) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4838801) q[0];
sx q[0];
rz(-1.9498511) q[0];
sx q[0];
rz(2.0188278) q[0];
x q[1];
rz(-0.27539092) q[2];
sx q[2];
rz(-0.73705929) q[2];
sx q[2];
rz(-3.0050957) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.43208968) q[1];
sx q[1];
rz(-2.5910834) q[1];
sx q[1];
rz(-3.1269323) q[1];
rz(1.9363251) q[3];
sx q[3];
rz(-2.6726279) q[3];
sx q[3];
rz(2.7597357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2510117) q[2];
sx q[2];
rz(-0.0958395) q[2];
sx q[2];
rz(-1.4061692) q[2];
rz(0.74875325) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(2.9296618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1521848) q[0];
sx q[0];
rz(-0.87225544) q[0];
sx q[0];
rz(-2.8023791) q[0];
rz(-2.5665414) q[1];
sx q[1];
rz(-1.1851858) q[1];
sx q[1];
rz(-2.7599879) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0038892) q[0];
sx q[0];
rz(-0.48332542) q[0];
sx q[0];
rz(2.0893731) q[0];
rz(0.26338844) q[2];
sx q[2];
rz(-1.2424505) q[2];
sx q[2];
rz(-1.245861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84880616) q[1];
sx q[1];
rz(-1.5804351) q[1];
sx q[1];
rz(-1.1718434) q[1];
x q[2];
rz(0.80581237) q[3];
sx q[3];
rz(-1.9776113) q[3];
sx q[3];
rz(0.47570566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40994689) q[2];
sx q[2];
rz(-0.68845981) q[2];
sx q[2];
rz(2.6551969) q[2];
rz(2.5977123) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(-0.16415183) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51089066) q[0];
sx q[0];
rz(-0.4158026) q[0];
sx q[0];
rz(2.1422332) q[0];
rz(0.73127812) q[1];
sx q[1];
rz(-0.46842289) q[1];
sx q[1];
rz(2.9670002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.664266) q[0];
sx q[0];
rz(-1.5683929) q[0];
sx q[0];
rz(1.1241358) q[0];
x q[1];
rz(2.9810437) q[2];
sx q[2];
rz(-1.7458651) q[2];
sx q[2];
rz(2.0714945) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5970351) q[1];
sx q[1];
rz(-1.6591151) q[1];
sx q[1];
rz(-3.0889325) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9932132) q[3];
sx q[3];
rz(-0.49634051) q[3];
sx q[3];
rz(-1.5143192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38873765) q[2];
sx q[2];
rz(-0.85192215) q[2];
sx q[2];
rz(1.8641776) q[2];
rz(-0.82469213) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(-2.9741014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.326062) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(-1.1482358) q[0];
rz(2.8094021) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(-2.0567599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6713632) q[0];
sx q[0];
rz(-1.7106904) q[0];
sx q[0];
rz(1.8663919) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81697322) q[2];
sx q[2];
rz(-1.8228056) q[2];
sx q[2];
rz(-1.9793881) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12382896) q[1];
sx q[1];
rz(-0.74018407) q[1];
sx q[1];
rz(0.077880903) q[1];
rz(-2.3938136) q[3];
sx q[3];
rz(-2.5430395) q[3];
sx q[3];
rz(0.35394305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3499902) q[2];
sx q[2];
rz(-2.377066) q[2];
sx q[2];
rz(3.1349658) q[2];
rz(-0.22348063) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(-2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4357736) q[0];
sx q[0];
rz(-2.3888102) q[0];
sx q[0];
rz(1.1821795) q[0];
rz(-1.8403107) q[1];
sx q[1];
rz(-0.92295206) q[1];
sx q[1];
rz(-2.9822947) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.015291) q[0];
sx q[0];
rz(-1.2833923) q[0];
sx q[0];
rz(-1.6096344) q[0];
x q[1];
rz(1.2550687) q[2];
sx q[2];
rz(-1.7509394) q[2];
sx q[2];
rz(2.2314928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0634126) q[1];
sx q[1];
rz(-2.0495119) q[1];
sx q[1];
rz(-0.53005752) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0924657) q[3];
sx q[3];
rz(-0.73724607) q[3];
sx q[3];
rz(-2.3741219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84676877) q[2];
sx q[2];
rz(-0.20192768) q[2];
sx q[2];
rz(1.3429886) q[2];
rz(-2.790847) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(-2.8158367) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419234) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(2.9747466) q[0];
rz(2.9150561) q[1];
sx q[1];
rz(-1.8140565) q[1];
sx q[1];
rz(-0.47795263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27875459) q[0];
sx q[0];
rz(-2.4669381) q[0];
sx q[0];
rz(1.007403) q[0];
x q[1];
rz(1.208838) q[2];
sx q[2];
rz(-2.2565292) q[2];
sx q[2];
rz(0.63206965) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30580995) q[1];
sx q[1];
rz(-2.1205582) q[1];
sx q[1];
rz(0.98819403) q[1];
rz(2.7109954) q[3];
sx q[3];
rz(-1.9816508) q[3];
sx q[3];
rz(1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53092521) q[2];
sx q[2];
rz(-1.3767367) q[2];
sx q[2];
rz(-1.2923856) q[2];
rz(0.90062201) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(1.5463411) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6677299) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(-1.6170172) q[0];
rz(1.4432817) q[1];
sx q[1];
rz(-0.89569211) q[1];
sx q[1];
rz(-0.5009833) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65376908) q[0];
sx q[0];
rz(-2.1697576) q[0];
sx q[0];
rz(0.45287153) q[0];
rz(-0.55866629) q[2];
sx q[2];
rz(-1.8353454) q[2];
sx q[2];
rz(3.0857744) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2051207) q[1];
sx q[1];
rz(-2.3354482) q[1];
sx q[1];
rz(-2.2938726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30876183) q[3];
sx q[3];
rz(-1.689286) q[3];
sx q[3];
rz(-1.0498966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3524126) q[2];
sx q[2];
rz(-3.0173306) q[2];
sx q[2];
rz(-2.6782356) q[2];
rz(3.0739259) q[3];
sx q[3];
rz(-1.8384408) q[3];
sx q[3];
rz(2.9786003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8891334) q[0];
sx q[0];
rz(-1.8661789) q[0];
sx q[0];
rz(2.6305991) q[0];
rz(2.4976318) q[1];
sx q[1];
rz(-1.124758) q[1];
sx q[1];
rz(-0.28800979) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8159008) q[0];
sx q[0];
rz(-1.3790695) q[0];
sx q[0];
rz(1.6098538) q[0];
x q[1];
rz(-2.5019849) q[2];
sx q[2];
rz(-2.158463) q[2];
sx q[2];
rz(-1.7220875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0916413) q[1];
sx q[1];
rz(-0.63739751) q[1];
sx q[1];
rz(2.6693488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8657416) q[3];
sx q[3];
rz(-2.7497661) q[3];
sx q[3];
rz(-1.0759169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94706941) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(0.21491773) q[2];
rz(1.8042709) q[3];
sx q[3];
rz(-0.53842068) q[3];
sx q[3];
rz(0.18068331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27338481) q[0];
sx q[0];
rz(-2.2053563) q[0];
sx q[0];
rz(-1.0816164) q[0];
rz(0.99010211) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(2.8066011) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2401857) q[0];
sx q[0];
rz(-1.6044582) q[0];
sx q[0];
rz(1.89374) q[0];
x q[1];
rz(1.9926979) q[2];
sx q[2];
rz(-1.5765911) q[2];
sx q[2];
rz(-0.96265477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8674842) q[1];
sx q[1];
rz(-1.47164) q[1];
sx q[1];
rz(-1.4851557) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5906628) q[3];
sx q[3];
rz(-0.40106138) q[3];
sx q[3];
rz(3.0326518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.058502402) q[2];
sx q[2];
rz(-0.52512705) q[2];
sx q[2];
rz(-0.38880175) q[2];
rz(0.36924103) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19287547) q[0];
sx q[0];
rz(-0.95272869) q[0];
sx q[0];
rz(0.70190758) q[0];
rz(-1.6096055) q[1];
sx q[1];
rz(-0.67370266) q[1];
sx q[1];
rz(2.7598377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8030539) q[0];
sx q[0];
rz(-0.10333867) q[0];
sx q[0];
rz(2.4843012) q[0];
rz(-pi) q[1];
rz(1.6959305) q[2];
sx q[2];
rz(-1.9788673) q[2];
sx q[2];
rz(0.073485188) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81667852) q[1];
sx q[1];
rz(-2.7902998) q[1];
sx q[1];
rz(2.3354946) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66017229) q[3];
sx q[3];
rz(-2.2260465) q[3];
sx q[3];
rz(-1.8546263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4664885) q[2];
sx q[2];
rz(-1.0305104) q[2];
sx q[2];
rz(-2.4934736) q[2];
rz(1.0068007) q[3];
sx q[3];
rz(-1.1454134) q[3];
sx q[3];
rz(0.072331585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61160144) q[0];
sx q[0];
rz(-1.7399104) q[0];
sx q[0];
rz(-2.4763784) q[0];
rz(2.8499659) q[1];
sx q[1];
rz(-1.7886152) q[1];
sx q[1];
rz(2.0033966) q[1];
rz(-2.6219528) q[2];
sx q[2];
rz(-1.8815132) q[2];
sx q[2];
rz(2.5757679) q[2];
rz(-1.9129561) q[3];
sx q[3];
rz(-0.4426601) q[3];
sx q[3];
rz(1.9667251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
