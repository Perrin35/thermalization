OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(3.5482121) q[0];
sx q[0];
rz(9.1756048) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43139797) q[0];
sx q[0];
rz(-0.22127998) q[0];
sx q[0];
rz(-1.6943323) q[0];
x q[1];
rz(-2.2416441) q[2];
sx q[2];
rz(-0.64621011) q[2];
sx q[2];
rz(2.3044569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33221153) q[1];
sx q[1];
rz(-0.84377938) q[1];
sx q[1];
rz(1.0512645) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95991858) q[3];
sx q[3];
rz(-1.6392518) q[3];
sx q[3];
rz(1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41574079) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7249811) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(0.38831098) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8578051) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(-2.6736205) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1158893) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(2.2504836) q[1];
x q[2];
rz(-0.87725957) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(0.045543268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2018532) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2597044) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(2.846068) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(2.3957516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5999818) q[0];
sx q[0];
rz(-0.2395425) q[0];
sx q[0];
rz(-1.177686) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52650555) q[2];
sx q[2];
rz(-1.7130934) q[2];
sx q[2];
rz(-1.7775747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6927166) q[1];
sx q[1];
rz(-1.2926896) q[1];
sx q[1];
rz(2.6443308) q[1];
x q[2];
rz(1.8196705) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(-0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(-0.27329683) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.435794) q[0];
sx q[0];
rz(-1.6596284) q[0];
sx q[0];
rz(-1.8949132) q[0];
rz(2.4934713) q[2];
sx q[2];
rz(-1.1037877) q[2];
sx q[2];
rz(1.311122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8400152) q[1];
sx q[1];
rz(-1.8043648) q[1];
sx q[1];
rz(0.15109269) q[1];
rz(0.64951879) q[3];
sx q[3];
rz(-1.7139385) q[3];
sx q[3];
rz(2.3424145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(3.058847) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(2.1544429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50901978) q[0];
sx q[0];
rz(-2.1633254) q[0];
sx q[0];
rz(0.49125262) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5802025) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(0.26088342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8122711) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(-0.66004628) q[1];
rz(1.1779285) q[3];
sx q[3];
rz(-1.6097027) q[3];
sx q[3];
rz(1.2272569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(-0.42090297) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6565276) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(-2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.8621559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56023843) q[0];
sx q[0];
rz(-1.2890153) q[0];
sx q[0];
rz(3.1262585) q[0];
rz(1.860414) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(2.9786125) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93950677) q[1];
sx q[1];
rz(-2.0923951) q[1];
sx q[1];
rz(-3.107198) q[1];
rz(1.7518696) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(0.74336038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(-2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8188748) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-0.0035704426) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17655003) q[0];
sx q[0];
rz(-1.1561484) q[0];
sx q[0];
rz(1.167017) q[0];
x q[1];
rz(0.75242075) q[2];
sx q[2];
rz(-1.8634897) q[2];
sx q[2];
rz(-0.21586403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7230941) q[1];
sx q[1];
rz(-1.2730036) q[1];
sx q[1];
rz(-0.94350068) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3725029) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(1.2411181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(0.92010951) q[2];
rz(-0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(2.7238817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(2.4628941) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6060651) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(0.33959629) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6762814) q[2];
sx q[2];
rz(-2.2530472) q[2];
sx q[2];
rz(-1.0516143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.855367) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(-1.0463868) q[1];
x q[2];
rz(-1.7292761) q[3];
sx q[3];
rz(-1.9915808) q[3];
sx q[3];
rz(-2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(0.91782451) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(2.881799) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(-1.2277463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(2.3234141) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7582015) q[1];
sx q[1];
rz(-1.2609298) q[1];
sx q[1];
rz(-1.8499225) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6960877) q[3];
sx q[3];
rz(-1.8337436) q[3];
sx q[3];
rz(1.6642237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-0.79088598) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(2.1561484) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26786131) q[0];
sx q[0];
rz(-0.16927969) q[0];
sx q[0];
rz(-1.4219567) q[0];
x q[1];
rz(0.12014328) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(0.10211589) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4611778) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(-1.3780891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4264614) q[3];
sx q[3];
rz(-0.22278856) q[3];
sx q[3];
rz(-0.95105329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(1.7760989) q[2];
sx q[2];
rz(-1.4751954) q[2];
sx q[2];
rz(1.8196646) q[2];
rz(2.2684569) q[3];
sx q[3];
rz(-2.1921039) q[3];
sx q[3];
rz(2.3154142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];