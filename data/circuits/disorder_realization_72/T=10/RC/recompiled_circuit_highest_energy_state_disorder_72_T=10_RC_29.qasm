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
rz(-1.8347972) q[0];
sx q[0];
rz(-2.2936294) q[0];
sx q[0];
rz(0.037394878) q[0];
rz(-1.0132064) q[1];
sx q[1];
rz(-0.4816882) q[1];
sx q[1];
rz(-1.6964635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28222154) q[0];
sx q[0];
rz(-1.9992775) q[0];
sx q[0];
rz(1.3708049) q[0];
rz(-pi) q[1];
x q[1];
rz(0.070568512) q[2];
sx q[2];
rz(-1.6733992) q[2];
sx q[2];
rz(-3.0377394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5132222) q[1];
sx q[1];
rz(-2.0525825) q[1];
sx q[1];
rz(-2.1867153) q[1];
rz(-pi) q[2];
rz(-1.0057428) q[3];
sx q[3];
rz(-1.4742719) q[3];
sx q[3];
rz(2.9526476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41363132) q[2];
sx q[2];
rz(-0.093158826) q[2];
sx q[2];
rz(1.5067345) q[2];
rz(-0.19351752) q[3];
sx q[3];
rz(-0.7842803) q[3];
sx q[3];
rz(-2.1957652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0993318) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(0.17901626) q[0];
rz(1.5366813) q[1];
sx q[1];
rz(-0.34114006) q[1];
sx q[1];
rz(-1.5515597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66288469) q[0];
sx q[0];
rz(-0.30618069) q[0];
sx q[0];
rz(-1.2423489) q[0];
rz(-pi) q[1];
rz(0.29921542) q[2];
sx q[2];
rz(-0.73198527) q[2];
sx q[2];
rz(-0.74886403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1074321) q[1];
sx q[1];
rz(-1.1723659) q[1];
sx q[1];
rz(-2.855019) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47081866) q[3];
sx q[3];
rz(-2.0662466) q[3];
sx q[3];
rz(-1.9268056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67148525) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(3.1191077) q[2];
rz(-0.89618987) q[3];
sx q[3];
rz(-0.78138566) q[3];
sx q[3];
rz(1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80980587) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(1.5326387) q[0];
rz(2.0017852) q[1];
sx q[1];
rz(-0.31875691) q[1];
sx q[1];
rz(-2.2679813) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667359) q[0];
sx q[0];
rz(-2.2169161) q[0];
sx q[0];
rz(2.2227051) q[0];
x q[1];
rz(-2.1223091) q[2];
sx q[2];
rz(-2.1423369) q[2];
sx q[2];
rz(-0.71046605) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.552678) q[1];
sx q[1];
rz(-2.6323018) q[1];
sx q[1];
rz(-0.31917787) q[1];
rz(-1.4954044) q[3];
sx q[3];
rz(-1.2986548) q[3];
sx q[3];
rz(2.4484602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56796271) q[2];
sx q[2];
rz(-2.6831388) q[2];
sx q[2];
rz(-2.6652794) q[2];
rz(-2.1161946) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(0.29388139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3697701) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(-2.5452132) q[0];
rz(0.75309938) q[1];
sx q[1];
rz(-1.785708) q[1];
sx q[1];
rz(0.3515884) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542587) q[0];
sx q[0];
rz(-2.0947347) q[0];
sx q[0];
rz(-1.41904) q[0];
rz(0.48248191) q[2];
sx q[2];
rz(-0.94581476) q[2];
sx q[2];
rz(-0.04908726) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6476485) q[1];
sx q[1];
rz(-2.679849) q[1];
sx q[1];
rz(2.1069202) q[1];
x q[2];
rz(-0.80644108) q[3];
sx q[3];
rz(-2.3262352) q[3];
sx q[3];
rz(-1.7585332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4190462) q[2];
sx q[2];
rz(-1.498035) q[2];
sx q[2];
rz(2.2124115) q[2];
rz(1.9123214) q[3];
sx q[3];
rz(-0.036673948) q[3];
sx q[3];
rz(1.8721972) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1790328) q[0];
sx q[0];
rz(-0.61229175) q[0];
sx q[0];
rz(-0.52781934) q[0];
rz(-2.5714696) q[1];
sx q[1];
rz(-2.5716883) q[1];
sx q[1];
rz(2.747587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4657426) q[0];
sx q[0];
rz(-2.619474) q[0];
sx q[0];
rz(2.7444849) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2038241) q[2];
sx q[2];
rz(-1.0093371) q[2];
sx q[2];
rz(0.16753627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3574419) q[1];
sx q[1];
rz(-1.7926182) q[1];
sx q[1];
rz(-0.12089575) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3218003) q[3];
sx q[3];
rz(-0.91320437) q[3];
sx q[3];
rz(-2.4321604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7492619) q[2];
sx q[2];
rz(-1.2886084) q[2];
sx q[2];
rz(1.1345081) q[2];
rz(-0.025731651) q[3];
sx q[3];
rz(-1.0545571) q[3];
sx q[3];
rz(-2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57539097) q[0];
sx q[0];
rz(-1.3302777) q[0];
sx q[0];
rz(-2.6824685) q[0];
rz(-2.3205914) q[1];
sx q[1];
rz(-0.88290015) q[1];
sx q[1];
rz(0.51148907) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82152459) q[0];
sx q[0];
rz(-1.2787158) q[0];
sx q[0];
rz(-2.065737) q[0];
x q[1];
rz(0.66368565) q[2];
sx q[2];
rz(-1.384871) q[2];
sx q[2];
rz(2.3472361) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.59907179) q[1];
sx q[1];
rz(-0.59172809) q[1];
sx q[1];
rz(-0.27246957) q[1];
x q[2];
rz(1.4019764) q[3];
sx q[3];
rz(-0.69069117) q[3];
sx q[3];
rz(2.9258779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5795827) q[2];
sx q[2];
rz(-1.1082114) q[2];
sx q[2];
rz(0.75135922) q[2];
rz(-2.5747418) q[3];
sx q[3];
rz(-1.6040498) q[3];
sx q[3];
rz(-1.3273299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71996561) q[0];
sx q[0];
rz(-2.9599074) q[0];
sx q[0];
rz(3.1340461) q[0];
rz(-0.94002062) q[1];
sx q[1];
rz(-0.66497856) q[1];
sx q[1];
rz(3.0090581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7791393) q[0];
sx q[0];
rz(-0.7374239) q[0];
sx q[0];
rz(2.2064184) q[0];
x q[1];
rz(-1.7039677) q[2];
sx q[2];
rz(-1.9308763) q[2];
sx q[2];
rz(-0.11225004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30250994) q[1];
sx q[1];
rz(-1.1801475) q[1];
sx q[1];
rz(-1.7724228) q[1];
x q[2];
rz(1.9136732) q[3];
sx q[3];
rz(-1.7352805) q[3];
sx q[3];
rz(3.0066688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14564766) q[2];
sx q[2];
rz(-1.6463248) q[2];
sx q[2];
rz(0.10124595) q[2];
rz(0.64949399) q[3];
sx q[3];
rz(-0.5972623) q[3];
sx q[3];
rz(-0.35972843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7262909) q[0];
sx q[0];
rz(-1.2969718) q[0];
sx q[0];
rz(1.8743961) q[0];
rz(1.2665117) q[1];
sx q[1];
rz(-2.9837065) q[1];
sx q[1];
rz(-2.3983541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6544391) q[0];
sx q[0];
rz(-1.7726573) q[0];
sx q[0];
rz(2.5623723) q[0];
rz(-pi) q[1];
rz(1.2879804) q[2];
sx q[2];
rz(-2.3415945) q[2];
sx q[2];
rz(-1.204042) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.69276492) q[1];
sx q[1];
rz(-1.3346486) q[1];
sx q[1];
rz(-2.9797878) q[1];
rz(2.4489347) q[3];
sx q[3];
rz(-2.1941278) q[3];
sx q[3];
rz(0.75875635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0172952) q[2];
sx q[2];
rz(-2.9399019) q[2];
sx q[2];
rz(1.4264433) q[2];
rz(2.1258449) q[3];
sx q[3];
rz(-1.7044715) q[3];
sx q[3];
rz(-0.82971853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5468686) q[0];
sx q[0];
rz(-2.4877553) q[0];
sx q[0];
rz(3.1256909) q[0];
rz(-1.6390027) q[1];
sx q[1];
rz(-2.465261) q[1];
sx q[1];
rz(-2.1471088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.573365) q[0];
sx q[0];
rz(-1.957743) q[0];
sx q[0];
rz(-0.9181482) q[0];
rz(1.5076958) q[2];
sx q[2];
rz(-1.1931538) q[2];
sx q[2];
rz(-1.3444195) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2331824) q[1];
sx q[1];
rz(-0.841978) q[1];
sx q[1];
rz(-0.89836095) q[1];
rz(1.0570378) q[3];
sx q[3];
rz(-1.4389068) q[3];
sx q[3];
rz(-1.7943602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9859621) q[2];
sx q[2];
rz(-1.7113643) q[2];
sx q[2];
rz(-0.17781167) q[2];
rz(1.6832247) q[3];
sx q[3];
rz(-2.7198313) q[3];
sx q[3];
rz(-1.873707) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47852248) q[0];
sx q[0];
rz(-1.8926184) q[0];
sx q[0];
rz(1.2836237) q[0];
rz(-0.78370699) q[1];
sx q[1];
rz(-0.77373928) q[1];
sx q[1];
rz(1.8342038) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2268902) q[0];
sx q[0];
rz(-2.693334) q[0];
sx q[0];
rz(-0.21313711) q[0];
rz(-pi) q[1];
rz(3.1112017) q[2];
sx q[2];
rz(-1.1830529) q[2];
sx q[2];
rz(-3.0770242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5518903) q[1];
sx q[1];
rz(-1.74472) q[1];
sx q[1];
rz(-3.0822192) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94759192) q[3];
sx q[3];
rz(-2.4952808) q[3];
sx q[3];
rz(1.0912967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6002097) q[2];
sx q[2];
rz(-1.4098488) q[2];
sx q[2];
rz(-0.48995885) q[2];
rz(2.0269035) q[3];
sx q[3];
rz(-1.9652818) q[3];
sx q[3];
rz(-2.8118242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(3.1212696) q[0];
sx q[0];
rz(-1.2280432) q[0];
sx q[0];
rz(-1.8228774) q[0];
rz(1.2976788) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(-1.8280468) q[2];
sx q[2];
rz(-1.8218818) q[2];
sx q[2];
rz(2.8043361) q[2];
rz(-0.86434563) q[3];
sx q[3];
rz(-1.4180687) q[3];
sx q[3];
rz(-1.1579469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
