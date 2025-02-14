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
rz(0.12953144) q[0];
sx q[0];
rz(3.01053) q[0];
sx q[0];
rz(13.170903) q[0];
rz(-0.52855748) q[1];
sx q[1];
rz(2.8837535) q[1];
sx q[1];
rz(16.937994) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.226868) q[0];
sx q[0];
rz(-2.719923) q[0];
sx q[0];
rz(2.492621) q[0];
rz(2.5276353) q[2];
sx q[2];
rz(-1.7579305) q[2];
sx q[2];
rz(2.2452315) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9046272) q[1];
sx q[1];
rz(-0.89348999) q[1];
sx q[1];
rz(-0.64579247) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0139066) q[3];
sx q[3];
rz(-2.0632072) q[3];
sx q[3];
rz(-3.1000053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0872385) q[2];
sx q[2];
rz(-1.1776935) q[2];
sx q[2];
rz(-3.0053906) q[2];
rz(2.6916091) q[3];
sx q[3];
rz(-0.68322244) q[3];
sx q[3];
rz(-2.3581678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1012652) q[0];
sx q[0];
rz(-2.3430921) q[0];
sx q[0];
rz(-0.26327565) q[0];
rz(2.1030262) q[1];
sx q[1];
rz(-0.34663215) q[1];
sx q[1];
rz(-0.77478772) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8717125) q[0];
sx q[0];
rz(-0.96202606) q[0];
sx q[0];
rz(-1.8222429) q[0];
rz(2.4500174) q[2];
sx q[2];
rz(-2.1082768) q[2];
sx q[2];
rz(0.36851685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0475743) q[1];
sx q[1];
rz(-1.1915922) q[1];
sx q[1];
rz(3.1017041) q[1];
x q[2];
rz(0.56258308) q[3];
sx q[3];
rz(-2.354971) q[3];
sx q[3];
rz(-2.7275769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.420555) q[2];
sx q[2];
rz(-1.6246395) q[2];
sx q[2];
rz(-2.5431385) q[2];
rz(1.362644) q[3];
sx q[3];
rz(-2.2612031) q[3];
sx q[3];
rz(2.8708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8552928) q[0];
sx q[0];
rz(-0.4628276) q[0];
sx q[0];
rz(1.6337974) q[0];
rz(0.15448054) q[1];
sx q[1];
rz(-1.721761) q[1];
sx q[1];
rz(-0.78937626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3837483) q[0];
sx q[0];
rz(-2.0998635) q[0];
sx q[0];
rz(-1.6711414) q[0];
x q[1];
rz(2.8604498) q[2];
sx q[2];
rz(-1.4340377) q[2];
sx q[2];
rz(-3.1355372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7676413) q[1];
sx q[1];
rz(-2.1726296) q[1];
sx q[1];
rz(-0.37202073) q[1];
rz(-2.3422935) q[3];
sx q[3];
rz(-1.2033312) q[3];
sx q[3];
rz(0.41900837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7827451) q[2];
sx q[2];
rz(-2.2300356) q[2];
sx q[2];
rz(1.4770799) q[2];
rz(-1.2137671) q[3];
sx q[3];
rz(-1.8002847) q[3];
sx q[3];
rz(1.7267797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.7724991) q[0];
sx q[0];
rz(-2.0857683) q[0];
sx q[0];
rz(2.5308727) q[0];
rz(-2.4702813) q[1];
sx q[1];
rz(-1.5076312) q[1];
sx q[1];
rz(1.7866561) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1917065) q[0];
sx q[0];
rz(-1.0903795) q[0];
sx q[0];
rz(-2.0818162) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2569095) q[2];
sx q[2];
rz(-2.7531024) q[2];
sx q[2];
rz(1.9635278) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0555041) q[1];
sx q[1];
rz(-1.633606) q[1];
sx q[1];
rz(-1.8163866) q[1];
rz(-0.98396222) q[3];
sx q[3];
rz(-1.4915183) q[3];
sx q[3];
rz(2.6266499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2063724) q[2];
sx q[2];
rz(-2.361203) q[2];
sx q[2];
rz(-2.8727403) q[2];
rz(2.3937461) q[3];
sx q[3];
rz(-0.81955376) q[3];
sx q[3];
rz(1.6929172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.185323) q[0];
sx q[0];
rz(-1.7522426) q[0];
sx q[0];
rz(-0.67778936) q[0];
rz(-2.9970844) q[1];
sx q[1];
rz(-0.78737193) q[1];
sx q[1];
rz(1.6544624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9155898) q[0];
sx q[0];
rz(-1.5446413) q[0];
sx q[0];
rz(1.7067065) q[0];
rz(1.9727835) q[2];
sx q[2];
rz(-1.5176519) q[2];
sx q[2];
rz(1.0099908) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37461149) q[1];
sx q[1];
rz(-2.0882147) q[1];
sx q[1];
rz(-2.6042038) q[1];
rz(-pi) q[2];
rz(-1.6923589) q[3];
sx q[3];
rz(-0.37827493) q[3];
sx q[3];
rz(0.88751436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3378478) q[2];
sx q[2];
rz(-2.8779112) q[2];
sx q[2];
rz(0.35550508) q[2];
rz(-1.0454987) q[3];
sx q[3];
rz(-1.239536) q[3];
sx q[3];
rz(2.579328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1657555) q[0];
sx q[0];
rz(-0.58335692) q[0];
sx q[0];
rz(-0.088706644) q[0];
rz(1.6070131) q[1];
sx q[1];
rz(-1.3101703) q[1];
sx q[1];
rz(-3.065899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.158094) q[0];
sx q[0];
rz(-0.69639054) q[0];
sx q[0];
rz(-2.8386001) q[0];
x q[1];
rz(0.050932542) q[2];
sx q[2];
rz(-0.83980745) q[2];
sx q[2];
rz(-2.5185846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20670375) q[1];
sx q[1];
rz(-2.5674501) q[1];
sx q[1];
rz(-1.8341792) q[1];
rz(-pi) q[2];
rz(0.47786062) q[3];
sx q[3];
rz(-2.5091672) q[3];
sx q[3];
rz(0.53148182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82464108) q[2];
sx q[2];
rz(-1.757916) q[2];
sx q[2];
rz(1.0392044) q[2];
rz(0.21229395) q[3];
sx q[3];
rz(-2.0391235) q[3];
sx q[3];
rz(1.4958517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742663) q[0];
sx q[0];
rz(-2.3889611) q[0];
sx q[0];
rz(2.0972032) q[0];
rz(-2.3692865) q[1];
sx q[1];
rz(-2.4817395) q[1];
sx q[1];
rz(-0.73371249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6873467) q[0];
sx q[0];
rz(-2.2139619) q[0];
sx q[0];
rz(2.1976794) q[0];
rz(1.1999424) q[2];
sx q[2];
rz(-2.7463253) q[2];
sx q[2];
rz(2.7140009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0391239) q[1];
sx q[1];
rz(-0.3141292) q[1];
sx q[1];
rz(0.38508319) q[1];
rz(-pi) q[2];
rz(0.10802631) q[3];
sx q[3];
rz(-2.3026534) q[3];
sx q[3];
rz(0.15209231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5250728) q[2];
sx q[2];
rz(-0.51023444) q[2];
sx q[2];
rz(0.60866848) q[2];
rz(-2.0742553) q[3];
sx q[3];
rz(-2.267024) q[3];
sx q[3];
rz(-0.55060351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50255018) q[0];
sx q[0];
rz(-2.783343) q[0];
sx q[0];
rz(1.6453561) q[0];
rz(-1.7773588) q[1];
sx q[1];
rz(-0.61449209) q[1];
sx q[1];
rz(-2.7552557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2992512) q[0];
sx q[0];
rz(-2.4066397) q[0];
sx q[0];
rz(-0.39836653) q[0];
x q[1];
rz(-0.069938439) q[2];
sx q[2];
rz(-1.3856263) q[2];
sx q[2];
rz(-2.1500146) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.98294241) q[1];
sx q[1];
rz(-1.6931931) q[1];
sx q[1];
rz(0.33315181) q[1];
rz(1.6101077) q[3];
sx q[3];
rz(-2.1609801) q[3];
sx q[3];
rz(-0.055824669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6284457) q[2];
sx q[2];
rz(-1.3725504) q[2];
sx q[2];
rz(-2.9885542) q[2];
rz(-2.2506574) q[3];
sx q[3];
rz(-0.95518437) q[3];
sx q[3];
rz(2.4002767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.9855758) q[0];
sx q[0];
rz(-2.9122536) q[0];
sx q[0];
rz(-0.34737059) q[0];
rz(-0.037847606) q[1];
sx q[1];
rz(-1.1696576) q[1];
sx q[1];
rz(0.52245021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.396401) q[0];
sx q[0];
rz(-2.7184882) q[0];
sx q[0];
rz(1.3518831) q[0];
rz(-pi) q[1];
rz(3.1407479) q[2];
sx q[2];
rz(-2.7760996) q[2];
sx q[2];
rz(1.5216684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6647937) q[1];
sx q[1];
rz(-2.0719686) q[1];
sx q[1];
rz(-1.0887926) q[1];
x q[2];
rz(2.7344875) q[3];
sx q[3];
rz(-1.0488172) q[3];
sx q[3];
rz(-0.84195053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5106421) q[2];
sx q[2];
rz(-2.6783671) q[2];
sx q[2];
rz(1.9412712) q[2];
rz(2.5837894) q[3];
sx q[3];
rz(-1.6226945) q[3];
sx q[3];
rz(0.38448486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2931622) q[0];
sx q[0];
rz(-2.1577305) q[0];
sx q[0];
rz(-1.4137319) q[0];
rz(0.14104715) q[1];
sx q[1];
rz(-1.3755362) q[1];
sx q[1];
rz(-2.486855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4439693) q[0];
sx q[0];
rz(-1.807741) q[0];
sx q[0];
rz(0.21392845) q[0];
x q[1];
rz(-2.4528772) q[2];
sx q[2];
rz(-2.42197) q[2];
sx q[2];
rz(1.0154427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6605986) q[1];
sx q[1];
rz(-0.29876041) q[1];
sx q[1];
rz(-2.1887652) q[1];
x q[2];
rz(0.61956866) q[3];
sx q[3];
rz(-1.1244785) q[3];
sx q[3];
rz(-2.5321162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35141382) q[2];
sx q[2];
rz(-1.979579) q[2];
sx q[2];
rz(0.72511017) q[2];
rz(2.2509947) q[3];
sx q[3];
rz(-2.2913439) q[3];
sx q[3];
rz(-0.58722535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17726041) q[0];
sx q[0];
rz(-1.5252508) q[0];
sx q[0];
rz(-1.9510212) q[0];
rz(1.1217077) q[1];
sx q[1];
rz(-1.5737166) q[1];
sx q[1];
rz(2.166688) q[1];
rz(0.45050889) q[2];
sx q[2];
rz(-2.4263739) q[2];
sx q[2];
rz(0.77745773) q[2];
rz(-1.7305456) q[3];
sx q[3];
rz(-0.86426576) q[3];
sx q[3];
rz(-0.15440253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
