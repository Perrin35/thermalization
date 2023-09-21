OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5117447) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(-1.758979) q[0];
x q[1];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(1.7125318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.558555) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(1.0298883) q[1];
rz(-pi) q[2];
rz(-1.3052985) q[3];
sx q[3];
rz(-1.7146829) q[3];
sx q[3];
rz(-0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-0.077117292) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.6024626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3829271) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(1.4325607) q[0];
rz(-2.402926) q[2];
sx q[2];
rz(-1.4091485) q[2];
sx q[2];
rz(2.6612298) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1748845) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(-2.944988) q[1];
x q[2];
rz(-2.7178571) q[3];
sx q[3];
rz(-1.3556004) q[3];
sx q[3];
rz(-2.7620897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(-0.69491274) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-2.7084896) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(2.5862397) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1486737) q[0];
sx q[0];
rz(-2.2745471) q[0];
sx q[0];
rz(2.2253195) q[0];
rz(-pi) q[1];
rz(2.136134) q[2];
sx q[2];
rz(-1.3034504) q[2];
sx q[2];
rz(-0.29758673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8086116) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(1.8379184) q[1];
rz(-pi) q[2];
rz(1.2461353) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(-0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.2505442) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53084757) q[0];
sx q[0];
rz(-1.8662211) q[0];
sx q[0];
rz(2.7644964) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43222506) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(1.150711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7123588) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(0.952094) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0905686) q[3];
sx q[3];
rz(-2.0912366) q[3];
sx q[3];
rz(-0.40363064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.3758855) q[0];
rz(1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8431906) q[0];
sx q[0];
rz(-1.0023596) q[0];
sx q[0];
rz(-2.0060904) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34727879) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(-1.4594644) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5406815) q[1];
sx q[1];
rz(-1.2680506) q[1];
sx q[1];
rz(1.5555698) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0104996) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(1.8492941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4650824) q[0];
sx q[0];
rz(-2.0970779) q[0];
sx q[0];
rz(0.17980534) q[0];
x q[1];
rz(0.94021057) q[2];
sx q[2];
rz(-1.5988837) q[2];
sx q[2];
rz(1.1379776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8050025) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(0.57979433) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8618705) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(2.1978956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(0.61002237) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726495) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(2.6555496) q[0];
rz(-0.91673135) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(0.23840657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94177946) q[1];
sx q[1];
rz(-2.4738414) q[1];
sx q[1];
rz(1.6112531) q[1];
rz(-2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-2.4050074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026222762) q[0];
sx q[0];
rz(-1.6008458) q[0];
sx q[0];
rz(-3.1307334) q[0];
rz(-pi) q[1];
rz(-1.893115) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(-1.581574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0524307) q[1];
sx q[1];
rz(-0.4757291) q[1];
sx q[1];
rz(-2.9176941) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8183476) q[3];
sx q[3];
rz(-1.9976227) q[3];
sx q[3];
rz(-2.1904898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(-1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(2.396778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841524) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(1.1293344) q[0];
rz(-0.3955598) q[2];
sx q[2];
rz(-0.8330847) q[2];
sx q[2];
rz(-1.2566483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71036584) q[1];
sx q[1];
rz(-1.2214298) q[1];
sx q[1];
rz(-2.9037895) q[1];
rz(-0.78299384) q[3];
sx q[3];
rz(-1.8564312) q[3];
sx q[3];
rz(-2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(0.40400305) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.9706479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5924581) q[0];
sx q[0];
rz(-0.42692703) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8659586) q[2];
sx q[2];
rz(-2.3466957) q[2];
sx q[2];
rz(2.8385712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78632894) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(1.5481871) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5578299) q[3];
sx q[3];
rz(-2.6423892) q[3];
sx q[3];
rz(-1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(-0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-2.4377433) q[2];
sx q[2];
rz(-0.8740295) q[2];
sx q[2];
rz(-0.340273) q[2];
rz(-3.1254461) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];