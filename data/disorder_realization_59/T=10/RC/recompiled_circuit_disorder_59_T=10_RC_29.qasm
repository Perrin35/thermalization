OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(-2.864569) q[1];
sx q[1];
rz(-2.6695873) q[1];
sx q[1];
rz(-0.0013874887) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.815925) q[0];
sx q[0];
rz(-0.40580931) q[0];
sx q[0];
rz(-1.9947467) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6127365) q[2];
sx q[2];
rz(-2.0176) q[2];
sx q[2];
rz(2.1697793) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5056155) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(-0.14640267) q[1];
x q[2];
rz(-1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(0.97066561) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(-0.81545365) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0464697) q[0];
sx q[0];
rz(-2.2182584) q[0];
sx q[0];
rz(1.6460653) q[0];
x q[1];
rz(-0.30351992) q[2];
sx q[2];
rz(-2.3887861) q[2];
sx q[2];
rz(0.34740651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7533469) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(1.6080329) q[1];
rz(-2.4957982) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(-1.4327232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(-2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(-1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-2.1420746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8529352) q[0];
sx q[0];
rz(-1.7203727) q[0];
sx q[0];
rz(-1.268671) q[0];
rz(2.524316) q[2];
sx q[2];
rz(-1.0059788) q[2];
sx q[2];
rz(1.9895983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1332902) q[1];
sx q[1];
rz(-1.8255207) q[1];
sx q[1];
rz(-0.50526527) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0619377) q[3];
sx q[3];
rz(-1.0632535) q[3];
sx q[3];
rz(1.5910651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(-2.3245658) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(0.91039175) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(2.8667563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1823605) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(1.4819281) q[0];
rz(1.574013) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(-0.38052961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7819314) q[1];
sx q[1];
rz(-1.1210103) q[1];
sx q[1];
rz(0.91314544) q[1];
x q[2];
rz(-0.12675385) q[3];
sx q[3];
rz(-2.2634014) q[3];
sx q[3];
rz(-0.17514378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.794902) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.516974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18407962) q[0];
sx q[0];
rz(-0.74680579) q[0];
sx q[0];
rz(1.0190796) q[0];
x q[1];
rz(-1.9618271) q[2];
sx q[2];
rz(-2.070825) q[2];
sx q[2];
rz(0.21991877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72488898) q[1];
sx q[1];
rz(-0.54425889) q[1];
sx q[1];
rz(-2.4328028) q[1];
rz(0.98832163) q[3];
sx q[3];
rz(-2.0867996) q[3];
sx q[3];
rz(0.9048681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8020442) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(-2.8175763) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76535392) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(2.8009159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6402123) q[0];
sx q[0];
rz(-3.0684154) q[0];
sx q[0];
rz(1.120938) q[0];
rz(-pi) q[1];
rz(0.31008115) q[2];
sx q[2];
rz(-2.3611464) q[2];
sx q[2];
rz(-2.9030637) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.016547116) q[1];
sx q[1];
rz(-0.53495896) q[1];
sx q[1];
rz(2.6782481) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7192781) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(-0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(0.34564885) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3195254) q[0];
sx q[0];
rz(-1.7626581) q[0];
sx q[0];
rz(1.096154) q[0];
rz(2.6480688) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(1.7871737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3073687) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(2.6878396) q[1];
x q[2];
rz(1.0602337) q[3];
sx q[3];
rz(-0.052882346) q[3];
sx q[3];
rz(-0.75170654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1376301) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(0.078358738) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5057482) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(-pi) q[1];
rz(1.3252844) q[2];
sx q[2];
rz(-1.8010555) q[2];
sx q[2];
rz(-0.85207176) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52787493) q[1];
sx q[1];
rz(-0.67786874) q[1];
sx q[1];
rz(-1.2916958) q[1];
x q[2];
rz(-1.3571635) q[3];
sx q[3];
rz(-0.16031081) q[3];
sx q[3];
rz(-0.3233288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-0.051368512) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-0.87402469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0747781) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(0.76807036) q[0];
rz(-pi) q[1];
rz(-0.36058493) q[2];
sx q[2];
rz(-1.1869831) q[2];
sx q[2];
rz(-1.4521445) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45174949) q[1];
sx q[1];
rz(-2.6362231) q[1];
sx q[1];
rz(-1.4228574) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5328818) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(-1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-0.61974636) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(-1.1788517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26319474) q[0];
sx q[0];
rz(-2.7976755) q[0];
sx q[0];
rz(-0.92349903) q[0];
x q[1];
rz(-2.5095021) q[2];
sx q[2];
rz(-0.1569911) q[2];
sx q[2];
rz(-1.1319515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19503838) q[1];
sx q[1];
rz(-0.95609162) q[1];
sx q[1];
rz(2.9123995) q[1];
rz(-pi) q[2];
rz(0.52272777) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(-2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.1307217) q[2];
sx q[2];
rz(-1.5390839) q[2];
sx q[2];
rz(2.5577953) q[2];
rz(0.27579565) q[3];
sx q[3];
rz(-1.4281359) q[3];
sx q[3];
rz(1.7208163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];