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
rz(1.037984) q[0];
sx q[0];
rz(4.8259566) q[0];
sx q[0];
rz(11.2136) q[0];
rz(-5.0985131) q[1];
sx q[1];
rz(2.4689622) q[1];
sx q[1];
rz(11.050635) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2555849) q[0];
sx q[0];
rz(-2.8681462) q[0];
sx q[0];
rz(-2.3257564) q[0];
rz(-pi) q[1];
rz(-2.0667384) q[2];
sx q[2];
rz(-2.3124394) q[2];
sx q[2];
rz(3.0665891) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7616854) q[1];
sx q[1];
rz(-0.60098472) q[1];
sx q[1];
rz(-2.4883449) q[1];
rz(-pi) q[2];
rz(1.2836841) q[3];
sx q[3];
rz(-1.1653882) q[3];
sx q[3];
rz(1.4472345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1379913) q[2];
sx q[2];
rz(-3.0333952) q[2];
sx q[2];
rz(2.9555964) q[2];
rz(3.1229535) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(2.0712461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828534) q[0];
sx q[0];
rz(-0.56621972) q[0];
sx q[0];
rz(-0.3183611) q[0];
rz(2.5045577) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(0.46924082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56457389) q[0];
sx q[0];
rz(-1.3718318) q[0];
sx q[0];
rz(1.3964723) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3606684) q[2];
sx q[2];
rz(-2.7581425) q[2];
sx q[2];
rz(2.0166778) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9904203) q[1];
sx q[1];
rz(-1.797533) q[1];
sx q[1];
rz(-2.8137035) q[1];
rz(-pi) q[2];
rz(-0.0057159609) q[3];
sx q[3];
rz(-2.2607231) q[3];
sx q[3];
rz(-1.6315414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8239173) q[2];
sx q[2];
rz(-2.9313512) q[2];
sx q[2];
rz(-2.4483185) q[2];
rz(-1.8215826) q[3];
sx q[3];
rz(-1.8157248) q[3];
sx q[3];
rz(2.0461953) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005118) q[0];
sx q[0];
rz(-2.5040099) q[0];
sx q[0];
rz(-2.6620423) q[0];
rz(-0.50411049) q[1];
sx q[1];
rz(-1.1481552) q[1];
sx q[1];
rz(0.058301059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8970015) q[0];
sx q[0];
rz(-1.0718597) q[0];
sx q[0];
rz(1.8656436) q[0];
x q[1];
rz(2.451926) q[2];
sx q[2];
rz(-1.153077) q[2];
sx q[2];
rz(2.6725685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1830777) q[1];
sx q[1];
rz(-0.56499186) q[1];
sx q[1];
rz(-1.3426258) q[1];
rz(-pi) q[2];
rz(-0.4543484) q[3];
sx q[3];
rz(-0.55255167) q[3];
sx q[3];
rz(2.1228028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7551859) q[2];
sx q[2];
rz(-1.9599954) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(1.5197598) q[3];
sx q[3];
rz(-2.683679) q[3];
sx q[3];
rz(0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195049) q[0];
sx q[0];
rz(-2.6730972) q[0];
sx q[0];
rz(3.0750437) q[0];
rz(-1.5244124) q[1];
sx q[1];
rz(-1.2434554) q[1];
sx q[1];
rz(0.7349416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733058) q[0];
sx q[0];
rz(-2.5096008) q[0];
sx q[0];
rz(-2.7422264) q[0];
rz(-pi) q[1];
rz(-2.4114716) q[2];
sx q[2];
rz(-1.3956304) q[2];
sx q[2];
rz(-2.1736682) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7801223) q[1];
sx q[1];
rz(-1.2166942) q[1];
sx q[1];
rz(1.2722871) q[1];
x q[2];
rz(1.4052584) q[3];
sx q[3];
rz(-0.72742262) q[3];
sx q[3];
rz(-0.35044985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8526326) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(-3.0328879) q[2];
rz(0.92681256) q[3];
sx q[3];
rz(-2.1709397) q[3];
sx q[3];
rz(-1.577781) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0573334) q[0];
sx q[0];
rz(-0.63705343) q[0];
sx q[0];
rz(-1.0908701) q[0];
rz(0.57251656) q[1];
sx q[1];
rz(-1.9520091) q[1];
sx q[1];
rz(2.3172839) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59412724) q[0];
sx q[0];
rz(-2.6931433) q[0];
sx q[0];
rz(-2.1350103) q[0];
rz(-pi) q[1];
rz(-0.90472533) q[2];
sx q[2];
rz(-0.76132971) q[2];
sx q[2];
rz(0.29354087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.854519) q[1];
sx q[1];
rz(-1.9654915) q[1];
sx q[1];
rz(-0.97868195) q[1];
rz(-pi) q[2];
rz(-0.40559988) q[3];
sx q[3];
rz(-1.1685017) q[3];
sx q[3];
rz(1.9674439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.480964) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(2.8503897) q[2];
rz(-0.63699841) q[3];
sx q[3];
rz(-1.4642508) q[3];
sx q[3];
rz(-1.252482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7630735) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(3.0105403) q[0];
rz(1.9746926) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(-0.022620591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2289091) q[0];
sx q[0];
rz(-2.1118409) q[0];
sx q[0];
rz(-1.156927) q[0];
rz(0.51949595) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(1.0102538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80161051) q[1];
sx q[1];
rz(-1.2622854) q[1];
sx q[1];
rz(-2.1478081) q[1];
rz(-pi) q[2];
rz(0.87225391) q[3];
sx q[3];
rz(-2.8917851) q[3];
sx q[3];
rz(1.2605309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8733946) q[2];
sx q[2];
rz(-2.4384273) q[2];
sx q[2];
rz(-2.0568636) q[2];
rz(-1.1449413) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(2.1196608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5522083) q[0];
sx q[0];
rz(-1.6227868) q[0];
sx q[0];
rz(-2.6680706) q[0];
rz(-1.2208968) q[1];
sx q[1];
rz(-2.7291606) q[1];
sx q[1];
rz(1.268505) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.761577) q[0];
sx q[0];
rz(-2.3748739) q[0];
sx q[0];
rz(0.79464998) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73641554) q[2];
sx q[2];
rz(-2.066095) q[2];
sx q[2];
rz(0.30799949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.699325) q[1];
sx q[1];
rz(-1.066477) q[1];
sx q[1];
rz(2.0513351) q[1];
rz(1.7470737) q[3];
sx q[3];
rz(-2.699614) q[3];
sx q[3];
rz(-2.7689742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1641757) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(0.24641985) q[2];
rz(2.1988846) q[3];
sx q[3];
rz(-1.0859414) q[3];
sx q[3];
rz(-0.76343083) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134269) q[0];
sx q[0];
rz(-1.270371) q[0];
sx q[0];
rz(-3.0810007) q[0];
rz(1.5024441) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(1.5477808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029744) q[0];
sx q[0];
rz(-1.5718223) q[0];
sx q[0];
rz(-0.0038504168) q[0];
x q[1];
rz(-0.28569371) q[2];
sx q[2];
rz(-2.5473352) q[2];
sx q[2];
rz(0.93462925) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8503849) q[1];
sx q[1];
rz(-0.7813079) q[1];
sx q[1];
rz(-2.2154097) q[1];
x q[2];
rz(0.5129359) q[3];
sx q[3];
rz(-0.78425927) q[3];
sx q[3];
rz(1.5363786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26555201) q[2];
sx q[2];
rz(-0.92090845) q[2];
sx q[2];
rz(-1.457224) q[2];
rz(3.0217116) q[3];
sx q[3];
rz(-2.9370152) q[3];
sx q[3];
rz(2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(-2.0062398) q[0];
rz(3.0382233) q[1];
sx q[1];
rz(-1.4048978) q[1];
sx q[1];
rz(0.9476544) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68610588) q[0];
sx q[0];
rz(-2.4999833) q[0];
sx q[0];
rz(-3.0296586) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88236188) q[2];
sx q[2];
rz(-1.244592) q[2];
sx q[2];
rz(0.38089124) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3909437) q[1];
sx q[1];
rz(-2.237197) q[1];
sx q[1];
rz(2.8007068) q[1];
rz(0.89449785) q[3];
sx q[3];
rz(-0.94493659) q[3];
sx q[3];
rz(-1.2813527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4244298) q[2];
sx q[2];
rz(-0.62817502) q[2];
sx q[2];
rz(0.90432811) q[2];
rz(1.8782015) q[3];
sx q[3];
rz(-1.9046013) q[3];
sx q[3];
rz(-2.625107) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(-0.98536056) q[0];
rz(2.789978) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(-2.7947289) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37736402) q[0];
sx q[0];
rz(-0.94616156) q[0];
sx q[0];
rz(-2.8714116) q[0];
rz(1.3386334) q[2];
sx q[2];
rz(-0.45660161) q[2];
sx q[2];
rz(-0.62397623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.092791768) q[1];
sx q[1];
rz(-1.4818703) q[1];
sx q[1];
rz(-0.6912749) q[1];
rz(-pi) q[2];
rz(0.68315398) q[3];
sx q[3];
rz(-2.5954163) q[3];
sx q[3];
rz(2.7462296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2615307) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(3.0066709) q[2];
rz(-3.0123582) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(1.038704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1226817) q[0];
sx q[0];
rz(-1.8862579) q[0];
sx q[0];
rz(0.12611623) q[0];
rz(2.6799754) q[1];
sx q[1];
rz(-1.7414265) q[1];
sx q[1];
rz(0.19207676) q[1];
rz(2.1264524) q[2];
sx q[2];
rz(-0.75426741) q[2];
sx q[2];
rz(-1.5744899) q[2];
rz(-2.8830388) q[3];
sx q[3];
rz(-1.7333442) q[3];
sx q[3];
rz(2.6556591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
