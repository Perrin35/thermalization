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
rz(-1.4909622) q[0];
sx q[0];
rz(-1.7234252) q[0];
sx q[0];
rz(-1.6562847) q[0];
rz(1.0506884) q[1];
sx q[1];
rz(-1.7424072) q[1];
sx q[1];
rz(2.4032226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0169058) q[0];
sx q[0];
rz(-1.7217741) q[0];
sx q[0];
rz(-0.24589234) q[0];
x q[1];
rz(-0.15282571) q[2];
sx q[2];
rz(-1.9597561) q[2];
sx q[2];
rz(-0.14755421) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0289606) q[1];
sx q[1];
rz(-2.6765209) q[1];
sx q[1];
rz(-2.8264224) q[1];
rz(-pi) q[2];
rz(1.2655696) q[3];
sx q[3];
rz(-0.97346993) q[3];
sx q[3];
rz(-2.0862153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38306132) q[2];
sx q[2];
rz(-0.28306511) q[2];
sx q[2];
rz(-2.7281249) q[2];
rz(0.56420285) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(1.1845425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7269932) q[0];
sx q[0];
rz(-2.0818384) q[0];
sx q[0];
rz(0.10398908) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-2.7242421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39825059) q[0];
sx q[0];
rz(-2.5210327) q[0];
sx q[0];
rz(-2.1055566) q[0];
x q[1];
rz(0.012243791) q[2];
sx q[2];
rz(-1.2070388) q[2];
sx q[2];
rz(-3.0573483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5896259) q[1];
sx q[1];
rz(-1.6124649) q[1];
sx q[1];
rz(-1.2677626) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1363229) q[3];
sx q[3];
rz(-1.1124753) q[3];
sx q[3];
rz(1.895973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61791164) q[2];
sx q[2];
rz(-1.0289501) q[2];
sx q[2];
rz(2.9376612) q[2];
rz(-1.6265053) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6025036) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(-2.0643945) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(-2.5881252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3368126) q[0];
sx q[0];
rz(-0.83641988) q[0];
sx q[0];
rz(-1.5848716) q[0];
rz(-pi) q[1];
rz(-1.7856423) q[2];
sx q[2];
rz(-0.87085491) q[2];
sx q[2];
rz(2.6756949) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3543713) q[1];
sx q[1];
rz(-1.0007035) q[1];
sx q[1];
rz(1.9891504) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19016805) q[3];
sx q[3];
rz(-1.3983418) q[3];
sx q[3];
rz(1.8784539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0674949) q[2];
sx q[2];
rz(-2.7390538) q[2];
sx q[2];
rz(1.925776) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(-1.8892939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2860586) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(1.4372987) q[0];
rz(-1.8057436) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(-2.6864973) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44667654) q[0];
sx q[0];
rz(-1.4908264) q[0];
sx q[0];
rz(-2.5237501) q[0];
rz(-2.0680769) q[2];
sx q[2];
rz(-0.4182564) q[2];
sx q[2];
rz(0.69415316) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0008853) q[1];
sx q[1];
rz(-1.326084) q[1];
sx q[1];
rz(0.34414704) q[1];
rz(-pi) q[2];
rz(-0.20858553) q[3];
sx q[3];
rz(-1.6634395) q[3];
sx q[3];
rz(-1.1535597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2745634) q[2];
sx q[2];
rz(-2.6919591) q[2];
sx q[2];
rz(-2.3095798) q[2];
rz(-2.4882107) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(-2.466989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117383) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(2.4526556) q[0];
rz(-0.2116994) q[1];
sx q[1];
rz(-1.7023106) q[1];
sx q[1];
rz(2.6874218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70063299) q[0];
sx q[0];
rz(-2.0122177) q[0];
sx q[0];
rz(2.174413) q[0];
rz(0.93438678) q[2];
sx q[2];
rz(-1.2678522) q[2];
sx q[2];
rz(-2.7222939) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3836174) q[1];
sx q[1];
rz(-1.0838008) q[1];
sx q[1];
rz(2.3112554) q[1];
rz(2.1319785) q[3];
sx q[3];
rz(-2.8343763) q[3];
sx q[3];
rz(3.0209783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6058558) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(-0.5109171) q[2];
rz(-1.2153252) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(-2.5069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257492) q[0];
sx q[0];
rz(-0.31084335) q[0];
sx q[0];
rz(-2.6281443) q[0];
rz(-2.2250941) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(1.7880012) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5149263) q[0];
sx q[0];
rz(-2.0139222) q[0];
sx q[0];
rz(-1.4791545) q[0];
rz(-pi) q[1];
x q[1];
rz(2.774418) q[2];
sx q[2];
rz(-1.9556139) q[2];
sx q[2];
rz(3.0462616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.533355) q[1];
sx q[1];
rz(-1.5441193) q[1];
sx q[1];
rz(-2.6538626) q[1];
x q[2];
rz(2.2013389) q[3];
sx q[3];
rz(-3.0315786) q[3];
sx q[3];
rz(2.2785254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4765656) q[2];
sx q[2];
rz(-0.55321425) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(0.50470662) q[3];
sx q[3];
rz(-1.7028156) q[3];
sx q[3];
rz(-0.75468841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1228834) q[0];
sx q[0];
rz(-1.9174734) q[0];
sx q[0];
rz(-1.1705742) q[0];
rz(0.19078828) q[1];
sx q[1];
rz(-2.6759594) q[1];
sx q[1];
rz(-0.64291397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2633154) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(-3.0356867) q[0];
x q[1];
rz(0.91251164) q[2];
sx q[2];
rz(-0.40900074) q[2];
sx q[2];
rz(2.8308979) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2126522) q[1];
sx q[1];
rz(-1.3319322) q[1];
sx q[1];
rz(-2.8293306) q[1];
x q[2];
rz(-1.4516109) q[3];
sx q[3];
rz(-1.9728807) q[3];
sx q[3];
rz(-1.5086255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(1.9996803) q[2];
rz(-0.029684639) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(-0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90758816) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(1.3080066) q[0];
rz(2.0246778) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(-2.9294779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40529272) q[0];
sx q[0];
rz(-2.2515902) q[0];
sx q[0];
rz(-0.41053562) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8372357) q[2];
sx q[2];
rz(-1.7946912) q[2];
sx q[2];
rz(-0.18098132) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37112889) q[1];
sx q[1];
rz(-1.1327289) q[1];
sx q[1];
rz(0.99645946) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2905825) q[3];
sx q[3];
rz(-0.86736263) q[3];
sx q[3];
rz(1.855576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9461296) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(-1.8606261) q[2];
rz(1.7383176) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(-2.4173071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4487149) q[0];
sx q[0];
rz(-2.4249478) q[0];
sx q[0];
rz(-2.2591059) q[0];
rz(2.8485883) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(-1.1262456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6544271) q[0];
sx q[0];
rz(-1.9357271) q[0];
sx q[0];
rz(0.69661822) q[0];
x q[1];
rz(-1.5746905) q[2];
sx q[2];
rz(-0.97087395) q[2];
sx q[2];
rz(-1.1614161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3346904) q[1];
sx q[1];
rz(-0.77065361) q[1];
sx q[1];
rz(-1.4937964) q[1];
x q[2];
rz(-0.4057986) q[3];
sx q[3];
rz(-1.278459) q[3];
sx q[3];
rz(-2.1617011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5249411) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(-0.16723995) q[2];
rz(-0.031115726) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(2.4955366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50877082) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(2.3311145) q[0];
rz(-1.8327911) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(-0.097361758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3980078) q[0];
sx q[0];
rz(-0.79774374) q[0];
sx q[0];
rz(1.3378851) q[0];
x q[1];
rz(-1.6426769) q[2];
sx q[2];
rz(-0.38205636) q[2];
sx q[2];
rz(-2.9660513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88579455) q[1];
sx q[1];
rz(-0.58572873) q[1];
sx q[1];
rz(-0.08759193) q[1];
rz(-pi) q[2];
rz(-1.9769659) q[3];
sx q[3];
rz(-0.28986922) q[3];
sx q[3];
rz(2.4932662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8436766) q[2];
sx q[2];
rz(-2.127779) q[2];
sx q[2];
rz(1.8664912) q[2];
rz(-0.44025931) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(-2.8662203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.1210099) q[0];
sx q[0];
rz(-0.6211716) q[0];
sx q[0];
rz(2.4564263) q[0];
rz(-2.5195925) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(-2.1828281) q[2];
sx q[2];
rz(-3.0059881) q[2];
sx q[2];
rz(-0.81632951) q[2];
rz(-2.9138184) q[3];
sx q[3];
rz(-0.65613272) q[3];
sx q[3];
rz(-1.6348742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
