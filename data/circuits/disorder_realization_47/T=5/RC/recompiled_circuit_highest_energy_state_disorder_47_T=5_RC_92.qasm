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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(1.6562847) q[0];
rz(-2.0909042) q[1];
sx q[1];
rz(-1.3991855) q[1];
sx q[1];
rz(-2.4032226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6577539) q[0];
sx q[0];
rz(-1.8138348) q[0];
sx q[0];
rz(1.7263822) q[0];
rz(-0.15282571) q[2];
sx q[2];
rz(-1.9597561) q[2];
sx q[2];
rz(2.9940384) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0289606) q[1];
sx q[1];
rz(-2.6765209) q[1];
sx q[1];
rz(0.3151703) q[1];
rz(-0.41599993) q[3];
sx q[3];
rz(-0.66222727) q[3];
sx q[3];
rz(1.5755788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(0.4134678) q[2];
rz(-0.56420285) q[3];
sx q[3];
rz(-1.4671624) q[3];
sx q[3];
rz(1.1845425) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7269932) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(-3.0376036) q[0];
rz(0.14101401) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-0.41735059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0275729) q[0];
sx q[0];
rz(-2.0947523) q[0];
sx q[0];
rz(-0.34932507) q[0];
rz(-pi) q[1];
rz(3.1293489) q[2];
sx q[2];
rz(-1.2070388) q[2];
sx q[2];
rz(-0.084244339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88637832) q[1];
sx q[1];
rz(-2.8357949) q[1];
sx q[1];
rz(1.4319819) q[1];
rz(-0.82666918) q[3];
sx q[3];
rz(-2.429768) q[3];
sx q[3];
rz(-0.28361187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.523681) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(0.20393142) q[2];
rz(1.5150874) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(0.55346742) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7837778) q[0];
sx q[0];
rz(-0.73448616) q[0];
sx q[0];
rz(0.015588394) q[0];
x q[1];
rz(1.7856423) q[2];
sx q[2];
rz(-0.87085491) q[2];
sx q[2];
rz(-2.6756949) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7872214) q[1];
sx q[1];
rz(-2.1408892) q[1];
sx q[1];
rz(1.1524423) q[1];
rz(0.74456498) q[3];
sx q[3];
rz(-0.25601632) q[3];
sx q[3];
rz(0.42041966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.074097721) q[2];
sx q[2];
rz(-2.7390538) q[2];
sx q[2];
rz(1.2158166) q[2];
rz(-2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(-1.2522987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2860586) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(-1.704294) q[0];
rz(1.8057436) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(2.6864973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2360596) q[0];
sx q[0];
rz(-0.62232557) q[0];
sx q[0];
rz(-0.13747352) q[0];
rz(-pi) q[1];
rz(-1.1983775) q[2];
sx q[2];
rz(-1.3758059) q[2];
sx q[2];
rz(0.41620987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0008853) q[1];
sx q[1];
rz(-1.8155087) q[1];
sx q[1];
rz(-0.34414704) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7198445) q[3];
sx q[3];
rz(-2.9136326) q[3];
sx q[3];
rz(2.3123119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(0.8320128) q[2];
rz(2.4882107) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(-0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117383) q[0];
sx q[0];
rz(-1.8269704) q[0];
sx q[0];
rz(-2.4526556) q[0];
rz(-0.2116994) q[1];
sx q[1];
rz(-1.7023106) q[1];
sx q[1];
rz(2.6874218) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70063299) q[0];
sx q[0];
rz(-1.129375) q[0];
sx q[0];
rz(-0.9671797) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0549704) q[2];
sx q[2];
rz(-0.6956898) q[2];
sx q[2];
rz(1.5351968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8555357) q[1];
sx q[1];
rz(-2.2814676) q[1];
sx q[1];
rz(0.90531207) q[1];
x q[2];
rz(-1.0096142) q[3];
sx q[3];
rz(-2.8343763) q[3];
sx q[3];
rz(3.0209783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6058558) q[2];
sx q[2];
rz(-1.7822632) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(1.9262675) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(-2.5069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257492) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(-0.51344839) q[0];
rz(-2.2250941) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(1.7880012) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90474709) q[0];
sx q[0];
rz(-1.6535656) q[0];
sx q[0];
rz(-2.6968357) q[0];
rz(-2.2960179) q[2];
sx q[2];
rz(-2.6161199) q[2];
sx q[2];
rz(0.89286823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1648851) q[1];
sx q[1];
rz(-1.083255) q[1];
sx q[1];
rz(-1.6009925) q[1];
x q[2];
rz(0.94025375) q[3];
sx q[3];
rz(-0.11001405) q[3];
sx q[3];
rz(2.2785254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6650271) q[2];
sx q[2];
rz(-2.5883784) q[2];
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
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018709239) q[0];
sx q[0];
rz(-1.9174734) q[0];
sx q[0];
rz(-1.9710185) q[0];
rz(2.9508044) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(2.4986787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87827728) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(-0.10590597) q[0];
x q[1];
rz(0.25919886) q[2];
sx q[2];
rz(-1.890823) q[2];
sx q[2];
rz(2.7519403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2126522) q[1];
sx q[1];
rz(-1.8096605) q[1];
sx q[1];
rz(-2.8293306) q[1];
rz(-pi) q[2];
rz(-2.7369381) q[3];
sx q[3];
rz(-1.6804365) q[3];
sx q[3];
rz(3.0325923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.99870318) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(-1.9996803) q[2];
rz(-3.111908) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(-2.3619385) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90758816) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(1.3080066) q[0];
rz(-2.0246778) q[1];
sx q[1];
rz(-1.4968548) q[1];
sx q[1];
rz(0.21211472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19973913) q[0];
sx q[0];
rz(-0.77773577) q[0];
sx q[0];
rz(1.1129802) q[0];
x q[1];
rz(1.3364974) q[2];
sx q[2];
rz(-1.8673225) q[2];
sx q[2];
rz(1.4594452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37112889) q[1];
sx q[1];
rz(-2.0088638) q[1];
sx q[1];
rz(-0.99645946) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1954631) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(1.2809666) q[2];
rz(-1.403275) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(0.72428552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.69287777) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(-0.88248673) q[0];
rz(-2.8485883) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(-2.015347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4868797) q[0];
sx q[0];
rz(-2.3695787) q[0];
sx q[0];
rz(-2.604542) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59992591) q[2];
sx q[2];
rz(-1.5675822) q[2];
sx q[2];
rz(-0.41157882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4418151) q[1];
sx q[1];
rz(-2.338577) q[1];
sx q[1];
rz(-3.0670428) q[1];
rz(-pi) q[2];
rz(-2.4902017) q[3];
sx q[3];
rz(-0.49534251) q[3];
sx q[3];
rz(1.1817396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5249411) q[2];
sx q[2];
rz(-2.5280819) q[2];
sx q[2];
rz(-0.16723995) q[2];
rz(-3.1104769) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(0.64605609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6328218) q[0];
sx q[0];
rz(-1.333586) q[0];
sx q[0];
rz(-2.3311145) q[0];
rz(-1.8327911) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(3.0442309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705436) q[0];
sx q[0];
rz(-0.80035058) q[0];
sx q[0];
rz(-0.23231028) q[0];
rz(3.112744) q[2];
sx q[2];
rz(-1.1897773) q[2];
sx q[2];
rz(2.888607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2557981) q[1];
sx q[1];
rz(-0.58572873) q[1];
sx q[1];
rz(-0.08759193) q[1];
rz(-pi) q[2];
rz(-1.8382357) q[3];
sx q[3];
rz(-1.6839661) q[3];
sx q[3];
rz(-1.3133996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29791609) q[2];
sx q[2];
rz(-2.127779) q[2];
sx q[2];
rz(-1.2751014) q[2];
rz(-2.7013333) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(2.8662203) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(2.5195925) q[1];
sx q[1];
rz(-1.2983464) q[1];
sx q[1];
rz(1.3141528) q[1];
rz(0.95876454) q[2];
sx q[2];
rz(-3.0059881) q[2];
sx q[2];
rz(-0.81632951) q[2];
rz(-1.7429327) q[3];
sx q[3];
rz(-0.93440104) q[3];
sx q[3];
rz(1.222119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
