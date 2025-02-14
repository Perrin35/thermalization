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
rz(1.485308) q[0];
rz(-2.0909042) q[1];
sx q[1];
rz(-1.3991855) q[1];
sx q[1];
rz(-2.4032226) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12468689) q[0];
sx q[0];
rz(-1.4198185) q[0];
sx q[0];
rz(-0.24589234) q[0];
rz(-pi) q[1];
rz(-1.2151488) q[2];
sx q[2];
rz(-2.7251149) q[2];
sx q[2];
rz(2.6082325) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.112632) q[1];
sx q[1];
rz(-0.46507177) q[1];
sx q[1];
rz(2.8264224) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61951903) q[3];
sx q[3];
rz(-1.3196919) q[3];
sx q[3];
rz(-0.34003231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7585313) q[2];
sx q[2];
rz(-0.28306511) q[2];
sx q[2];
rz(-0.4134678) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41459945) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(-0.10398908) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-2.7242421) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7234897) q[0];
sx q[0];
rz(-1.8716629) q[0];
sx q[0];
rz(2.1221493) q[0];
rz(-pi) q[1];
rz(-1.9345788) q[2];
sx q[2];
rz(-1.5822389) q[2];
sx q[2];
rz(1.4909084) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0318533) q[1];
sx q[1];
rz(-1.8735587) q[1];
sx q[1];
rz(-0.043655386) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52882282) q[3];
sx q[3];
rz(-1.0695056) q[3];
sx q[3];
rz(-0.59893417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.523681) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(-2.9376612) q[2];
rz(1.6265053) q[3];
sx q[3];
rz(-2.6650059) q[3];
sx q[3];
rz(-1.509607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5390891) q[0];
sx q[0];
rz(-1.6747549) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(2.0643945) q[1];
sx q[1];
rz(-2.4927683) q[1];
sx q[1];
rz(-2.5881252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3368126) q[0];
sx q[0];
rz(-2.3051728) q[0];
sx q[0];
rz(1.5848716) q[0];
rz(-pi) q[1];
rz(-2.8936549) q[2];
sx q[2];
rz(-2.4147779) q[2];
sx q[2];
rz(-3.0023129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1225153) q[1];
sx q[1];
rz(-1.919849) q[1];
sx q[1];
rz(-0.61180361) q[1];
rz(2.9514246) q[3];
sx q[3];
rz(-1.7432508) q[3];
sx q[3];
rz(-1.8784539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0674949) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(1.925776) q[2];
rz(1.0007693) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(-1.2522987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2860586) q[0];
sx q[0];
rz(-1.0581886) q[0];
sx q[0];
rz(-1.4372987) q[0];
rz(-1.3358491) q[1];
sx q[1];
rz(-2.5156486) q[1];
sx q[1];
rz(2.6864973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.06741) q[0];
sx q[0];
rz(-2.1863696) q[0];
sx q[0];
rz(-1.4727946) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1983775) q[2];
sx q[2];
rz(-1.3758059) q[2];
sx q[2];
rz(0.41620987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3434717) q[1];
sx q[1];
rz(-1.9042854) q[1];
sx q[1];
rz(-1.8300959) q[1];
rz(-1.6654799) q[3];
sx q[3];
rz(-1.3631184) q[3];
sx q[3];
rz(-2.7439347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8670292) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(0.8320128) q[2];
rz(-2.4882107) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(-2.466989) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117383) q[0];
sx q[0];
rz(-1.8269704) q[0];
sx q[0];
rz(-2.4526556) q[0];
rz(0.2116994) q[1];
sx q[1];
rz(-1.7023106) q[1];
sx q[1];
rz(0.45417085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9849761) q[0];
sx q[0];
rz(-1.0318705) q[0];
sx q[0];
rz(-2.6205553) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0866223) q[2];
sx q[2];
rz(-0.6956898) q[2];
sx q[2];
rz(1.6063959) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.28605697) q[1];
sx q[1];
rz(-0.86012506) q[1];
sx q[1];
rz(2.2362806) q[1];
rz(-1.0096142) q[3];
sx q[3];
rz(-0.30721634) q[3];
sx q[3];
rz(0.12061435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6058558) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(-1.9262675) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(-0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257492) q[0];
sx q[0];
rz(-0.31084335) q[0];
sx q[0];
rz(0.51344839) q[0];
rz(-0.91649857) q[1];
sx q[1];
rz(-1.1767574) q[1];
sx q[1];
rz(-1.3535915) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83780629) q[0];
sx q[0];
rz(-0.45188658) q[0];
sx q[0];
rz(0.19048283) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2960179) q[2];
sx q[2];
rz(-0.52547272) q[2];
sx q[2];
rz(-0.89286823) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60823764) q[1];
sx q[1];
rz(-1.5441193) q[1];
sx q[1];
rz(0.4877301) q[1];
rz(-pi) q[2];
rz(-2.2013389) q[3];
sx q[3];
rz(-3.0315786) q[3];
sx q[3];
rz(0.86306727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6650271) q[2];
sx q[2];
rz(-0.55321425) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(2.636886) q[3];
sx q[3];
rz(-1.4387771) q[3];
sx q[3];
rz(-0.75468841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.1228834) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(1.9710185) q[0];
rz(2.9508044) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(2.4986787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64775) q[0];
sx q[0];
rz(-0.28573418) q[0];
sx q[0];
rz(1.9410551) q[0];
rz(-pi) q[1];
rz(0.25919886) q[2];
sx q[2];
rz(-1.2507696) q[2];
sx q[2];
rz(-2.7519403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8596803) q[1];
sx q[1];
rz(-1.2676867) q[1];
sx q[1];
rz(-1.3202841) q[1];
rz(-0.40465458) q[3];
sx q[3];
rz(-1.4611562) q[3];
sx q[3];
rz(3.0325923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(-1.1419123) q[2];
rz(-3.111908) q[3];
sx q[3];
rz(-1.6073062) q[3];
sx q[3];
rz(-0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90758816) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(1.3080066) q[0];
rz(1.1169149) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(2.9294779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19973913) q[0];
sx q[0];
rz(-2.3638569) q[0];
sx q[0];
rz(-2.0286125) q[0];
rz(-pi) q[1];
rz(-1.3364974) q[2];
sx q[2];
rz(-1.2742701) q[2];
sx q[2];
rz(1.4594452) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5221716) q[1];
sx q[1];
rz(-2.4345349) q[1];
sx q[1];
rz(-0.85925428) q[1];
rz(0.85101012) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(-1.2860166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1954631) q[2];
sx q[2];
rz(-1.3529494) q[2];
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
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69287777) q[0];
sx q[0];
rz(-2.4249478) q[0];
sx q[0];
rz(0.88248673) q[0];
rz(2.8485883) q[1];
sx q[1];
rz(-2.2182783) q[1];
sx q[1];
rz(-2.015347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348081) q[0];
sx q[0];
rz(-0.92807584) q[0];
sx q[0];
rz(-2.0329143) q[0];
x q[1];
rz(0.59992591) q[2];
sx q[2];
rz(-1.5675822) q[2];
sx q[2];
rz(0.41157882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3346904) q[1];
sx q[1];
rz(-2.370939) q[1];
sx q[1];
rz(1.6477963) q[1];
rz(-1.8873439) q[3];
sx q[3];
rz(-1.9584362) q[3];
sx q[3];
rz(2.6738808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61665159) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(0.16723995) q[2];
rz(-0.031115726) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(2.4955366) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50877082) q[0];
sx q[0];
rz(-1.333586) q[0];
sx q[0];
rz(0.81047812) q[0];
rz(1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(3.0442309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0705436) q[0];
sx q[0];
rz(-2.3412421) q[0];
sx q[0];
rz(0.23231028) q[0];
rz(-pi) q[1];
rz(-3.112744) q[2];
sx q[2];
rz(-1.9518153) q[2];
sx q[2];
rz(2.888607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75805) q[1];
sx q[1];
rz(-1.6191747) q[1];
sx q[1];
rz(0.58396062) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3033569) q[3];
sx q[3];
rz(-1.6839661) q[3];
sx q[3];
rz(-1.3133996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8436766) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(-1.2751014) q[2];
rz(0.44025931) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(2.8662203) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(0.62200017) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(0.078230187) q[2];
sx q[2];
rz(-1.6816734) q[2];
sx q[2];
rz(1.7088877) q[2];
rz(-0.6435271) q[3];
sx q[3];
rz(-1.7089927) q[3];
sx q[3];
rz(-0.24571936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
