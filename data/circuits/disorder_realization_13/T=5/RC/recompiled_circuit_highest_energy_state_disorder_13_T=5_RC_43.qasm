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
rz(0.19652551) q[0];
sx q[0];
rz(-0.82634574) q[0];
sx q[0];
rz(0.91808051) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(-1.5674051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.471431) q[0];
sx q[0];
rz(-1.5613371) q[0];
sx q[0];
rz(0.0034250101) q[0];
rz(-1.9419045) q[2];
sx q[2];
rz(-2.1016309) q[2];
sx q[2];
rz(-2.4715854) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.543055) q[1];
sx q[1];
rz(-2.8208113) q[1];
sx q[1];
rz(-1.2307274) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2229177) q[3];
sx q[3];
rz(-2.6379728) q[3];
sx q[3];
rz(-0.38583392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.427318) q[2];
sx q[2];
rz(-1.486342) q[2];
sx q[2];
rz(-1.2662668) q[2];
rz(-2.0990939) q[3];
sx q[3];
rz(-1.6008585) q[3];
sx q[3];
rz(-1.3955759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0365486) q[0];
sx q[0];
rz(-0.62750134) q[0];
sx q[0];
rz(0.096916048) q[0];
rz(-2.3333683) q[1];
sx q[1];
rz(-0.40882912) q[1];
sx q[1];
rz(-1.2867297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1306123) q[0];
sx q[0];
rz(-1.6063147) q[0];
sx q[0];
rz(1.1208486) q[0];
rz(-pi) q[1];
rz(-2.419728) q[2];
sx q[2];
rz(-2.0316191) q[2];
sx q[2];
rz(-0.45022717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0285312) q[1];
sx q[1];
rz(-0.82958991) q[1];
sx q[1];
rz(2.100551) q[1];
rz(-2.8032885) q[3];
sx q[3];
rz(-1.2354992) q[3];
sx q[3];
rz(1.1962868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6404932) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(-2.7596149) q[2];
rz(2.6169418) q[3];
sx q[3];
rz(-1.4068406) q[3];
sx q[3];
rz(-2.9978571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76661888) q[0];
sx q[0];
rz(-1.5559649) q[0];
sx q[0];
rz(2.4554456) q[0];
rz(1.067591) q[1];
sx q[1];
rz(-0.99718863) q[1];
sx q[1];
rz(3.0608665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7537621) q[0];
sx q[0];
rz(-1.7942611) q[0];
sx q[0];
rz(0.039964393) q[0];
rz(2.4915495) q[2];
sx q[2];
rz(-2.1116684) q[2];
sx q[2];
rz(-2.385187) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0559433) q[1];
sx q[1];
rz(-2.9102444) q[1];
sx q[1];
rz(-3.0401708) q[1];
rz(-pi) q[2];
rz(2.3537797) q[3];
sx q[3];
rz(-1.9160877) q[3];
sx q[3];
rz(-2.6097176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2637691) q[2];
sx q[2];
rz(-2.4407385) q[2];
sx q[2];
rz(2.7027255) q[2];
rz(1.2980488) q[3];
sx q[3];
rz(-2.7545007) q[3];
sx q[3];
rz(2.7264061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2271659) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(0.67657226) q[0];
rz(3.0034972) q[1];
sx q[1];
rz(-1.0905677) q[1];
sx q[1];
rz(0.20060435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56708589) q[0];
sx q[0];
rz(-0.73154035) q[0];
sx q[0];
rz(2.125243) q[0];
rz(-pi) q[1];
rz(2.5121763) q[2];
sx q[2];
rz(-1.417629) q[2];
sx q[2];
rz(1.5063733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.323206) q[1];
sx q[1];
rz(-0.96508677) q[1];
sx q[1];
rz(-1.731864) q[1];
rz(-pi) q[2];
rz(-2.58039) q[3];
sx q[3];
rz(-1.4491014) q[3];
sx q[3];
rz(2.7336371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5033919) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(-1.7207883) q[2];
rz(-0.11451379) q[3];
sx q[3];
rz(-1.4736466) q[3];
sx q[3];
rz(0.99944559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223709) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(2.1694699) q[0];
rz(2.8179893) q[1];
sx q[1];
rz(-1.690003) q[1];
sx q[1];
rz(1.9591029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2401449) q[0];
sx q[0];
rz(-1.6109544) q[0];
sx q[0];
rz(-1.3754649) q[0];
rz(-pi) q[1];
rz(-1.576831) q[2];
sx q[2];
rz(-2.1409935) q[2];
sx q[2];
rz(-1.6369199) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3283493) q[1];
sx q[1];
rz(-2.9651838) q[1];
sx q[1];
rz(-1.117068) q[1];
rz(-pi) q[2];
rz(2.0911522) q[3];
sx q[3];
rz(-1.0624059) q[3];
sx q[3];
rz(-0.45724487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9779196) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(-2.6046806) q[2];
rz(-2.4162857) q[3];
sx q[3];
rz(-3.091843) q[3];
sx q[3];
rz(2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.86730114) q[0];
sx q[0];
rz(-1.0649571) q[0];
sx q[0];
rz(-1.4428447) q[0];
rz(-2.9310215) q[1];
sx q[1];
rz(-1.4434394) q[1];
sx q[1];
rz(2.494536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7919901) q[0];
sx q[0];
rz(-0.46030435) q[0];
sx q[0];
rz(-0.60636284) q[0];
rz(-2.5850614) q[2];
sx q[2];
rz(-1.663563) q[2];
sx q[2];
rz(-1.1586729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14110672) q[1];
sx q[1];
rz(-0.55505156) q[1];
sx q[1];
rz(-2.4480113) q[1];
rz(1.296455) q[3];
sx q[3];
rz(-1.2452494) q[3];
sx q[3];
rz(-0.59301585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1232542) q[2];
sx q[2];
rz(-1.0945357) q[2];
sx q[2];
rz(-1.9617762) q[2];
rz(1.2849464) q[3];
sx q[3];
rz(-2.732087) q[3];
sx q[3];
rz(3.0783317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061851) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(0.89757288) q[0];
rz(1.3342185) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(-0.99475494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332116) q[0];
sx q[0];
rz(-1.8803839) q[0];
sx q[0];
rz(0.82010834) q[0];
rz(-pi) q[1];
rz(0.9681692) q[2];
sx q[2];
rz(-1.954284) q[2];
sx q[2];
rz(2.9619975) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23235944) q[1];
sx q[1];
rz(-0.5482175) q[1];
sx q[1];
rz(2.1562804) q[1];
x q[2];
rz(-3.1174421) q[3];
sx q[3];
rz(-1.8514048) q[3];
sx q[3];
rz(2.5606972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.083954088) q[2];
sx q[2];
rz(-1.8123764) q[2];
sx q[2];
rz(-2.8482021) q[2];
rz(-0.053226274) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(-1.7085541) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5695802) q[0];
sx q[0];
rz(-1.3895637) q[0];
sx q[0];
rz(-2.8386175) q[0];
rz(1.6527893) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(0.66910076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55941033) q[0];
sx q[0];
rz(-0.56920496) q[0];
sx q[0];
rz(-2.4962382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9225408) q[2];
sx q[2];
rz(-1.7985849) q[2];
sx q[2];
rz(0.57727376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7380306) q[1];
sx q[1];
rz(-2.4564017) q[1];
sx q[1];
rz(0.86802796) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8101997) q[3];
sx q[3];
rz(-1.9839915) q[3];
sx q[3];
rz(-1.5904782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9992493) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(0.86714253) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(1.1423906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.860723) q[0];
sx q[0];
rz(-2.6779802) q[0];
sx q[0];
rz(0.11775693) q[0];
rz(-0.87751687) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(-2.1868475) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4539826) q[0];
sx q[0];
rz(-1.9508525) q[0];
sx q[0];
rz(-0.65681547) q[0];
x q[1];
rz(2.0659201) q[2];
sx q[2];
rz(-1.9946163) q[2];
sx q[2];
rz(-2.0002805) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92312557) q[1];
sx q[1];
rz(-2.4187741) q[1];
sx q[1];
rz(-1.9002401) q[1];
rz(0.68902632) q[3];
sx q[3];
rz(-0.9456767) q[3];
sx q[3];
rz(0.90009119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8942326) q[2];
sx q[2];
rz(-0.33984137) q[2];
sx q[2];
rz(1.4840688) q[2];
rz(-1.5225211) q[3];
sx q[3];
rz(-1.3831474) q[3];
sx q[3];
rz(1.1744261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5840983) q[0];
sx q[0];
rz(-1.431594) q[0];
sx q[0];
rz(0.75310055) q[0];
rz(-1.151471) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(-1.6023191) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879047) q[0];
sx q[0];
rz(-0.48948727) q[0];
sx q[0];
rz(1.1734084) q[0];
rz(2.1588155) q[2];
sx q[2];
rz(-0.78938198) q[2];
sx q[2];
rz(-0.090580926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6021612) q[1];
sx q[1];
rz(-1.185158) q[1];
sx q[1];
rz(-1.7785317) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33892858) q[3];
sx q[3];
rz(-2.458161) q[3];
sx q[3];
rz(-0.11250699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3573542) q[2];
sx q[2];
rz(-2.9604762) q[2];
sx q[2];
rz(0.46585807) q[2];
rz(-1.0957796) q[3];
sx q[3];
rz(-2.1491137) q[3];
sx q[3];
rz(2.5994658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096238484) q[0];
sx q[0];
rz(-1.5662554) q[0];
sx q[0];
rz(1.5731496) q[0];
rz(1.0678328) q[1];
sx q[1];
rz(-0.44034958) q[1];
sx q[1];
rz(-0.82028295) q[1];
rz(2.8779262) q[2];
sx q[2];
rz(-1.34524) q[2];
sx q[2];
rz(-3.0312579) q[2];
rz(-2.3798306) q[3];
sx q[3];
rz(-2.4488425) q[3];
sx q[3];
rz(-0.52656534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
