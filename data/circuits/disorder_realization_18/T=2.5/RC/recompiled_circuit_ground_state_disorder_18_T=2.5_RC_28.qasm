OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48489269) q[0];
sx q[0];
rz(-3.0527053) q[0];
sx q[0];
rz(0.77699295) q[0];
rz(-2.4465893) q[1];
sx q[1];
rz(-0.13091317) q[1];
sx q[1];
rz(-0.34832365) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3366813) q[0];
sx q[0];
rz(-2.3279193) q[0];
sx q[0];
rz(-1.8437852) q[0];
rz(-pi) q[1];
rz(0.25599249) q[2];
sx q[2];
rz(-0.29689327) q[2];
sx q[2];
rz(-0.96742899) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6887853) q[1];
sx q[1];
rz(-1.0579351) q[1];
sx q[1];
rz(0.60223438) q[1];
rz(-pi) q[2];
rz(1.5218842) q[3];
sx q[3];
rz(-1.9659316) q[3];
sx q[3];
rz(1.5103769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6473039) q[2];
sx q[2];
rz(-2.5755197) q[2];
sx q[2];
rz(-1.0757793) q[2];
rz(-3.0268269) q[3];
sx q[3];
rz(-1.4817295) q[3];
sx q[3];
rz(0.65576321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455416) q[0];
sx q[0];
rz(-0.93463722) q[0];
sx q[0];
rz(-0.70282394) q[0];
rz(-1.5952236) q[1];
sx q[1];
rz(-2.399235) q[1];
sx q[1];
rz(-2.5667618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0718577) q[0];
sx q[0];
rz(-1.1848579) q[0];
sx q[0];
rz(1.2161847) q[0];
rz(-pi) q[1];
rz(-2.0971336) q[2];
sx q[2];
rz(-2.1169726) q[2];
sx q[2];
rz(2.1072497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.489255) q[1];
sx q[1];
rz(-1.4397419) q[1];
sx q[1];
rz(-2.5136957) q[1];
x q[2];
rz(2.094926) q[3];
sx q[3];
rz(-1.0051749) q[3];
sx q[3];
rz(-3.0622967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90291643) q[2];
sx q[2];
rz(-1.5006337) q[2];
sx q[2];
rz(2.404876) q[2];
rz(-0.77218562) q[3];
sx q[3];
rz(-2.9530647) q[3];
sx q[3];
rz(0.34672117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18948983) q[0];
sx q[0];
rz(-2.0661856) q[0];
sx q[0];
rz(-1.9814251) q[0];
rz(0.53933764) q[1];
sx q[1];
rz(-0.62349206) q[1];
sx q[1];
rz(-3.0711874) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3905593) q[0];
sx q[0];
rz(-0.91827276) q[0];
sx q[0];
rz(-2.5120048) q[0];
x q[1];
rz(2.3023476) q[2];
sx q[2];
rz(-0.23942023) q[2];
sx q[2];
rz(0.024062238) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7097096) q[1];
sx q[1];
rz(-1.6358915) q[1];
sx q[1];
rz(0.38159926) q[1];
rz(0.60172867) q[3];
sx q[3];
rz(-2.3961107) q[3];
sx q[3];
rz(-1.1423542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77014852) q[2];
sx q[2];
rz(-1.4263209) q[2];
sx q[2];
rz(-0.16177978) q[2];
rz(-2.3700304) q[3];
sx q[3];
rz(-3.0050889) q[3];
sx q[3];
rz(1.0426883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55018392) q[0];
sx q[0];
rz(-0.78224459) q[0];
sx q[0];
rz(2.8915306) q[0];
rz(-2.1024044) q[1];
sx q[1];
rz(-2.3999374) q[1];
sx q[1];
rz(-2.340462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68576061) q[0];
sx q[0];
rz(-1.5859817) q[0];
sx q[0];
rz(-0.21160971) q[0];
x q[1];
rz(-2.4014339) q[2];
sx q[2];
rz(-2.7974786) q[2];
sx q[2];
rz(-1.0482074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2755807) q[1];
sx q[1];
rz(-0.51803127) q[1];
sx q[1];
rz(0.14054246) q[1];
x q[2];
rz(-1.2053554) q[3];
sx q[3];
rz(-1.1574928) q[3];
sx q[3];
rz(-0.45729056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73418009) q[2];
sx q[2];
rz(-1.4449747) q[2];
sx q[2];
rz(0.59047353) q[2];
rz(-2.7967795) q[3];
sx q[3];
rz(-2.7761288) q[3];
sx q[3];
rz(-1.4709681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94575268) q[0];
sx q[0];
rz(-1.4386289) q[0];
sx q[0];
rz(-1.6666743) q[0];
rz(1.6189812) q[1];
sx q[1];
rz(-1.5432065) q[1];
sx q[1];
rz(-0.408907) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58137888) q[0];
sx q[0];
rz(-1.5239851) q[0];
sx q[0];
rz(-1.5623832) q[0];
rz(-pi) q[1];
rz(-0.40134238) q[2];
sx q[2];
rz(-2.231488) q[2];
sx q[2];
rz(1.5398538) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.77264841) q[1];
sx q[1];
rz(-1.1330422) q[1];
sx q[1];
rz(2.6681929) q[1];
x q[2];
rz(2.4675995) q[3];
sx q[3];
rz(-1.6445275) q[3];
sx q[3];
rz(-2.0115832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54395479) q[2];
sx q[2];
rz(-3.0089617) q[2];
sx q[2];
rz(0.65072101) q[2];
rz(-1.4607653) q[3];
sx q[3];
rz(-1.4131578) q[3];
sx q[3];
rz(0.4167324) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2862709) q[0];
sx q[0];
rz(-2.0725508) q[0];
sx q[0];
rz(1.8977813) q[0];
rz(2.6142201) q[1];
sx q[1];
rz(-1.6970044) q[1];
sx q[1];
rz(1.7605555) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6207989) q[0];
sx q[0];
rz(-2.8993911) q[0];
sx q[0];
rz(3.0907756) q[0];
x q[1];
rz(-1.4970793) q[2];
sx q[2];
rz(-1.6699162) q[2];
sx q[2];
rz(-0.69621554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0308847) q[1];
sx q[1];
rz(-2.0763016) q[1];
sx q[1];
rz(0.55191374) q[1];
x q[2];
rz(-2.1072371) q[3];
sx q[3];
rz(-1.1749246) q[3];
sx q[3];
rz(0.27908868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1018299) q[2];
sx q[2];
rz(-1.7888174) q[2];
sx q[2];
rz(-1.9977894) q[2];
rz(0.77491289) q[3];
sx q[3];
rz(-1.9655922) q[3];
sx q[3];
rz(-2.9162143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5483122) q[0];
sx q[0];
rz(-0.53638023) q[0];
sx q[0];
rz(0.29658741) q[0];
rz(2.2071154) q[1];
sx q[1];
rz(-2.2521033) q[1];
sx q[1];
rz(-2.8737822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6952301) q[0];
sx q[0];
rz(-2.8382549) q[0];
sx q[0];
rz(0.78064755) q[0];
rz(1.3828692) q[2];
sx q[2];
rz(-0.96804995) q[2];
sx q[2];
rz(1.2056269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0632656) q[1];
sx q[1];
rz(-2.4111679) q[1];
sx q[1];
rz(2.6759381) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7095836) q[3];
sx q[3];
rz(-1.7273081) q[3];
sx q[3];
rz(1.9902347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3725793) q[2];
sx q[2];
rz(-1.0472426) q[2];
sx q[2];
rz(-0.11865842) q[2];
rz(-2.3246121) q[3];
sx q[3];
rz(-1.9575565) q[3];
sx q[3];
rz(-2.584804) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34767041) q[0];
sx q[0];
rz(-2.6642647) q[0];
sx q[0];
rz(-2.0678066) q[0];
rz(2.1444164) q[1];
sx q[1];
rz(-1.247783) q[1];
sx q[1];
rz(-2.1352077) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7770588) q[0];
sx q[0];
rz(-1.5617592) q[0];
sx q[0];
rz(-1.6994411) q[0];
rz(-0.7147737) q[2];
sx q[2];
rz(-2.6539485) q[2];
sx q[2];
rz(2.2478916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80015228) q[1];
sx q[1];
rz(-1.7625215) q[1];
sx q[1];
rz(-0.33445332) q[1];
rz(-pi) q[2];
rz(0.48000042) q[3];
sx q[3];
rz(-0.98491231) q[3];
sx q[3];
rz(-0.65174499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5239079) q[2];
sx q[2];
rz(-1.074581) q[2];
sx q[2];
rz(2.1802444) q[2];
rz(0.62054595) q[3];
sx q[3];
rz(-0.70182645) q[3];
sx q[3];
rz(2.3294241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0945702) q[0];
sx q[0];
rz(-0.71806878) q[0];
sx q[0];
rz(-2.4156003) q[0];
rz(0.30966169) q[1];
sx q[1];
rz(-2.1084712) q[1];
sx q[1];
rz(2.9594701) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0326906) q[0];
sx q[0];
rz(-1.4497524) q[0];
sx q[0];
rz(-2.3897479) q[0];
rz(-pi) q[1];
rz(-1.0605889) q[2];
sx q[2];
rz(-2.0296718) q[2];
sx q[2];
rz(-1.4993678) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9023064) q[1];
sx q[1];
rz(-0.42879802) q[1];
sx q[1];
rz(-2.8132755) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5043132) q[3];
sx q[3];
rz(-1.2361119) q[3];
sx q[3];
rz(-0.94663564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2187659) q[2];
sx q[2];
rz(-1.4309692) q[2];
sx q[2];
rz(-1.9975086) q[2];
rz(-2.4036582) q[3];
sx q[3];
rz(-0.72665557) q[3];
sx q[3];
rz(1.1206333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69796193) q[0];
sx q[0];
rz(-1.3314891) q[0];
sx q[0];
rz(0.6829845) q[0];
rz(1.5525275) q[1];
sx q[1];
rz(-0.88337675) q[1];
sx q[1];
rz(-0.89673269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67499298) q[0];
sx q[0];
rz(-1.3764125) q[0];
sx q[0];
rz(2.7036576) q[0];
x q[1];
rz(0.14210578) q[2];
sx q[2];
rz(-2.1143205) q[2];
sx q[2];
rz(-2.6212111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1410213) q[1];
sx q[1];
rz(-1.7011756) q[1];
sx q[1];
rz(-1.4818373) q[1];
rz(1.1223144) q[3];
sx q[3];
rz(-2.0297673) q[3];
sx q[3];
rz(2.36623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5985742) q[2];
sx q[2];
rz(-1.2335346) q[2];
sx q[2];
rz(2.7198071) q[2];
rz(-2.5770523) q[3];
sx q[3];
rz(-2.3081503) q[3];
sx q[3];
rz(-0.88258266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17557872) q[0];
sx q[0];
rz(-1.5400374) q[0];
sx q[0];
rz(-1.2651428) q[0];
rz(1.9059146) q[1];
sx q[1];
rz(-1.0503678) q[1];
sx q[1];
rz(-2.8138524) q[1];
rz(0.50516456) q[2];
sx q[2];
rz(-2.724833) q[2];
sx q[2];
rz(-0.45062296) q[2];
rz(-0.47673419) q[3];
sx q[3];
rz(-1.1179964) q[3];
sx q[3];
rz(-2.6159058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
