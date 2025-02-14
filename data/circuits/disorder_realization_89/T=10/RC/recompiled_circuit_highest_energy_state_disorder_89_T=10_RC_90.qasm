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
rz(3.2193174) q[0];
sx q[0];
rz(5.2585703) q[0];
sx q[0];
rz(10.724714) q[0];
rz(-0.44252244) q[1];
sx q[1];
rz(-1.1392925) q[1];
sx q[1];
rz(0.43599573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47600043) q[0];
sx q[0];
rz(-1.9622335) q[0];
sx q[0];
rz(-3.0276171) q[0];
rz(-pi) q[1];
rz(0.23171111) q[2];
sx q[2];
rz(-0.75318906) q[2];
sx q[2];
rz(1.9875658) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30937815) q[1];
sx q[1];
rz(-1.9736027) q[1];
sx q[1];
rz(-1.554053) q[1];
rz(2.1840582) q[3];
sx q[3];
rz(-1.4671578) q[3];
sx q[3];
rz(2.9985425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0787597) q[2];
sx q[2];
rz(-1.4818413) q[2];
sx q[2];
rz(3.1256342) q[2];
rz(1.4912841) q[3];
sx q[3];
rz(-1.8989547) q[3];
sx q[3];
rz(-2.7119467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780739) q[0];
sx q[0];
rz(-1.0915382) q[0];
sx q[0];
rz(1.7953405) q[0];
rz(-1.9740055) q[1];
sx q[1];
rz(-2.0474032) q[1];
sx q[1];
rz(-1.8928554) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48544504) q[0];
sx q[0];
rz(-2.1753575) q[0];
sx q[0];
rz(-2.6660835) q[0];
x q[1];
rz(-1.6925519) q[2];
sx q[2];
rz(-1.7261862) q[2];
sx q[2];
rz(-0.64823417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.84184835) q[1];
sx q[1];
rz(-1.9711291) q[1];
sx q[1];
rz(-3.0818066) q[1];
rz(-pi) q[2];
rz(-0.98442673) q[3];
sx q[3];
rz(-1.274144) q[3];
sx q[3];
rz(0.082060952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0060129082) q[2];
sx q[2];
rz(-2.4485782) q[2];
sx q[2];
rz(-0.47323027) q[2];
rz(3.0114975) q[3];
sx q[3];
rz(-1.7546763) q[3];
sx q[3];
rz(0.010802833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8148282) q[0];
sx q[0];
rz(-0.15693754) q[0];
sx q[0];
rz(2.9275295) q[0];
rz(-1.7474984) q[1];
sx q[1];
rz(-0.97624818) q[1];
sx q[1];
rz(-0.086437978) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.466216) q[0];
sx q[0];
rz(-1.4450184) q[0];
sx q[0];
rz(-1.2745538) q[0];
rz(-2.2485224) q[2];
sx q[2];
rz(-1.1249591) q[2];
sx q[2];
rz(2.0634212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1972034) q[1];
sx q[1];
rz(-1.3968588) q[1];
sx q[1];
rz(-2.7157213) q[1];
x q[2];
rz(-0.53303252) q[3];
sx q[3];
rz(-2.1500476) q[3];
sx q[3];
rz(-1.90383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4702686) q[2];
sx q[2];
rz(-0.9321804) q[2];
sx q[2];
rz(1.1364802) q[2];
rz(1.7911576) q[3];
sx q[3];
rz(-1.330749) q[3];
sx q[3];
rz(2.7854846) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7478624) q[0];
sx q[0];
rz(-0.50734729) q[0];
sx q[0];
rz(-0.90079975) q[0];
rz(-1.9081217) q[1];
sx q[1];
rz(-0.8546468) q[1];
sx q[1];
rz(-2.5305117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.184) q[0];
sx q[0];
rz(-1.5636356) q[0];
sx q[0];
rz(2.0462799) q[0];
rz(0.15005269) q[2];
sx q[2];
rz(-0.23153472) q[2];
sx q[2];
rz(0.48047149) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6590699) q[1];
sx q[1];
rz(-2.540641) q[1];
sx q[1];
rz(2.3769578) q[1];
x q[2];
rz(-2.7902725) q[3];
sx q[3];
rz(-1.4610944) q[3];
sx q[3];
rz(0.20449311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5078807) q[2];
sx q[2];
rz(-2.4780126) q[2];
sx q[2];
rz(0.12082417) q[2];
rz(3.1145596) q[3];
sx q[3];
rz(-0.20172541) q[3];
sx q[3];
rz(-2.2656608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19288572) q[0];
sx q[0];
rz(-2.3159733) q[0];
sx q[0];
rz(1.3407619) q[0];
rz(-1.6138529) q[1];
sx q[1];
rz(-1.8283045) q[1];
sx q[1];
rz(-2.5921879) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2717239) q[0];
sx q[0];
rz(-2.3921674) q[0];
sx q[0];
rz(-0.11336993) q[0];
x q[1];
rz(-1.5169296) q[2];
sx q[2];
rz(-2.6096811) q[2];
sx q[2];
rz(-0.96964449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.96619256) q[1];
sx q[1];
rz(-2.0560802) q[1];
sx q[1];
rz(3.0851689) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0321376) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(-2.1738659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.01650979) q[2];
sx q[2];
rz(-1.5089401) q[2];
sx q[2];
rz(2.5588918) q[2];
rz(-1.6216283) q[3];
sx q[3];
rz(-0.71526066) q[3];
sx q[3];
rz(-1.8250072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3358066) q[0];
sx q[0];
rz(-2.057322) q[0];
sx q[0];
rz(1.8871319) q[0];
rz(1.9460024) q[1];
sx q[1];
rz(-1.4072199) q[1];
sx q[1];
rz(-1.5171299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74745893) q[0];
sx q[0];
rz(-1.582495) q[0];
sx q[0];
rz(1.5040843) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5810691) q[2];
sx q[2];
rz(-0.88489489) q[2];
sx q[2];
rz(1.0087412) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75116762) q[1];
sx q[1];
rz(-1.8590305) q[1];
sx q[1];
rz(-0.61528968) q[1];
rz(-pi) q[2];
rz(-1.5593525) q[3];
sx q[3];
rz(-2.0858039) q[3];
sx q[3];
rz(0.29667618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17455165) q[2];
sx q[2];
rz(-1.6074564) q[2];
sx q[2];
rz(0.045844585) q[2];
rz(0.52715078) q[3];
sx q[3];
rz(-1.0322626) q[3];
sx q[3];
rz(-2.3003858) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884278) q[0];
sx q[0];
rz(-2.4808352) q[0];
sx q[0];
rz(0.38594693) q[0];
rz(1.7653607) q[1];
sx q[1];
rz(-2.0538581) q[1];
sx q[1];
rz(-1.2765346) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71300426) q[0];
sx q[0];
rz(-2.3516293) q[0];
sx q[0];
rz(2.6053564) q[0];
rz(0.047983147) q[2];
sx q[2];
rz(-1.9429824) q[2];
sx q[2];
rz(-1.7593149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0110198) q[1];
sx q[1];
rz(-1.5351194) q[1];
sx q[1];
rz(-3.1336354) q[1];
rz(-1.5064729) q[3];
sx q[3];
rz(-1.5350047) q[3];
sx q[3];
rz(2.2944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71629268) q[2];
sx q[2];
rz(-1.7159117) q[2];
sx q[2];
rz(-1.3670134) q[2];
rz(-3.0573209) q[3];
sx q[3];
rz(-1.1140946) q[3];
sx q[3];
rz(0.082898609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080634236) q[0];
sx q[0];
rz(-0.11140379) q[0];
sx q[0];
rz(-1.2782619) q[0];
rz(0.06079611) q[1];
sx q[1];
rz(-2.1863329) q[1];
sx q[1];
rz(-1.0999058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9644224) q[0];
sx q[0];
rz(-1.9510498) q[0];
sx q[0];
rz(-2.9833262) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75662778) q[2];
sx q[2];
rz(-0.95284684) q[2];
sx q[2];
rz(-0.47060395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5085735) q[1];
sx q[1];
rz(-0.46357511) q[1];
sx q[1];
rz(-2.4615088) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6693484) q[3];
sx q[3];
rz(-2.409777) q[3];
sx q[3];
rz(-1.2199618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3595769) q[2];
sx q[2];
rz(-2.4961553) q[2];
sx q[2];
rz(0.47425708) q[2];
rz(0.54840243) q[3];
sx q[3];
rz(-0.67265284) q[3];
sx q[3];
rz(-1.7614346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7062374) q[0];
sx q[0];
rz(-2.7147003) q[0];
sx q[0];
rz(-1.2109582) q[0];
rz(1.8636761) q[1];
sx q[1];
rz(-1.5457109) q[1];
sx q[1];
rz(-0.011215297) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8497711) q[0];
sx q[0];
rz(-2.1309843) q[0];
sx q[0];
rz(-1.4683506) q[0];
x q[1];
rz(-2.5915543) q[2];
sx q[2];
rz(-0.42114741) q[2];
sx q[2];
rz(2.5530346) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4947191) q[1];
sx q[1];
rz(-0.74999627) q[1];
sx q[1];
rz(1.4916625) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2645006) q[3];
sx q[3];
rz(-0.3031177) q[3];
sx q[3];
rz(2.1949286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8212905) q[2];
sx q[2];
rz(-1.236329) q[2];
sx q[2];
rz(-2.383929) q[2];
rz(2.441794) q[3];
sx q[3];
rz(-2.564513) q[3];
sx q[3];
rz(-2.2363766) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59117544) q[0];
sx q[0];
rz(-1.2400405) q[0];
sx q[0];
rz(-2.8072939) q[0];
rz(0.39407691) q[1];
sx q[1];
rz(-2.1912992) q[1];
sx q[1];
rz(-1.2324415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1533547) q[0];
sx q[0];
rz(-1.7890837) q[0];
sx q[0];
rz(1.2984896) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7712461) q[2];
sx q[2];
rz(-2.2589189) q[2];
sx q[2];
rz(-1.8375979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10127549) q[1];
sx q[1];
rz(-0.35666944) q[1];
sx q[1];
rz(2.2471395) q[1];
rz(-pi) q[2];
rz(-2.8673792) q[3];
sx q[3];
rz(-1.0875487) q[3];
sx q[3];
rz(0.26457126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1831827) q[2];
sx q[2];
rz(-2.5739058) q[2];
sx q[2];
rz(-3.0779823) q[2];
rz(-0.76572865) q[3];
sx q[3];
rz(-1.1252517) q[3];
sx q[3];
rz(3.0797899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80339377) q[0];
sx q[0];
rz(-0.82427187) q[0];
sx q[0];
rz(-1.6765539) q[0];
rz(-1.6784531) q[1];
sx q[1];
rz(-1.2668162) q[1];
sx q[1];
rz(-0.75513671) q[1];
rz(-2.4974291) q[2];
sx q[2];
rz(-1.7048057) q[2];
sx q[2];
rz(-2.6028462) q[2];
rz(0.16514292) q[3];
sx q[3];
rz(-0.84652918) q[3];
sx q[3];
rz(-2.3565945) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
