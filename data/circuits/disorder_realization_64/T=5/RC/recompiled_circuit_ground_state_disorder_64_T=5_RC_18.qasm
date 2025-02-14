OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7944827) q[0];
sx q[0];
rz(-1.4248983) q[0];
sx q[0];
rz(0.27118924) q[0];
rz(3.5085161) q[1];
sx q[1];
rz(2.38382) q[1];
sx q[1];
rz(7.741306) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9957033) q[0];
sx q[0];
rz(-1.567739) q[0];
sx q[0];
rz(-1.5778592) q[0];
rz(-pi) q[1];
rz(2.053399) q[2];
sx q[2];
rz(-1.4652243) q[2];
sx q[2];
rz(-2.2029049) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7768727) q[1];
sx q[1];
rz(-2.128771) q[1];
sx q[1];
rz(0.67529894) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28370195) q[3];
sx q[3];
rz(-1.5531504) q[3];
sx q[3];
rz(-2.4775503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1272507) q[2];
sx q[2];
rz(-0.6659826) q[2];
sx q[2];
rz(1.2464397) q[2];
rz(1.0154826) q[3];
sx q[3];
rz(-2.6452439) q[3];
sx q[3];
rz(2.523876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.604973) q[0];
sx q[0];
rz(-0.07285694) q[0];
sx q[0];
rz(0.68742043) q[0];
rz(2.5303326) q[1];
sx q[1];
rz(-1.5115073) q[1];
sx q[1];
rz(2.2329109) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1686104) q[0];
sx q[0];
rz(-1.4797204) q[0];
sx q[0];
rz(2.220972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15678997) q[2];
sx q[2];
rz(-2.1794381) q[2];
sx q[2];
rz(1.551924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7327888) q[1];
sx q[1];
rz(-2.286274) q[1];
sx q[1];
rz(-0.73231952) q[1];
x q[2];
rz(2.9396003) q[3];
sx q[3];
rz(-0.56623161) q[3];
sx q[3];
rz(-0.30893886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0841792) q[2];
sx q[2];
rz(-2.1647191) q[2];
sx q[2];
rz(-1.9834391) q[2];
rz(-0.17364764) q[3];
sx q[3];
rz(-1.6896788) q[3];
sx q[3];
rz(-2.4524073) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77370206) q[0];
sx q[0];
rz(-2.9359718) q[0];
sx q[0];
rz(-0.50262991) q[0];
rz(-1.3357119) q[1];
sx q[1];
rz(-3.0323961) q[1];
sx q[1];
rz(1.4044382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96539474) q[0];
sx q[0];
rz(-2.0098916) q[0];
sx q[0];
rz(-2.0140854) q[0];
rz(-pi) q[1];
rz(0.37356202) q[2];
sx q[2];
rz(-1.7508011) q[2];
sx q[2];
rz(1.9587868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3223649) q[1];
sx q[1];
rz(-0.33999264) q[1];
sx q[1];
rz(-0.26588666) q[1];
x q[2];
rz(2.1039906) q[3];
sx q[3];
rz(-2.0948862) q[3];
sx q[3];
rz(2.1253824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3001083) q[2];
sx q[2];
rz(-0.69977641) q[2];
sx q[2];
rz(-2.4669199) q[2];
rz(-2.789433) q[3];
sx q[3];
rz(-1.9904741) q[3];
sx q[3];
rz(1.5153511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.801055) q[0];
sx q[0];
rz(-3.099649) q[0];
sx q[0];
rz(1.1658143) q[0];
rz(0.55533987) q[1];
sx q[1];
rz(-0.97709877) q[1];
sx q[1];
rz(-1.7305444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87303091) q[0];
sx q[0];
rz(-1.8202204) q[0];
sx q[0];
rz(-1.5536352) q[0];
rz(-pi) q[1];
rz(1.0434112) q[2];
sx q[2];
rz(-1.4031271) q[2];
sx q[2];
rz(-2.4737918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.725961) q[1];
sx q[1];
rz(-1.9187732) q[1];
sx q[1];
rz(2.1044974) q[1];
rz(-0.81713809) q[3];
sx q[3];
rz(-1.3970831) q[3];
sx q[3];
rz(-0.44718633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0393151) q[2];
sx q[2];
rz(-0.5216051) q[2];
sx q[2];
rz(-2.6585141) q[2];
rz(-0.77110243) q[3];
sx q[3];
rz(-1.5811788) q[3];
sx q[3];
rz(-1.7980827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.022472) q[0];
sx q[0];
rz(-2.6799057) q[0];
sx q[0];
rz(-2.9126677) q[0];
rz(-1.4643033) q[1];
sx q[1];
rz(-0.56956446) q[1];
sx q[1];
rz(-0.48008188) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0374566) q[0];
sx q[0];
rz(-2.7711282) q[0];
sx q[0];
rz(-0.78471009) q[0];
rz(-2.6297683) q[2];
sx q[2];
rz(-2.6445342) q[2];
sx q[2];
rz(1.1443421) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11949524) q[1];
sx q[1];
rz(-1.9429617) q[1];
sx q[1];
rz(1.6254025) q[1];
x q[2];
rz(-2.1958417) q[3];
sx q[3];
rz(-1.4212928) q[3];
sx q[3];
rz(0.76555071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0698801) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(-1.9055535) q[2];
rz(-1.49336) q[3];
sx q[3];
rz(-0.83087102) q[3];
sx q[3];
rz(-1.6911471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1278095) q[0];
sx q[0];
rz(-0.74087983) q[0];
sx q[0];
rz(-2.2770449) q[0];
rz(-2.3132482) q[1];
sx q[1];
rz(-1.8048077) q[1];
sx q[1];
rz(-2.4116662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975271) q[0];
sx q[0];
rz(-2.0407009) q[0];
sx q[0];
rz(-0.92961981) q[0];
x q[1];
rz(-2.4721844) q[2];
sx q[2];
rz(-1.4140437) q[2];
sx q[2];
rz(-2.0956958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.266842) q[1];
sx q[1];
rz(-2.0952941) q[1];
sx q[1];
rz(2.8598818) q[1];
x q[2];
rz(3.0804068) q[3];
sx q[3];
rz(-2.4172999) q[3];
sx q[3];
rz(-2.7887087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.56968969) q[2];
sx q[2];
rz(-2.9986585) q[2];
sx q[2];
rz(1.8611056) q[2];
rz(-1.5754383) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(-2.2787826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5363619) q[0];
sx q[0];
rz(-2.1385758) q[0];
sx q[0];
rz(-0.59148106) q[0];
rz(-0.65762562) q[1];
sx q[1];
rz(-1.769519) q[1];
sx q[1];
rz(3.0234911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0765083) q[0];
sx q[0];
rz(-1.8002836) q[0];
sx q[0];
rz(0.35603948) q[0];
x q[1];
rz(0.51609765) q[2];
sx q[2];
rz(-1.5071259) q[2];
sx q[2];
rz(-0.89833591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.98886314) q[1];
sx q[1];
rz(-2.0901457) q[1];
sx q[1];
rz(-1.4593655) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.348446) q[3];
sx q[3];
rz(-1.7365713) q[3];
sx q[3];
rz(-0.83525688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9401231) q[2];
sx q[2];
rz(-0.38429364) q[2];
sx q[2];
rz(1.2598134) q[2];
rz(1.2093557) q[3];
sx q[3];
rz(-1.7989379) q[3];
sx q[3];
rz(-2.2945837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2131293) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(1.553836) q[0];
rz(1.8153048) q[1];
sx q[1];
rz(-2.1554048) q[1];
sx q[1];
rz(-0.80002588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5309842) q[0];
sx q[0];
rz(-1.2982695) q[0];
sx q[0];
rz(0.39604183) q[0];
rz(2.3800092) q[2];
sx q[2];
rz(-1.2084476) q[2];
sx q[2];
rz(-1.5630388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8064373) q[1];
sx q[1];
rz(-1.0943606) q[1];
sx q[1];
rz(2.128674) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0142087) q[3];
sx q[3];
rz(-1.9032904) q[3];
sx q[3];
rz(-2.1115164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.942261) q[2];
sx q[2];
rz(-2.3776725) q[2];
sx q[2];
rz(1.2660816) q[2];
rz(2.074504) q[3];
sx q[3];
rz(-1.7231924) q[3];
sx q[3];
rz(-3.0838695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1570194) q[0];
sx q[0];
rz(-2.3047801) q[0];
sx q[0];
rz(-2.9465604) q[0];
rz(-0.86206478) q[1];
sx q[1];
rz(-1.74086) q[1];
sx q[1];
rz(-2.92235) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57306141) q[0];
sx q[0];
rz(-0.6210621) q[0];
sx q[0];
rz(-2.0173418) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17941189) q[2];
sx q[2];
rz(-2.0127014) q[2];
sx q[2];
rz(-2.4561858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5089607) q[1];
sx q[1];
rz(-1.0380942) q[1];
sx q[1];
rz(-1.6091634) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52893838) q[3];
sx q[3];
rz(-2.248507) q[3];
sx q[3];
rz(-2.4916745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7361043) q[2];
sx q[2];
rz(-0.89724237) q[2];
sx q[2];
rz(-1.415095) q[2];
rz(1.803558) q[3];
sx q[3];
rz(-1.8357364) q[3];
sx q[3];
rz(-3.0726748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5445484) q[0];
sx q[0];
rz(-0.91240779) q[0];
sx q[0];
rz(1.190881) q[0];
rz(-1.4471588) q[1];
sx q[1];
rz(-1.2971327) q[1];
sx q[1];
rz(1.7097293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3668684) q[0];
sx q[0];
rz(-2.4503243) q[0];
sx q[0];
rz(-0.33095215) q[0];
x q[1];
rz(0.58379057) q[2];
sx q[2];
rz(-2.4707831) q[2];
sx q[2];
rz(-1.5354615) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14395844) q[1];
sx q[1];
rz(-1.9737509) q[1];
sx q[1];
rz(-2.2657299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3031381) q[3];
sx q[3];
rz(-1.5946232) q[3];
sx q[3];
rz(1.4122653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65149629) q[2];
sx q[2];
rz(-0.33612529) q[2];
sx q[2];
rz(2.2641342) q[2];
rz(1.6089926) q[3];
sx q[3];
rz(-1.420615) q[3];
sx q[3];
rz(1.0271614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.132591) q[0];
sx q[0];
rz(-1.4132211) q[0];
sx q[0];
rz(2.6268517) q[0];
rz(2.4280836) q[1];
sx q[1];
rz(-0.82095725) q[1];
sx q[1];
rz(1.531442) q[1];
rz(-1.8590676) q[2];
sx q[2];
rz(-1.2119515) q[2];
sx q[2];
rz(-1.8744946) q[2];
rz(-0.26127451) q[3];
sx q[3];
rz(-1.3521104) q[3];
sx q[3];
rz(2.1043652) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
