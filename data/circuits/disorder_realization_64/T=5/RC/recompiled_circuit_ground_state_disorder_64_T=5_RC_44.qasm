OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34710994) q[0];
sx q[0];
rz(-1.7166944) q[0];
sx q[0];
rz(-0.27118924) q[0];
rz(-2.7746692) q[1];
sx q[1];
rz(-0.75777268) q[1];
sx q[1];
rz(-1.4581207) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0163907) q[0];
sx q[0];
rz(-0.0076961829) q[0];
sx q[0];
rz(1.1622692) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0225296) q[2];
sx q[2];
rz(-2.0504842) q[2];
sx q[2];
rz(-2.4543311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36472) q[1];
sx q[1];
rz(-1.0128216) q[1];
sx q[1];
rz(-2.4662937) q[1];
x q[2];
rz(0.062964418) q[3];
sx q[3];
rz(-2.8573572) q[3];
sx q[3];
rz(2.1743944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.014341982) q[2];
sx q[2];
rz(-0.6659826) q[2];
sx q[2];
rz(1.8951529) q[2];
rz(1.0154826) q[3];
sx q[3];
rz(-0.4963488) q[3];
sx q[3];
rz(-2.523876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366197) q[0];
sx q[0];
rz(-3.0687357) q[0];
sx q[0];
rz(-2.4541722) q[0];
rz(-0.61126002) q[1];
sx q[1];
rz(-1.5115073) q[1];
sx q[1];
rz(-0.90868178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8584367) q[0];
sx q[0];
rz(-0.65560616) q[0];
sx q[0];
rz(1.7205419) q[0];
rz(-1.7912175) q[2];
sx q[2];
rz(-0.62602717) q[2];
sx q[2];
rz(-1.859425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7327888) q[1];
sx q[1];
rz(-0.85531863) q[1];
sx q[1];
rz(2.4092731) q[1];
rz(-0.20199235) q[3];
sx q[3];
rz(-2.575361) q[3];
sx q[3];
rz(-2.8326538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0841792) q[2];
sx q[2];
rz(-0.97687352) q[2];
sx q[2];
rz(-1.1581536) q[2];
rz(0.17364764) q[3];
sx q[3];
rz(-1.4519139) q[3];
sx q[3];
rz(0.68918532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77370206) q[0];
sx q[0];
rz(-2.9359718) q[0];
sx q[0];
rz(2.6389627) q[0];
rz(-1.3357119) q[1];
sx q[1];
rz(-3.0323961) q[1];
sx q[1];
rz(-1.7371545) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40622845) q[0];
sx q[0];
rz(-1.9695008) q[0];
sx q[0];
rz(-2.6621292) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3777755) q[2];
sx q[2];
rz(-1.9380331) q[2];
sx q[2];
rz(-2.6835416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3223649) q[1];
sx q[1];
rz(-2.8016) q[1];
sx q[1];
rz(2.875706) q[1];
x q[2];
rz(-0.59111528) q[3];
sx q[3];
rz(-1.1151259) q[3];
sx q[3];
rz(2.299813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8414843) q[2];
sx q[2];
rz(-0.69977641) q[2];
sx q[2];
rz(2.4669199) q[2];
rz(2.789433) q[3];
sx q[3];
rz(-1.1511185) q[3];
sx q[3];
rz(1.5153511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053764) q[0];
sx q[0];
rz(-0.041943701) q[0];
sx q[0];
rz(1.1658143) q[0];
rz(2.5862528) q[1];
sx q[1];
rz(-2.1644939) q[1];
sx q[1];
rz(1.4110483) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4395907) q[0];
sx q[0];
rz(-1.5874264) q[0];
sx q[0];
rz(2.8921333) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.246366) q[2];
sx q[2];
rz(-0.55098767) q[2];
sx q[2];
rz(-1.9595326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4630314) q[1];
sx q[1];
rz(-0.62778607) q[1];
sx q[1];
rz(-2.1902172) q[1];
x q[2];
rz(1.8218231) q[3];
sx q[3];
rz(-0.76956144) q[3];
sx q[3];
rz(-1.8358474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0393151) q[2];
sx q[2];
rz(-0.5216051) q[2];
sx q[2];
rz(0.48307854) q[2];
rz(-0.77110243) q[3];
sx q[3];
rz(-1.5604138) q[3];
sx q[3];
rz(1.7980827) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1191206) q[0];
sx q[0];
rz(-2.6799057) q[0];
sx q[0];
rz(-0.22892496) q[0];
rz(-1.4643033) q[1];
sx q[1];
rz(-0.56956446) q[1];
sx q[1];
rz(-0.48008188) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92394476) q[0];
sx q[0];
rz(-1.311725) q[0];
sx q[0];
rz(-1.3029419) q[0];
rz(2.6998015) q[2];
sx q[2];
rz(-1.3350772) q[2];
sx q[2];
rz(-3.1093895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.029705392) q[1];
sx q[1];
rz(-0.3759653) q[1];
sx q[1];
rz(0.13891797) q[1];
rz(2.1958417) q[3];
sx q[3];
rz(-1.4212928) q[3];
sx q[3];
rz(2.3760419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.071712581) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(-1.9055535) q[2];
rz(1.49336) q[3];
sx q[3];
rz(-0.83087102) q[3];
sx q[3];
rz(1.6911471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013783197) q[0];
sx q[0];
rz(-2.4007128) q[0];
sx q[0];
rz(-2.2770449) q[0];
rz(2.3132482) q[1];
sx q[1];
rz(-1.8048077) q[1];
sx q[1];
rz(2.4116662) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406554) q[0];
sx q[0];
rz(-1.1008917) q[0];
sx q[0];
rz(2.2119728) q[0];
rz(-1.3719158) q[2];
sx q[2];
rz(-0.91105295) q[2];
sx q[2];
rz(-2.4937862) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3507281) q[1];
sx q[1];
rz(-0.58907382) q[1];
sx q[1];
rz(2.0187316) q[1];
rz(-1.5167522) q[3];
sx q[3];
rz(-0.84815787) q[3];
sx q[3];
rz(-0.27127008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56968969) q[2];
sx q[2];
rz(-2.9986585) q[2];
sx q[2];
rz(-1.8611056) q[2];
rz(1.5754383) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(-0.86281002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(2.5363619) q[0];
sx q[0];
rz(-1.0030168) q[0];
sx q[0];
rz(2.5501116) q[0];
rz(-0.65762562) q[1];
sx q[1];
rz(-1.769519) q[1];
sx q[1];
rz(-0.11810158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0869136) q[0];
sx q[0];
rz(-2.7206693) q[0];
sx q[0];
rz(0.59043365) q[0];
x q[1];
rz(-2.625495) q[2];
sx q[2];
rz(-1.6344667) q[2];
sx q[2];
rz(0.89833591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76712045) q[1];
sx q[1];
rz(-0.53009696) q[1];
sx q[1];
rz(2.9494826) q[1];
rz(-0.16987993) q[3];
sx q[3];
rz(-1.7900482) q[3];
sx q[3];
rz(-2.3687621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9401231) q[2];
sx q[2];
rz(-0.38429364) q[2];
sx q[2];
rz(1.2598134) q[2];
rz(-1.9322369) q[3];
sx q[3];
rz(-1.7989379) q[3];
sx q[3];
rz(0.847009) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92846337) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(-1.553836) q[0];
rz(-1.3262879) q[1];
sx q[1];
rz(-0.98618788) q[1];
sx q[1];
rz(0.80002588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6094506) q[0];
sx q[0];
rz(-2.6649619) q[0];
sx q[0];
rz(2.5146288) q[0];
rz(-0.7615835) q[2];
sx q[2];
rz(-1.2084476) q[2];
sx q[2];
rz(-1.5630388) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8064373) q[1];
sx q[1];
rz(-1.0943606) q[1];
sx q[1];
rz(1.0129187) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1497277) q[3];
sx q[3];
rz(-0.63922113) q[3];
sx q[3];
rz(1.0238436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19933166) q[2];
sx q[2];
rz(-0.76392019) q[2];
sx q[2];
rz(1.2660816) q[2];
rz(1.0670886) q[3];
sx q[3];
rz(-1.7231924) q[3];
sx q[3];
rz(-0.057723109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1570194) q[0];
sx q[0];
rz(-2.3047801) q[0];
sx q[0];
rz(-0.19503221) q[0];
rz(-0.86206478) q[1];
sx q[1];
rz(-1.4007327) q[1];
sx q[1];
rz(2.92235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1006323) q[0];
sx q[0];
rz(-1.0182683) q[0];
sx q[0];
rz(-2.841903) q[0];
rz(1.2101096) q[2];
sx q[2];
rz(-2.6668913) q[2];
sx q[2];
rz(1.086496) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4335208) q[1];
sx q[1];
rz(-0.5339491) q[1];
sx q[1];
rz(-3.0766218) q[1];
rz(-pi) q[2];
rz(-0.52893838) q[3];
sx q[3];
rz(-0.89308561) q[3];
sx q[3];
rz(2.4916745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7361043) q[2];
sx q[2];
rz(-2.2443503) q[2];
sx q[2];
rz(1.415095) q[2];
rz(-1.803558) q[3];
sx q[3];
rz(-1.8357364) q[3];
sx q[3];
rz(-0.068917902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5445484) q[0];
sx q[0];
rz(-2.2291849) q[0];
sx q[0];
rz(1.9507116) q[0];
rz(-1.4471588) q[1];
sx q[1];
rz(-1.2971327) q[1];
sx q[1];
rz(1.7097293) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3668684) q[0];
sx q[0];
rz(-2.4503243) q[0];
sx q[0];
rz(0.33095215) q[0];
rz(1.1584615) q[2];
sx q[2];
rz(-2.1160876) q[2];
sx q[2];
rz(2.3067428) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.398795) q[1];
sx q[1];
rz(-0.94091641) q[1];
sx q[1];
rz(0.50666084) q[1];
rz(1.5351684) q[3];
sx q[3];
rz(-0.7326574) q[3];
sx q[3];
rz(-3.0095524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4900964) q[2];
sx q[2];
rz(-2.8054674) q[2];
sx q[2];
rz(0.87745848) q[2];
rz(1.5326001) q[3];
sx q[3];
rz(-1.420615) q[3];
sx q[3];
rz(-1.0271614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.132591) q[0];
sx q[0];
rz(-1.7283716) q[0];
sx q[0];
rz(-0.51474095) q[0];
rz(-2.4280836) q[1];
sx q[1];
rz(-2.3206354) q[1];
sx q[1];
rz(-1.6101507) q[1];
rz(0.3729214) q[2];
sx q[2];
rz(-1.3013617) q[2];
sx q[2];
rz(-0.19993275) q[2];
rz(1.3446829) q[3];
sx q[3];
rz(-1.825708) q[3];
sx q[3];
rz(0.59151266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
