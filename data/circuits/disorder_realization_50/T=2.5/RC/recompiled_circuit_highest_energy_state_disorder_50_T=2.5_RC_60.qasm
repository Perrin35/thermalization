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
rz(-1.4122352) q[0];
sx q[0];
rz(-0.5923624) q[0];
sx q[0];
rz(1.9368197) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(1.9538716) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.710085) q[0];
sx q[0];
rz(-1.0241226) q[0];
sx q[0];
rz(1.3923079) q[0];
rz(2.4815791) q[2];
sx q[2];
rz(-2.4610991) q[2];
sx q[2];
rz(-1.6689491) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0443327) q[1];
sx q[1];
rz(-2.9625434) q[1];
sx q[1];
rz(1.3023085) q[1];
rz(-pi) q[2];
rz(1.9378565) q[3];
sx q[3];
rz(-1.4407183) q[3];
sx q[3];
rz(-2.6405622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54186934) q[2];
sx q[2];
rz(-2.5399127) q[2];
sx q[2];
rz(-1.0428693) q[2];
rz(-0.045182191) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(-1.6056304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586108) q[0];
sx q[0];
rz(-0.11317145) q[0];
sx q[0];
rz(-0.97852069) q[0];
rz(-1.9784031) q[1];
sx q[1];
rz(-2.1470943) q[1];
sx q[1];
rz(1.3226343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037905692) q[0];
sx q[0];
rz(-0.9390489) q[0];
sx q[0];
rz(2.5292255) q[0];
x q[1];
rz(2.959526) q[2];
sx q[2];
rz(-1.8497502) q[2];
sx q[2];
rz(-1.830991) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.757431) q[1];
sx q[1];
rz(-2.4520284) q[1];
sx q[1];
rz(1.227725) q[1];
x q[2];
rz(2.9632872) q[3];
sx q[3];
rz(-0.70209356) q[3];
sx q[3];
rz(-2.3012379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80403745) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(3.1031754) q[2];
rz(1.5086959) q[3];
sx q[3];
rz(-1.326606) q[3];
sx q[3];
rz(-1.4940777) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22817336) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(2.9056554) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.6940247) q[1];
sx q[1];
rz(-0.49740121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.936393) q[0];
sx q[0];
rz(-0.25435624) q[0];
sx q[0];
rz(-2.3533871) q[0];
rz(1.986662) q[2];
sx q[2];
rz(-1.9344887) q[2];
sx q[2];
rz(1.9726582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.062089109) q[1];
sx q[1];
rz(-1.5487284) q[1];
sx q[1];
rz(0.39366054) q[1];
x q[2];
rz(-0.80422359) q[3];
sx q[3];
rz(-2.7905593) q[3];
sx q[3];
rz(-1.1277175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95911372) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(-1.4500827) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(-0.85795295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813783) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(-0.1928992) q[0];
rz(-1.3554205) q[1];
sx q[1];
rz(-1.5813634) q[1];
sx q[1];
rz(0.17098175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6714013) q[0];
sx q[0];
rz(-0.51958749) q[0];
sx q[0];
rz(1.6891102) q[0];
x q[1];
rz(2.1095721) q[2];
sx q[2];
rz(-2.2316885) q[2];
sx q[2];
rz(-0.62884841) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2747171) q[1];
sx q[1];
rz(-1.9993385) q[1];
sx q[1];
rz(3.0264167) q[1];
x q[2];
rz(2.234262) q[3];
sx q[3];
rz(-1.3795092) q[3];
sx q[3];
rz(-0.80463791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1811447) q[2];
sx q[2];
rz(-1.170271) q[2];
sx q[2];
rz(0.32285264) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-2.1026473) q[3];
sx q[3];
rz(-0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5906931) q[0];
sx q[0];
rz(-2.1192079) q[0];
sx q[0];
rz(0.92556959) q[0];
rz(0.12807056) q[1];
sx q[1];
rz(-1.4984727) q[1];
sx q[1];
rz(0.099460348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29716897) q[0];
sx q[0];
rz(-0.71650973) q[0];
sx q[0];
rz(1.5708718) q[0];
x q[1];
rz(1.9173309) q[2];
sx q[2];
rz(-1.9111897) q[2];
sx q[2];
rz(0.16763359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9156218) q[1];
sx q[1];
rz(-0.56949896) q[1];
sx q[1];
rz(-1.2036588) q[1];
rz(-pi) q[2];
rz(0.080187967) q[3];
sx q[3];
rz(-1.8155451) q[3];
sx q[3];
rz(-0.65746869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9596935) q[2];
sx q[2];
rz(-1.5450031) q[2];
sx q[2];
rz(-2.7679475) q[2];
rz(-2.2419808) q[3];
sx q[3];
rz(-2.3989232) q[3];
sx q[3];
rz(-2.0067748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(0.59189558) q[0];
rz(-0.20467155) q[1];
sx q[1];
rz(-0.99212956) q[1];
sx q[1];
rz(-0.051102292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.013553) q[0];
sx q[0];
rz(-0.2464412) q[0];
sx q[0];
rz(0.67009573) q[0];
x q[1];
rz(-1.2418397) q[2];
sx q[2];
rz(-2.4566413) q[2];
sx q[2];
rz(2.9732413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8569736) q[1];
sx q[1];
rz(-1.2220807) q[1];
sx q[1];
rz(-0.32964175) q[1];
rz(1.072364) q[3];
sx q[3];
rz(-2.0494771) q[3];
sx q[3];
rz(2.5414027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0245612) q[2];
sx q[2];
rz(-1.1144964) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(2.7545605) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(-2.835137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36560202) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(0.24755092) q[0];
rz(1.4546825) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(-0.036458485) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2941847) q[0];
sx q[0];
rz(-0.80889055) q[0];
sx q[0];
rz(-2.7436867) q[0];
rz(-pi) q[1];
rz(2.9159732) q[2];
sx q[2];
rz(-0.7157514) q[2];
sx q[2];
rz(0.72474397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6853906) q[1];
sx q[1];
rz(-1.0674607) q[1];
sx q[1];
rz(-0.53629843) q[1];
x q[2];
rz(-3.0667449) q[3];
sx q[3];
rz(-1.4818496) q[3];
sx q[3];
rz(1.68881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.054691943) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(-1.3521693) q[2];
rz(1.2482721) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(-1.1917535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9098814) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(1.6491718) q[0];
rz(-2.9630648) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(0.55007225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15353841) q[0];
sx q[0];
rz(-0.037153553) q[0];
sx q[0];
rz(-1.0038934) q[0];
x q[1];
rz(-1.8185316) q[2];
sx q[2];
rz(-1.6322517) q[2];
sx q[2];
rz(0.8601735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8892563) q[1];
sx q[1];
rz(-2.1057711) q[1];
sx q[1];
rz(1.3843797) q[1];
rz(0.40591235) q[3];
sx q[3];
rz(-2.0024869) q[3];
sx q[3];
rz(1.622705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(2.3883635) q[2];
rz(-0.29427648) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(-0.068597138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2733961) q[0];
sx q[0];
rz(-2.8559339) q[0];
sx q[0];
rz(-1.6896601) q[0];
rz(-2.946335) q[1];
sx q[1];
rz(-1.0804907) q[1];
sx q[1];
rz(1.4917699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7719472) q[0];
sx q[0];
rz(-0.5777244) q[0];
sx q[0];
rz(0.58322241) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0259616) q[2];
sx q[2];
rz(-2.3732269) q[2];
sx q[2];
rz(-0.40711428) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0780219) q[1];
sx q[1];
rz(-1.6352904) q[1];
sx q[1];
rz(-1.4149354) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1141485) q[3];
sx q[3];
rz(-1.5641091) q[3];
sx q[3];
rz(2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0407654) q[2];
sx q[2];
rz(-2.6746174) q[2];
sx q[2];
rz(0.85774285) q[2];
rz(1.6929251) q[3];
sx q[3];
rz(-1.492615) q[3];
sx q[3];
rz(-2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345074) q[0];
sx q[0];
rz(-2.9869098) q[0];
sx q[0];
rz(0.78306985) q[0];
rz(1.1231517) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(-0.060308594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9674112) q[0];
sx q[0];
rz(-2.7740712) q[0];
sx q[0];
rz(-1.307748) q[0];
rz(-3.0357828) q[2];
sx q[2];
rz(-2.1544837) q[2];
sx q[2];
rz(2.2367144) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3871413) q[1];
sx q[1];
rz(-1.6384513) q[1];
sx q[1];
rz(-0.67004244) q[1];
x q[2];
rz(-1.4571244) q[3];
sx q[3];
rz(-1.9316564) q[3];
sx q[3];
rz(1.833503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15410885) q[2];
sx q[2];
rz(-0.86468148) q[2];
sx q[2];
rz(1.3101428) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(1.3069299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780728) q[0];
sx q[0];
rz(-2.0250043) q[0];
sx q[0];
rz(2.9153839) q[0];
rz(-2.3851867) q[1];
sx q[1];
rz(-2.6466128) q[1];
sx q[1];
rz(-0.80642798) q[1];
rz(2.7112167) q[2];
sx q[2];
rz(-1.5622361) q[2];
sx q[2];
rz(2.1282276) q[2];
rz(-2.0404048) q[3];
sx q[3];
rz(-2.7536177) q[3];
sx q[3];
rz(-2.9256647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
