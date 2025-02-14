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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76554322) q[0];
sx q[0];
rz(-0.57224579) q[0];
sx q[0];
rz(0.28387286) q[0];
x q[1];
rz(1.1101025) q[2];
sx q[2];
rz(-2.0909799) q[2];
sx q[2];
rz(-0.68797639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.097259911) q[1];
sx q[1];
rz(-9/(16*pi)) q[1];
sx q[1];
rz(1.8392842) q[1];
x q[2];
rz(-0.13924573) q[3];
sx q[3];
rz(-1.2069824) q[3];
sx q[3];
rz(2.0219959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-0.60167998) q[2];
sx q[2];
rz(1.0428693) q[2];
rz(3.0964105) q[3];
sx q[3];
rz(-0.16862814) q[3];
sx q[3];
rz(-1.5359623) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586108) q[0];
sx q[0];
rz(-0.11317145) q[0];
sx q[0];
rz(-0.97852069) q[0];
rz(1.1631896) q[1];
sx q[1];
rz(-0.99449831) q[1];
sx q[1];
rz(1.8189583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0019313) q[0];
sx q[0];
rz(-2.0531512) q[0];
sx q[0];
rz(2.300452) q[0];
rz(-pi) q[1];
x q[1];
rz(1.00707) q[2];
sx q[2];
rz(-0.33180922) q[2];
sx q[2];
rz(-2.4203468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45554494) q[1];
sx q[1];
rz(-1.7864704) q[1];
sx q[1];
rz(-0.91051813) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9632872) q[3];
sx q[3];
rz(-2.4394991) q[3];
sx q[3];
rz(2.3012379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3375552) q[2];
sx q[2];
rz(-2.940371) q[2];
sx q[2];
rz(3.1031754) q[2];
rz(1.6328968) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(-1.4940777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22817336) q[0];
sx q[0];
rz(-0.67263022) q[0];
sx q[0];
rz(-0.2359373) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(-2.6441914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59374124) q[0];
sx q[0];
rz(-1.7501795) q[0];
sx q[0];
rz(0.18130882) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1549306) q[2];
sx q[2];
rz(-1.9344887) q[2];
sx q[2];
rz(-1.9726582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5798074) q[1];
sx q[1];
rz(-0.39424636) q[1];
sx q[1];
rz(-0.057478776) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24876066) q[3];
sx q[3];
rz(-1.3205055) q[3];
sx q[3];
rz(-0.32969013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1824789) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(-1.2314931) q[2];
rz(-1.4500827) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(-0.85795295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602144) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(0.1928992) q[0];
rz(1.7861722) q[1];
sx q[1];
rz(-1.5602292) q[1];
sx q[1];
rz(-0.17098175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3340958) q[0];
sx q[0];
rz(-2.0863895) q[0];
sx q[0];
rz(-3.0741755) q[0];
rz(-pi) q[1];
rz(2.1095721) q[2];
sx q[2];
rz(-2.2316885) q[2];
sx q[2];
rz(-0.62884841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0031922) q[1];
sx q[1];
rz(-0.44282162) q[1];
sx q[1];
rz(1.3243616) q[1];
rz(-0.90733068) q[3];
sx q[3];
rz(-1.3795092) q[3];
sx q[3];
rz(-0.80463791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1811447) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(-2.81874) q[2];
rz(-1.7668096) q[3];
sx q[3];
rz(-1.0389453) q[3];
sx q[3];
rz(0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5906931) q[0];
sx q[0];
rz(-2.1192079) q[0];
sx q[0];
rz(-0.92556959) q[0];
rz(-3.0135221) q[1];
sx q[1];
rz(-1.4984727) q[1];
sx q[1];
rz(0.099460348) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29726899) q[0];
sx q[0];
rz(-2.2873061) q[0];
sx q[0];
rz(3.141527) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9173309) q[2];
sx q[2];
rz(-1.230403) q[2];
sx q[2];
rz(0.16763359) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0316098) q[1];
sx q[1];
rz(-1.3760202) q[1];
sx q[1];
rz(-1.0321478) q[1];
rz(0.080187967) q[3];
sx q[3];
rz(-1.8155451) q[3];
sx q[3];
rz(2.484124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1818992) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(2.7679475) q[2];
rz(0.89961189) q[3];
sx q[3];
rz(-2.3989232) q[3];
sx q[3];
rz(1.1348178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9099092) q[0];
sx q[0];
rz(-0.21655701) q[0];
sx q[0];
rz(-0.59189558) q[0];
rz(2.9369211) q[1];
sx q[1];
rz(-0.99212956) q[1];
sx q[1];
rz(3.0904904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0436193) q[0];
sx q[0];
rz(-1.7228925) q[0];
sx q[0];
rz(-2.9469304) q[0];
rz(1.899753) q[2];
sx q[2];
rz(-0.68495132) q[2];
sx q[2];
rz(0.16835131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8569736) q[1];
sx q[1];
rz(-1.9195119) q[1];
sx q[1];
rz(0.32964175) q[1];
x q[2];
rz(0.53364086) q[3];
sx q[3];
rz(-1.1325877) q[3];
sx q[3];
rz(-0.72497382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11703141) q[2];
sx q[2];
rz(-1.1144964) q[2];
sx q[2];
rz(2.047211) q[2];
rz(0.38703212) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(-0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7759906) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(2.8940417) q[0];
rz(-1.6869102) q[1];
sx q[1];
rz(-1.8984112) q[1];
sx q[1];
rz(0.036458485) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55900967) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(-2.3733632) q[0];
rz(-pi) q[1];
rz(-1.7629303) q[2];
sx q[2];
rz(-2.2647144) q[2];
sx q[2];
rz(-2.1215699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4562021) q[1];
sx q[1];
rz(-1.0674607) q[1];
sx q[1];
rz(0.53629843) q[1];
rz(-pi) q[2];
rz(-3.0667449) q[3];
sx q[3];
rz(-1.4818496) q[3];
sx q[3];
rz(1.68881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.054691943) q[2];
sx q[2];
rz(-2.1754103) q[2];
sx q[2];
rz(1.7894233) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5779481) q[3];
sx q[3];
rz(-1.9498391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9098814) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(-1.6491718) q[0];
rz(-0.17852783) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(2.5915204) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4208385) q[0];
sx q[0];
rz(-1.6021358) q[0];
sx q[0];
rz(0.019958812) q[0];
rz(-pi) q[1];
rz(1.8166601) q[2];
sx q[2];
rz(-0.255092) q[2];
sx q[2];
rz(2.6691797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7272721) q[1];
sx q[1];
rz(-1.7309233) q[1];
sx q[1];
rz(0.54267197) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86217238) q[3];
sx q[3];
rz(-0.58364999) q[3];
sx q[3];
rz(-2.421511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.7358235) q[2];
sx q[2];
rz(2.3883635) q[2];
rz(-2.8473162) q[3];
sx q[3];
rz(-3.0615276) q[3];
sx q[3];
rz(3.0729955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2733961) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(-1.4519325) q[0];
rz(0.1952576) q[1];
sx q[1];
rz(-2.0611019) q[1];
sx q[1];
rz(1.6498227) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30376745) q[0];
sx q[0];
rz(-1.8762824) q[0];
sx q[0];
rz(-0.49834337) q[0];
x q[1];
rz(-0.765018) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(1.8945872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0635707) q[1];
sx q[1];
rz(-1.6352904) q[1];
sx q[1];
rz(-1.4149354) q[1];
x q[2];
rz(-3.1141485) q[3];
sx q[3];
rz(-1.5774836) q[3];
sx q[3];
rz(-0.96379507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1008272) q[2];
sx q[2];
rz(-2.6746174) q[2];
sx q[2];
rz(-0.85774285) q[2];
rz(-1.4486676) q[3];
sx q[3];
rz(-1.6489776) q[3];
sx q[3];
rz(2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0070852) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(-0.78306985) q[0];
rz(-2.018441) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(3.0812841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9674112) q[0];
sx q[0];
rz(-0.36752146) q[0];
sx q[0];
rz(-1.307748) q[0];
rz(1.4122293) q[2];
sx q[2];
rz(-2.5494908) q[2];
sx q[2];
rz(2.427096) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8728208) q[1];
sx q[1];
rz(-2.4686681) q[1];
sx q[1];
rz(-0.10867837) q[1];
x q[2];
rz(2.7785886) q[3];
sx q[3];
rz(-1.6771183) q[3];
sx q[3];
rz(-2.8385988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15410885) q[2];
sx q[2];
rz(-0.86468148) q[2];
sx q[2];
rz(1.3101428) q[2];
rz(-1.8080669) q[3];
sx q[3];
rz(-1.3430877) q[3];
sx q[3];
rz(-1.3069299) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780728) q[0];
sx q[0];
rz(-2.0250043) q[0];
sx q[0];
rz(2.9153839) q[0];
rz(-0.75640596) q[1];
sx q[1];
rz(-0.49497985) q[1];
sx q[1];
rz(2.3351647) q[1];
rz(2.7112167) q[2];
sx q[2];
rz(-1.5622361) q[2];
sx q[2];
rz(2.1282276) q[2];
rz(1.1011878) q[3];
sx q[3];
rz(-2.7536177) q[3];
sx q[3];
rz(-2.9256647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
