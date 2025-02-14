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
rz(-2.2235121) q[0];
rz(3.0300568) q[1];
sx q[1];
rz(-1.7938951) q[1];
sx q[1];
rz(1.5741875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2409256) q[0];
sx q[0];
rz(-1.5742212) q[0];
sx q[0];
rz(1.5802556) q[0];
x q[1];
rz(-1.9419045) q[2];
sx q[2];
rz(-2.1016309) q[2];
sx q[2];
rz(0.67000721) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5985377) q[1];
sx q[1];
rz(-2.8208113) q[1];
sx q[1];
rz(1.2307274) q[1];
rz(2.918675) q[3];
sx q[3];
rz(-2.6379728) q[3];
sx q[3];
rz(0.38583392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.427318) q[2];
sx q[2];
rz(-1.486342) q[2];
sx q[2];
rz(1.8753258) q[2];
rz(-1.0424987) q[3];
sx q[3];
rz(-1.5407341) q[3];
sx q[3];
rz(-1.3955759) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.105044) q[0];
sx q[0];
rz(-0.62750134) q[0];
sx q[0];
rz(3.0446766) q[0];
rz(-2.3333683) q[1];
sx q[1];
rz(-2.7327635) q[1];
sx q[1];
rz(-1.854863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185583) q[0];
sx q[0];
rz(-2.0204394) q[0];
sx q[0];
rz(0.03944035) q[0];
x q[1];
rz(-0.72186462) q[2];
sx q[2];
rz(-2.0316191) q[2];
sx q[2];
rz(0.45022717) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0285312) q[1];
sx q[1];
rz(-2.3120027) q[1];
sx q[1];
rz(2.100551) q[1];
x q[2];
rz(0.33830418) q[3];
sx q[3];
rz(-1.9060935) q[3];
sx q[3];
rz(1.9453059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6404932) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(-0.38197771) q[2];
rz(0.52465087) q[3];
sx q[3];
rz(-1.4068406) q[3];
sx q[3];
rz(2.9978571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.3749738) q[0];
sx q[0];
rz(-1.5856278) q[0];
sx q[0];
rz(-2.4554456) q[0];
rz(1.067591) q[1];
sx q[1];
rz(-0.99718863) q[1];
sx q[1];
rz(3.0608665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3878306) q[0];
sx q[0];
rz(-1.7942611) q[0];
sx q[0];
rz(3.1016283) q[0];
x q[1];
rz(-2.3600134) q[2];
sx q[2];
rz(-2.321876) q[2];
sx q[2];
rz(0.21910659) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0856493) q[1];
sx q[1];
rz(-0.23134821) q[1];
sx q[1];
rz(-0.10142188) q[1];
rz(-pi) q[2];
rz(2.0423546) q[3];
sx q[3];
rz(-2.3009389) q[3];
sx q[3];
rz(1.7748347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8778235) q[2];
sx q[2];
rz(-0.70085415) q[2];
sx q[2];
rz(-0.43886718) q[2];
rz(-1.8435439) q[3];
sx q[3];
rz(-2.7545007) q[3];
sx q[3];
rz(-0.41518655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91442672) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(0.67657226) q[0];
rz(3.0034972) q[1];
sx q[1];
rz(-2.0510249) q[1];
sx q[1];
rz(2.9409883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.706079) q[0];
sx q[0];
rz(-1.2114176) q[0];
sx q[0];
rz(-2.2227004) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5121763) q[2];
sx q[2];
rz(-1.417629) q[2];
sx q[2];
rz(1.6352194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8183867) q[1];
sx q[1];
rz(-2.1765059) q[1];
sx q[1];
rz(-1.731864) q[1];
rz(1.7142606) q[3];
sx q[3];
rz(-1.0142361) q[3];
sx q[3];
rz(-1.0866764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5033919) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(-1.7207883) q[2];
rz(-3.0270789) q[3];
sx q[3];
rz(-1.6679461) q[3];
sx q[3];
rz(-2.1421471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81922174) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(0.97212273) q[0];
rz(2.8179893) q[1];
sx q[1];
rz(-1.690003) q[1];
sx q[1];
rz(1.9591029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33859461) q[0];
sx q[0];
rz(-1.3756244) q[0];
sx q[0];
rz(3.1006569) q[0];
x q[1];
rz(-0.57020541) q[2];
sx q[2];
rz(-1.5758762) q[2];
sx q[2];
rz(-3.0722116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8132434) q[1];
sx q[1];
rz(-2.9651838) q[1];
sx q[1];
rz(-1.117068) q[1];
rz(-0.72850169) q[3];
sx q[3];
rz(-2.430901) q[3];
sx q[3];
rz(-1.817734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16367308) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(0.53691205) q[2];
rz(-0.72530693) q[3];
sx q[3];
rz(-0.0497497) q[3];
sx q[3];
rz(2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86730114) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(-1.6987479) q[0];
rz(-2.9310215) q[1];
sx q[1];
rz(-1.6981533) q[1];
sx q[1];
rz(0.64705667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7919901) q[0];
sx q[0];
rz(-2.6812883) q[0];
sx q[0];
rz(-2.5352298) q[0];
rz(2.9672616) q[2];
sx q[2];
rz(-2.5781879) q[2];
sx q[2];
rz(2.8772815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9155026) q[1];
sx q[1];
rz(-1.9880297) q[1];
sx q[1];
rz(-1.1933865) q[1];
rz(0.67630597) q[3];
sx q[3];
rz(-0.42255536) q[3];
sx q[3];
rz(-3.0129715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1232542) q[2];
sx q[2];
rz(-1.0945357) q[2];
sx q[2];
rz(-1.1798165) q[2];
rz(-1.8566462) q[3];
sx q[3];
rz(-2.732087) q[3];
sx q[3];
rz(3.0783317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0354075) q[0];
sx q[0];
rz(-1.6949061) q[0];
sx q[0];
rz(-2.2440198) q[0];
rz(-1.8073742) q[1];
sx q[1];
rz(-1.370627) q[1];
sx q[1];
rz(2.1468377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16066954) q[0];
sx q[0];
rz(-0.86359016) q[0];
sx q[0];
rz(-0.41236931) q[0];
x q[1];
rz(-2.1894073) q[2];
sx q[2];
rz(-2.4403009) q[2];
sx q[2];
rz(-1.2522956) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89284507) q[1];
sx q[1];
rz(-1.1214646) q[1];
sx q[1];
rz(0.32545089) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8514824) q[3];
sx q[3];
rz(-1.5475905) q[3];
sx q[3];
rz(0.9965903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0576386) q[2];
sx q[2];
rz(-1.8123764) q[2];
sx q[2];
rz(0.29339054) q[2];
rz(3.0883664) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(1.4330385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5720125) q[0];
sx q[0];
rz(-1.3895637) q[0];
sx q[0];
rz(0.30297512) q[0];
rz(-1.4888034) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(-2.4724919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5821823) q[0];
sx q[0];
rz(-0.56920496) q[0];
sx q[0];
rz(0.64535443) q[0];
x q[1];
rz(2.1630493) q[2];
sx q[2];
rz(-2.7251232) q[2];
sx q[2];
rz(1.5964674) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2337906) q[1];
sx q[1];
rz(-1.0668584) q[1];
sx q[1];
rz(-0.48598098) q[1];
rz(-pi) q[2];
rz(-0.49578285) q[3];
sx q[3];
rz(-0.47405973) q[3];
sx q[3];
rz(-2.1366675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1423433) q[2];
sx q[2];
rz(-1.3331022) q[2];
sx q[2];
rz(-0.86714253) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(-1.999202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2808696) q[0];
sx q[0];
rz(-2.6779802) q[0];
sx q[0];
rz(3.0238357) q[0];
rz(0.87751687) q[1];
sx q[1];
rz(-2.120647) q[1];
sx q[1];
rz(-0.95474517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60459899) q[0];
sx q[0];
rz(-2.1736896) q[0];
sx q[0];
rz(2.0379809) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81127848) q[2];
sx q[2];
rz(-2.5015321) q[2];
sx q[2];
rz(-2.0612353) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92312557) q[1];
sx q[1];
rz(-0.72281853) q[1];
sx q[1];
rz(1.9002401) q[1];
x q[2];
rz(-2.3225962) q[3];
sx q[3];
rz(-1.0291417) q[3];
sx q[3];
rz(-1.1198695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8942326) q[2];
sx q[2];
rz(-0.33984137) q[2];
sx q[2];
rz(-1.6575238) q[2];
rz(-1.5225211) q[3];
sx q[3];
rz(-1.3831474) q[3];
sx q[3];
rz(-1.9671666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5574944) q[0];
sx q[0];
rz(-1.7099986) q[0];
sx q[0];
rz(-0.75310055) q[0];
rz(-1.9901216) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(1.6023191) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70645151) q[0];
sx q[0];
rz(-1.1223842) q[0];
sx q[0];
rz(0.20332341) q[0];
rz(-2.2686917) q[2];
sx q[2];
rz(-1.1660327) q[2];
sx q[2];
rz(-1.0412316) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5394314) q[1];
sx q[1];
rz(-1.185158) q[1];
sx q[1];
rz(1.3630609) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3063823) q[3];
sx q[3];
rz(-0.9328649) q[3];
sx q[3];
rz(2.8273945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7842385) q[2];
sx q[2];
rz(-2.9604762) q[2];
sx q[2];
rz(0.46585807) q[2];
rz(-2.0458131) q[3];
sx q[3];
rz(-0.992479) q[3];
sx q[3];
rz(-0.54212681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.096238484) q[0];
sx q[0];
rz(-1.5753373) q[0];
sx q[0];
rz(-1.5684431) q[0];
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
rz(1.0505843) q[3];
sx q[3];
rz(-2.0512085) q[3];
sx q[3];
rz(-2.7762085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
