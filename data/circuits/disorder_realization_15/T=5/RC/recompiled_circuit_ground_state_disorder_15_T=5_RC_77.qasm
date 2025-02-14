OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48159596) q[0];
sx q[0];
rz(3.5591535) q[0];
sx q[0];
rz(9.4661718) q[0];
rz(-2.5490835) q[1];
sx q[1];
rz(-2.9046287) q[1];
sx q[1];
rz(-2.2466329) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0314908) q[0];
sx q[0];
rz(-0.43927017) q[0];
sx q[0];
rz(1.8861559) q[0];
rz(2.1092806) q[2];
sx q[2];
rz(-0.56864029) q[2];
sx q[2];
rz(1.2775354) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0022566) q[1];
sx q[1];
rz(-1.8426241) q[1];
sx q[1];
rz(2.8546531) q[1];
x q[2];
rz(-1.6144606) q[3];
sx q[3];
rz(-1.6796675) q[3];
sx q[3];
rz(-0.085069503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.71901739) q[2];
sx q[2];
rz(-0.97192478) q[2];
sx q[2];
rz(1.0389339) q[2];
rz(-0.59764189) q[3];
sx q[3];
rz(-0.20331764) q[3];
sx q[3];
rz(1.982127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0036302) q[0];
sx q[0];
rz(-2.2901386) q[0];
sx q[0];
rz(0.33388579) q[0];
rz(-2.4489898) q[1];
sx q[1];
rz(-1.875016) q[1];
sx q[1];
rz(-1.0612706) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740299) q[0];
sx q[0];
rz(-2.2252796) q[0];
sx q[0];
rz(0.45463134) q[0];
rz(-pi) q[1];
rz(0.40410903) q[2];
sx q[2];
rz(-1.3349018) q[2];
sx q[2];
rz(2.1020232) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088527352) q[1];
sx q[1];
rz(-0.47855166) q[1];
sx q[1];
rz(0.38691945) q[1];
rz(-pi) q[2];
rz(0.96978404) q[3];
sx q[3];
rz(-0.24479391) q[3];
sx q[3];
rz(0.56772619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65427238) q[2];
sx q[2];
rz(-1.6795748) q[2];
sx q[2];
rz(-0.1839323) q[2];
rz(-0.0025302689) q[3];
sx q[3];
rz(-0.80357426) q[3];
sx q[3];
rz(2.7185503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.84394395) q[0];
sx q[0];
rz(-0.14744814) q[0];
sx q[0];
rz(2.3609128) q[0];
rz(0.88790226) q[1];
sx q[1];
rz(-1.0941411) q[1];
sx q[1];
rz(-0.81781864) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.458815) q[0];
sx q[0];
rz(-1.11894) q[0];
sx q[0];
rz(-2.8917612) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70728126) q[2];
sx q[2];
rz(-0.92628252) q[2];
sx q[2];
rz(-1.4799581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0101658) q[1];
sx q[1];
rz(-1.5701811) q[1];
sx q[1];
rz(-0.58227957) q[1];
rz(0.45210101) q[3];
sx q[3];
rz(-1.107405) q[3];
sx q[3];
rz(2.2832561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4320977) q[2];
sx q[2];
rz(-0.12505394) q[2];
sx q[2];
rz(2.1091667) q[2];
rz(0.89140511) q[3];
sx q[3];
rz(-0.67474198) q[3];
sx q[3];
rz(-2.3064822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29435232) q[0];
sx q[0];
rz(-2.8645741) q[0];
sx q[0];
rz(2.5936122) q[0];
rz(3.0849988) q[1];
sx q[1];
rz(-1.2013925) q[1];
sx q[1];
rz(0.024554575) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5875754) q[0];
sx q[0];
rz(-1.2535106) q[0];
sx q[0];
rz(1.1486828) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33523021) q[2];
sx q[2];
rz(-1.9465228) q[2];
sx q[2];
rz(-0.59240985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78950846) q[1];
sx q[1];
rz(-2.2768436) q[1];
sx q[1];
rz(-0.45070453) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4931498) q[3];
sx q[3];
rz(-2.7973865) q[3];
sx q[3];
rz(1.6370809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6294488) q[2];
sx q[2];
rz(-0.63261837) q[2];
sx q[2];
rz(2.4150685) q[2];
rz(1.4350545) q[3];
sx q[3];
rz(-1.8972998) q[3];
sx q[3];
rz(-1.6632891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5479945) q[0];
sx q[0];
rz(-1.7578121) q[0];
sx q[0];
rz(-1.7543678) q[0];
rz(2.6902426) q[1];
sx q[1];
rz(-0.74140048) q[1];
sx q[1];
rz(-0.50419921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58012843) q[0];
sx q[0];
rz(-0.21004349) q[0];
sx q[0];
rz(2.5145636) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3875972) q[2];
sx q[2];
rz(-2.0903433) q[2];
sx q[2];
rz(-2.4343672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54795116) q[1];
sx q[1];
rz(-0.91386139) q[1];
sx q[1];
rz(2.2398578) q[1];
rz(-pi) q[2];
rz(-0.80662722) q[3];
sx q[3];
rz(-0.99116814) q[3];
sx q[3];
rz(-1.5416073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7826358) q[2];
sx q[2];
rz(-0.88023305) q[2];
sx q[2];
rz(3.1202988) q[2];
rz(3.0535789) q[3];
sx q[3];
rz(-0.26282495) q[3];
sx q[3];
rz(2.9737441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39040318) q[0];
sx q[0];
rz(-0.57325554) q[0];
sx q[0];
rz(2.7943352) q[0];
rz(-1.685453) q[1];
sx q[1];
rz(-0.29734722) q[1];
sx q[1];
rz(-0.28486326) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69377855) q[0];
sx q[0];
rz(-1.8577357) q[0];
sx q[0];
rz(1.8836431) q[0];
rz(-1.6256394) q[2];
sx q[2];
rz(-0.82717878) q[2];
sx q[2];
rz(1.0393927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9544475) q[1];
sx q[1];
rz(-1.3055077) q[1];
sx q[1];
rz(0.89098661) q[1];
x q[2];
rz(0.29355704) q[3];
sx q[3];
rz(-1.9133074) q[3];
sx q[3];
rz(-2.5409043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.138729) q[2];
sx q[2];
rz(-1.8314654) q[2];
sx q[2];
rz(2.1076473) q[2];
rz(-0.73505861) q[3];
sx q[3];
rz(-1.6481954) q[3];
sx q[3];
rz(-2.471931) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1351778) q[0];
sx q[0];
rz(-0.66295755) q[0];
sx q[0];
rz(-2.2354777) q[0];
rz(0.16618973) q[1];
sx q[1];
rz(-0.48887417) q[1];
sx q[1];
rz(-0.30950549) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85649895) q[0];
sx q[0];
rz(-1.4011607) q[0];
sx q[0];
rz(-1.7986675) q[0];
x q[1];
rz(-0.79321547) q[2];
sx q[2];
rz(-1.4275107) q[2];
sx q[2];
rz(-2.6806841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55927235) q[1];
sx q[1];
rz(-2.4666335) q[1];
sx q[1];
rz(1.9155424) q[1];
rz(-pi) q[2];
rz(2.0887435) q[3];
sx q[3];
rz(-0.098909698) q[3];
sx q[3];
rz(-0.78341752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6918148) q[2];
sx q[2];
rz(-0.66902995) q[2];
sx q[2];
rz(-1.4867894) q[2];
rz(-0.4977704) q[3];
sx q[3];
rz(-0.83511746) q[3];
sx q[3];
rz(-0.2275137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787702) q[0];
sx q[0];
rz(-2.4202122) q[0];
sx q[0];
rz(2.873514) q[0];
rz(-2.2396741) q[1];
sx q[1];
rz(-1.0682769) q[1];
sx q[1];
rz(1.8394151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3288181) q[0];
sx q[0];
rz(-0.98085603) q[0];
sx q[0];
rz(1.3716212) q[0];
rz(1.5728358) q[2];
sx q[2];
rz(-1.1241978) q[2];
sx q[2];
rz(2.0799321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5561192) q[1];
sx q[1];
rz(-0.98656174) q[1];
sx q[1];
rz(0.16884905) q[1];
rz(-pi) q[2];
rz(-1.4623649) q[3];
sx q[3];
rz(-2.83395) q[3];
sx q[3];
rz(2.3189312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7997718) q[2];
sx q[2];
rz(-1.2398961) q[2];
sx q[2];
rz(-0.58490252) q[2];
rz(-1.2029485) q[3];
sx q[3];
rz(-1.7479618) q[3];
sx q[3];
rz(1.5929619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0840266) q[0];
sx q[0];
rz(-0.49087048) q[0];
sx q[0];
rz(0.67114818) q[0];
rz(-2.8533543) q[1];
sx q[1];
rz(-0.75575525) q[1];
sx q[1];
rz(2.5075358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6995981) q[0];
sx q[0];
rz(-1.6091248) q[0];
sx q[0];
rz(-0.015873578) q[0];
rz(1.4052275) q[2];
sx q[2];
rz(-1.6201303) q[2];
sx q[2];
rz(-2.0967576) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.523484) q[1];
sx q[1];
rz(-1.5963703) q[1];
sx q[1];
rz(0.32483883) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67231549) q[3];
sx q[3];
rz(-1.3642715) q[3];
sx q[3];
rz(0.12609161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59342283) q[2];
sx q[2];
rz(-1.5735184) q[2];
sx q[2];
rz(-0.28369743) q[2];
rz(0.13188322) q[3];
sx q[3];
rz(-2.7851084) q[3];
sx q[3];
rz(-0.62294817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.236096) q[0];
sx q[0];
rz(-0.14362366) q[0];
sx q[0];
rz(0.96963257) q[0];
rz(0.49599221) q[1];
sx q[1];
rz(-0.6876567) q[1];
sx q[1];
rz(-2.2775441) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0893893) q[0];
sx q[0];
rz(-1.6411575) q[0];
sx q[0];
rz(2.2785827) q[0];
rz(2.8882083) q[2];
sx q[2];
rz(-1.14865) q[2];
sx q[2];
rz(-2.6146023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3083777) q[1];
sx q[1];
rz(-1.0755441) q[1];
sx q[1];
rz(2.2332195) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0122369) q[3];
sx q[3];
rz(-0.73775916) q[3];
sx q[3];
rz(-0.078928909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8833017) q[2];
sx q[2];
rz(-2.624888) q[2];
sx q[2];
rz(-1.0374163) q[2];
rz(0.20710219) q[3];
sx q[3];
rz(-2.243302) q[3];
sx q[3];
rz(-2.5521539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0265738) q[0];
sx q[0];
rz(-1.6828231) q[0];
sx q[0];
rz(-1.9376301) q[0];
rz(-0.013982458) q[1];
sx q[1];
rz(-1.7054066) q[1];
sx q[1];
rz(-1.1048497) q[1];
rz(1.7183279) q[2];
sx q[2];
rz(-2.9083283) q[2];
sx q[2];
rz(1.8988594) q[2];
rz(-1.1963853) q[3];
sx q[3];
rz(-0.76210124) q[3];
sx q[3];
rz(-2.6740554) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
