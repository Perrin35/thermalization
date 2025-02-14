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
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(-2.8859477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50798049) q[0];
sx q[0];
rz(-2.2574212) q[0];
sx q[0];
rz(0.79721862) q[0];
x q[1];
rz(-2.9886888) q[2];
sx q[2];
rz(-0.39275482) q[2];
sx q[2];
rz(-0.037234779) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5140443) q[1];
sx q[1];
rz(-1.6395634) q[1];
sx q[1];
rz(-1.7538944) q[1];
rz(1.3575451) q[3];
sx q[3];
rz(-1.0825601) q[3];
sx q[3];
rz(-2.8247716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6306182) q[2];
sx q[2];
rz(-1.3773842) q[2];
sx q[2];
rz(-1.1154491) q[2];
rz(-3.1274146) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9344591) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(0.44218749) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(0.8173379) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38430518) q[0];
sx q[0];
rz(-1.7624859) q[0];
sx q[0];
rz(1.0377645) q[0];
rz(2.8645682) q[2];
sx q[2];
rz(-1.1454586) q[2];
sx q[2];
rz(-2.9556731) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0603588) q[1];
sx q[1];
rz(-2.4801755) q[1];
sx q[1];
rz(2.8702186) q[1];
rz(-pi) q[2];
rz(-2.3807008) q[3];
sx q[3];
rz(-2.0756654) q[3];
sx q[3];
rz(-0.99772108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77357972) q[2];
sx q[2];
rz(-0.37675884) q[2];
sx q[2];
rz(0.22029857) q[2];
rz(1.1822654) q[3];
sx q[3];
rz(-1.9672829) q[3];
sx q[3];
rz(-3.0784472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050506266) q[0];
sx q[0];
rz(-1.8244705) q[0];
sx q[0];
rz(2.0029946) q[0];
rz(-2.7298722) q[1];
sx q[1];
rz(-2.3717334) q[1];
sx q[1];
rz(1.0134816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6665812) q[0];
sx q[0];
rz(-1.5856167) q[0];
sx q[0];
rz(0.25304742) q[0];
x q[1];
rz(3.0092952) q[2];
sx q[2];
rz(-1.1331285) q[2];
sx q[2];
rz(-0.79298151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1057844) q[1];
sx q[1];
rz(-0.86584751) q[1];
sx q[1];
rz(0.99143274) q[1];
x q[2];
rz(-1.7168379) q[3];
sx q[3];
rz(-2.4096294) q[3];
sx q[3];
rz(-0.41106658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1608405) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(1.8864924) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(0.14061558) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18144064) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(2.5312359) q[0];
rz(-0.70862526) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(-1.5708539) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7371668) q[0];
sx q[0];
rz(-0.20658399) q[0];
sx q[0];
rz(-0.73295672) q[0];
rz(-1.6897292) q[2];
sx q[2];
rz(-1.2858675) q[2];
sx q[2];
rz(-0.97078427) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29250568) q[1];
sx q[1];
rz(-1.6189112) q[1];
sx q[1];
rz(-2.8699257) q[1];
x q[2];
rz(-0.80067486) q[3];
sx q[3];
rz(-1.8599556) q[3];
sx q[3];
rz(-1.9137933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9307956) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(2.5861758) q[2];
rz(2.4173315) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-2.485062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.637218) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(-2.8952307) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(-1.4341644) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61489366) q[0];
sx q[0];
rz(-1.5759567) q[0];
sx q[0];
rz(-0.75684383) q[0];
rz(-pi) q[1];
rz(0.44539659) q[2];
sx q[2];
rz(-1.4120308) q[2];
sx q[2];
rz(-0.36728718) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.757548) q[1];
sx q[1];
rz(-1.9725807) q[1];
sx q[1];
rz(-1.5753788) q[1];
rz(-pi) q[2];
rz(-1.9798109) q[3];
sx q[3];
rz(-2.8944765) q[3];
sx q[3];
rz(0.97774001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(-2.0595713) q[2];
rz(-0.79484445) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.865888) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(0.24359447) q[0];
rz(2.1547735) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(-2.2056244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2631131) q[0];
sx q[0];
rz(-1.8917483) q[0];
sx q[0];
rz(1.1522989) q[0];
rz(-pi) q[1];
rz(-2.317993) q[2];
sx q[2];
rz(-1.8317501) q[2];
sx q[2];
rz(1.5023155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.080612506) q[1];
sx q[1];
rz(-0.83773936) q[1];
sx q[1];
rz(2.5886646) q[1];
rz(0.021401568) q[3];
sx q[3];
rz(-0.97107065) q[3];
sx q[3];
rz(-2.3779526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.662107) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(2.7916419) q[2];
rz(-2.5838666) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(-0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.4418942) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(2.8398474) q[0];
rz(0.31983495) q[1];
sx q[1];
rz(-1.6022857) q[1];
sx q[1];
rz(2.005827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0253191) q[0];
sx q[0];
rz(-0.92705446) q[0];
sx q[0];
rz(1.6380861) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93923969) q[2];
sx q[2];
rz(-1.1432836) q[2];
sx q[2];
rz(-1.7981426) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3150683) q[1];
sx q[1];
rz(-0.60501912) q[1];
sx q[1];
rz(1.4304763) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3507422) q[3];
sx q[3];
rz(-1.6501657) q[3];
sx q[3];
rz(-2.1805003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4096421) q[2];
sx q[2];
rz(-2.4739517) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(-1.7536633) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.1801572) q[0];
sx q[0];
rz(-3.1249983) q[0];
sx q[0];
rz(-2.5855682) q[0];
rz(0.22932886) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(-1.75846) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97998226) q[0];
sx q[0];
rz(-1.5357657) q[0];
sx q[0];
rz(-2.4039414) q[0];
rz(0.11923125) q[2];
sx q[2];
rz(-2.6802353) q[2];
sx q[2];
rz(-0.90864321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.759093) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(0.31633693) q[1];
rz(-pi) q[2];
rz(2.4501958) q[3];
sx q[3];
rz(-0.56400877) q[3];
sx q[3];
rz(0.15182748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51656276) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(-1.2410835) q[2];
rz(-0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(-0.82714287) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271707) q[0];
sx q[0];
rz(-1.4143455) q[0];
sx q[0];
rz(2.993809) q[0];
rz(2.6328909) q[1];
sx q[1];
rz(-1.3721507) q[1];
sx q[1];
rz(0.37857372) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62703778) q[0];
sx q[0];
rz(-0.48279027) q[0];
sx q[0];
rz(-0.52282368) q[0];
rz(1.0580334) q[2];
sx q[2];
rz(-0.65750018) q[2];
sx q[2];
rz(0.12557827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8151103) q[1];
sx q[1];
rz(-1.2604144) q[1];
sx q[1];
rz(-1.2182477) q[1];
rz(-pi) q[2];
rz(1.6949953) q[3];
sx q[3];
rz(-2.327965) q[3];
sx q[3];
rz(0.99909335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.96616894) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.8625205) q[2];
rz(-1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(0.14555791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3928423) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(0.064099126) q[0];
rz(0.52866689) q[1];
sx q[1];
rz(-1.8730947) q[1];
sx q[1];
rz(-0.97506964) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.753016) q[0];
sx q[0];
rz(-0.56682366) q[0];
sx q[0];
rz(0.50811572) q[0];
x q[1];
rz(-0.80143546) q[2];
sx q[2];
rz(-2.3972627) q[2];
sx q[2];
rz(2.3542885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8243557) q[1];
sx q[1];
rz(-1.1350766) q[1];
sx q[1];
rz(1.812326) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75210877) q[3];
sx q[3];
rz(-2.8127713) q[3];
sx q[3];
rz(-3.0260835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9958682) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(-2.5346942) q[2];
rz(-2.8912344) q[3];
sx q[3];
rz(-1.7253877) q[3];
sx q[3];
rz(2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470227) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(-2.6429214) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(-0.90169883) q[2];
sx q[2];
rz(-2.31291) q[2];
sx q[2];
rz(-0.55533021) q[2];
rz(-0.12228431) q[3];
sx q[3];
rz(-1.7752083) q[3];
sx q[3];
rz(1.0158252) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
