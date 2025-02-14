OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1801017) q[0];
sx q[0];
rz(-0.56668007) q[0];
sx q[0];
rz(2.3508747) q[0];
rz(-0.92791954) q[1];
sx q[1];
rz(-0.035049573) q[1];
sx q[1];
rz(0.30326453) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68214455) q[0];
sx q[0];
rz(-1.3431699) q[0];
sx q[0];
rz(-2.4683003) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.319807) q[2];
sx q[2];
rz(-2.6394834) q[2];
sx q[2];
rz(-0.443845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3825079) q[1];
sx q[1];
rz(-2.0469991) q[1];
sx q[1];
rz(-2.166283) q[1];
x q[2];
rz(2.4062626) q[3];
sx q[3];
rz(-1.6323252) q[3];
sx q[3];
rz(1.3593591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7683893) q[2];
sx q[2];
rz(-2.781227) q[2];
sx q[2];
rz(-2.4599794) q[2];
rz(2.5810589) q[3];
sx q[3];
rz(-2.3571641) q[3];
sx q[3];
rz(2.4973629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.949837) q[0];
sx q[0];
rz(-2.2332709) q[0];
sx q[0];
rz(-0.35582304) q[0];
rz(-0.25335723) q[1];
sx q[1];
rz(-0.8496049) q[1];
sx q[1];
rz(-1.2454978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6424774) q[0];
sx q[0];
rz(-2.0428951) q[0];
sx q[0];
rz(1.319122) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63388692) q[2];
sx q[2];
rz(-1.1340464) q[2];
sx q[2];
rz(0.89981198) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51735866) q[1];
sx q[1];
rz(-1.1453724) q[1];
sx q[1];
rz(-2.5555771) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50470501) q[3];
sx q[3];
rz(-1.0796121) q[3];
sx q[3];
rz(2.1757656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3972828) q[2];
sx q[2];
rz(-0.77703589) q[2];
sx q[2];
rz(0.12766078) q[2];
rz(-2.4658261) q[3];
sx q[3];
rz(-2.6830169) q[3];
sx q[3];
rz(2.7202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965939) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(-1.1340207) q[0];
rz(-2.2484312) q[1];
sx q[1];
rz(-2.2375219) q[1];
sx q[1];
rz(1.2718511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29933481) q[0];
sx q[0];
rz(-2.133373) q[0];
sx q[0];
rz(-0.03118731) q[0];
rz(-1.3273622) q[2];
sx q[2];
rz(-0.31844246) q[2];
sx q[2];
rz(-1.8093579) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33128502) q[1];
sx q[1];
rz(-0.87281094) q[1];
sx q[1];
rz(0.26878243) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16678129) q[3];
sx q[3];
rz(-1.7875009) q[3];
sx q[3];
rz(-3.0971017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9280055) q[2];
sx q[2];
rz(-0.86557937) q[2];
sx q[2];
rz(-1.8013087) q[2];
rz(-1.853893) q[3];
sx q[3];
rz(-1.4475334) q[3];
sx q[3];
rz(-2.652216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5260148) q[0];
sx q[0];
rz(-0.57681334) q[0];
sx q[0];
rz(-0.73337972) q[0];
rz(-2.3461657) q[1];
sx q[1];
rz(-2.9454102) q[1];
sx q[1];
rz(-1.1682074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35287898) q[0];
sx q[0];
rz(-1.1601935) q[0];
sx q[0];
rz(-3.0263508) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2890786) q[2];
sx q[2];
rz(-2.3688683) q[2];
sx q[2];
rz(2.1022405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2760022) q[1];
sx q[1];
rz(-1.6011213) q[1];
sx q[1];
rz(-3.0986165) q[1];
x q[2];
rz(-2.4766604) q[3];
sx q[3];
rz(-2.0789609) q[3];
sx q[3];
rz(0.45979655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9613793) q[2];
sx q[2];
rz(-0.81575477) q[2];
sx q[2];
rz(-0.95480314) q[2];
rz(-1.0500326) q[3];
sx q[3];
rz(-2.3539216) q[3];
sx q[3];
rz(-0.55154705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.10057218) q[0];
sx q[0];
rz(-2.6028778) q[0];
sx q[0];
rz(-2.4464497) q[0];
rz(-1.9255385) q[1];
sx q[1];
rz(-1.0168409) q[1];
sx q[1];
rz(-0.020811828) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6489844) q[0];
sx q[0];
rz(-2.2243315) q[0];
sx q[0];
rz(0.18006353) q[0];
x q[1];
rz(-0.32692744) q[2];
sx q[2];
rz(-2.280664) q[2];
sx q[2];
rz(-0.29151379) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6823618) q[1];
sx q[1];
rz(-1.8098847) q[1];
sx q[1];
rz(0.011929529) q[1];
rz(0.5848123) q[3];
sx q[3];
rz(-3.0303114) q[3];
sx q[3];
rz(-2.5113784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1702599) q[2];
sx q[2];
rz(-2.8016165) q[2];
sx q[2];
rz(2.5374832) q[2];
rz(0.64770925) q[3];
sx q[3];
rz(-2.6205687) q[3];
sx q[3];
rz(0.46085301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.0172417) q[0];
sx q[0];
rz(-2.0483973) q[0];
sx q[0];
rz(-0.35853115) q[0];
rz(-0.56309807) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(0.30555746) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9193076) q[0];
sx q[0];
rz(-2.3940175) q[0];
sx q[0];
rz(-0.54493746) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9592909) q[2];
sx q[2];
rz(-1.4023702) q[2];
sx q[2];
rz(1.7163079) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2626052) q[1];
sx q[1];
rz(-2.0630347) q[1];
sx q[1];
rz(2.5237971) q[1];
x q[2];
rz(-2.4044068) q[3];
sx q[3];
rz(-3.0372826) q[3];
sx q[3];
rz(2.7953469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2629338) q[2];
sx q[2];
rz(-0.96076751) q[2];
sx q[2];
rz(0.052138694) q[2];
rz(-0.42929286) q[3];
sx q[3];
rz(-1.41058) q[3];
sx q[3];
rz(2.8894292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4572064) q[0];
sx q[0];
rz(-0.52427137) q[0];
sx q[0];
rz(-0.23505178) q[0];
rz(-0.60335195) q[1];
sx q[1];
rz(-0.82913202) q[1];
sx q[1];
rz(-2.5650909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7991195) q[0];
sx q[0];
rz(-0.89460374) q[0];
sx q[0];
rz(-0.66838092) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7965661) q[2];
sx q[2];
rz(-2.3357311) q[2];
sx q[2];
rz(-2.8553183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.723878) q[1];
sx q[1];
rz(-2.6195994) q[1];
sx q[1];
rz(0.88845171) q[1];
x q[2];
rz(2.0527538) q[3];
sx q[3];
rz(-1.4432943) q[3];
sx q[3];
rz(-1.3393928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.06360402) q[2];
sx q[2];
rz(-2.3006907) q[2];
sx q[2];
rz(-1.3464751) q[2];
rz(0.50955647) q[3];
sx q[3];
rz(-1.8815123) q[3];
sx q[3];
rz(-2.8349561) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52617306) q[0];
sx q[0];
rz(-2.058448) q[0];
sx q[0];
rz(0.35588595) q[0];
rz(-0.7016167) q[1];
sx q[1];
rz(-0.18467782) q[1];
sx q[1];
rz(-2.1772749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0052781) q[0];
sx q[0];
rz(-1.1169254) q[0];
sx q[0];
rz(-0.14880609) q[0];
rz(-pi) q[1];
rz(1.6379328) q[2];
sx q[2];
rz(-0.29342857) q[2];
sx q[2];
rz(-0.98240438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0605392) q[1];
sx q[1];
rz(-0.49526981) q[1];
sx q[1];
rz(-3.1002863) q[1];
rz(-pi) q[2];
rz(2.8921732) q[3];
sx q[3];
rz(-0.32085996) q[3];
sx q[3];
rz(1.3414931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4926766) q[2];
sx q[2];
rz(-2.4339088) q[2];
sx q[2];
rz(2.3259582) q[2];
rz(-0.91551578) q[3];
sx q[3];
rz(-2.2736277) q[3];
sx q[3];
rz(2.9149132) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30064073) q[0];
sx q[0];
rz(-0.17423593) q[0];
sx q[0];
rz(-1.9006282) q[0];
rz(3.0366483) q[1];
sx q[1];
rz(-2.8009156) q[1];
sx q[1];
rz(-2.6822283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7631885) q[0];
sx q[0];
rz(-1.8109461) q[0];
sx q[0];
rz(1.9394919) q[0];
rz(-1.2392339) q[2];
sx q[2];
rz(-1.825983) q[2];
sx q[2];
rz(2.55748) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85377598) q[1];
sx q[1];
rz(-0.48928919) q[1];
sx q[1];
rz(0.46238203) q[1];
rz(-pi) q[2];
rz(-0.21441238) q[3];
sx q[3];
rz(-0.58480703) q[3];
sx q[3];
rz(2.6406276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0910054) q[2];
sx q[2];
rz(-0.688474) q[2];
sx q[2];
rz(0.090593226) q[2];
rz(-0.76720864) q[3];
sx q[3];
rz(-1.4827261) q[3];
sx q[3];
rz(-1.4513133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7200658) q[0];
sx q[0];
rz(-0.92362112) q[0];
sx q[0];
rz(1.5371171) q[0];
rz(0.9556669) q[1];
sx q[1];
rz(-1.6803398) q[1];
sx q[1];
rz(-2.59424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065009251) q[0];
sx q[0];
rz(-3.0616548) q[0];
sx q[0];
rz(1.4968833) q[0];
x q[1];
rz(2.5432406) q[2];
sx q[2];
rz(-0.6171591) q[2];
sx q[2];
rz(-1.1197108) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9319084) q[1];
sx q[1];
rz(-0.80291498) q[1];
sx q[1];
rz(-1.9002537) q[1];
x q[2];
rz(0.71228446) q[3];
sx q[3];
rz(-1.1761208) q[3];
sx q[3];
rz(-2.8443732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0305816) q[2];
sx q[2];
rz(-2.9647398) q[2];
sx q[2];
rz(-2.492823) q[2];
rz(0.23160058) q[3];
sx q[3];
rz(-0.88445556) q[3];
sx q[3];
rz(3.054936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34407525) q[0];
sx q[0];
rz(-1.9301013) q[0];
sx q[0];
rz(1.7420266) q[0];
rz(-2.4823785) q[1];
sx q[1];
rz(-1.2792239) q[1];
sx q[1];
rz(-1.4469133) q[1];
rz(1.7805889) q[2];
sx q[2];
rz(-2.9083512) q[2];
sx q[2];
rz(-0.19610263) q[2];
rz(1.9271196) q[3];
sx q[3];
rz(-0.15542843) q[3];
sx q[3];
rz(-0.036416362) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
