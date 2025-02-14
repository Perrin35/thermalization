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
rz(2.0514367) q[0];
sx q[0];
rz(-3.0213254) q[0];
sx q[0];
rz(0.076844849) q[0];
rz(2.285217) q[1];
sx q[1];
rz(5.2353274) q[1];
sx q[1];
rz(8.6040001) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79362291) q[0];
sx q[0];
rz(-1.4948778) q[0];
sx q[0];
rz(2.6708219) q[0];
rz(1.3581543) q[2];
sx q[2];
rz(-0.66591381) q[2];
sx q[2];
rz(-2.7261811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32626881) q[1];
sx q[1];
rz(-2.2827049) q[1];
sx q[1];
rz(2.2253042) q[1];
x q[2];
rz(-1.3935116) q[3];
sx q[3];
rz(-1.0330135) q[3];
sx q[3];
rz(3.1238525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0277752) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(1.6119733) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(0.99942708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73794) q[0];
sx q[0];
rz(-2.0495179) q[0];
sx q[0];
rz(1.8988443) q[0];
rz(-2.2620762) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(1.8125199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14318109) q[0];
sx q[0];
rz(-0.77309009) q[0];
sx q[0];
rz(-2.4804658) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0602559) q[2];
sx q[2];
rz(-1.1266409) q[2];
sx q[2];
rz(-1.3655438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14790598) q[1];
sx q[1];
rz(-2.4588486) q[1];
sx q[1];
rz(0.77862378) q[1];
rz(-pi) q[2];
rz(-1.4089083) q[3];
sx q[3];
rz(-2.257785) q[3];
sx q[3];
rz(1.5713728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30895147) q[2];
sx q[2];
rz(-1.917104) q[2];
sx q[2];
rz(-1.2844757) q[2];
rz(2.4744611) q[3];
sx q[3];
rz(-2.7693373) q[3];
sx q[3];
rz(2.4767806) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852378) q[0];
sx q[0];
rz(-0.41362718) q[0];
sx q[0];
rz(1.941823) q[0];
rz(1.2780227) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(-0.75495458) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4729378) q[0];
sx q[0];
rz(-1.1140214) q[0];
sx q[0];
rz(-3.1371389) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9725482) q[2];
sx q[2];
rz(-0.72578207) q[2];
sx q[2];
rz(-0.58453416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7717786) q[1];
sx q[1];
rz(-1.6057456) q[1];
sx q[1];
rz(-0.82280226) q[1];
rz(-pi) q[2];
rz(-0.6471031) q[3];
sx q[3];
rz(-0.80296338) q[3];
sx q[3];
rz(2.3784007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7167012) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(1.971604) q[2];
rz(0.26015002) q[3];
sx q[3];
rz(-1.6157506) q[3];
sx q[3];
rz(-2.9366176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78885704) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(-0.66057551) q[0];
rz(-0.13064101) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(2.7052243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6141967) q[0];
sx q[0];
rz(-2.578311) q[0];
sx q[0];
rz(-2.0296996) q[0];
x q[1];
rz(0.5882259) q[2];
sx q[2];
rz(-0.67413051) q[2];
sx q[2];
rz(2.9720486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.032724) q[1];
sx q[1];
rz(-2.0744262) q[1];
sx q[1];
rz(-0.41028604) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.800501) q[3];
sx q[3];
rz(-1.6636831) q[3];
sx q[3];
rz(-0.22985458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.074162) q[2];
sx q[2];
rz(-2.6219411) q[2];
sx q[2];
rz(1.1379498) q[2];
rz(-0.84749311) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.725089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6013019) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(0.20703319) q[0];
rz(-1.3218468) q[1];
sx q[1];
rz(-1.9925947) q[1];
sx q[1];
rz(-2.102899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.032271) q[0];
sx q[0];
rz(-2.3242852) q[0];
sx q[0];
rz(-0.82839806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62049753) q[2];
sx q[2];
rz(-1.7417241) q[2];
sx q[2];
rz(-0.88384274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6398479) q[1];
sx q[1];
rz(-1.5875438) q[1];
sx q[1];
rz(2.5052951) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30666017) q[3];
sx q[3];
rz(-1.6741317) q[3];
sx q[3];
rz(0.54319004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98628226) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(0.60688621) q[2];
rz(-0.33489975) q[3];
sx q[3];
rz(-1.5891985) q[3];
sx q[3];
rz(1.5437532) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63810054) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(2.959429) q[0];
rz(-2.4496574) q[1];
sx q[1];
rz(-1.7057995) q[1];
sx q[1];
rz(-1.8686132) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385524) q[0];
sx q[0];
rz(-1.5219757) q[0];
sx q[0];
rz(3.131034) q[0];
rz(-1.7164549) q[2];
sx q[2];
rz(-1.7985059) q[2];
sx q[2];
rz(2.6109254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23964393) q[1];
sx q[1];
rz(-1.8297503) q[1];
sx q[1];
rz(-2.3560897) q[1];
x q[2];
rz(-2.7837378) q[3];
sx q[3];
rz(-1.5113514) q[3];
sx q[3];
rz(-0.72261963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89340297) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(-0.87598652) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(2.7178154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.9559602) q[0];
sx q[0];
rz(-2.1946588) q[0];
sx q[0];
rz(-0.53674269) q[0];
rz(2.5770889) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(-1.5435262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604328) q[0];
sx q[0];
rz(-0.33012046) q[0];
sx q[0];
rz(0.31417234) q[0];
x q[1];
rz(0.44223152) q[2];
sx q[2];
rz(-2.6870607) q[2];
sx q[2];
rz(1.0584128) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1607223) q[1];
sx q[1];
rz(-2.0288784) q[1];
sx q[1];
rz(-2.0715817) q[1];
rz(-pi) q[2];
rz(0.48637511) q[3];
sx q[3];
rz(-2.1280299) q[3];
sx q[3];
rz(-0.0096021113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87968612) q[2];
sx q[2];
rz(-1.8330816) q[2];
sx q[2];
rz(-1.3685263) q[2];
rz(-2.8241099) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(0.69664636) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9484061) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(1.0411221) q[0];
rz(0.55024534) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(0.1951898) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0597417) q[0];
sx q[0];
rz(-1.7072816) q[0];
sx q[0];
rz(-1.5166212) q[0];
rz(1.3615029) q[2];
sx q[2];
rz(-2.6511764) q[2];
sx q[2];
rz(1.121466) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.079551815) q[1];
sx q[1];
rz(-1.5925381) q[1];
sx q[1];
rz(1.8575536) q[1];
rz(-pi) q[2];
rz(0.82542586) q[3];
sx q[3];
rz(-1.22681) q[3];
sx q[3];
rz(2.0478562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0027085) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(1.8769334) q[2];
rz(-1.5876611) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(-1.9896265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3569262) q[0];
sx q[0];
rz(-2.9279828) q[0];
sx q[0];
rz(0.40823644) q[0];
rz(-1.4534072) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(-0.85809392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5018117) q[0];
sx q[0];
rz(-0.94450356) q[0];
sx q[0];
rz(-0.39299742) q[0];
rz(-pi) q[1];
rz(1.4420322) q[2];
sx q[2];
rz(-1.6673281) q[2];
sx q[2];
rz(-1.0573204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.063805997) q[1];
sx q[1];
rz(-1.1764472) q[1];
sx q[1];
rz(1.1358244) q[1];
rz(-pi) q[2];
rz(-1.7210049) q[3];
sx q[3];
rz(-2.4433063) q[3];
sx q[3];
rz(2.5774389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66576362) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(-2.8909454) q[2];
rz(0.88746873) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(2.7953776) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733521) q[0];
sx q[0];
rz(-2.1526985) q[0];
sx q[0];
rz(-2.2762779) q[0];
rz(-0.52262753) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(1.685198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9495372) q[0];
sx q[0];
rz(-2.9350781) q[0];
sx q[0];
rz(-1.2386991) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7865042) q[2];
sx q[2];
rz(-0.68195365) q[2];
sx q[2];
rz(-2.8485659) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.870112) q[1];
sx q[1];
rz(-2.5462954) q[1];
sx q[1];
rz(-0.65956302) q[1];
rz(-pi) q[2];
rz(2.6271712) q[3];
sx q[3];
rz(-1.2687917) q[3];
sx q[3];
rz(0.3946886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7437848) q[2];
sx q[2];
rz(-2.2297771) q[2];
sx q[2];
rz(-1.816681) q[2];
rz(-1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(-2.3769784) q[3];
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
rz(-2.926173) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(2.678395) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(1.8824957) q[2];
sx q[2];
rz(-2.1915956) q[2];
sx q[2];
rz(-1.6804463) q[2];
rz(-0.15179612) q[3];
sx q[3];
rz(-0.67852466) q[3];
sx q[3];
rz(-1.8192374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
