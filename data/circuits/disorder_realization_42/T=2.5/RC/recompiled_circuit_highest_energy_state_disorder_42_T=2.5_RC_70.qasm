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
rz(3.2618599) q[0];
sx q[0];
rz(9.5016228) q[0];
rz(2.285217) q[1];
sx q[1];
rz(5.2353274) q[1];
sx q[1];
rz(8.6040001) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79362291) q[0];
sx q[0];
rz(-1.6467148) q[0];
sx q[0];
rz(0.47077076) q[0];
x q[1];
rz(-1.3581543) q[2];
sx q[2];
rz(-2.4756788) q[2];
sx q[2];
rz(-2.7261811) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3617435) q[1];
sx q[1];
rz(-2.0498154) q[1];
sx q[1];
rz(0.82734014) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5968339) q[3];
sx q[3];
rz(-1.4187471) q[3];
sx q[3];
rz(-1.6800547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0277752) q[2];
sx q[2];
rz(-0.57089388) q[2];
sx q[2];
rz(-1.6119733) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(-2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73794) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(1.8988443) q[0];
rz(-2.2620762) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(-1.3290728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984116) q[0];
sx q[0];
rz(-0.77309009) q[0];
sx q[0];
rz(0.66112681) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3622386) q[2];
sx q[2];
rz(-2.4930304) q[2];
sx q[2];
rz(-2.2575602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0893007) q[1];
sx q[1];
rz(-2.0366021) q[1];
sx q[1];
rz(-2.0897081) q[1];
rz(-pi) q[2];
rz(-2.4481421) q[3];
sx q[3];
rz(-1.4458522) q[3];
sx q[3];
rz(0.10378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30895147) q[2];
sx q[2];
rz(-1.917104) q[2];
sx q[2];
rz(-1.2844757) q[2];
rz(-0.66713157) q[3];
sx q[3];
rz(-0.37225538) q[3];
sx q[3];
rz(0.66481203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-0.27690241) q[1];
sx q[1];
rz(-2.3866381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4729378) q[0];
sx q[0];
rz(-2.0275712) q[0];
sx q[0];
rz(-3.1371389) q[0];
rz(0.16904449) q[2];
sx q[2];
rz(-0.72578207) q[2];
sx q[2];
rz(-2.5570585) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36981407) q[1];
sx q[1];
rz(-1.6057456) q[1];
sx q[1];
rz(0.82280226) q[1];
x q[2];
rz(2.4509804) q[3];
sx q[3];
rz(-1.1221782) q[3];
sx q[3];
rz(1.2909362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42489147) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(-1.1699886) q[2];
rz(-2.8814426) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(2.9366176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78885704) q[0];
sx q[0];
rz(-1.641888) q[0];
sx q[0];
rz(2.4810171) q[0];
rz(0.13064101) q[1];
sx q[1];
rz(-0.6159997) q[1];
sx q[1];
rz(2.7052243) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.527396) q[0];
sx q[0];
rz(-0.56328162) q[0];
sx q[0];
rz(2.0296996) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9881093) q[2];
sx q[2];
rz(-2.1168323) q[2];
sx q[2];
rz(-0.87620902) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.032724) q[1];
sx q[1];
rz(-2.0744262) q[1];
sx q[1];
rz(2.7313066) q[1];
rz(-pi) q[2];
rz(-2.8699974) q[3];
sx q[3];
rz(-2.7885572) q[3];
sx q[3];
rz(-2.0562382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.067430647) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(1.1379498) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(1.4165037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6013019) q[0];
sx q[0];
rz(-1.8922946) q[0];
sx q[0];
rz(-0.20703319) q[0];
rz(-1.8197458) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(1.0386937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0393676) q[0];
sx q[0];
rz(-1.0034585) q[0];
sx q[0];
rz(-0.62444425) q[0];
rz(0.62049753) q[2];
sx q[2];
rz(-1.7417241) q[2];
sx q[2];
rz(-0.88384274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0601692) q[1];
sx q[1];
rz(-2.2069903) q[1];
sx q[1];
rz(1.549975) q[1];
x q[2];
rz(-2.8349325) q[3];
sx q[3];
rz(-1.467461) q[3];
sx q[3];
rz(-2.5984026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.98628226) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(-2.5347064) q[2];
rz(-0.33489975) q[3];
sx q[3];
rz(-1.5891985) q[3];
sx q[3];
rz(-1.5978395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5034921) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(0.18216369) q[0];
rz(2.4496574) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(1.2729794) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733213) q[0];
sx q[0];
rz(-1.5813424) q[0];
sx q[0];
rz(-1.521973) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9115306) q[2];
sx q[2];
rz(-1.7126691) q[2];
sx q[2];
rz(1.0070247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.581785) q[1];
sx q[1];
rz(-2.323287) q[1];
sx q[1];
rz(-0.35840987) q[1];
rz(-pi) q[2];
rz(0.35785488) q[3];
sx q[3];
rz(-1.5113514) q[3];
sx q[3];
rz(-0.72261963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2481897) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(-2.2656061) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(2.7178154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9559602) q[0];
sx q[0];
rz(-2.1946588) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(2.5770889) q[1];
sx q[1];
rz(-0.91465488) q[1];
sx q[1];
rz(1.5435262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7537345) q[0];
sx q[0];
rz(-1.470454) q[0];
sx q[0];
rz(-2.8265586) q[0];
x q[1];
rz(0.41588628) q[2];
sx q[2];
rz(-1.759811) q[2];
sx q[2];
rz(3.031446) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98087039) q[1];
sx q[1];
rz(-2.0288784) q[1];
sx q[1];
rz(-1.0700109) q[1];
rz(-0.95682896) q[3];
sx q[3];
rz(-1.162863) q[3];
sx q[3];
rz(-1.3077426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87968612) q[2];
sx q[2];
rz(-1.3085111) q[2];
sx q[2];
rz(-1.7730664) q[2];
rz(-2.8241099) q[3];
sx q[3];
rz(-1.8902794) q[3];
sx q[3];
rz(-0.69664636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19318652) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(2.1004706) q[0];
rz(0.55024534) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(-2.9464029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6379162) q[0];
sx q[0];
rz(-1.5171255) q[0];
sx q[0];
rz(3.0049094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3615029) q[2];
sx q[2];
rz(-2.6511764) q[2];
sx q[2];
rz(-1.121466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4176538) q[1];
sx q[1];
rz(-0.28755763) q[1];
sx q[1];
rz(-1.647526) q[1];
rz(-2.0567254) q[3];
sx q[3];
rz(-0.80683364) q[3];
sx q[3];
rz(-3.0148447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13888415) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(-1.8769334) q[2];
rz(1.5539315) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78466648) q[0];
sx q[0];
rz(-2.9279828) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(-1.4534072) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(2.2834987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2553697) q[0];
sx q[0];
rz(-2.4165034) q[0];
sx q[0];
rz(2.057647) q[0];
rz(-2.2169154) q[2];
sx q[2];
rz(-0.16077006) q[2];
sx q[2];
rz(1.9882261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.063805997) q[1];
sx q[1];
rz(-1.9651455) q[1];
sx q[1];
rz(1.1358244) q[1];
rz(-pi) q[2];
x q[2];
rz(3.016641) q[3];
sx q[3];
rz(-2.2596684) q[3];
sx q[3];
rz(0.75923789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66576362) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(2.8909454) q[2];
rz(-0.88746873) q[3];
sx q[3];
rz(-1.9991425) q[3];
sx q[3];
rz(-0.34621507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733521) q[0];
sx q[0];
rz(-2.1526985) q[0];
sx q[0];
rz(2.2762779) q[0];
rz(2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(-1.4563947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9495372) q[0];
sx q[0];
rz(-2.9350781) q[0];
sx q[0];
rz(-1.9028936) q[0];
rz(-pi) q[1];
rz(0.17205949) q[2];
sx q[2];
rz(-0.90751782) q[2];
sx q[2];
rz(3.1236529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48113808) q[1];
sx q[1];
rz(-2.0298971) q[1];
sx q[1];
rz(-1.1774239) q[1];
rz(-pi) q[2];
rz(-1.2271475) q[3];
sx q[3];
rz(-2.0598186) q[3];
sx q[3];
rz(1.0095613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7437848) q[2];
sx q[2];
rz(-2.2297771) q[2];
sx q[2];
rz(-1.816681) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(2.3769784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.926173) q[0];
sx q[0];
rz(-1.8126491) q[0];
sx q[0];
rz(1.5480315) q[0];
rz(-2.678395) q[1];
sx q[1];
rz(-1.6013655) q[1];
sx q[1];
rz(-0.2688437) q[1];
rz(-1.8824957) q[2];
sx q[2];
rz(-0.94999708) q[2];
sx q[2];
rz(1.4611464) q[2];
rz(-2.4687121) q[3];
sx q[3];
rz(-1.665848) q[3];
sx q[3];
rz(3.0116826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
