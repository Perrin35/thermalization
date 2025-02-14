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
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(0.82077789) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92511237) q[0];
sx q[0];
rz(-0.47639936) q[0];
sx q[0];
rz(0.16615482) q[0];
x q[1];
rz(-2.2256741) q[2];
sx q[2];
rz(-1.7015463) q[2];
sx q[2];
rz(-2.1543825) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1917853) q[1];
sx q[1];
rz(-2.2151183) q[1];
sx q[1];
rz(0.61442805) q[1];
rz(2.5968339) q[3];
sx q[3];
rz(-1.4187471) q[3];
sx q[3];
rz(-1.6800547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1138175) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(-1.5296193) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(-2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40365264) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(1.8988443) q[0];
rz(2.2620762) q[1];
sx q[1];
rz(-0.95006919) q[1];
sx q[1];
rz(1.8125199) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97016818) q[0];
sx q[0];
rz(-0.98698436) q[0];
sx q[0];
rz(-2.1105355) q[0];
x q[1];
rz(-0.49449253) q[2];
sx q[2];
rz(-1.1323511) q[2];
sx q[2];
rz(3.1218253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.14790598) q[1];
sx q[1];
rz(-0.68274409) q[1];
sx q[1];
rz(2.3629689) q[1];
rz(-pi) q[2];
rz(2.9475698) q[3];
sx q[3];
rz(-0.7027773) q[3];
sx q[3];
rz(1.3181835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30895147) q[2];
sx q[2];
rz(-1.917104) q[2];
sx q[2];
rz(-1.857117) q[2];
rz(-2.4744611) q[3];
sx q[3];
rz(-2.7693373) q[3];
sx q[3];
rz(-2.4767806) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852378) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(1.941823) q[0];
rz(1.86357) q[1];
sx q[1];
rz(-0.27690241) q[1];
sx q[1];
rz(2.3866381) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2374868) q[0];
sx q[0];
rz(-1.5747935) q[0];
sx q[0];
rz(-1.1140175) q[0];
rz(0.16904449) q[2];
sx q[2];
rz(-0.72578207) q[2];
sx q[2];
rz(0.58453416) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9082005) q[1];
sx q[1];
rz(-2.3182239) q[1];
sx q[1];
rz(-0.047659454) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0126013) q[3];
sx q[3];
rz(-0.95944384) q[3];
sx q[3];
rz(-0.064289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42489147) q[2];
sx q[2];
rz(-2.5548013) q[2];
sx q[2];
rz(-1.1699886) q[2];
rz(2.8814426) q[3];
sx q[3];
rz(-1.6157506) q[3];
sx q[3];
rz(-0.20497504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3527356) q[0];
sx q[0];
rz(-1.641888) q[0];
sx q[0];
rz(-0.66057551) q[0];
rz(0.13064101) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(-2.7052243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35231464) q[0];
sx q[0];
rz(-1.8095865) q[0];
sx q[0];
rz(2.0859857) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9881093) q[2];
sx q[2];
rz(-1.0247604) q[2];
sx q[2];
rz(0.87620902) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37533942) q[1];
sx q[1];
rz(-0.63823593) q[1];
sx q[1];
rz(-2.1973647) q[1];
x q[2];
rz(-1.6693258) q[3];
sx q[3];
rz(-1.9103582) q[3];
sx q[3];
rz(-1.7677386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.074162) q[2];
sx q[2];
rz(-2.6219411) q[2];
sx q[2];
rz(-1.1379498) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.7646101) q[3];
sx q[3];
rz(-1.4165037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(1.5402907) q[0];
sx q[0];
rz(-1.8922946) q[0];
sx q[0];
rz(-0.20703319) q[0];
rz(1.3218468) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(1.0386937) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0393676) q[0];
sx q[0];
rz(-1.0034585) q[0];
sx q[0];
rz(2.5171484) q[0];
rz(0.62049753) q[2];
sx q[2];
rz(-1.7417241) q[2];
sx q[2];
rz(-0.88384274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50174471) q[1];
sx q[1];
rz(-1.5540489) q[1];
sx q[1];
rz(2.5052951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33090277) q[3];
sx q[3];
rz(-0.32308137) q[3];
sx q[3];
rz(-0.71268717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1553104) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(-2.5347064) q[2];
rz(2.8066929) q[3];
sx q[3];
rz(-1.5891985) q[3];
sx q[3];
rz(1.5437532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63810054) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(0.18216369) q[0];
rz(-2.4496574) q[1];
sx q[1];
rz(-1.7057995) q[1];
sx q[1];
rz(1.2729794) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86827134) q[0];
sx q[0];
rz(-1.5602503) q[0];
sx q[0];
rz(-1.6196197) q[0];
rz(-pi) q[1];
rz(-1.4251377) q[2];
sx q[2];
rz(-1.3430867) q[2];
sx q[2];
rz(2.6109254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5598076) q[1];
sx q[1];
rz(-0.81830561) q[1];
sx q[1];
rz(2.7831828) q[1];
rz(-2.7837378) q[3];
sx q[3];
rz(-1.6302413) q[3];
sx q[3];
rz(0.72261963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2481897) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(2.2656061) q[2];
rz(-2.857483) q[3];
sx q[3];
rz(-2.1967389) q[3];
sx q[3];
rz(2.7178154) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1856325) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(-2.5770889) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(-1.5980665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3878582) q[0];
sx q[0];
rz(-1.470454) q[0];
sx q[0];
rz(0.31503408) q[0];
rz(1.7769496) q[2];
sx q[2];
rz(-1.9788303) q[2];
sx q[2];
rz(1.598151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7891414) q[1];
sx q[1];
rz(-1.1256212) q[1];
sx q[1];
rz(0.51207249) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2143977) q[3];
sx q[3];
rz(-2.4193086) q[3];
sx q[3];
rz(0.77570206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2619065) q[2];
sx q[2];
rz(-1.3085111) q[2];
sx q[2];
rz(1.3685263) q[2];
rz(0.3174828) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(-2.4449463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9484061) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(-1.0411221) q[0];
rz(-0.55024534) q[1];
sx q[1];
rz(-2.7698275) q[1];
sx q[1];
rz(0.1951898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0597417) q[0];
sx q[0];
rz(-1.7072816) q[0];
sx q[0];
rz(-1.6249715) q[0];
x q[1];
rz(-0.11048079) q[2];
sx q[2];
rz(-2.049597) q[2];
sx q[2];
rz(2.2564121) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4976552) q[1];
sx q[1];
rz(-1.8574839) q[1];
sx q[1];
rz(-0.02266702) q[1];
x q[2];
rz(2.3161668) q[3];
sx q[3];
rz(-1.22681) q[3];
sx q[3];
rz(1.0937364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13888415) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(-1.8769334) q[2];
rz(1.5876611) q[3];
sx q[3];
rz(-1.2041661) q[3];
sx q[3];
rz(1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3569262) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(-2.7333562) q[0];
rz(-1.6881855) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(0.85809392) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30735874) q[0];
sx q[0];
rz(-1.8862794) q[0];
sx q[0];
rz(-0.90641109) q[0];
x q[1];
rz(-0.9246773) q[2];
sx q[2];
rz(-0.16077006) q[2];
sx q[2];
rz(-1.9882261) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9436904) q[1];
sx q[1];
rz(-2.5630782) q[1];
sx q[1];
rz(2.3499418) q[1];
rz(-pi) q[2];
rz(3.016641) q[3];
sx q[3];
rz(-0.88192421) q[3];
sx q[3];
rz(-0.75923789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66576362) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(-2.8909454) q[2];
rz(0.88746873) q[3];
sx q[3];
rz(-1.9991425) q[3];
sx q[3];
rz(0.34621507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86824054) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(-0.86531472) q[0];
rz(2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(1.685198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0532074) q[0];
sx q[0];
rz(-1.6376978) q[0];
sx q[0];
rz(1.3752723) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7865042) q[2];
sx q[2];
rz(-0.68195365) q[2];
sx q[2];
rz(-0.29302675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.870112) q[1];
sx q[1];
rz(-2.5462954) q[1];
sx q[1];
rz(0.65956302) q[1];
rz(-pi) q[2];
rz(-2.5771477) q[3];
sx q[3];
rz(-2.5520177) q[3];
sx q[3];
rz(-1.6605476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3978079) q[2];
sx q[2];
rz(-0.91181552) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.926173) q[0];
sx q[0];
rz(-1.8126491) q[0];
sx q[0];
rz(1.5480315) q[0];
rz(2.678395) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(-0.40512591) q[2];
sx q[2];
rz(-2.4562757) q[2];
sx q[2];
rz(1.9669878) q[2];
rz(1.4494827) q[3];
sx q[3];
rz(-2.2400845) q[3];
sx q[3];
rz(-1.6252124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
