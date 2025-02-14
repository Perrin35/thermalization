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
rz(-1.090156) q[0];
sx q[0];
rz(-0.12026726) q[0];
sx q[0];
rz(3.0647478) q[0];
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(-2.3208148) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92511237) q[0];
sx q[0];
rz(-2.6651933) q[0];
sx q[0];
rz(0.16615482) q[0];
x q[1];
rz(1.7834383) q[2];
sx q[2];
rz(-0.66591381) q[2];
sx q[2];
rz(2.7261811) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8153238) q[1];
sx q[1];
rz(-0.85888777) q[1];
sx q[1];
rz(2.2253042) q[1];
rz(0.54475875) q[3];
sx q[3];
rz(-1.7228456) q[3];
sx q[3];
rz(1.461538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0277752) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(-1.5296193) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.8703929) q[3];
sx q[3];
rz(-0.99942708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40365264) q[0];
sx q[0];
rz(-2.0495179) q[0];
sx q[0];
rz(-1.2427484) q[0];
rz(-0.87951648) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(-1.8125199) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2220228) q[0];
sx q[0];
rz(-1.1276414) q[0];
sx q[0];
rz(2.485347) q[0];
rz(-pi) q[1];
rz(-2.3622386) q[2];
sx q[2];
rz(-0.64856224) q[2];
sx q[2];
rz(2.2575602) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3719888) q[1];
sx q[1];
rz(-1.1117443) q[1];
sx q[1];
rz(0.52476669) q[1];
rz(-pi) q[2];
rz(-1.7326844) q[3];
sx q[3];
rz(-0.88380764) q[3];
sx q[3];
rz(1.5713728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30895147) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(1.857117) q[2];
rz(2.4744611) q[3];
sx q[3];
rz(-0.37225538) q[3];
sx q[3];
rz(-2.4767806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852378) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(-1.941823) q[0];
rz(-1.2780227) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(-2.3866381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4729378) q[0];
sx q[0];
rz(-1.1140214) q[0];
sx q[0];
rz(0.0044538023) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16904449) q[2];
sx q[2];
rz(-2.4158106) q[2];
sx q[2];
rz(2.5570585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1633411) q[1];
sx q[1];
rz(-0.74865197) q[1];
sx q[1];
rz(1.6221552) q[1];
rz(0.69061227) q[3];
sx q[3];
rz(-2.0194144) q[3];
sx q[3];
rz(1.2909362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42489147) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(-1.1699886) q[2];
rz(0.26015002) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(-0.20497504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3527356) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(-0.66057551) q[0];
rz(-3.0109516) q[1];
sx q[1];
rz(-0.6159997) q[1];
sx q[1];
rz(2.7052243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35231464) q[0];
sx q[0];
rz(-1.8095865) q[0];
sx q[0];
rz(-2.0859857) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9881093) q[2];
sx q[2];
rz(-2.1168323) q[2];
sx q[2];
rz(-2.2653836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.032724) q[1];
sx q[1];
rz(-1.0671664) q[1];
sx q[1];
rz(2.7313066) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27159528) q[3];
sx q[3];
rz(-0.35303548) q[3];
sx q[3];
rz(-2.0562382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.074162) q[2];
sx q[2];
rz(-2.6219411) q[2];
sx q[2];
rz(-1.1379498) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.7646101) q[3];
sx q[3];
rz(1.725089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5402907) q[0];
sx q[0];
rz(-1.8922946) q[0];
sx q[0];
rz(-2.9345595) q[0];
rz(1.8197458) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(-1.0386937) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0393676) q[0];
sx q[0];
rz(-1.0034585) q[0];
sx q[0];
rz(-0.62444425) q[0];
rz(-pi) q[1];
rz(2.8530099) q[2];
sx q[2];
rz(-2.5009857) q[2];
sx q[2];
rz(-0.45329061) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50174471) q[1];
sx q[1];
rz(-1.5540489) q[1];
sx q[1];
rz(-0.63629758) q[1];
rz(-pi) q[2];
rz(-0.33090277) q[3];
sx q[3];
rz(-2.8185113) q[3];
sx q[3];
rz(2.4289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98628226) q[2];
sx q[2];
rz(-1.0460514) q[2];
sx q[2];
rz(2.5347064) q[2];
rz(-0.33489975) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(-1.5437532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63810054) q[0];
sx q[0];
rz(-1.426921) q[0];
sx q[0];
rz(2.959429) q[0];
rz(2.4496574) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(1.2729794) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70304027) q[0];
sx q[0];
rz(-1.5219757) q[0];
sx q[0];
rz(3.131034) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9115306) q[2];
sx q[2];
rz(-1.4289235) q[2];
sx q[2];
rz(-2.134568) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.581785) q[1];
sx q[1];
rz(-0.81830561) q[1];
sx q[1];
rz(0.35840987) q[1];
rz(1.6342513) q[3];
sx q[3];
rz(-1.2136019) q[3];
sx q[3];
rz(-0.82596367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89340297) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(2.2656061) q[2];
rz(-2.857483) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(0.42377728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1856325) q[0];
sx q[0];
rz(-2.1946588) q[0];
sx q[0];
rz(2.60485) q[0];
rz(2.5770889) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(-1.5435262) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604328) q[0];
sx q[0];
rz(-2.8114722) q[0];
sx q[0];
rz(2.8274203) q[0];
rz(-0.41588628) q[2];
sx q[2];
rz(-1.3817817) q[2];
sx q[2];
rz(-0.11014665) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8720756) q[1];
sx q[1];
rz(-2.4763514) q[1];
sx q[1];
rz(2.3694982) q[1];
rz(2.1847637) q[3];
sx q[3];
rz(-1.162863) q[3];
sx q[3];
rz(-1.3077426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-2.4449463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9484061) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(-2.1004706) q[0];
rz(0.55024534) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(-2.9464029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.081851) q[0];
sx q[0];
rz(-1.7072816) q[0];
sx q[0];
rz(-1.5166212) q[0];
x q[1];
rz(1.3615029) q[2];
sx q[2];
rz(-0.49041623) q[2];
sx q[2];
rz(2.0201266) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.079551815) q[1];
sx q[1];
rz(-1.5925381) q[1];
sx q[1];
rz(-1.2840391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82542586) q[3];
sx q[3];
rz(-1.9147827) q[3];
sx q[3];
rz(-1.0937364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13888415) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(1.8769334) q[2];
rz(1.5539315) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(1.4534072) q[1];
sx q[1];
rz(-1.6994349) q[1];
sx q[1];
rz(-0.85809392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6397809) q[0];
sx q[0];
rz(-2.1970891) q[0];
sx q[0];
rz(-0.39299742) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0442601) q[2];
sx q[2];
rz(-1.4426349) q[2];
sx q[2];
rz(-2.6405957) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1979023) q[1];
sx q[1];
rz(-2.5630782) q[1];
sx q[1];
rz(-2.3499418) q[1];
x q[2];
rz(1.7210049) q[3];
sx q[3];
rz(-2.4433063) q[3];
sx q[3];
rz(0.56415375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.475829) q[2];
sx q[2];
rz(-2.2593311) q[2];
sx q[2];
rz(2.8909454) q[2];
rz(2.2541239) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(-2.7953776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733521) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(-0.86531472) q[0];
rz(-2.6189651) q[1];
sx q[1];
rz(-2.3326383) q[1];
sx q[1];
rz(-1.4563947) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6107643) q[0];
sx q[0];
rz(-1.3757154) q[0];
sx q[0];
rz(3.0733956) q[0];
rz(1.3550884) q[2];
sx q[2];
rz(-0.68195365) q[2];
sx q[2];
rz(2.8485659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6604546) q[1];
sx q[1];
rz(-1.1116955) q[1];
sx q[1];
rz(-1.9641687) q[1];
x q[2];
rz(-1.9144451) q[3];
sx q[3];
rz(-1.081774) q[3];
sx q[3];
rz(-2.1320313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7437848) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(-1.816681) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(-0.76461422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21541967) q[0];
sx q[0];
rz(-1.8126491) q[0];
sx q[0];
rz(1.5480315) q[0];
rz(2.678395) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(-0.64434509) q[2];
sx q[2];
rz(-1.8229137) q[2];
sx q[2];
rz(0.075621123) q[2];
rz(0.67288053) q[3];
sx q[3];
rz(-1.665848) q[3];
sx q[3];
rz(3.0116826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
