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
rz(-0.076844849) q[0];
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(-2.3208148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73859204) q[0];
sx q[0];
rz(-2.0401017) q[0];
sx q[0];
rz(1.4856536) q[0];
x q[1];
rz(-2.2256741) q[2];
sx q[2];
rz(-1.7015463) q[2];
sx q[2];
rz(0.98721013) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1917853) q[1];
sx q[1];
rz(-2.2151183) q[1];
sx q[1];
rz(-2.5271646) q[1];
rz(-0.54475875) q[3];
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
rz(-0.57089388) q[2];
sx q[2];
rz(-1.6119733) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.8703929) q[3];
sx q[3];
rz(2.1421656) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73794) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(-1.8988443) q[0];
rz(0.87951648) q[1];
sx q[1];
rz(-0.95006919) q[1];
sx q[1];
rz(-1.8125199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14318109) q[0];
sx q[0];
rz(-0.77309009) q[0];
sx q[0];
rz(0.66112681) q[0];
x q[1];
rz(2.3622386) q[2];
sx q[2];
rz(-2.4930304) q[2];
sx q[2];
rz(-0.88403242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0893007) q[1];
sx q[1];
rz(-1.1049905) q[1];
sx q[1];
rz(2.0897081) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69345052) q[3];
sx q[3];
rz(-1.4458522) q[3];
sx q[3];
rz(-3.0378064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8326412) q[2];
sx q[2];
rz(-1.917104) q[2];
sx q[2];
rz(1.857117) q[2];
rz(-0.66713157) q[3];
sx q[3];
rz(-0.37225538) q[3];
sx q[3];
rz(-2.4767806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852378) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(1.941823) q[0];
rz(1.2780227) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(2.3866381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2374868) q[0];
sx q[0];
rz(-1.5747935) q[0];
sx q[0];
rz(2.0275751) q[0];
x q[1];
rz(2.9725482) q[2];
sx q[2];
rz(-2.4158106) q[2];
sx q[2];
rz(-2.5570585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1633411) q[1];
sx q[1];
rz(-2.3929407) q[1];
sx q[1];
rz(-1.6221552) q[1];
x q[2];
rz(-0.69061227) q[3];
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
rz(0.26015002) q[3];
sx q[3];
rz(-1.6157506) q[3];
sx q[3];
rz(-2.9366176) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78885704) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(-2.4810171) q[0];
rz(-0.13064101) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(-0.43636838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0562387) q[0];
sx q[0];
rz(-1.071601) q[0];
sx q[0];
rz(0.27277314) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5533668) q[2];
sx q[2];
rz(-2.4674621) q[2];
sx q[2];
rz(0.16954409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1088686) q[1];
sx q[1];
rz(-2.0744262) q[1];
sx q[1];
rz(-2.7313066) q[1];
rz(-pi) q[2];
rz(1.4722669) q[3];
sx q[3];
rz(-1.9103582) q[3];
sx q[3];
rz(1.3738541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.074162) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(-1.1379498) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.725089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013019) q[0];
sx q[0];
rz(-1.8922946) q[0];
sx q[0];
rz(-2.9345595) q[0];
rz(1.8197458) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(-1.0386937) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90090758) q[0];
sx q[0];
rz(-1.0552013) q[0];
sx q[0];
rz(2.2365966) q[0];
x q[1];
rz(-2.8530099) q[2];
sx q[2];
rz(-2.5009857) q[2];
sx q[2];
rz(-2.688302) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0463883) q[1];
sx q[1];
rz(-0.63648736) q[1];
sx q[1];
rz(3.1134136) q[1];
x q[2];
rz(0.33090277) q[3];
sx q[3];
rz(-2.8185113) q[3];
sx q[3];
rz(-2.4289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1553104) q[2];
sx q[2];
rz(-1.0460514) q[2];
sx q[2];
rz(-2.5347064) q[2];
rz(0.33489975) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(1.5437532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5034921) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(-2.959429) q[0];
rz(-0.69193524) q[1];
sx q[1];
rz(-1.7057995) q[1];
sx q[1];
rz(1.8686132) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48995724) q[0];
sx q[0];
rz(-3.0916442) q[0];
sx q[0];
rz(-1.7836216) q[0];
rz(-0.23006205) q[2];
sx q[2];
rz(-1.4289235) q[2];
sx q[2];
rz(1.0070247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.581785) q[1];
sx q[1];
rz(-0.81830561) q[1];
sx q[1];
rz(0.35840987) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5073413) q[3];
sx q[3];
rz(-1.9279908) q[3];
sx q[3];
rz(-0.82596367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2481897) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(-2.2656061) q[2];
rz(0.28410965) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(-2.7178154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9559602) q[0];
sx q[0];
rz(-2.1946588) q[0];
sx q[0];
rz(2.60485) q[0];
rz(-2.5770889) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(1.5435262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3878582) q[0];
sx q[0];
rz(-1.6711387) q[0];
sx q[0];
rz(-2.8265586) q[0];
rz(-pi) q[1];
rz(2.6993611) q[2];
sx q[2];
rz(-2.6870607) q[2];
sx q[2];
rz(-1.0584128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8720756) q[1];
sx q[1];
rz(-2.4763514) q[1];
sx q[1];
rz(-0.77209446) q[1];
rz(0.95682896) q[3];
sx q[3];
rz(-1.9787297) q[3];
sx q[3];
rz(1.8338501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87968612) q[2];
sx q[2];
rz(-1.8330816) q[2];
sx q[2];
rz(-1.7730664) q[2];
rz(0.3174828) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(0.69664636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.19318652) q[0];
sx q[0];
rz(-0.53006154) q[0];
sx q[0];
rz(-2.1004706) q[0];
rz(2.5913473) q[1];
sx q[1];
rz(-2.7698275) q[1];
sx q[1];
rz(0.1951898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.081851) q[0];
sx q[0];
rz(-1.4343111) q[0];
sx q[0];
rz(-1.6249715) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7800897) q[2];
sx q[2];
rz(-2.6511764) q[2];
sx q[2];
rz(-2.0201266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.079551815) q[1];
sx q[1];
rz(-1.5490546) q[1];
sx q[1];
rz(-1.2840391) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45359277) q[3];
sx q[3];
rz(-2.2633584) q[3];
sx q[3];
rz(-0.778824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0027085) q[2];
sx q[2];
rz(-2.0033074) q[2];
sx q[2];
rz(-1.8769334) q[2];
rz(-1.5539315) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(1.9896265) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(-2.7333562) q[0];
rz(-1.4534072) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(-0.85809392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6397809) q[0];
sx q[0];
rz(-2.1970891) q[0];
sx q[0];
rz(0.39299742) q[0];
rz(-pi) q[1];
rz(1.4420322) q[2];
sx q[2];
rz(-1.4742645) q[2];
sx q[2];
rz(-2.0842722) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0777867) q[1];
sx q[1];
rz(-1.1764472) q[1];
sx q[1];
rz(2.0057682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4205877) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(-0.56415375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.475829) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(0.25064722) q[2];
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
rz(-pi/2) q[3];
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
rz(-0.86824054) q[0];
sx q[0];
rz(-2.1526985) q[0];
sx q[0];
rz(2.2762779) q[0];
rz(2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(1.685198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0532074) q[0];
sx q[0];
rz(-1.6376978) q[0];
sx q[0];
rz(-1.7663203) q[0];
x q[1];
rz(0.9002879) q[2];
sx q[2];
rz(-1.7061201) q[2];
sx q[2];
rz(-1.4462665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2715312) q[1];
sx q[1];
rz(-1.9215596) q[1];
sx q[1];
rz(2.6501772) q[1];
rz(-0.56444498) q[3];
sx q[3];
rz(-0.58957499) q[3];
sx q[3];
rz(1.481045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.3978079) q[2];
sx q[2];
rz(-2.2297771) q[2];
sx q[2];
rz(1.3249116) q[2];
rz(1.2815453) q[3];
sx q[3];
rz(-0.89373389) q[3];
sx q[3];
rz(2.3769784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21541967) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(0.46319766) q[1];
sx q[1];
rz(-1.6013655) q[1];
sx q[1];
rz(-0.2688437) q[1];
rz(2.4972476) q[2];
sx q[2];
rz(-1.8229137) q[2];
sx q[2];
rz(0.075621123) q[2];
rz(0.15179612) q[3];
sx q[3];
rz(-2.463068) q[3];
sx q[3];
rz(1.3223552) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
