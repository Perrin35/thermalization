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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92511237) q[0];
sx q[0];
rz(-0.47639936) q[0];
sx q[0];
rz(2.9754378) q[0];
x q[1];
rz(-2.9772867) q[2];
sx q[2];
rz(-0.92245692) q[2];
sx q[2];
rz(0.68337464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3617435) q[1];
sx q[1];
rz(-1.0917772) q[1];
sx q[1];
rz(-0.82734014) q[1];
rz(-pi) q[2];
rz(-1.7480811) q[3];
sx q[3];
rz(-2.1085792) q[3];
sx q[3];
rz(3.1238525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1138175) q[2];
sx q[2];
rz(-0.57089388) q[2];
sx q[2];
rz(1.6119733) q[2];
rz(2.2714603) q[3];
sx q[3];
rz(-1.8703929) q[3];
sx q[3];
rz(2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73794) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(1.8988443) q[0];
rz(-0.87951648) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(1.3290728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97016818) q[0];
sx q[0];
rz(-2.1546083) q[0];
sx q[0];
rz(1.0310571) q[0];
rz(-pi) q[1];
rz(1.0813367) q[2];
sx q[2];
rz(-1.1266409) q[2];
sx q[2];
rz(1.3655438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0893007) q[1];
sx q[1];
rz(-2.0366021) q[1];
sx q[1];
rz(-1.0518846) q[1];
x q[2];
rz(-2.4481421) q[3];
sx q[3];
rz(-1.6957404) q[3];
sx q[3];
rz(-0.10378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8326412) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(1.857117) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2892147) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(-1.1997696) q[0];
rz(-1.2780227) q[1];
sx q[1];
rz(-0.27690241) q[1];
sx q[1];
rz(2.3866381) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65855712) q[0];
sx q[0];
rz(-0.45679507) q[0];
sx q[0];
rz(-1.5798588) q[0];
rz(0.71866106) q[2];
sx q[2];
rz(-1.6826944) q[2];
sx q[2];
rz(1.1132357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2333922) q[1];
sx q[1];
rz(-0.82336872) q[1];
sx q[1];
rz(3.0939332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4944896) q[3];
sx q[3];
rz(-0.80296338) q[3];
sx q[3];
rz(0.76319198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42489147) q[2];
sx q[2];
rz(-2.5548013) q[2];
sx q[2];
rz(1.1699886) q[2];
rz(2.8814426) q[3];
sx q[3];
rz(-1.6157506) q[3];
sx q[3];
rz(-0.20497504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3527356) q[0];
sx q[0];
rz(-1.641888) q[0];
sx q[0];
rz(2.4810171) q[0];
rz(-3.0109516) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(0.43636838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789278) q[0];
sx q[0];
rz(-1.8095865) q[0];
sx q[0];
rz(-1.055607) q[0];
x q[1];
rz(1.1534833) q[2];
sx q[2];
rz(-2.1168323) q[2];
sx q[2];
rz(-2.2653836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37533942) q[1];
sx q[1];
rz(-2.5033567) q[1];
sx q[1];
rz(2.1973647) q[1];
x q[2];
rz(-1.6693258) q[3];
sx q[3];
rz(-1.9103582) q[3];
sx q[3];
rz(1.3738541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.074162) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(1.1379498) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.7646101) q[3];
sx q[3];
rz(1.725089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-1.6013019) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(0.20703319) q[0];
rz(1.3218468) q[1];
sx q[1];
rz(-1.9925947) q[1];
sx q[1];
rz(-1.0386937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90090758) q[0];
sx q[0];
rz(-1.0552013) q[0];
sx q[0];
rz(2.2365966) q[0];
rz(-pi) q[1];
rz(-1.3617351) q[2];
sx q[2];
rz(-0.96067498) q[2];
sx q[2];
rz(0.80792141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0952044) q[1];
sx q[1];
rz(-0.63648736) q[1];
sx q[1];
rz(0.028179006) q[1];
rz(-0.30666017) q[3];
sx q[3];
rz(-1.6741317) q[3];
sx q[3];
rz(0.54319004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1553104) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(-0.60688621) q[2];
rz(2.8066929) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(-1.5437532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63810054) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(2.959429) q[0];
rz(-0.69193524) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(1.2729794) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48995724) q[0];
sx q[0];
rz(-0.049948461) q[0];
sx q[0];
rz(-1.3579711) q[0];
x q[1];
rz(-1.7164549) q[2];
sx q[2];
rz(-1.7985059) q[2];
sx q[2];
rz(2.6109254) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9019487) q[1];
sx q[1];
rz(-1.8297503) q[1];
sx q[1];
rz(-0.78550291) q[1];
rz(0.16830651) q[3];
sx q[3];
rz(-0.3625488) q[3];
sx q[3];
rz(1.0057288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2481897) q[2];
sx q[2];
rz(-0.71257633) q[2];
sx q[2];
rz(0.87598652) q[2];
rz(2.857483) q[3];
sx q[3];
rz(-2.1967389) q[3];
sx q[3];
rz(0.42377728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1856325) q[0];
sx q[0];
rz(-2.1946588) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(0.56450379) q[1];
sx q[1];
rz(-0.91465488) q[1];
sx q[1];
rz(1.5980665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48115981) q[0];
sx q[0];
rz(-0.33012046) q[0];
sx q[0];
rz(2.8274203) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41588628) q[2];
sx q[2];
rz(-1.759811) q[2];
sx q[2];
rz(-0.11014665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8720756) q[1];
sx q[1];
rz(-0.66524129) q[1];
sx q[1];
rz(2.3694982) q[1];
rz(0.48637511) q[3];
sx q[3];
rz(-2.1280299) q[3];
sx q[3];
rz(3.1319905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87968612) q[2];
sx q[2];
rz(-1.3085111) q[2];
sx q[2];
rz(1.3685263) q[2];
rz(-2.8241099) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(-2.4449463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9484061) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(2.1004706) q[0];
rz(2.5913473) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(-0.1951898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6379162) q[0];
sx q[0];
rz(-1.6244672) q[0];
sx q[0];
rz(-0.1366833) q[0];
x q[1];
rz(0.11048079) q[2];
sx q[2];
rz(-2.049597) q[2];
sx q[2];
rz(0.88518054) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7239388) q[1];
sx q[1];
rz(-0.28755763) q[1];
sx q[1];
rz(1.4940667) q[1];
x q[2];
rz(0.82542586) q[3];
sx q[3];
rz(-1.22681) q[3];
sx q[3];
rz(2.0478562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0027085) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(1.2646593) q[2];
rz(1.5539315) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(-1.9896265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(-1.6881855) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(-2.2834987) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88622296) q[0];
sx q[0];
rz(-0.72508922) q[0];
sx q[0];
rz(-1.0839456) q[0];
rz(-pi) q[1];
rz(3.0442601) q[2];
sx q[2];
rz(-1.6989577) q[2];
sx q[2];
rz(-0.50099696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1979023) q[1];
sx q[1];
rz(-2.5630782) q[1];
sx q[1];
rz(2.3499418) q[1];
x q[2];
rz(-1.4205877) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(2.5774389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.475829) q[2];
sx q[2];
rz(-2.2593311) q[2];
sx q[2];
rz(-0.25064722) q[2];
rz(0.88746873) q[3];
sx q[3];
rz(-1.9991425) q[3];
sx q[3];
rz(0.34621507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.2733521) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(0.86531472) q[0];
rz(-0.52262753) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(-1.4563947) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0532074) q[0];
sx q[0];
rz(-1.6376978) q[0];
sx q[0];
rz(-1.3752723) q[0];
rz(-pi) q[1];
rz(-1.7865042) q[2];
sx q[2];
rz(-2.459639) q[2];
sx q[2];
rz(-2.8485659) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48113808) q[1];
sx q[1];
rz(-2.0298971) q[1];
sx q[1];
rz(1.1774239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2271475) q[3];
sx q[3];
rz(-2.0598186) q[3];
sx q[3];
rz(1.0095613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3978079) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(1.816681) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(-0.76461422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.926173) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(-0.46319766) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(1.2590969) q[2];
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
