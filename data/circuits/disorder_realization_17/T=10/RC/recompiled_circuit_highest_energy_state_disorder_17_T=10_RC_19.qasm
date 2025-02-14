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
rz(1.849527) q[0];
sx q[0];
rz(-0.82395616) q[0];
sx q[0];
rz(-0.92744654) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(-2.3271022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789291) q[0];
sx q[0];
rz(-1.9534847) q[0];
sx q[0];
rz(2.0576059) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1998946) q[2];
sx q[2];
rz(-0.88028958) q[2];
sx q[2];
rz(-0.59213537) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7128752) q[1];
sx q[1];
rz(-2.4230885) q[1];
sx q[1];
rz(-2.3796622) q[1];
x q[2];
rz(-0.67072224) q[3];
sx q[3];
rz(-2.1938087) q[3];
sx q[3];
rz(-0.58341208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3689975) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(-2.7102846) q[2];
rz(1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5876708) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(-1.824463) q[0];
rz(2.092284) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(0.61872331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4275682) q[0];
sx q[0];
rz(-0.62494266) q[0];
sx q[0];
rz(-0.24863909) q[0];
rz(-2.2631386) q[2];
sx q[2];
rz(-1.8900423) q[2];
sx q[2];
rz(0.4934267) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8730388) q[1];
sx q[1];
rz(-1.7015) q[1];
sx q[1];
rz(1.9221582) q[1];
rz(-0.5388362) q[3];
sx q[3];
rz(-1.2528786) q[3];
sx q[3];
rz(2.2054118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72545663) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(-0.49187342) q[2];
rz(1.5623931) q[3];
sx q[3];
rz(-1.648936) q[3];
sx q[3];
rz(0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014515011) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(0.28011093) q[0];
rz(-3.0666871) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(1.3704971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2477596) q[0];
sx q[0];
rz(-0.77447184) q[0];
sx q[0];
rz(-1.2694556) q[0];
x q[1];
rz(2.3241256) q[2];
sx q[2];
rz(-2.3139179) q[2];
sx q[2];
rz(-0.65929123) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0947511) q[1];
sx q[1];
rz(-1.0221204) q[1];
sx q[1];
rz(1.8803174) q[1];
rz(-pi) q[2];
rz(-0.13223572) q[3];
sx q[3];
rz(-1.5894798) q[3];
sx q[3];
rz(1.736369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8568153) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(-0.80186296) q[2];
rz(2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(0.97889939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(-2.5208852) q[0];
rz(2.8630818) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(-2.1381569) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706948) q[0];
sx q[0];
rz(-1.8866294) q[0];
sx q[0];
rz(1.210808) q[0];
rz(-1.9915646) q[2];
sx q[2];
rz(-1.6347957) q[2];
sx q[2];
rz(0.45129946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6764073) q[1];
sx q[1];
rz(-0.27449755) q[1];
sx q[1];
rz(-0.20250116) q[1];
rz(-pi) q[2];
rz(1.3347606) q[3];
sx q[3];
rz(-2.1563362) q[3];
sx q[3];
rz(-1.745831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3044546) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(-1.1262013) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-3.059721) q[0];
sx q[0];
rz(0.78731147) q[0];
rz(-2.2832504) q[1];
sx q[1];
rz(-0.97802496) q[1];
sx q[1];
rz(-1.0816921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7855378) q[0];
sx q[0];
rz(-2.303146) q[0];
sx q[0];
rz(0.69489702) q[0];
rz(2.2901847) q[2];
sx q[2];
rz(-2.5485196) q[2];
sx q[2];
rz(-2.3874758) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9138152) q[1];
sx q[1];
rz(-0.28474131) q[1];
sx q[1];
rz(0.022241994) q[1];
rz(-pi) q[2];
rz(-0.95329625) q[3];
sx q[3];
rz(-0.82103182) q[3];
sx q[3];
rz(-1.2987069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16296884) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(0.84257379) q[2];
rz(-2.8652371) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.855298) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(1.9516161) q[0];
rz(-1.2225993) q[1];
sx q[1];
rz(-1.2346376) q[1];
sx q[1];
rz(-2.55866) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.615743) q[0];
sx q[0];
rz(-1.2229563) q[0];
sx q[0];
rz(0.49973947) q[0];
x q[1];
rz(-1.5924443) q[2];
sx q[2];
rz(-0.57367838) q[2];
sx q[2];
rz(-1.4230185) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1057824) q[1];
sx q[1];
rz(-1.3329734) q[1];
sx q[1];
rz(1.4877122) q[1];
rz(-pi) q[2];
rz(0.56598466) q[3];
sx q[3];
rz(-2.4275587) q[3];
sx q[3];
rz(-0.66431724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0827834) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(1.0397376) q[2];
rz(2.5522363) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(2.1350433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37650087) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(-1.804922) q[0];
rz(-2.916015) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(2.371686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8827646) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(0.33035853) q[0];
rz(-pi) q[1];
rz(-2.3138397) q[2];
sx q[2];
rz(-1.079664) q[2];
sx q[2];
rz(0.14358768) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6837782) q[1];
sx q[1];
rz(-1.646161) q[1];
sx q[1];
rz(-0.1981258) q[1];
x q[2];
rz(1.6490593) q[3];
sx q[3];
rz(-2.1746082) q[3];
sx q[3];
rz(-0.1426689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(-2.1262271) q[2];
rz(-0.20396248) q[3];
sx q[3];
rz(-1.5488307) q[3];
sx q[3];
rz(1.3913733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400951) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(-2.748306) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(-3.111305) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6323551) q[0];
sx q[0];
rz(-1.3233174) q[0];
sx q[0];
rz(-1.6028274) q[0];
rz(-pi) q[1];
x q[1];
rz(1.302839) q[2];
sx q[2];
rz(-0.7795802) q[2];
sx q[2];
rz(-1.5991581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4612164) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(3.097018) q[1];
rz(-2.6859517) q[3];
sx q[3];
rz(-0.29248387) q[3];
sx q[3];
rz(-0.71575056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0344737) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(1.1220453) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17290641) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(0.77447844) q[0];
rz(-2.5075746) q[1];
sx q[1];
rz(-1.0301544) q[1];
sx q[1];
rz(-1.0672306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7726583) q[0];
sx q[0];
rz(-2.308929) q[0];
sx q[0];
rz(-1.5475818) q[0];
rz(-1.2634694) q[2];
sx q[2];
rz(-1.5241703) q[2];
sx q[2];
rz(1.7152552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53505234) q[1];
sx q[1];
rz(-0.4036223) q[1];
sx q[1];
rz(-0.89151793) q[1];
x q[2];
rz(-1.2379856) q[3];
sx q[3];
rz(-2.6389671) q[3];
sx q[3];
rz(0.42325936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5344703) q[2];
sx q[2];
rz(-3.0202713) q[2];
sx q[2];
rz(1.8758476) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6604615) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.8817425) q[0];
rz(-1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(3.1390417) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5887506) q[0];
sx q[0];
rz(-1.6903983) q[0];
sx q[0];
rz(-1.4445067) q[0];
rz(-pi) q[1];
rz(1.506729) q[2];
sx q[2];
rz(-1.5612215) q[2];
sx q[2];
rz(0.14450444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59460179) q[1];
sx q[1];
rz(-0.78557116) q[1];
sx q[1];
rz(-2.5681096) q[1];
rz(2.0052913) q[3];
sx q[3];
rz(-1.5948203) q[3];
sx q[3];
rz(2.0498118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(2.0241578) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085070327) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(2.7936735) q[1];
sx q[1];
rz(-1.2384474) q[1];
sx q[1];
rz(1.7927982) q[1];
rz(1.7835708) q[2];
sx q[2];
rz(-1.909303) q[2];
sx q[2];
rz(2.252966) q[2];
rz(1.4244637) q[3];
sx q[3];
rz(-1.5884593) q[3];
sx q[3];
rz(-2.9149169) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
