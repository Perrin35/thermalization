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
rz(5.4592291) q[0];
sx q[0];
rz(8.4973314) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(-2.3271022) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0448582) q[0];
sx q[0];
rz(-1.1219026) q[0];
sx q[0];
rz(-2.7142224) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1998946) q[2];
sx q[2];
rz(-2.2613031) q[2];
sx q[2];
rz(-2.5494573) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42871745) q[1];
sx q[1];
rz(-0.71850417) q[1];
sx q[1];
rz(0.76193049) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85759576) q[3];
sx q[3];
rz(-2.2603545) q[3];
sx q[3];
rz(1.5201207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3689975) q[2];
sx q[2];
rz(-1.3694222) q[2];
sx q[2];
rz(2.7102846) q[2];
rz(-1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-0.78540426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5539219) q[0];
sx q[0];
rz(-0.19667721) q[0];
sx q[0];
rz(1.824463) q[0];
rz(1.0493086) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(-2.5228693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878889) q[0];
sx q[0];
rz(-1.426322) q[0];
sx q[0];
rz(0.61019759) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0929865) q[2];
sx q[2];
rz(-0.75121643) q[2];
sx q[2];
rz(0.71556811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2545027) q[1];
sx q[1];
rz(-1.2225593) q[1];
sx q[1];
rz(-3.00249) q[1];
x q[2];
rz(-0.57024184) q[3];
sx q[3];
rz(-2.5240353) q[3];
sx q[3];
rz(2.9887426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.416136) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(0.49187342) q[2];
rz(1.5791996) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.1270776) q[0];
sx q[0];
rz(-1.0579695) q[0];
sx q[0];
rz(-2.8614817) q[0];
rz(0.07490553) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(1.3704971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89383306) q[0];
sx q[0];
rz(-2.3671208) q[0];
sx q[0];
rz(-1.872137) q[0];
x q[1];
rz(2.2417775) q[2];
sx q[2];
rz(-2.0986989) q[2];
sx q[2];
rz(1.4769276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.83082) q[1];
sx q[1];
rz(-1.3078863) q[1];
sx q[1];
rz(2.5710158) q[1];
x q[2];
rz(1.5519484) q[3];
sx q[3];
rz(-1.4385838) q[3];
sx q[3];
rz(0.16308768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8568153) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(0.80186296) q[2];
rz(2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(-2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
rz(2.3367679) q[0];
sx q[0];
rz(-0.66773361) q[0];
sx q[0];
rz(-2.5208852) q[0];
rz(-2.8630818) q[1];
sx q[1];
rz(-1.4418437) q[1];
sx q[1];
rz(2.1381569) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8580756) q[0];
sx q[0];
rz(-1.2293613) q[0];
sx q[0];
rz(0.33591875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7264257) q[2];
sx q[2];
rz(-0.4253201) q[2];
sx q[2];
rz(2.164054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4662557) q[1];
sx q[1];
rz(-1.3020483) q[1];
sx q[1];
rz(1.6273725) q[1];
rz(-pi) q[2];
rz(-1.8068321) q[3];
sx q[3];
rz(-0.98525648) q[3];
sx q[3];
rz(1.745831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8371381) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(3.0755074) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-0.76032138) q[3];
sx q[3];
rz(0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5758301) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(-2.3542812) q[0];
rz(-0.85834223) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(-1.0816921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4612686) q[0];
sx q[0];
rz(-2.1787365) q[0];
sx q[0];
rz(-2.1896088) q[0];
rz(-pi) q[1];
x q[1];
rz(1.101564) q[2];
sx q[2];
rz(-1.9479556) q[2];
sx q[2];
rz(0.18845972) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9138152) q[1];
sx q[1];
rz(-0.28474131) q[1];
sx q[1];
rz(-3.1193507) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55629913) q[3];
sx q[3];
rz(-2.2101758) q[3];
sx q[3];
rz(-2.104708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(0.84257379) q[2];
rz(0.27635559) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(1.3493376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-1.1899765) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(-2.55866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806429) q[0];
sx q[0];
rz(-1.1034729) q[0];
sx q[0];
rz(1.1790465) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.013986258) q[2];
sx q[2];
rz(-0.99726935) q[2];
sx q[2];
rz(1.397246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5546023) q[1];
sx q[1];
rz(-1.6515367) q[1];
sx q[1];
rz(0.23861532) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63151367) q[3];
sx q[3];
rz(-1.9296292) q[3];
sx q[3];
rz(1.7874908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.8913816) q[2];
sx q[2];
rz(1.0397376) q[2];
rz(0.58935634) q[3];
sx q[3];
rz(-1.6683234) q[3];
sx q[3];
rz(-1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650918) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(-1.804922) q[0];
rz(-0.22557766) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(0.76990661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25882803) q[0];
sx q[0];
rz(-2.6443091) q[0];
sx q[0];
rz(-0.33035853) q[0];
rz(-2.3138397) q[2];
sx q[2];
rz(-1.079664) q[2];
sx q[2];
rz(0.14358768) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0437255) q[1];
sx q[1];
rz(-1.3732404) q[1];
sx q[1];
rz(1.6476589) q[1];
x q[2];
rz(3.0287241) q[3];
sx q[3];
rz(-0.60823554) q[3];
sx q[3];
rz(0.0054159482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(2.1262271) q[2];
rz(-2.9376302) q[3];
sx q[3];
rz(-1.5488307) q[3];
sx q[3];
rz(1.7502194) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400951) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(0.39328662) q[0];
rz(-2.6914864) q[1];
sx q[1];
rz(-1.422912) q[1];
sx q[1];
rz(-0.030287655) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6393042) q[0];
sx q[0];
rz(-0.24950108) q[0];
sx q[0];
rz(-0.12608237) q[0];
rz(2.8856336) q[2];
sx q[2];
rz(-0.82590196) q[2];
sx q[2];
rz(1.2307375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4612164) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(0.044574634) q[1];
rz(-2.6859517) q[3];
sx q[3];
rz(-2.8491088) q[3];
sx q[3];
rz(-2.4258421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0344737) q[2];
sx q[2];
rz(-1.5566885) q[2];
sx q[2];
rz(2.0195473) q[2];
rz(-1.0551039) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(-1.8900227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1114383) q[1];
sx q[1];
rz(-2.074362) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36893435) q[0];
sx q[0];
rz(-2.308929) q[0];
sx q[0];
rz(-1.5940109) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7238317) q[2];
sx q[2];
rz(-2.8308597) q[2];
sx q[2];
rz(2.8513214) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6744819) q[1];
sx q[1];
rz(-1.8201105) q[1];
sx q[1];
rz(1.8915908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9638941) q[3];
sx q[3];
rz(-2.0434994) q[3];
sx q[3];
rz(0.79897579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.60712236) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(-1.2657451) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-0.37467343) q[3];
sx q[3];
rz(-3.1025187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(1.2598502) q[0];
rz(1.6549567) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(-3.1390417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5887506) q[0];
sx q[0];
rz(-1.4511943) q[0];
sx q[0];
rz(1.6970859) q[0];
rz(-1.4223406) q[2];
sx q[2];
rz(-3.0768148) q[2];
sx q[2];
rz(-1.5744408) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8067183) q[1];
sx q[1];
rz(-2.2069227) q[1];
sx q[1];
rz(-2.0680554) q[1];
rz(-pi) q[2];
rz(-2.0052913) q[3];
sx q[3];
rz(-1.5467724) q[3];
sx q[3];
rz(2.0498118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(0.67443887) q[2];
rz(-1.1174348) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-0.31291541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0565223) q[0];
sx q[0];
rz(-2.0105579) q[0];
sx q[0];
rz(2.6149268) q[0];
rz(-2.7936735) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(2.7958776) q[2];
sx q[2];
rz(-1.7713265) q[2];
sx q[2];
rz(-2.3878018) q[2];
rz(-1.4244637) q[3];
sx q[3];
rz(-1.5531333) q[3];
sx q[3];
rz(0.2266758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
