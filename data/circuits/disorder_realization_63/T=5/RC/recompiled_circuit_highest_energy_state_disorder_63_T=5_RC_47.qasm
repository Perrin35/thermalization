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
rz(2.0923848) q[0];
sx q[0];
rz(-0.75924358) q[0];
sx q[0];
rz(2.7651751) q[0];
rz(1.5578101) q[1];
sx q[1];
rz(-1.0695142) q[1];
sx q[1];
rz(0.95626107) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.555684) q[0];
sx q[0];
rz(-0.95147419) q[0];
sx q[0];
rz(-2.7509718) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7189288) q[2];
sx q[2];
rz(-1.3939121) q[2];
sx q[2];
rz(1.6852578) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9701582) q[1];
sx q[1];
rz(-1.6925319) q[1];
sx q[1];
rz(2.3129169) q[1];
rz(-2.8703635) q[3];
sx q[3];
rz(-1.800897) q[3];
sx q[3];
rz(-2.7903737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3236397) q[2];
sx q[2];
rz(-0.91922593) q[2];
sx q[2];
rz(1.9383355) q[2];
rz(1.0407) q[3];
sx q[3];
rz(-0.87472707) q[3];
sx q[3];
rz(-3.021595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284718) q[0];
sx q[0];
rz(-0.57336837) q[0];
sx q[0];
rz(2.0290802) q[0];
rz(-1.4828697) q[1];
sx q[1];
rz(-0.590938) q[1];
sx q[1];
rz(-2.9842751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1293662) q[0];
sx q[0];
rz(-2.2528045) q[0];
sx q[0];
rz(1.7971044) q[0];
rz(-pi) q[1];
rz(2.7685086) q[2];
sx q[2];
rz(-1.6519128) q[2];
sx q[2];
rz(-0.099992601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78029481) q[1];
sx q[1];
rz(-2.2289702) q[1];
sx q[1];
rz(2.8315413) q[1];
x q[2];
rz(2.9187536) q[3];
sx q[3];
rz(-1.6442219) q[3];
sx q[3];
rz(-1.8334695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6436254) q[2];
sx q[2];
rz(-0.47933856) q[2];
sx q[2];
rz(-0.41110006) q[2];
rz(2.5655365) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(-2.1435553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.75329798) q[0];
sx q[0];
rz(-2.1385215) q[0];
sx q[0];
rz(2.4004747) q[0];
rz(-1.6916212) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(2.5695739) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8328089) q[0];
sx q[0];
rz(-0.35825142) q[0];
sx q[0];
rz(-2.1705296) q[0];
rz(-pi) q[1];
rz(1.3640312) q[2];
sx q[2];
rz(-1.8305147) q[2];
sx q[2];
rz(-2.5842359) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1831844) q[1];
sx q[1];
rz(-2.5377877) q[1];
sx q[1];
rz(0.81012695) q[1];
rz(2.9111262) q[3];
sx q[3];
rz(-1.2433882) q[3];
sx q[3];
rz(-0.62124204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5570306) q[2];
sx q[2];
rz(-1.1921554) q[2];
sx q[2];
rz(0.7977879) q[2];
rz(1.8320463) q[3];
sx q[3];
rz(-1.8833501) q[3];
sx q[3];
rz(1.4057188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13084594) q[0];
sx q[0];
rz(-2.2444785) q[0];
sx q[0];
rz(2.77453) q[0];
rz(-1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(2.9201115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9629134) q[0];
sx q[0];
rz(-1.3604465) q[0];
sx q[0];
rz(-1.6400997) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7738557) q[2];
sx q[2];
rz(-1.443111) q[2];
sx q[2];
rz(0.62386306) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96326107) q[1];
sx q[1];
rz(-1.8823922) q[1];
sx q[1];
rz(2.7072886) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7189923) q[3];
sx q[3];
rz(-0.56767332) q[3];
sx q[3];
rz(1.6740396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9941142) q[2];
sx q[2];
rz(-1.1722112) q[2];
sx q[2];
rz(1.8032903) q[2];
rz(-0.5021247) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(2.8978735) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5153656) q[0];
sx q[0];
rz(-2.5318662) q[0];
sx q[0];
rz(1.5355661) q[0];
rz(-3.0961127) q[1];
sx q[1];
rz(-1.3810424) q[1];
sx q[1];
rz(0.46612003) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2855466) q[0];
sx q[0];
rz(-1.4338126) q[0];
sx q[0];
rz(-1.3187506) q[0];
rz(2.8117198) q[2];
sx q[2];
rz(-2.1876642) q[2];
sx q[2];
rz(-0.62451476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.22039761) q[1];
sx q[1];
rz(-0.092893727) q[1];
sx q[1];
rz(0.18735207) q[1];
rz(0.97753559) q[3];
sx q[3];
rz(-2.1338507) q[3];
sx q[3];
rz(-0.42773358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8972299) q[2];
sx q[2];
rz(-1.4655317) q[2];
sx q[2];
rz(-2.3587295) q[2];
rz(-0.016294567) q[3];
sx q[3];
rz(-1.8330845) q[3];
sx q[3];
rz(-2.079594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29430729) q[0];
sx q[0];
rz(-0.48536244) q[0];
sx q[0];
rz(2.7622727) q[0];
rz(-0.30917057) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(-0.22383037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1334134) q[0];
sx q[0];
rz(-2.2531765) q[0];
sx q[0];
rz(-0.40485415) q[0];
x q[1];
rz(1.0439453) q[2];
sx q[2];
rz(-0.2204403) q[2];
sx q[2];
rz(1.8048087) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1081738) q[1];
sx q[1];
rz(-1.464132) q[1];
sx q[1];
rz(-1.4640385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70586127) q[3];
sx q[3];
rz(-2.0076723) q[3];
sx q[3];
rz(2.049751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0763187) q[2];
sx q[2];
rz(-1.9741917) q[2];
sx q[2];
rz(2.8889612) q[2];
rz(-2.1624055) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(2.7948936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.613649) q[0];
sx q[0];
rz(-2.376463) q[0];
sx q[0];
rz(-1.8494404) q[0];
rz(2.6076803) q[1];
sx q[1];
rz(-1.6894059) q[1];
sx q[1];
rz(0.36010489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394516) q[0];
sx q[0];
rz(-1.2107673) q[0];
sx q[0];
rz(0.074010297) q[0];
x q[1];
rz(2.656237) q[2];
sx q[2];
rz(-1.2872641) q[2];
sx q[2];
rz(-0.6142217) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8812027) q[1];
sx q[1];
rz(-2.5298928) q[1];
sx q[1];
rz(-0.42914294) q[1];
rz(-pi) q[2];
rz(-1.4680181) q[3];
sx q[3];
rz(-1.53781) q[3];
sx q[3];
rz(1.7793293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99122938) q[2];
sx q[2];
rz(-1.6232619) q[2];
sx q[2];
rz(0.72597996) q[2];
rz(-2.7109072) q[3];
sx q[3];
rz(-2.3613598) q[3];
sx q[3];
rz(0.45346692) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3579269) q[0];
sx q[0];
rz(-1.4877321) q[0];
sx q[0];
rz(0.7861535) q[0];
rz(0.17732009) q[1];
sx q[1];
rz(-1.0847849) q[1];
sx q[1];
rz(-1.0160944) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1209512) q[0];
sx q[0];
rz(-0.30123152) q[0];
sx q[0];
rz(-1.2094638) q[0];
x q[1];
rz(-0.54722007) q[2];
sx q[2];
rz(-2.8946218) q[2];
sx q[2];
rz(-0.17683593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7238771) q[1];
sx q[1];
rz(-2.2347365) q[1];
sx q[1];
rz(1.592403) q[1];
rz(-2.3887064) q[3];
sx q[3];
rz(-1.2883923) q[3];
sx q[3];
rz(-1.8771859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5440386) q[2];
sx q[2];
rz(-2.1860213) q[2];
sx q[2];
rz(-2.3717144) q[2];
rz(-2.9652413) q[3];
sx q[3];
rz(-1.4498815) q[3];
sx q[3];
rz(-0.31340733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4835994) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(-1.0461079) q[0];
rz(-1.5931891) q[1];
sx q[1];
rz(-1.4175339) q[1];
sx q[1];
rz(-1.9200602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9402091) q[0];
sx q[0];
rz(-2.7701158) q[0];
sx q[0];
rz(1.8626271) q[0];
x q[1];
rz(-1.8893355) q[2];
sx q[2];
rz(-1.7143068) q[2];
sx q[2];
rz(-0.39538639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0455148) q[1];
sx q[1];
rz(-1.7723262) q[1];
sx q[1];
rz(-1.3526826) q[1];
rz(-pi) q[2];
rz(-1.5658978) q[3];
sx q[3];
rz(-1.5272445) q[3];
sx q[3];
rz(-2.5033308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1239502) q[2];
sx q[2];
rz(-1.4175748) q[2];
sx q[2];
rz(-1.9564015) q[2];
rz(-2.8017398) q[3];
sx q[3];
rz(-0.9797107) q[3];
sx q[3];
rz(-0.11782304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3621984) q[0];
sx q[0];
rz(-1.3267936) q[0];
sx q[0];
rz(2.4438044) q[0];
rz(0.70480529) q[1];
sx q[1];
rz(-2.0185399) q[1];
sx q[1];
rz(2.8885081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9390209) q[0];
sx q[0];
rz(-2.8968985) q[0];
sx q[0];
rz(-1.8807202) q[0];
rz(-pi) q[1];
rz(-1.0519876) q[2];
sx q[2];
rz(-0.98271433) q[2];
sx q[2];
rz(1.9668417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53348161) q[1];
sx q[1];
rz(-2.2072756) q[1];
sx q[1];
rz(2.5627506) q[1];
rz(-pi) q[2];
rz(2.4109957) q[3];
sx q[3];
rz(-1.7640778) q[3];
sx q[3];
rz(2.766942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42085984) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(2.4534658) q[2];
rz(-2.1963035) q[3];
sx q[3];
rz(-1.1752081) q[3];
sx q[3];
rz(2.2139886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5505493) q[0];
sx q[0];
rz(-1.3591546) q[0];
sx q[0];
rz(-1.2426283) q[0];
rz(0.91724829) q[1];
sx q[1];
rz(-1.4731673) q[1];
sx q[1];
rz(0.57631667) q[1];
rz(-2.6131931) q[2];
sx q[2];
rz(-2.4965299) q[2];
sx q[2];
rz(-1.0685558) q[2];
rz(0.096002738) q[3];
sx q[3];
rz(-1.9954372) q[3];
sx q[3];
rz(0.33653997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
