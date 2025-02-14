OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9659757) q[0];
sx q[0];
rz(-1.5953925) q[0];
sx q[0];
rz(-1.5104798) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(4.4463867) q[1];
sx q[1];
rz(10.211853) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.118282) q[0];
sx q[0];
rz(-2.2459789) q[0];
sx q[0];
rz(-1.7143296) q[0];
x q[1];
rz(-3.0120968) q[2];
sx q[2];
rz(-0.39948764) q[2];
sx q[2];
rz(-1.1286038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2752645) q[1];
sx q[1];
rz(-1.6073078) q[1];
sx q[1];
rz(2.9300772) q[1];
x q[2];
rz(-2.5236457) q[3];
sx q[3];
rz(-1.071458) q[3];
sx q[3];
rz(-1.7955347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8481019) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(-0.21609406) q[2];
rz(0.50348336) q[3];
sx q[3];
rz(-1.8899906) q[3];
sx q[3];
rz(3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1450495) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-0.37522069) q[0];
rz(2.5253865) q[1];
sx q[1];
rz(-1.0169949) q[1];
sx q[1];
rz(-0.18888758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4746162) q[0];
sx q[0];
rz(-2.24581) q[0];
sx q[0];
rz(0.53860177) q[0];
rz(-pi) q[1];
rz(-2.6134256) q[2];
sx q[2];
rz(-2.2303827) q[2];
sx q[2];
rz(-0.41925493) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5423373) q[1];
sx q[1];
rz(-2.2399678) q[1];
sx q[1];
rz(-2.8821844) q[1];
rz(2.5152605) q[3];
sx q[3];
rz(-0.62343684) q[3];
sx q[3];
rz(1.3960736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6836493) q[2];
sx q[2];
rz(-1.318149) q[2];
sx q[2];
rz(1.1434435) q[2];
rz(0.6360561) q[3];
sx q[3];
rz(-3.0885185) q[3];
sx q[3];
rz(-1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627593) q[0];
sx q[0];
rz(-2.9017359) q[0];
sx q[0];
rz(1.1676316) q[0];
rz(3.0150343) q[1];
sx q[1];
rz(-1.626222) q[1];
sx q[1];
rz(-2.5680241) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180998) q[0];
sx q[0];
rz(-2.3947236) q[0];
sx q[0];
rz(1.7043714) q[0];
rz(-pi) q[1];
rz(-1.838348) q[2];
sx q[2];
rz(-1.9256136) q[2];
sx q[2];
rz(1.4894384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3237649) q[1];
sx q[1];
rz(-2.7080302) q[1];
sx q[1];
rz(1.6852412) q[1];
rz(-2.6010547) q[3];
sx q[3];
rz(-2.4063568) q[3];
sx q[3];
rz(-2.102004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48587376) q[2];
sx q[2];
rz(-1.3081009) q[2];
sx q[2];
rz(1.6386848) q[2];
rz(-1.2949299) q[3];
sx q[3];
rz(-2.86627) q[3];
sx q[3];
rz(0.82726971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5824222) q[0];
sx q[0];
rz(-0.14635135) q[0];
sx q[0];
rz(-0.91648066) q[0];
rz(2.9811663) q[1];
sx q[1];
rz(-1.8564936) q[1];
sx q[1];
rz(-0.697335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6354016) q[0];
sx q[0];
rz(-2.7317193) q[0];
sx q[0];
rz(0.91110595) q[0];
rz(1.2658046) q[2];
sx q[2];
rz(-2.2100497) q[2];
sx q[2];
rz(1.6219307) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1707499) q[1];
sx q[1];
rz(-0.80180556) q[1];
sx q[1];
rz(2.6796653) q[1];
rz(-pi) q[2];
rz(3.011441) q[3];
sx q[3];
rz(-1.5846545) q[3];
sx q[3];
rz(-2.7133301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15630284) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(1.1227597) q[2];
rz(1.5924234) q[3];
sx q[3];
rz(-1.4607818) q[3];
sx q[3];
rz(2.7482225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321946) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(1.1291809) q[0];
rz(0.87942266) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(0.62854356) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.382269) q[0];
sx q[0];
rz(-0.17062561) q[0];
sx q[0];
rz(-2.3294446) q[0];
x q[1];
rz(-2.4714575) q[2];
sx q[2];
rz(-2.3628855) q[2];
sx q[2];
rz(-2.7846732) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0776808) q[1];
sx q[1];
rz(-1.6110629) q[1];
sx q[1];
rz(-1.0185373) q[1];
rz(2.3102949) q[3];
sx q[3];
rz(-1.4834838) q[3];
sx q[3];
rz(-1.0333453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2200372) q[2];
sx q[2];
rz(-0.76801378) q[2];
sx q[2];
rz(1.3735695) q[2];
rz(-0.31275648) q[3];
sx q[3];
rz(-1.0765358) q[3];
sx q[3];
rz(-3.055174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81109989) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(-2.9789441) q[0];
rz(1.2341518) q[1];
sx q[1];
rz(-1.1470497) q[1];
sx q[1];
rz(-0.92076921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61495691) q[0];
sx q[0];
rz(-1.4702073) q[0];
sx q[0];
rz(1.7593764) q[0];
rz(-0.8779432) q[2];
sx q[2];
rz(-0.77550807) q[2];
sx q[2];
rz(0.42693769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15359503) q[1];
sx q[1];
rz(-2.3210166) q[1];
sx q[1];
rz(-1.4819891) q[1];
x q[2];
rz(0.7821857) q[3];
sx q[3];
rz(-1.4504004) q[3];
sx q[3];
rz(2.6162755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2938701) q[2];
sx q[2];
rz(-2.4203478) q[2];
sx q[2];
rz(3.0976307) q[2];
rz(1.4219159) q[3];
sx q[3];
rz(-2.7094262) q[3];
sx q[3];
rz(-0.034984263) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1304355) q[0];
sx q[0];
rz(-0.80428094) q[0];
sx q[0];
rz(-0.093611896) q[0];
rz(1.0253819) q[1];
sx q[1];
rz(-1.5461812) q[1];
sx q[1];
rz(0.78790087) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073613361) q[0];
sx q[0];
rz(-2.1755621) q[0];
sx q[0];
rz(0.66901274) q[0];
x q[1];
rz(0.56886832) q[2];
sx q[2];
rz(-1.427257) q[2];
sx q[2];
rz(1.6678866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3026784) q[1];
sx q[1];
rz(-2.5414349) q[1];
sx q[1];
rz(-3.0693574) q[1];
x q[2];
rz(-0.8696733) q[3];
sx q[3];
rz(-0.18571412) q[3];
sx q[3];
rz(2.7839422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0755783) q[2];
sx q[2];
rz(-2.5983073) q[2];
sx q[2];
rz(0.96599609) q[2];
rz(0.14723369) q[3];
sx q[3];
rz(-1.6817663) q[3];
sx q[3];
rz(2.537263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886993) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(-1.8889486) q[0];
rz(0.057627536) q[1];
sx q[1];
rz(-1.7689972) q[1];
sx q[1];
rz(2.7431814) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3301901) q[0];
sx q[0];
rz(-0.63985642) q[0];
sx q[0];
rz(0.29307239) q[0];
x q[1];
rz(-0.49968991) q[2];
sx q[2];
rz(-1.4213287) q[2];
sx q[2];
rz(0.88525822) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1507453) q[1];
sx q[1];
rz(-1.7707033) q[1];
sx q[1];
rz(-3.1374817) q[1];
rz(-1.6774954) q[3];
sx q[3];
rz(-2.1114285) q[3];
sx q[3];
rz(1.4832527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56732279) q[2];
sx q[2];
rz(-1.8696573) q[2];
sx q[2];
rz(-0.92180139) q[2];
rz(2.4343893) q[3];
sx q[3];
rz(-2.0965529) q[3];
sx q[3];
rz(-0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(-2.2344672) q[0];
rz(-0.47997296) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(2.5352246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7146695) q[0];
sx q[0];
rz(-1.6712609) q[0];
sx q[0];
rz(2.8313374) q[0];
rz(-pi) q[1];
rz(1.0125514) q[2];
sx q[2];
rz(-1.1079567) q[2];
sx q[2];
rz(-1.3430581) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7545112) q[1];
sx q[1];
rz(-0.84085995) q[1];
sx q[1];
rz(2.74387) q[1];
rz(-pi) q[2];
rz(0.14056345) q[3];
sx q[3];
rz(-2.478001) q[3];
sx q[3];
rz(0.47534761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85764641) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(0.66961163) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-1.0457057) q[3];
sx q[3];
rz(1.4001575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46212101) q[0];
sx q[0];
rz(-1.7993878) q[0];
sx q[0];
rz(-0.8959499) q[0];
rz(-0.040291928) q[1];
sx q[1];
rz(-0.94172421) q[1];
sx q[1];
rz(1.5641854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5524254) q[0];
sx q[0];
rz(-1.5755113) q[0];
sx q[0];
rz(-1.5472359) q[0];
rz(-pi) q[1];
rz(-0.68564685) q[2];
sx q[2];
rz(-1.2276936) q[2];
sx q[2];
rz(-0.94742638) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7388707) q[1];
sx q[1];
rz(-2.4634857) q[1];
sx q[1];
rz(0.88749591) q[1];
rz(-pi) q[2];
rz(-0.96503105) q[3];
sx q[3];
rz(-1.26767) q[3];
sx q[3];
rz(0.81699634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3080052) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(-1.8316815) q[2];
rz(-1.9561907) q[3];
sx q[3];
rz(-0.88806051) q[3];
sx q[3];
rz(1.6185224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0014521) q[0];
sx q[0];
rz(-1.8403213) q[0];
sx q[0];
rz(2.1885827) q[0];
rz(3.0630655) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-1.8917812) q[2];
sx q[2];
rz(-2.4413296) q[2];
sx q[2];
rz(1.6632947) q[2];
rz(-0.20542145) q[3];
sx q[3];
rz(-1.4170437) q[3];
sx q[3];
rz(0.99544256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
