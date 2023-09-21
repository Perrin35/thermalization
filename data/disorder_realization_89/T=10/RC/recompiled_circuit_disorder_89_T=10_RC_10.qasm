OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(4.8737704) q[0];
sx q[0];
rz(9.1302172) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298676) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(3.1153326) q[0];
rz(0.047702567) q[2];
sx q[2];
rz(-0.68714266) q[2];
sx q[2];
rz(1.7575761) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4946343) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(-1.3691982) q[1];
x q[2];
rz(-2.5476417) q[3];
sx q[3];
rz(-1.6868601) q[3];
sx q[3];
rz(-0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0191779) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0533957) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(-0.85900599) q[0];
rz(0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.4650311) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61629907) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(0.3152245) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6071762) q[2];
sx q[2];
rz(-0.60306163) q[2];
sx q[2];
rz(0.39819983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4620004) q[1];
sx q[1];
rz(-2.3604217) q[1];
sx q[1];
rz(-2.6189234) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6724085) q[3];
sx q[3];
rz(-0.57933925) q[3];
sx q[3];
rz(-2.7964696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64840245) q[0];
sx q[0];
rz(-1.5957386) q[0];
sx q[0];
rz(-0.52669749) q[0];
x q[1];
rz(-1.221237) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(-1.0457872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.219017) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(-2.4465357) q[1];
rz(-pi) q[2];
rz(0.80188607) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0107161) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(2.1219357) q[2];
rz(-2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-2.5286634) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(-1.5377195) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(0.61027169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71826868) q[0];
sx q[0];
rz(-1.4615131) q[0];
sx q[0];
rz(0.011358326) q[0];
x q[1];
rz(-0.35595603) q[2];
sx q[2];
rz(-1.2095371) q[2];
sx q[2];
rz(1.4796096) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8295633) q[1];
sx q[1];
rz(-2.2514572) q[1];
sx q[1];
rz(-2.8042253) q[1];
rz(-pi) q[2];
rz(-2.652076) q[3];
sx q[3];
rz(-2.094141) q[3];
sx q[3];
rz(-0.19260064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6326555) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(0.73348796) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.3132494) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35996138) q[0];
sx q[0];
rz(-0.93802035) q[0];
sx q[0];
rz(-0.40681337) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0079185) q[2];
sx q[2];
rz(-0.067069947) q[2];
sx q[2];
rz(-1.9474533) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6254127) q[1];
sx q[1];
rz(-1.1986294) q[1];
sx q[1];
rz(0.401293) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89574121) q[3];
sx q[3];
rz(-2.3268019) q[3];
sx q[3];
rz(2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-2.4772947) q[2];
rz(-0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(0.48318133) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7532363) q[0];
sx q[0];
rz(-1.5090794) q[0];
sx q[0];
rz(-1.7524936) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1088164) q[2];
sx q[2];
rz(-1.1083229) q[2];
sx q[2];
rz(-2.0605161) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.253787) q[1];
sx q[1];
rz(-2.0704381) q[1];
sx q[1];
rz(-2.3274724) q[1];
rz(1.2504962) q[3];
sx q[3];
rz(-0.34892504) q[3];
sx q[3];
rz(-2.2826209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(1.9412458) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(2.8894997) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82089059) q[0];
sx q[0];
rz(-1.1517236) q[0];
sx q[0];
rz(1.4120031) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(2.1512254) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7455709) q[1];
sx q[1];
rz(-2.3704297) q[1];
sx q[1];
rz(1.039617) q[1];
rz(0.89594706) q[3];
sx q[3];
rz(-0.5657256) q[3];
sx q[3];
rz(1.56324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(2.2371116) q[2];
rz(2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.3716912) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6721281) q[0];
sx q[0];
rz(-2.1882595) q[0];
sx q[0];
rz(-3.0332964) q[0];
x q[1];
rz(-0.76099446) q[2];
sx q[2];
rz(-1.7241663) q[2];
sx q[2];
rz(2.209687) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9904069) q[1];
sx q[1];
rz(-1.2745665) q[1];
sx q[1];
rz(2.5469261) q[1];
rz(1.4189818) q[3];
sx q[3];
rz(-1.2137128) q[3];
sx q[3];
rz(-2.9861772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(2.7436942) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(2.871002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4280677) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(0.2268976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19240304) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(2.2601688) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8447664) q[1];
sx q[1];
rz(-1.7080194) q[1];
sx q[1];
rz(1.574135) q[1];
x q[2];
rz(1.3769334) q[3];
sx q[3];
rz(-0.48741515) q[3];
sx q[3];
rz(2.2262115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(-1.489893) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.9809013) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437317) q[0];
sx q[0];
rz(-1.568856) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(0.74238366) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091632) q[0];
sx q[0];
rz(-2.2262648) q[0];
sx q[0];
rz(-0.79188939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97787751) q[2];
sx q[2];
rz(-1.2061384) q[2];
sx q[2];
rz(-2.8713771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87634516) q[1];
sx q[1];
rz(-2.7956388) q[1];
sx q[1];
rz(-0.032851263) q[1];
x q[2];
rz(0.19974444) q[3];
sx q[3];
rz(-2.9962857) q[3];
sx q[3];
rz(2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(-0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-2.8200091) q[2];
sx q[2];
rz(-2.2939773) q[2];
sx q[2];
rz(1.1800223) q[2];
rz(-1.7066163) q[3];
sx q[3];
rz(-1.4321696) q[3];
sx q[3];
rz(-1.1925478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
