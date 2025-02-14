OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.175617) q[0];
sx q[0];
rz(4.7369851) q[0];
sx q[0];
rz(10.935258) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(-1.8367986) q[1];
sx q[1];
rz(0.78707492) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6841976) q[0];
sx q[0];
rz(-1.4589063) q[0];
sx q[0];
rz(0.68023139) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12949582) q[2];
sx q[2];
rz(-0.39948764) q[2];
sx q[2];
rz(-1.1286038) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2752645) q[1];
sx q[1];
rz(-1.5342848) q[1];
sx q[1];
rz(-2.9300772) q[1];
rz(2.3863411) q[3];
sx q[3];
rz(-0.77338615) q[3];
sx q[3];
rz(-0.36808792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2934908) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(2.9254986) q[2];
rz(2.6381093) q[3];
sx q[3];
rz(-1.8899906) q[3];
sx q[3];
rz(-3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9965432) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-2.766372) q[0];
rz(2.5253865) q[1];
sx q[1];
rz(-1.0169949) q[1];
sx q[1];
rz(2.9527051) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66697648) q[0];
sx q[0];
rz(-2.24581) q[0];
sx q[0];
rz(2.6029909) q[0];
rz(-0.52816708) q[2];
sx q[2];
rz(-0.91120992) q[2];
sx q[2];
rz(-0.41925493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5423373) q[1];
sx q[1];
rz(-2.2399678) q[1];
sx q[1];
rz(0.25940829) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6140574) q[3];
sx q[3];
rz(-1.9200824) q[3];
sx q[3];
rz(0.35638242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45794332) q[2];
sx q[2];
rz(-1.318149) q[2];
sx q[2];
rz(1.9981492) q[2];
rz(0.6360561) q[3];
sx q[3];
rz(-3.0885185) q[3];
sx q[3];
rz(-1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77883333) q[0];
sx q[0];
rz(-0.23985671) q[0];
sx q[0];
rz(1.973961) q[0];
rz(-3.0150343) q[1];
sx q[1];
rz(-1.626222) q[1];
sx q[1];
rz(-0.57356858) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180998) q[0];
sx q[0];
rz(-0.74686909) q[0];
sx q[0];
rz(-1.7043714) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3032446) q[2];
sx q[2];
rz(-1.2159791) q[2];
sx q[2];
rz(1.4894384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4497609) q[1];
sx q[1];
rz(-2.0013325) q[1];
sx q[1];
rz(3.0887763) q[1];
rz(-pi) q[2];
rz(-2.0063806) q[3];
sx q[3];
rz(-0.95802973) q[3];
sx q[3];
rz(-1.719914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.48587376) q[2];
sx q[2];
rz(-1.3081009) q[2];
sx q[2];
rz(-1.6386848) q[2];
rz(-1.8466628) q[3];
sx q[3];
rz(-0.27532268) q[3];
sx q[3];
rz(-2.3143229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5591705) q[0];
sx q[0];
rz(-2.9952413) q[0];
sx q[0];
rz(-0.91648066) q[0];
rz(-0.16042635) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(-2.4442576) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19577414) q[0];
sx q[0];
rz(-1.8911288) q[0];
sx q[0];
rz(-0.26024241) q[0];
x q[1];
rz(-1.2658046) q[2];
sx q[2];
rz(-2.2100497) q[2];
sx q[2];
rz(-1.6219307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5494255) q[1];
sx q[1];
rz(-0.87199713) q[1];
sx q[1];
rz(1.1392085) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13015161) q[3];
sx q[3];
rz(-1.5846545) q[3];
sx q[3];
rz(-0.42826251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9852898) q[2];
sx q[2];
rz(-0.52729815) q[2];
sx q[2];
rz(-2.018833) q[2];
rz(1.5491693) q[3];
sx q[3];
rz(-1.6808108) q[3];
sx q[3];
rz(-0.39337015) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321946) q[0];
sx q[0];
rz(-0.10646146) q[0];
sx q[0];
rz(-1.1291809) q[0];
rz(-0.87942266) q[1];
sx q[1];
rz(-1.3112473) q[1];
sx q[1];
rz(0.62854356) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.382269) q[0];
sx q[0];
rz(-2.970967) q[0];
sx q[0];
rz(0.81214805) q[0];
rz(-2.4833335) q[2];
sx q[2];
rz(-1.1193917) q[2];
sx q[2];
rz(-2.4413204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5694738) q[1];
sx q[1];
rz(-0.55357305) q[1];
sx q[1];
rz(1.4941503) q[1];
x q[2];
rz(-2.3102949) q[3];
sx q[3];
rz(-1.4834838) q[3];
sx q[3];
rz(-2.1082474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2200372) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(1.3735695) q[2];
rz(-2.8288362) q[3];
sx q[3];
rz(-1.0765358) q[3];
sx q[3];
rz(3.055174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3304928) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(2.9789441) q[0];
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
rz(-1.4402078) q[0];
sx q[0];
rz(-0.21344859) q[0];
sx q[0];
rz(-2.0646854) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92451325) q[2];
sx q[2];
rz(-2.0343896) q[2];
sx q[2];
rz(-1.6789503) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28345767) q[1];
sx q[1];
rz(-2.3871536) q[1];
sx q[1];
rz(0.094875022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.359407) q[3];
sx q[3];
rz(-1.4504004) q[3];
sx q[3];
rz(0.52531717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2938701) q[2];
sx q[2];
rz(-0.72124481) q[2];
sx q[2];
rz(3.0976307) q[2];
rz(1.7196767) q[3];
sx q[3];
rz(-2.7094262) q[3];
sx q[3];
rz(0.034984263) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1304355) q[0];
sx q[0];
rz(-2.3373117) q[0];
sx q[0];
rz(3.0479808) q[0];
rz(-1.0253819) q[1];
sx q[1];
rz(-1.5461812) q[1];
sx q[1];
rz(-0.78790087) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073613361) q[0];
sx q[0];
rz(-0.96603051) q[0];
sx q[0];
rz(-2.4725799) q[0];
rz(-pi) q[1];
rz(2.879456) q[2];
sx q[2];
rz(-0.58475625) q[2];
sx q[2];
rz(2.824397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3026784) q[1];
sx q[1];
rz(-0.60015772) q[1];
sx q[1];
rz(0.072235302) q[1];
rz(-3.0209846) q[3];
sx q[3];
rz(-1.4292307) q[3];
sx q[3];
rz(1.0673351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0660144) q[2];
sx q[2];
rz(-0.54328537) q[2];
sx q[2];
rz(-2.1755966) q[2];
rz(2.994359) q[3];
sx q[3];
rz(-1.6817663) q[3];
sx q[3];
rz(-2.537263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.052893355) q[0];
sx q[0];
rz(-1.7558782) q[0];
sx q[0];
rz(1.2526441) q[0];
rz(0.057627536) q[1];
sx q[1];
rz(-1.3725955) q[1];
sx q[1];
rz(-2.7431814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47810546) q[0];
sx q[0];
rz(-1.7441569) q[0];
sx q[0];
rz(-2.5224713) q[0];
x q[1];
rz(0.49968991) q[2];
sx q[2];
rz(-1.4213287) q[2];
sx q[2];
rz(-0.88525822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7207253) q[1];
sx q[1];
rz(-1.5667672) q[1];
sx q[1];
rz(1.770705) q[1];
rz(-0.1755821) q[3];
sx q[3];
rz(-2.5915603) q[3];
sx q[3];
rz(-1.2780785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5742699) q[2];
sx q[2];
rz(-1.8696573) q[2];
sx q[2];
rz(2.2197913) q[2];
rz(-2.4343893) q[3];
sx q[3];
rz(-2.0965529) q[3];
sx q[3];
rz(0.31064492) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-0.12117584) q[0];
sx q[0];
rz(-0.90712547) q[0];
rz(0.47997296) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(0.60636806) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42692318) q[0];
sx q[0];
rz(-1.6712609) q[0];
sx q[0];
rz(2.8313374) q[0];
x q[1];
rz(-2.6098183) q[2];
sx q[2];
rz(-2.0645879) q[2];
sx q[2];
rz(-0.49966787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94928259) q[1];
sx q[1];
rz(-2.3282166) q[1];
sx q[1];
rz(1.9793012) q[1];
rz(-pi) q[2];
rz(1.6799029) q[3];
sx q[3];
rz(-0.91489313) q[3];
sx q[3];
rz(-2.4885268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2839462) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(-0.66961163) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-2.0958869) q[3];
sx q[3];
rz(1.7414352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46212101) q[0];
sx q[0];
rz(-1.7993878) q[0];
sx q[0];
rz(2.2456428) q[0];
rz(-3.1013007) q[1];
sx q[1];
rz(-2.1998684) q[1];
sx q[1];
rz(1.5641854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1600747) q[0];
sx q[0];
rz(-1.5943564) q[0];
sx q[0];
rz(3.1368763) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0032004) q[2];
sx q[2];
rz(-2.2096425) q[2];
sx q[2];
rz(-2.2496719) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5385754) q[1];
sx q[1];
rz(-1.9780206) q[1];
sx q[1];
rz(2.1291705) q[1];
rz(-pi) q[2];
rz(2.0730998) q[3];
sx q[3];
rz(-2.4728259) q[3];
sx q[3];
rz(-2.7947102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3080052) q[2];
sx q[2];
rz(-0.31978017) q[2];
sx q[2];
rz(-1.8316815) q[2];
rz(1.9561907) q[3];
sx q[3];
rz(-2.2535321) q[3];
sx q[3];
rz(1.6185224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.0014521) q[0];
sx q[0];
rz(-1.8403213) q[0];
sx q[0];
rz(2.1885827) q[0];
rz(-3.0630655) q[1];
sx q[1];
rz(-1.086906) q[1];
sx q[1];
rz(-0.79798098) q[1];
rz(-1.2498114) q[2];
sx q[2];
rz(-0.70026308) q[2];
sx q[2];
rz(-1.478298) q[2];
rz(-1.4137951) q[3];
sx q[3];
rz(-1.3678322) q[3];
sx q[3];
rz(2.598138) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
