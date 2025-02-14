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
rz(1.6311128) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(4.4463867) q[1];
sx q[1];
rz(10.211853) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.118282) q[0];
sx q[0];
rz(-0.89561373) q[0];
sx q[0];
rz(-1.4272631) q[0];
x q[1];
rz(-1.6252615) q[2];
sx q[2];
rz(-1.174841) q[2];
sx q[2];
rz(-1.2690282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2752645) q[1];
sx q[1];
rz(-1.5342848) q[1];
sx q[1];
rz(-2.9300772) q[1];
rz(-pi) q[2];
rz(-2.5236457) q[3];
sx q[3];
rz(-2.0701346) q[3];
sx q[3];
rz(1.7955347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2934908) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(-0.21609406) q[2];
rz(-0.50348336) q[3];
sx q[3];
rz(-1.8899906) q[3];
sx q[3];
rz(-3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9965432) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-0.37522069) q[0];
rz(-2.5253865) q[1];
sx q[1];
rz(-2.1245978) q[1];
sx q[1];
rz(2.9527051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4746162) q[0];
sx q[0];
rz(-2.24581) q[0];
sx q[0];
rz(-2.6029909) q[0];
rz(-pi) q[1];
rz(0.52816708) q[2];
sx q[2];
rz(-0.91120992) q[2];
sx q[2];
rz(0.41925493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13470727) q[1];
sx q[1];
rz(-1.7733679) q[1];
sx q[1];
rz(-0.88500513) q[1];
rz(-0.62633212) q[3];
sx q[3];
rz(-0.62343684) q[3];
sx q[3];
rz(-1.745519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6836493) q[2];
sx q[2];
rz(-1.8234437) q[2];
sx q[2];
rz(-1.9981492) q[2];
rz(0.6360561) q[3];
sx q[3];
rz(-0.053074107) q[3];
sx q[3];
rz(1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3627593) q[0];
sx q[0];
rz(-2.9017359) q[0];
sx q[0];
rz(-1.1676316) q[0];
rz(-3.0150343) q[1];
sx q[1];
rz(-1.5153706) q[1];
sx q[1];
rz(0.57356858) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6535049) q[0];
sx q[0];
rz(-1.4801985) q[0];
sx q[0];
rz(2.3132064) q[0];
rz(-1.838348) q[2];
sx q[2];
rz(-1.9256136) q[2];
sx q[2];
rz(1.4894384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4497609) q[1];
sx q[1];
rz(-1.1402601) q[1];
sx q[1];
rz(-3.0887763) q[1];
x q[2];
rz(2.6010547) q[3];
sx q[3];
rz(-2.4063568) q[3];
sx q[3];
rz(2.102004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6557189) q[2];
sx q[2];
rz(-1.8334917) q[2];
sx q[2];
rz(-1.6386848) q[2];
rz(-1.8466628) q[3];
sx q[3];
rz(-0.27532268) q[3];
sx q[3];
rz(0.82726971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5591705) q[0];
sx q[0];
rz(-0.14635135) q[0];
sx q[0];
rz(-0.91648066) q[0];
rz(-0.16042635) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(0.697335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4586737) q[0];
sx q[0];
rz(-1.3240755) q[0];
sx q[0];
rz(1.2400859) q[0];
x q[1];
rz(-2.4795697) q[2];
sx q[2];
rz(-1.3274091) q[2];
sx q[2];
rz(0.13451974) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97084272) q[1];
sx q[1];
rz(-0.80180556) q[1];
sx q[1];
rz(0.46192731) q[1];
rz(-pi) q[2];
x q[2];
rz(3.011441) q[3];
sx q[3];
rz(-1.5569382) q[3];
sx q[3];
rz(2.7133301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15630284) q[2];
sx q[2];
rz(-0.52729815) q[2];
sx q[2];
rz(-2.018833) q[2];
rz(-1.5924234) q[3];
sx q[3];
rz(-1.6808108) q[3];
sx q[3];
rz(-0.39337015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3093981) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(-2.0124117) q[0];
rz(2.26217) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(2.5130491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7593236) q[0];
sx q[0];
rz(-2.970967) q[0];
sx q[0];
rz(0.81214805) q[0];
x q[1];
rz(2.4833335) q[2];
sx q[2];
rz(-2.022201) q[2];
sx q[2];
rz(0.70027225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6595093) q[1];
sx q[1];
rz(-2.1225559) q[1];
sx q[1];
rz(-0.047288069) q[1];
x q[2];
rz(-1.699963) q[3];
sx q[3];
rz(-0.74366513) q[3];
sx q[3];
rz(-0.63275933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2200372) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(-1.7680232) q[2];
rz(-2.8288362) q[3];
sx q[3];
rz(-2.0650568) q[3];
sx q[3];
rz(-3.055174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81109989) q[0];
sx q[0];
rz(-0.70527768) q[0];
sx q[0];
rz(-0.16264859) q[0];
rz(1.2341518) q[1];
sx q[1];
rz(-1.994543) q[1];
sx q[1];
rz(-2.2208234) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2049157) q[0];
sx q[0];
rz(-1.7584118) q[0];
sx q[0];
rz(-3.0392007) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5821136) q[2];
sx q[2];
rz(-1.0019511) q[2];
sx q[2];
rz(-0.43348962) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9879976) q[1];
sx q[1];
rz(-2.3210166) q[1];
sx q[1];
rz(1.6596036) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1699888) q[3];
sx q[3];
rz(-2.3521479) q[3];
sx q[3];
rz(1.1657749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84772253) q[2];
sx q[2];
rz(-2.4203478) q[2];
sx q[2];
rz(3.0976307) q[2];
rz(1.4219159) q[3];
sx q[3];
rz(-0.43216643) q[3];
sx q[3];
rz(0.034984263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111572) q[0];
sx q[0];
rz(-0.80428094) q[0];
sx q[0];
rz(-0.093611896) q[0];
rz(-2.1162107) q[1];
sx q[1];
rz(-1.5461812) q[1];
sx q[1];
rz(0.78790087) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9196531) q[0];
sx q[0];
rz(-2.1062231) q[0];
sx q[0];
rz(0.8485283) q[0];
rz(-0.2621367) q[2];
sx q[2];
rz(-2.5568364) q[2];
sx q[2];
rz(-2.824397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20847296) q[1];
sx q[1];
rz(-1.6115687) q[1];
sx q[1];
rz(2.5426514) q[1];
x q[2];
rz(-3.0209846) q[3];
sx q[3];
rz(-1.7123619) q[3];
sx q[3];
rz(-1.0673351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0660144) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886993) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(-1.8889486) q[0];
rz(3.0839651) q[1];
sx q[1];
rz(-1.3725955) q[1];
sx q[1];
rz(-0.39841121) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8114026) q[0];
sx q[0];
rz(-2.5017362) q[0];
sx q[0];
rz(-0.29307239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7407102) q[2];
sx q[2];
rz(-2.0644098) q[2];
sx q[2];
rz(2.5371671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0115464) q[1];
sx q[1];
rz(-0.1999487) q[1];
sx q[1];
rz(1.5910831) q[1];
rz(1.4640973) q[3];
sx q[3];
rz(-1.0301642) q[3];
sx q[3];
rz(-1.4832527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56732279) q[2];
sx q[2];
rz(-1.2719354) q[2];
sx q[2];
rz(-0.92180139) q[2];
rz(0.70720339) q[3];
sx q[3];
rz(-2.0965529) q[3];
sx q[3];
rz(0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-0.12117584) q[0];
sx q[0];
rz(2.2344672) q[0];
rz(-2.6616197) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(0.60636806) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3009038) q[0];
sx q[0];
rz(-0.32561526) q[0];
sx q[0];
rz(-2.8226844) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81524088) q[2];
sx q[2];
rz(-0.70910197) q[2];
sx q[2];
rz(-2.7486424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94928259) q[1];
sx q[1];
rz(-2.3282166) q[1];
sx q[1];
rz(1.9793012) q[1];
x q[2];
rz(1.6799029) q[3];
sx q[3];
rz(-0.91489313) q[3];
sx q[3];
rz(-2.4885268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.85764641) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(-2.471981) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-2.0958869) q[3];
sx q[3];
rz(-1.4001575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46212101) q[0];
sx q[0];
rz(-1.7993878) q[0];
sx q[0];
rz(-2.2456428) q[0];
rz(3.1013007) q[1];
sx q[1];
rz(-2.1998684) q[1];
sx q[1];
rz(-1.5641854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1791081) q[0];
sx q[0];
rz(-0.024027457) q[0];
sx q[0];
rz(1.3732617) q[0];
rz(-2.4559458) q[2];
sx q[2];
rz(-1.2276936) q[2];
sx q[2];
rz(-2.1941663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40272199) q[1];
sx q[1];
rz(-2.4634857) q[1];
sx q[1];
rz(-2.2540967) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0684929) q[3];
sx q[3];
rz(-0.66876679) q[3];
sx q[3];
rz(2.7947102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3080052) q[2];
sx q[2];
rz(-0.31978017) q[2];
sx q[2];
rz(1.8316815) q[2];
rz(1.1854019) q[3];
sx q[3];
rz(-0.88806051) q[3];
sx q[3];
rz(-1.5230702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1401405) q[0];
sx q[0];
rz(-1.8403213) q[0];
sx q[0];
rz(2.1885827) q[0];
rz(-0.078527191) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-1.8917812) q[2];
sx q[2];
rz(-2.4413296) q[2];
sx q[2];
rz(1.6632947) q[2];
rz(2.9361712) q[3];
sx q[3];
rz(-1.4170437) q[3];
sx q[3];
rz(0.99544256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
