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
rz(-1.2019914) q[0];
sx q[0];
rz(-2.658598) q[0];
sx q[0];
rz(1.510409) q[0];
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(-2.411627) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18099526) q[0];
sx q[0];
rz(-0.9708403) q[0];
sx q[0];
rz(2.9636895) q[0];
rz(1.3454622) q[2];
sx q[2];
rz(-2.1538434) q[2];
sx q[2];
rz(3.0449113) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.12355676) q[1];
sx q[1];
rz(-0.58096051) q[1];
sx q[1];
rz(0.058079795) q[1];
x q[2];
rz(0.9040622) q[3];
sx q[3];
rz(-2.781467) q[3];
sx q[3];
rz(-1.874543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1175179) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(2.3226341) q[2];
rz(-3.135318) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(-0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5903951) q[0];
sx q[0];
rz(-2.5421802) q[0];
rz(2.2741611) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(-2.3960466) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760147) q[0];
sx q[0];
rz(-1.6981372) q[0];
sx q[0];
rz(-2.8498184) q[0];
rz(-pi) q[1];
rz(-1.8944055) q[2];
sx q[2];
rz(-2.3043046) q[2];
sx q[2];
rz(2.3937283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6709969) q[1];
sx q[1];
rz(-1.4116916) q[1];
sx q[1];
rz(2.5039423) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2621731) q[3];
sx q[3];
rz(-1.169765) q[3];
sx q[3];
rz(1.5123617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-2.8440639) q[2];
sx q[2];
rz(-1.4073184) q[2];
rz(-0.84434858) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(2.1048529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5832962) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(-1.012828) q[0];
rz(2.5335675) q[1];
sx q[1];
rz(-1.5589747) q[1];
sx q[1];
rz(-1.2695405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657261) q[0];
sx q[0];
rz(-1.1320496) q[0];
sx q[0];
rz(-0.42515305) q[0];
rz(-pi) q[1];
rz(-0.082398675) q[2];
sx q[2];
rz(-0.90148704) q[2];
sx q[2];
rz(-1.5853887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3976058) q[1];
sx q[1];
rz(-0.74938801) q[1];
sx q[1];
rz(2.8515062) q[1];
x q[2];
rz(-0.57165159) q[3];
sx q[3];
rz(-2.5749675) q[3];
sx q[3];
rz(2.1376615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6262007) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(1.8499648) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(-2.9279809) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030647) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(1.7648765) q[0];
rz(-1.5273013) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4462858) q[0];
sx q[0];
rz(-0.55342544) q[0];
sx q[0];
rz(-2.2880716) q[0];
x q[1];
rz(-2.7815656) q[2];
sx q[2];
rz(-2.0930778) q[2];
sx q[2];
rz(1.6619267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5508073) q[1];
sx q[1];
rz(-1.1414764) q[1];
sx q[1];
rz(-2.5270259) q[1];
rz(2.2305829) q[3];
sx q[3];
rz(-1.1413594) q[3];
sx q[3];
rz(-0.24243203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.7899803) q[2];
rz(1.0194408) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(-1.9701689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(3.0426262) q[0];
rz(2.4348266) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(-2.6720572) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75824518) q[0];
sx q[0];
rz(-1.4329785) q[0];
sx q[0];
rz(0.024017781) q[0];
x q[1];
rz(2.1624998) q[2];
sx q[2];
rz(-2.5212477) q[2];
sx q[2];
rz(1.9949335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8626735) q[1];
sx q[1];
rz(-1.9544744) q[1];
sx q[1];
rz(2.5088599) q[1];
x q[2];
rz(-1.4245701) q[3];
sx q[3];
rz(-1.0195882) q[3];
sx q[3];
rz(2.220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(1.8355231) q[2];
rz(-2.7169054) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(-0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223618) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(1.3223883) q[0];
rz(2.0139096) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(1.3816396) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9499672) q[0];
sx q[0];
rz(-2.3268496) q[0];
sx q[0];
rz(3.0227468) q[0];
rz(-pi) q[1];
rz(0.83397978) q[2];
sx q[2];
rz(-0.72089855) q[2];
sx q[2];
rz(-0.79254675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7834251) q[1];
sx q[1];
rz(-1.3853711) q[1];
sx q[1];
rz(-0.76134759) q[1];
x q[2];
rz(1.6448037) q[3];
sx q[3];
rz(-1.8714219) q[3];
sx q[3];
rz(0.12412503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7606925) q[2];
sx q[2];
rz(-1.5194632) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(-0.5091269) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5803489) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(-0.089381889) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(2.1481029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.908113) q[0];
sx q[0];
rz(-2.3589239) q[0];
sx q[0];
rz(0.054305768) q[0];
x q[1];
rz(-2.2682796) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(-0.26604929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74718432) q[1];
sx q[1];
rz(-1.7176873) q[1];
sx q[1];
rz(1.9373075) q[1];
rz(-pi) q[2];
rz(-0.66331373) q[3];
sx q[3];
rz(-0.74771008) q[3];
sx q[3];
rz(-1.472689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79954687) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(1.0953995) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(-2.5636165) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47356975) q[0];
sx q[0];
rz(-0.80900017) q[0];
sx q[0];
rz(-1.8079669) q[0];
rz(-0.85583055) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(-2.2241101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8569023) q[0];
sx q[0];
rz(-1.7698341) q[0];
sx q[0];
rz(-3.0375098) q[0];
rz(-pi) q[1];
rz(2.1945004) q[2];
sx q[2];
rz(-1.6674041) q[2];
sx q[2];
rz(-1.4132702) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2189797) q[1];
sx q[1];
rz(-0.41480468) q[1];
sx q[1];
rz(-2.1475683) q[1];
rz(-pi) q[2];
rz(-2.4948984) q[3];
sx q[3];
rz(-1.9894674) q[3];
sx q[3];
rz(0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(-1.290192) q[2];
rz(-2.2318132) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(-0.040987404) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(0.27715096) q[0];
rz(1.2241036) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.3714429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8104102) q[0];
sx q[0];
rz(-0.91282377) q[0];
sx q[0];
rz(1.816279) q[0];
rz(-pi) q[1];
rz(-2.5941237) q[2];
sx q[2];
rz(-1.1927422) q[2];
sx q[2];
rz(-2.7965656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38199785) q[1];
sx q[1];
rz(-1.3642259) q[1];
sx q[1];
rz(-0.91464154) q[1];
x q[2];
rz(-0.17454608) q[3];
sx q[3];
rz(-0.8564328) q[3];
sx q[3];
rz(2.3256231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93854967) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(-0.68823632) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7670583) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(-0.64819711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72361904) q[0];
sx q[0];
rz(-3.1146568) q[0];
sx q[0];
rz(-2.1108146) q[0];
rz(-pi) q[1];
rz(-2.2806281) q[2];
sx q[2];
rz(-0.87536821) q[2];
sx q[2];
rz(2.8784213) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1163016) q[1];
sx q[1];
rz(-0.87193438) q[1];
sx q[1];
rz(-1.541733) q[1];
x q[2];
rz(-2.8043824) q[3];
sx q[3];
rz(-2.5993532) q[3];
sx q[3];
rz(-0.48619871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28006831) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(1.466922) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8661154) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(0.38446174) q[1];
sx q[1];
rz(-1.5059595) q[1];
sx q[1];
rz(-0.62677871) q[1];
rz(-1.2550884) q[2];
sx q[2];
rz(-1.9049302) q[2];
sx q[2];
rz(1.9096712) q[2];
rz(-0.81955688) q[3];
sx q[3];
rz(-0.71376505) q[3];
sx q[3];
rz(-2.5160088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
