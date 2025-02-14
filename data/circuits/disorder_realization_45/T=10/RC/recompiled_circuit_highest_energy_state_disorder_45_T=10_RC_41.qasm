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
rz(3.6245873) q[0];
sx q[0];
rz(10.935187) q[0];
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(0.72996563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6522927) q[0];
sx q[0];
rz(-2.518939) q[0];
sx q[0];
rz(1.3176509) q[0];
x q[1];
rz(-1.3454622) q[2];
sx q[2];
rz(-0.98774922) q[2];
sx q[2];
rz(3.0449113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12355676) q[1];
sx q[1];
rz(-2.5606321) q[1];
sx q[1];
rz(0.058079795) q[1];
x q[2];
rz(-0.22878756) q[3];
sx q[3];
rz(-1.8513894) q[3];
sx q[3];
rz(-0.56786637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0240747) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(-0.81895858) q[2];
rz(-0.0062746127) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5821424) q[0];
sx q[0];
rz(-1.5903951) q[0];
sx q[0];
rz(2.5421802) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(2.3960466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6744989) q[0];
sx q[0];
rz(-1.8601396) q[0];
sx q[0];
rz(1.4379005) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2471871) q[2];
sx q[2];
rz(-2.3043046) q[2];
sx q[2];
rz(2.3937283) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2520601) q[1];
sx q[1];
rz(-0.65450689) q[1];
sx q[1];
rz(2.8783074) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98394139) q[3];
sx q[3];
rz(-0.7823669) q[3];
sx q[3];
rz(0.38228211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45592371) q[2];
sx q[2];
rz(-2.8440639) q[2];
sx q[2];
rz(-1.4073184) q[2];
rz(-0.84434858) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(-1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5832962) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(2.1287647) q[0];
rz(0.60802513) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(1.8720522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8931173) q[0];
sx q[0];
rz(-0.60113827) q[0];
sx q[0];
rz(2.2918743) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.059194) q[2];
sx q[2];
rz(-2.2401056) q[2];
sx q[2];
rz(1.5562039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3976058) q[1];
sx q[1];
rz(-2.3922046) q[1];
sx q[1];
rz(2.8515062) q[1];
x q[2];
rz(-1.2392912) q[3];
sx q[3];
rz(-2.0391782) q[3];
sx q[3];
rz(-1.6553594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030647) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(1.7648765) q[0];
rz(1.5273013) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10258612) q[0];
sx q[0];
rz(-1.163536) q[0];
sx q[0];
rz(0.38577052) q[0];
x q[1];
rz(1.0215553) q[2];
sx q[2];
rz(-2.5168572) q[2];
sx q[2];
rz(-0.83323375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4473089) q[1];
sx q[1];
rz(-1.0189302) q[1];
sx q[1];
rz(-2.0815316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2305829) q[3];
sx q[3];
rz(-2.0002332) q[3];
sx q[3];
rz(-2.8991606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.7899803) q[2];
rz(2.1221519) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(-1.9701689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919507) q[0];
sx q[0];
rz(-1.4627946) q[0];
sx q[0];
rz(3.0426262) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(2.6720572) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5851368) q[0];
sx q[0];
rz(-3.0017108) q[0];
sx q[0];
rz(-1.742247) q[0];
rz(-pi) q[1];
rz(0.97909285) q[2];
sx q[2];
rz(-0.62034494) q[2];
sx q[2];
rz(-1.1466591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5818122) q[1];
sx q[1];
rz(-0.99039927) q[1];
sx q[1];
rz(-2.0349166) q[1];
x q[2];
rz(-0.23271493) q[3];
sx q[3];
rz(-2.5732627) q[3];
sx q[3];
rz(-1.9459023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88064757) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(1.3060695) q[2];
rz(2.7169054) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(-0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0223618) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(1.3223883) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(-1.7599531) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8441601) q[0];
sx q[0];
rz(-1.6571665) q[0];
sx q[0];
rz(2.3303836) q[0];
rz(-0.99389771) q[2];
sx q[2];
rz(-2.0303147) q[2];
sx q[2];
rz(1.3764868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7834251) q[1];
sx q[1];
rz(-1.7562215) q[1];
sx q[1];
rz(-2.3802451) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9074677) q[3];
sx q[3];
rz(-0.30933274) q[3];
sx q[3];
rz(-0.36946378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3809001) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(-2.6324658) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(-2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(3.0522108) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(2.1481029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3099965) q[0];
sx q[0];
rz(-0.78959268) q[0];
sx q[0];
rz(1.6247276) q[0];
rz(0.69665945) q[2];
sx q[2];
rz(-2.2242745) q[2];
sx q[2];
rz(-0.6763263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3741039) q[1];
sx q[1];
rz(-1.2084157) q[1];
sx q[1];
rz(-2.9844173) q[1];
x q[2];
rz(-0.63103038) q[3];
sx q[3];
rz(-2.0027805) q[3];
sx q[3];
rz(-0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3420458) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(-0.40204027) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(-0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.6680229) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(-1.8079669) q[0];
rz(2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(2.2241101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30675754) q[0];
sx q[0];
rz(-1.4687755) q[0];
sx q[0];
rz(-1.7708885) q[0];
x q[1];
rz(-1.7352261) q[2];
sx q[2];
rz(-2.5114369) q[2];
sx q[2];
rz(3.117331) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5405824) q[1];
sx q[1];
rz(-1.9153908) q[1];
sx q[1];
rz(0.2356727) q[1];
x q[2];
rz(2.5054974) q[3];
sx q[3];
rz(-2.3878752) q[3];
sx q[3];
rz(-1.5223283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.63460073) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(1.8514006) q[2];
rz(2.2318132) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(-0.040987404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(-0.27715096) q[0];
rz(-1.9174891) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(1.3714429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200206) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(-2.8369342) q[0];
rz(-pi) q[1];
rz(-1.1354228) q[2];
sx q[2];
rz(-1.0658385) q[2];
sx q[2];
rz(1.4471042) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4491357) q[1];
sx q[1];
rz(-2.4583011) q[1];
sx q[1];
rz(-1.9016674) q[1];
x q[2];
rz(-0.17454608) q[3];
sx q[3];
rz(-0.8564328) q[3];
sx q[3];
rz(-0.81596953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93854967) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(0.68823632) q[2];
rz(2.8271683) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7670583) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(0.64819711) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2637973) q[0];
sx q[0];
rz(-1.5476942) q[0];
sx q[0];
rz(-3.1277411) q[0];
rz(-0.83309116) q[2];
sx q[2];
rz(-1.0467741) q[2];
sx q[2];
rz(1.3307216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98013377) q[1];
sx q[1];
rz(-0.69936434) q[1];
sx q[1];
rz(3.1070263) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8043824) q[3];
sx q[3];
rz(-0.54223947) q[3];
sx q[3];
rz(0.48619871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28006831) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(2.9649949) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(-1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8661154) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(0.38446174) q[1];
sx q[1];
rz(-1.5059595) q[1];
sx q[1];
rz(-0.62677871) q[1];
rz(-1.8865042) q[2];
sx q[2];
rz(-1.2366625) q[2];
sx q[2];
rz(-1.2319215) q[2];
rz(2.6077059) q[3];
sx q[3];
rz(-1.0718828) q[3];
sx q[3];
rz(-0.26477118) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
