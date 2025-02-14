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
rz(-1.6753766) q[0];
sx q[0];
rz(-2.197062) q[0];
sx q[0];
rz(2.9401927) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(1.4210757) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3823377) q[0];
sx q[0];
rz(-1.6256285) q[0];
sx q[0];
rz(0.95983975) q[0];
rz(-pi) q[1];
rz(2.243648) q[2];
sx q[2];
rz(-1.9899758) q[2];
sx q[2];
rz(-1.3443332) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0396467) q[1];
sx q[1];
rz(-1.1803487) q[1];
sx q[1];
rz(-0.70266188) q[1];
rz(0.87405494) q[3];
sx q[3];
rz(-0.49421453) q[3];
sx q[3];
rz(0.53411247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(2.8625281) q[2];
rz(1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(-0.15160027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(-0.24638677) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.6349207) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4621861) q[0];
sx q[0];
rz(-2.6194318) q[0];
sx q[0];
rz(-2.0798111) q[0];
rz(-pi) q[1];
rz(-1.2678688) q[2];
sx q[2];
rz(-1.611735) q[2];
sx q[2];
rz(2.8589279) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7333201) q[1];
sx q[1];
rz(-2.773914) q[1];
sx q[1];
rz(-1.4094844) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.250397) q[3];
sx q[3];
rz(-0.63624428) q[3];
sx q[3];
rz(1.4435857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40111497) q[2];
sx q[2];
rz(-2.1666708) q[2];
sx q[2];
rz(1.9096036) q[2];
rz(-2.039382) q[3];
sx q[3];
rz(-3.1372742) q[3];
sx q[3];
rz(-1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44322893) q[0];
sx q[0];
rz(-0.6969499) q[0];
sx q[0];
rz(-2.2337636) q[0];
rz(-0.56697956) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(0.10279113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55926052) q[0];
sx q[0];
rz(-2.1301422) q[0];
sx q[0];
rz(1.3224959) q[0];
x q[1];
rz(1.4665098) q[2];
sx q[2];
rz(-1.6513746) q[2];
sx q[2];
rz(1.8957418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6425848) q[1];
sx q[1];
rz(-0.73621589) q[1];
sx q[1];
rz(-0.6249439) q[1];
rz(-pi) q[2];
rz(-2.0538883) q[3];
sx q[3];
rz(-2.5317041) q[3];
sx q[3];
rz(-2.1263378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28321442) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(1.4374479) q[2];
rz(-2.7490859) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0933541) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(-1.1887953) q[0];
rz(-0.28611723) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(-1.5123222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8522569) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(-2.1136978) q[0];
x q[1];
rz(2.7949906) q[2];
sx q[2];
rz(-0.79814974) q[2];
sx q[2];
rz(2.5426898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69763598) q[1];
sx q[1];
rz(-2.5165565) q[1];
sx q[1];
rz(1.2500863) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0779987) q[3];
sx q[3];
rz(-1.385687) q[3];
sx q[3];
rz(2.0349658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35147038) q[2];
sx q[2];
rz(-1.0834379) q[2];
sx q[2];
rz(-0.67698395) q[2];
rz(1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(-1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.93382728) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(0.024913464) q[0];
rz(1.0554396) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(0.28344646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993857) q[0];
sx q[0];
rz(-1.8309347) q[0];
sx q[0];
rz(-0.26630638) q[0];
x q[1];
rz(2.1589958) q[2];
sx q[2];
rz(-2.9074906) q[2];
sx q[2];
rz(-0.38123075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8185618) q[1];
sx q[1];
rz(-0.50172443) q[1];
sx q[1];
rz(-1.5072359) q[1];
rz(-2.8018762) q[3];
sx q[3];
rz(-0.92631683) q[3];
sx q[3];
rz(3.0655053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8289566) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(2.8157595) q[2];
rz(-1.2827986) q[3];
sx q[3];
rz(-1.9841586) q[3];
sx q[3];
rz(1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7376937) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(2.7666336) q[0];
rz(2.6629579) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(2.7714444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8106812) q[0];
sx q[0];
rz(-3.1157012) q[0];
sx q[0];
rz(2.6736892) q[0];
rz(-pi) q[1];
rz(-1.3856851) q[2];
sx q[2];
rz(-1.4338257) q[2];
sx q[2];
rz(2.1632407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3016175) q[1];
sx q[1];
rz(-2.2711922) q[1];
sx q[1];
rz(-2.5216334) q[1];
rz(-pi) q[2];
rz(3.1305976) q[3];
sx q[3];
rz(-2.7137626) q[3];
sx q[3];
rz(1.7181991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39773539) q[2];
sx q[2];
rz(-0.19146679) q[2];
sx q[2];
rz(-1.3810623) q[2];
rz(2.0928275) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(-2.5019808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.8884856) q[0];
sx q[0];
rz(-0.83919224) q[0];
sx q[0];
rz(-0.92415586) q[0];
rz(2.2249075) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(-1.4013269) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9144672) q[0];
sx q[0];
rz(-1.0773398) q[0];
sx q[0];
rz(-1.6363793) q[0];
x q[1];
rz(-1.3989053) q[2];
sx q[2];
rz(-2.2187833) q[2];
sx q[2];
rz(1.5251708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0192593) q[1];
sx q[1];
rz(-2.386555) q[1];
sx q[1];
rz(0.23438677) q[1];
rz(-pi) q[2];
rz(-0.21715607) q[3];
sx q[3];
rz(-1.0921061) q[3];
sx q[3];
rz(-1.1622805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.4618123) q[2];
sx q[2];
rz(-1.224996) q[2];
rz(3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(-2.2744961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55815721) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(2.661327) q[0];
rz(2.9971314) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(2.2772148) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1594161) q[0];
sx q[0];
rz(-1.5242531) q[0];
sx q[0];
rz(1.4070562) q[0];
rz(1.6046674) q[2];
sx q[2];
rz(-1.4761756) q[2];
sx q[2];
rz(0.44924863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0740114) q[1];
sx q[1];
rz(-1.0339104) q[1];
sx q[1];
rz(2.8606961) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79473991) q[3];
sx q[3];
rz(-2.7594372) q[3];
sx q[3];
rz(-0.69821815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.20624837) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(0.63507357) q[2];
rz(0.57279974) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(0.87535453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315345) q[0];
sx q[0];
rz(-2.3681971) q[0];
sx q[0];
rz(-0.12426678) q[0];
rz(2.7763413) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(-1.431538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54473684) q[0];
sx q[0];
rz(-1.7628306) q[0];
sx q[0];
rz(-1.7512683) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65406873) q[2];
sx q[2];
rz(-2.0383012) q[2];
sx q[2];
rz(-1.9273749) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2775949) q[1];
sx q[1];
rz(-0.91121948) q[1];
sx q[1];
rz(-1.9602174) q[1];
x q[2];
rz(-0.52732738) q[3];
sx q[3];
rz(-1.9550712) q[3];
sx q[3];
rz(-1.7155855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26620904) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.3318136) q[3];
sx q[3];
rz(-1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.59239546) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(-2.5469575) q[0];
rz(2.1356964) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(2.2524021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85814637) q[0];
sx q[0];
rz(-0.22617561) q[0];
sx q[0];
rz(1.7986138) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0311622) q[2];
sx q[2];
rz(-1.5074919) q[2];
sx q[2];
rz(-0.6748131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2584302) q[1];
sx q[1];
rz(-3.0297802) q[1];
sx q[1];
rz(2.1360141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6252006) q[3];
sx q[3];
rz(-1.5377196) q[3];
sx q[3];
rz(0.75504485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1903926) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(0.54538837) q[2];
rz(0.96327463) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.4639328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496898) q[0];
sx q[0];
rz(-0.90072537) q[0];
sx q[0];
rz(-1.8186722) q[0];
rz(2.2979965) q[1];
sx q[1];
rz(-1.0782764) q[1];
sx q[1];
rz(-1.2596399) q[1];
rz(-1.1651404) q[2];
sx q[2];
rz(-2.6688447) q[2];
sx q[2];
rz(1.3225854) q[2];
rz(2.5582377) q[3];
sx q[3];
rz(-1.1370549) q[3];
sx q[3];
rz(2.1255253) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
