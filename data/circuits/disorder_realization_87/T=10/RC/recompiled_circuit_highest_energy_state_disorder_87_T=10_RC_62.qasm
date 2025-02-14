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
rz(1.4662161) q[0];
sx q[0];
rz(5.3386547) q[0];
sx q[0];
rz(9.6261779) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(1.4210757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147707) q[0];
sx q[0];
rz(-2.1807007) q[0];
sx q[0];
rz(3.0746835) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1915053) q[2];
sx q[2];
rz(-0.7751152) q[2];
sx q[2];
rz(-2.895854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.101946) q[1];
sx q[1];
rz(-1.9612439) q[1];
sx q[1];
rz(-2.4389308) q[1];
x q[2];
rz(-0.33290036) q[3];
sx q[3];
rz(-1.9431291) q[3];
sx q[3];
rz(-1.8474735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1066771) q[2];
sx q[2];
rz(-2.3372529) q[2];
sx q[2];
rz(2.8625281) q[2];
rz(1.8883102) q[3];
sx q[3];
rz(-2.8918355) q[3];
sx q[3];
rz(0.15160027) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6107553) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(0.24638677) q[0];
rz(1.1950182) q[1];
sx q[1];
rz(-2.2476826) q[1];
sx q[1];
rz(1.5066719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55913299) q[0];
sx q[0];
rz(-1.3252859) q[0];
sx q[0];
rz(-1.1051635) q[0];
rz(-pi) q[1];
rz(1.8737239) q[2];
sx q[2];
rz(-1.611735) q[2];
sx q[2];
rz(-0.28266476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3132224) q[1];
sx q[1];
rz(-1.5130318) q[1];
sx q[1];
rz(1.934113) q[1];
rz(0.43466062) q[3];
sx q[3];
rz(-1.0903597) q[3];
sx q[3];
rz(2.2312589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7404777) q[2];
sx q[2];
rz(-0.97492188) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44322893) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(-2.2337636) q[0];
rz(-0.56697956) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(3.0388015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5823321) q[0];
sx q[0];
rz(-1.0114504) q[0];
sx q[0];
rz(1.3224959) q[0];
rz(-3.0605761) q[2];
sx q[2];
rz(-1.4668494) q[2];
sx q[2];
rz(2.8250717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4990078) q[1];
sx q[1];
rz(-0.73621589) q[1];
sx q[1];
rz(0.6249439) q[1];
rz(2.1249173) q[3];
sx q[3];
rz(-1.8401056) q[3];
sx q[3];
rz(-2.179972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8583782) q[2];
sx q[2];
rz(-0.80867043) q[2];
sx q[2];
rz(-1.4374479) q[2];
rz(-0.39250675) q[3];
sx q[3];
rz(-2.0065353) q[3];
sx q[3];
rz(-0.2984305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0482386) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(1.9527973) q[0];
rz(2.8554754) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(-1.5123222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8522569) q[0];
sx q[0];
rz(-1.2931543) q[0];
sx q[0];
rz(2.1136978) q[0];
x q[1];
rz(-1.9061162) q[2];
sx q[2];
rz(-2.309707) q[2];
sx q[2];
rz(-0.12140935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4439567) q[1];
sx q[1];
rz(-2.5165565) q[1];
sx q[1];
rz(-1.8915063) q[1];
rz(-0.21103255) q[3];
sx q[3];
rz(-1.0730626) q[3];
sx q[3];
rz(-0.56609234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35147038) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(0.67698395) q[2];
rz(-1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(-1.6251132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93382728) q[0];
sx q[0];
rz(-0.25660577) q[0];
sx q[0];
rz(0.024913464) q[0];
rz(2.086153) q[1];
sx q[1];
rz(-2.1565304) q[1];
sx q[1];
rz(-0.28344646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2569766) q[0];
sx q[0];
rz(-2.7714892) q[0];
sx q[0];
rz(0.79106541) q[0];
rz(-pi) q[1];
rz(-2.1589958) q[2];
sx q[2];
rz(-2.9074906) q[2];
sx q[2];
rz(-2.7603619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39549669) q[1];
sx q[1];
rz(-1.0701792) q[1];
sx q[1];
rz(0.034828111) q[1];
rz(-pi) q[2];
rz(-1.153454) q[3];
sx q[3];
rz(-0.71708369) q[3];
sx q[3];
rz(-2.5337608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8289566) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(-0.32583315) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.9841586) q[3];
sx q[3];
rz(-1.5475984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40389898) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(2.7666336) q[0];
rz(2.6629579) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(-0.37014827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2787196) q[0];
sx q[0];
rz(-1.5476883) q[0];
sx q[0];
rz(-1.5824759) q[0];
x q[1];
rz(-0.92801969) q[2];
sx q[2];
rz(-2.9117839) q[2];
sx q[2];
rz(1.9191051) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1380085) q[1];
sx q[1];
rz(-0.89911997) q[1];
sx q[1];
rz(-0.9673311) q[1];
rz(-pi) q[2];
rz(0.010995098) q[3];
sx q[3];
rz(-0.42783005) q[3];
sx q[3];
rz(-1.4233936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39773539) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.3810623) q[2];
rz(-1.0487652) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(2.5019808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.25310707) q[0];
sx q[0];
rz(-2.3024004) q[0];
sx q[0];
rz(0.92415586) q[0];
rz(-2.2249075) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(-1.7402657) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9144672) q[0];
sx q[0];
rz(-1.0773398) q[0];
sx q[0];
rz(-1.6363793) q[0];
rz(-1.7426874) q[2];
sx q[2];
rz(-2.2187833) q[2];
sx q[2];
rz(1.6164219) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8055011) q[1];
sx q[1];
rz(-0.84118836) q[1];
sx q[1];
rz(1.355624) q[1];
x q[2];
rz(2.9244366) q[3];
sx q[3];
rz(-1.0921061) q[3];
sx q[3];
rz(-1.1622805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2283198) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.224996) q[2];
rz(3.1332968) q[3];
sx q[3];
rz(-0.010992916) q[3];
sx q[3];
rz(-0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55815721) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(2.661327) q[0];
rz(-2.9971314) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(2.2772148) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68587559) q[0];
sx q[0];
rz(-2.9714231) q[0];
sx q[0];
rz(1.8491114) q[0];
rz(-pi) q[1];
rz(-0.34274613) q[2];
sx q[2];
rz(-3.0411093) q[2];
sx q[2];
rz(3.036694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56039366) q[1];
sx q[1];
rz(-2.542109) q[1];
sx q[1];
rz(-1.1349212) q[1];
x q[2];
rz(2.3468527) q[3];
sx q[3];
rz(-0.38215548) q[3];
sx q[3];
rz(-0.69821815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(-0.63507357) q[2];
rz(-2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(-2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315345) q[0];
sx q[0];
rz(-2.3681971) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(-0.36525137) q[1];
sx q[1];
rz(-1.7144014) q[1];
sx q[1];
rz(1.7100547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5968558) q[0];
sx q[0];
rz(-1.7628306) q[0];
sx q[0];
rz(-1.3903244) q[0];
rz(-pi) q[1];
rz(0.65406873) q[2];
sx q[2];
rz(-1.1032915) q[2];
sx q[2];
rz(-1.2142177) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86399779) q[1];
sx q[1];
rz(-0.91121948) q[1];
sx q[1];
rz(1.1813753) q[1];
x q[2];
rz(-2.6142653) q[3];
sx q[3];
rz(-1.1865215) q[3];
sx q[3];
rz(-1.7155855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26620904) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(1.4069936) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(1.9269358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491972) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(-0.59463516) q[0];
rz(-1.0058962) q[1];
sx q[1];
rz(-1.2024095) q[1];
sx q[1];
rz(2.2524021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.206736) q[0];
sx q[0];
rz(-1.5201269) q[0];
sx q[0];
rz(1.3502747) q[0];
rz(-pi) q[1];
rz(1.1104305) q[2];
sx q[2];
rz(-1.5074919) q[2];
sx q[2];
rz(-2.4667796) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2584302) q[1];
sx q[1];
rz(-3.0297802) q[1];
sx q[1];
rz(2.1360141) q[1];
rz(1.5327644) q[3];
sx q[3];
rz(-1.0547148) q[3];
sx q[3];
rz(-2.3446159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95120007) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(0.54538837) q[2];
rz(-0.96327463) q[3];
sx q[3];
rz(-1.1759718) q[3];
sx q[3];
rz(1.4639328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.9764523) q[2];
sx q[2];
rz(-0.47274796) q[2];
sx q[2];
rz(-1.8190073) q[2];
rz(0.69915184) q[3];
sx q[3];
rz(-2.4300162) q[3];
sx q[3];
rz(3.1292849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
