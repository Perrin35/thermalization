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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(-1.3036417) q[1];
sx q[1];
rz(1.5703896) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1364098) q[0];
sx q[0];
rz(-1.9290975) q[0];
sx q[0];
rz(0.65922064) q[0];
rz(-pi) q[1];
rz(-0.067367359) q[2];
sx q[2];
rz(-1.9695896) q[2];
sx q[2];
rz(2.452564) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6627854) q[1];
sx q[1];
rz(-2.7764455) q[1];
sx q[1];
rz(-1.7313596) q[1];
rz(-1.8800354) q[3];
sx q[3];
rz(-1.5133563) q[3];
sx q[3];
rz(0.32306898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6115173) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(-2.7008936) q[2];
rz(0.079744451) q[3];
sx q[3];
rz(-3.1414882) q[3];
sx q[3];
rz(-2.0174446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663064) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(0.16800198) q[0];
rz(-3.1213144) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(-1.5365938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.080606) q[0];
sx q[0];
rz(-1.3362552) q[0];
sx q[0];
rz(1.3479665) q[0];
x q[1];
rz(-1.5809143) q[2];
sx q[2];
rz(-1.5732906) q[2];
sx q[2];
rz(-0.025257142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78129301) q[1];
sx q[1];
rz(-1.5661245) q[1];
sx q[1];
rz(0.0010990573) q[1];
x q[2];
rz(-2.9320967) q[3];
sx q[3];
rz(-3.0983546) q[3];
sx q[3];
rz(-0.87856495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7961879) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(-1.4076642) q[2];
rz(2.0965072) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(-2.869587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8318091) q[0];
sx q[0];
rz(-0.97667664) q[0];
sx q[0];
rz(-2.580544) q[0];
rz(-0.27684119) q[1];
sx q[1];
rz(-0.012877348) q[1];
sx q[1];
rz(1.3078088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1501859) q[0];
sx q[0];
rz(-1.5298109) q[0];
sx q[0];
rz(-1.2738373) q[0];
x q[1];
rz(3.1342034) q[2];
sx q[2];
rz(-1.5708013) q[2];
sx q[2];
rz(0.50732343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0289291) q[1];
sx q[1];
rz(-1.512292) q[1];
sx q[1];
rz(-0.99716352) q[1];
x q[2];
rz(2.2482292) q[3];
sx q[3];
rz(-1.2033389) q[3];
sx q[3];
rz(1.8830255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7967367) q[2];
sx q[2];
rz(-3.1414746) q[2];
sx q[2];
rz(0.5893839) q[2];
rz(-1.2423337) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(1.7813659) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1347443) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(1.7929329) q[0];
rz(3.1349365) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(-0.033500813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2804256) q[0];
sx q[0];
rz(-0.67641947) q[0];
sx q[0];
rz(-1.2926213) q[0];
rz(-pi) q[1];
x q[1];
rz(0.037390402) q[2];
sx q[2];
rz(-0.11971902) q[2];
sx q[2];
rz(2.7604288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0211612) q[1];
sx q[1];
rz(-0.26909262) q[1];
sx q[1];
rz(-1.6090367) q[1];
rz(-1.7351522) q[3];
sx q[3];
rz(-1.4504315) q[3];
sx q[3];
rz(-3.0422378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1091619) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(3.0884009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.771516) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(0.44684967) q[0];
rz(-0.18613786) q[1];
sx q[1];
rz(-3.0797854) q[1];
sx q[1];
rz(1.7284547) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22767775) q[0];
sx q[0];
rz(-2.9083038) q[0];
sx q[0];
rz(1.4248217) q[0];
x q[1];
rz(-1.1485708) q[2];
sx q[2];
rz(-1.7686426) q[2];
sx q[2];
rz(2.5289218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91292101) q[1];
sx q[1];
rz(-0.067021772) q[1];
sx q[1];
rz(-2.3093501) q[1];
rz(-1.8522315) q[3];
sx q[3];
rz(-1.3556983) q[3];
sx q[3];
rz(-0.18751442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83429217) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(-2.6429122) q[2];
rz(-2.5636766) q[3];
sx q[3];
rz(-0.48360616) q[3];
sx q[3];
rz(2.5448866) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805098) q[0];
sx q[0];
rz(-2.0208277) q[0];
sx q[0];
rz(2.7951796) q[0];
rz(2.5397884) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(-2.3880889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0929209) q[0];
sx q[0];
rz(-1.3792097) q[0];
sx q[0];
rz(-1.3956153) q[0];
rz(-pi) q[1];
rz(-1.4592811) q[2];
sx q[2];
rz(-1.4961622) q[2];
sx q[2];
rz(-2.1837051) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1931455) q[1];
sx q[1];
rz(-2.3420534) q[1];
sx q[1];
rz(0.97481291) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0493204) q[3];
sx q[3];
rz(-1.408395) q[3];
sx q[3];
rz(0.55264651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5715013) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(-1.5241148) q[2];
rz(0.12024719) q[3];
sx q[3];
rz(-3.1382939) q[3];
sx q[3];
rz(-2.6038468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26871249) q[0];
sx q[0];
rz(-0.92395067) q[0];
sx q[0];
rz(-2.9203316) q[0];
rz(1.6809173) q[1];
sx q[1];
rz(-2.2082081) q[1];
sx q[1];
rz(3.0642919) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.59931) q[0];
sx q[0];
rz(-3.0998383) q[0];
sx q[0];
rz(1.5357273) q[0];
rz(-pi) q[1];
rz(2.0650778) q[2];
sx q[2];
rz(-3.1315127) q[2];
sx q[2];
rz(-0.28257698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.059346) q[1];
sx q[1];
rz(-1.4993854) q[1];
sx q[1];
rz(-0.17544984) q[1];
x q[2];
rz(-3.0139913) q[3];
sx q[3];
rz(-1.7831047) q[3];
sx q[3];
rz(-1.177976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3495425) q[2];
sx q[2];
rz(-3.1304066) q[2];
sx q[2];
rz(-2.1816317) q[2];
rz(-0.32944426) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(-0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9462117) q[0];
sx q[0];
rz(-2.52849) q[0];
sx q[0];
rz(3.0323113) q[0];
rz(-0.37336135) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(-1.9111309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7261911) q[0];
sx q[0];
rz(-0.93351948) q[0];
sx q[0];
rz(0.7755875) q[0];
rz(0.19725843) q[2];
sx q[2];
rz(-1.7575118) q[2];
sx q[2];
rz(-0.017053617) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6972887) q[1];
sx q[1];
rz(-3.0737801) q[1];
sx q[1];
rz(1.4013002) q[1];
rz(1.1730186) q[3];
sx q[3];
rz(-0.69044411) q[3];
sx q[3];
rz(0.08716128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5750778) q[2];
sx q[2];
rz(-1.2357624) q[2];
sx q[2];
rz(-1.3285948) q[2];
rz(-1.396842) q[3];
sx q[3];
rz(-3.1378523) q[3];
sx q[3];
rz(-1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0777271) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(-0.57300895) q[0];
rz(-0.30814463) q[1];
sx q[1];
rz(-0.40987086) q[1];
sx q[1];
rz(1.0073957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44237374) q[0];
sx q[0];
rz(-1.6773407) q[0];
sx q[0];
rz(3.1359948) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0976866) q[2];
sx q[2];
rz(-1.708235) q[2];
sx q[2];
rz(3.1071747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.17698174) q[1];
sx q[1];
rz(-1.4876502) q[1];
sx q[1];
rz(3.1008197) q[1];
rz(-pi) q[2];
rz(-1.9898207) q[3];
sx q[3];
rz(-3.1085827) q[3];
sx q[3];
rz(-2.9275683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8241626) q[2];
sx q[2];
rz(-0.62676668) q[2];
sx q[2];
rz(-0.38995788) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-3.1324813) q[3];
sx q[3];
rz(-0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214509) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(0.87156975) q[1];
sx q[1];
rz(-1.8337245) q[1];
sx q[1];
rz(-1.647324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0594306) q[0];
sx q[0];
rz(-0.63157394) q[0];
sx q[0];
rz(2.6078014) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.062537161) q[2];
sx q[2];
rz(-2.5258469) q[2];
sx q[2];
rz(-3.1341022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.941266) q[1];
sx q[1];
rz(-1.2709874) q[1];
sx q[1];
rz(0.34383641) q[1];
x q[2];
rz(2.787519) q[3];
sx q[3];
rz(-3.0210417) q[3];
sx q[3];
rz(-0.15124409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5747052) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(3.1097143) q[2];
rz(-0.77830642) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42302172) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(0.12693916) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(-1.4311287) q[2];
sx q[2];
rz(-1.5763379) q[2];
sx q[2];
rz(-1.3272641) q[2];
rz(2.245458) q[3];
sx q[3];
rz(-1.5265161) q[3];
sx q[3];
rz(2.0731887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
