OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(-1.2712103) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35933094) q[0];
sx q[0];
rz(-1.2027272) q[0];
sx q[0];
rz(2.9666535) q[0];
rz(-pi) q[1];
rz(-1.6127365) q[2];
sx q[2];
rz(-2.0176) q[2];
sx q[2];
rz(2.1697793) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.672294) q[1];
sx q[1];
rz(-2.0672332) q[1];
sx q[1];
rz(-1.4908355) q[1];
x q[2];
rz(-2.1625159) q[3];
sx q[3];
rz(-2.7032529) q[3];
sx q[3];
rz(-2.0548267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-2.326139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6204651) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(-2.492766) q[0];
x q[1];
rz(2.8380727) q[2];
sx q[2];
rz(-2.3887861) q[2];
sx q[2];
rz(0.34740651) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7533469) q[1];
sx q[1];
rz(-1.2699632) q[1];
sx q[1];
rz(1.6080329) q[1];
rz(-pi) q[2];
rz(0.29173298) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(3.047903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(0.63278502) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(-0.99951807) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(0.15655984) q[0];
x q[1];
rz(0.91018422) q[2];
sx q[2];
rz(-1.0599469) q[2];
sx q[2];
rz(2.3597033) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-2.058299) q[1];
sx q[1];
rz(-1.8600149) q[1];
rz(1.7129094) q[3];
sx q[3];
rz(-0.51321533) q[3];
sx q[3];
rz(1.7538479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-2.8667563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9592322) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(1.6596646) q[0];
x q[1];
rz(1.574013) q[2];
sx q[2];
rz(-2.6811757) q[2];
sx q[2];
rz(-2.761063) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3596613) q[1];
sx q[1];
rz(-2.0205824) q[1];
sx q[1];
rz(-2.2284472) q[1];
x q[2];
rz(1.7219909) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(-0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.5083195) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(-0.99779469) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.6246187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51380101) q[0];
sx q[0];
rz(-0.9538981) q[0];
sx q[0];
rz(-0.4517201) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60913182) q[2];
sx q[2];
rz(-0.62437781) q[2];
sx q[2];
rz(0.49027157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5114054) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(-1.9460815) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.153271) q[3];
sx q[3];
rz(-1.054793) q[3];
sx q[3];
rz(-0.9048681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(-1.8185395) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(0.34067672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9522889) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(-3.1097263) q[0];
x q[1];
rz(2.3855626) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(2.0331969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9601945) q[1];
sx q[1];
rz(-1.3409233) q[1];
sx q[1];
rz(-0.48744907) q[1];
rz(-pi) q[2];
rz(-2.7192781) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(0.54780444) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1691549) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8220673) q[0];
sx q[0];
rz(-1.3789346) q[0];
sx q[0];
rz(-2.0454387) q[0];
x q[1];
rz(2.1452227) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(-3.0976354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17864922) q[1];
sx q[1];
rz(-1.4649676) q[1];
sx q[1];
rz(2.9220198) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5246478) q[3];
sx q[3];
rz(-1.5966291) q[3];
sx q[3];
rz(1.3290562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-2.8239992) q[2];
rz(-2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(0.078358738) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5057482) q[0];
sx q[0];
rz(-0.76857476) q[0];
sx q[0];
rz(1.6705253) q[0];
rz(-pi) q[1];
rz(2.3381091) q[2];
sx q[2];
rz(-0.33499559) q[2];
sx q[2];
rz(-1.6840881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2625418) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(-2.2294728) q[1];
x q[2];
rz(1.4140698) q[3];
sx q[3];
rz(-1.6046451) q[3];
sx q[3];
rz(-2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(2.890214) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(2.267568) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6459991) q[0];
sx q[0];
rz(-2.3388303) q[0];
sx q[0];
rz(-1.5587224) q[0];
x q[1];
rz(-1.9782412) q[2];
sx q[2];
rz(-1.2375087) q[2];
sx q[2];
rz(2.8826706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45174949) q[1];
sx q[1];
rz(-2.6362231) q[1];
sx q[1];
rz(1.7187353) q[1];
x q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.2864283) q[3];
sx q[3];
rz(-1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154685) q[0];
sx q[0];
rz(-1.7755531) q[0];
sx q[0];
rz(1.8490851) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12702282) q[2];
sx q[2];
rz(-1.6633031) q[2];
sx q[2];
rz(-2.9542343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5621592) q[1];
sx q[1];
rz(-0.65083083) q[1];
sx q[1];
rz(1.2594373) q[1];
rz(-1.7970656) q[3];
sx q[3];
rz(-2.0824277) q[3];
sx q[3];
rz(2.1567878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.05802352) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(-0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(1.6806867) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5205004) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(0.035049546) q[2];
sx q[2];
rz(-1.1309584) q[2];
sx q[2];
rz(-2.1396648) q[2];
rz(-0.48537985) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
