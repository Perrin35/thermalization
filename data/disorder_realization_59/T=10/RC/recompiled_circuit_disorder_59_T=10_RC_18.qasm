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
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(-3.1402052) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3256677) q[0];
sx q[0];
rz(-2.7357833) q[0];
sx q[0];
rz(-1.146846) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5288562) q[2];
sx q[2];
rz(-2.0176) q[2];
sx q[2];
rz(-2.1697793) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.063349799) q[1];
sx q[1];
rz(-1.500505) q[1];
sx q[1];
rz(-0.49777828) q[1];
x q[2];
rz(2.8858521) q[3];
sx q[3];
rz(-1.2107953) q[3];
sx q[3];
rz(0.4482625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(-2.3922065) q[2];
rz(2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(2.326139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52112752) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(0.64882664) q[0];
x q[1];
rz(-2.4120861) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(-1.4480928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38824575) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(1.6080329) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8498597) q[3];
sx q[3];
rz(-0.66666224) q[3];
sx q[3];
rz(-3.047903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6796391) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(-1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777269) q[0];
sx q[0];
rz(-0.3361055) q[0];
sx q[0];
rz(1.1019812) q[0];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(1.0083025) q[1];
sx q[1];
rz(-1.8255207) q[1];
sx q[1];
rz(-2.6363274) q[1];
rz(0.079654982) q[3];
sx q[3];
rz(-1.0632535) q[3];
sx q[3];
rz(-1.5910651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6040566) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(0.58004722) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.375305) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(0.91039175) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(2.8667563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66514689) q[0];
sx q[0];
rz(-3.0229212) q[0];
sx q[0];
rz(2.2975886) q[0];
x q[1];
rz(3.1399973) q[2];
sx q[2];
rz(-2.0312107) q[2];
sx q[2];
rz(0.38412016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7819314) q[1];
sx q[1];
rz(-1.1210103) q[1];
sx q[1];
rz(-2.2284472) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12675385) q[3];
sx q[3];
rz(-2.2634014) q[3];
sx q[3];
rz(2.9664489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.516974) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51380101) q[0];
sx q[0];
rz(-2.1876946) q[0];
sx q[0];
rz(-2.6898726) q[0];
rz(-0.60913182) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(-2.6513211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6301873) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(-1.1955111) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5449104) q[3];
sx q[3];
rz(-1.0718857) q[3];
sx q[3];
rz(0.98017207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9522889) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(-3.1097263) q[0];
x q[1];
rz(-2.3855626) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(-2.0331969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6319879) q[1];
sx q[1];
rz(-2.0443516) q[1];
sx q[1];
rz(-1.3118841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0814704) q[3];
sx q[3];
rz(-1.193207) q[3];
sx q[3];
rz(1.1946354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-2.7959438) q[0];
rz(3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3195254) q[0];
sx q[0];
rz(-1.3789346) q[0];
sx q[0];
rz(2.0454387) q[0];
rz(2.1452227) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(-3.0976354) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.773015) q[1];
sx q[1];
rz(-1.7891208) q[1];
sx q[1];
rz(-1.4623843) q[1];
rz(-pi) q[2];
rz(1.0602337) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(0.75170654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(0.28433329) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(0.078358738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0067622234) q[0];
sx q[0];
rz(-1.5015331) q[0];
sx q[0];
rz(0.80471054) q[0];
rz(-1.3252844) q[2];
sx q[2];
rz(-1.3405372) q[2];
sx q[2];
rz(-0.85207176) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52787493) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(1.2916958) q[1];
x q[2];
rz(-1.7844291) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(2.8182639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(0.87402469) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066814518) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(2.3735223) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7810077) q[2];
sx q[2];
rz(-1.1869831) q[2];
sx q[2];
rz(-1.4521445) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.45174949) q[1];
sx q[1];
rz(-2.6362231) q[1];
sx q[1];
rz(1.7187353) q[1];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0062362) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26319474) q[0];
sx q[0];
rz(-0.34391719) q[0];
sx q[0];
rz(0.92349903) q[0];
rz(1.4775425) q[2];
sx q[2];
rz(-1.4443195) q[2];
sx q[2];
rz(1.3716413) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9465543) q[1];
sx q[1];
rz(-2.185501) q[1];
sx q[1];
rz(-0.22919319) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3445271) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(-2.1567878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.05802352) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(-0.87456885) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-2.3836366) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.1307217) q[2];
sx q[2];
rz(-1.5390839) q[2];
sx q[2];
rz(2.5577953) q[2];
rz(1.7189797) q[3];
sx q[3];
rz(-1.8437181) q[3];
sx q[3];
rz(-3.0317882) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];