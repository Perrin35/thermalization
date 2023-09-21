OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(-3.0602732) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9570219) q[0];
sx q[0];
rz(-1.5091358) q[0];
sx q[0];
rz(1.8210568) q[0];
rz(-1.9185669) q[2];
sx q[2];
rz(-2.7070621) q[2];
sx q[2];
rz(0.44473106) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1374955) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(1.1346243) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31580117) q[3];
sx q[3];
rz(-1.485802) q[3];
sx q[3];
rz(2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.2540934) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(-2.9017752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1447906) q[0];
sx q[0];
rz(-1.4142493) q[0];
sx q[0];
rz(2.8926204) q[0];
x q[1];
rz(1.2698783) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(0.84504715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6310196) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(-1.4960947) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5544942) q[3];
sx q[3];
rz(-1.5849515) q[3];
sx q[3];
rz(2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2770237) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(-1.4913303) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-1.0659165) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50915584) q[0];
sx q[0];
rz(-1.3226489) q[0];
sx q[0];
rz(-0.042225348) q[0];
x q[1];
rz(1.714528) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(-1.3155754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4855811) q[1];
sx q[1];
rz(-1.0080358) q[1];
sx q[1];
rz(2.807711) q[1];
rz(-pi) q[2];
rz(1.4170253) q[3];
sx q[3];
rz(-2.4745998) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86768326) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(-1.9354405) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.6548086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0006205) q[0];
sx q[0];
rz(-1.2396493) q[0];
sx q[0];
rz(1.8676057) q[0];
x q[1];
rz(1.9183667) q[2];
sx q[2];
rz(-2.1516557) q[2];
sx q[2];
rz(0.52085224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1292124) q[1];
sx q[1];
rz(-2.817569) q[1];
sx q[1];
rz(-0.48735168) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1348226) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(-0.063746728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(-3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(0.016074093) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92656266) q[0];
sx q[0];
rz(-0.83850551) q[0];
sx q[0];
rz(-2.3123884) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0083382567) q[2];
sx q[2];
rz(-1.1401046) q[2];
sx q[2];
rz(-2.8991933) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9863723) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(0.95814725) q[1];
x q[2];
rz(3.1297242) q[3];
sx q[3];
rz(-0.36232811) q[3];
sx q[3];
rz(-0.12479898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(-2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(-1.2247359) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46366102) q[0];
sx q[0];
rz(-1.7860798) q[0];
sx q[0];
rz(0.93925516) q[0];
rz(-pi) q[1];
rz(1.1512418) q[2];
sx q[2];
rz(-2.662979) q[2];
sx q[2];
rz(1.7350369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3983706) q[1];
sx q[1];
rz(-0.9658769) q[1];
sx q[1];
rz(-0.066992316) q[1];
rz(-pi) q[2];
rz(2.8858658) q[3];
sx q[3];
rz(-2.7284107) q[3];
sx q[3];
rz(-2.7030088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(2.9283004) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(1.0549818) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(-2.4408128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7535671) q[0];
sx q[0];
rz(-1.1209079) q[0];
sx q[0];
rz(-0.1247503) q[0];
rz(-pi) q[1];
rz(-1.6875661) q[2];
sx q[2];
rz(-1.2532557) q[2];
sx q[2];
rz(1.5955398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3157981) q[1];
sx q[1];
rz(-1.0072395) q[1];
sx q[1];
rz(-3.0347996) q[1];
rz(2.3015162) q[3];
sx q[3];
rz(-1.5080161) q[3];
sx q[3];
rz(1.1508133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(-0.74470216) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.3547524) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69670024) q[0];
sx q[0];
rz(-2.3377044) q[0];
sx q[0];
rz(-2.9219887) q[0];
rz(2.5052091) q[2];
sx q[2];
rz(-2.4578265) q[2];
sx q[2];
rz(-0.46905876) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8430431) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(0.72960735) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8250152) q[3];
sx q[3];
rz(-2.205875) q[3];
sx q[3];
rz(-2.3659335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(2.419557) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(1.7766215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9997172) q[0];
sx q[0];
rz(-1.6122072) q[0];
sx q[0];
rz(-2.3977604) q[0];
rz(-1.3806254) q[2];
sx q[2];
rz(-0.41751465) q[2];
sx q[2];
rz(-1.4830358) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2024467) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(-1.1254315) q[1];
rz(-pi) q[2];
rz(-0.98718231) q[3];
sx q[3];
rz(-2.7899788) q[3];
sx q[3];
rz(2.3021063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(-1.9101248) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(1.6636794) q[0];
rz(2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(0.96819425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4614842) q[0];
sx q[0];
rz(-2.619954) q[0];
sx q[0];
rz(-0.84043829) q[0];
rz(-pi) q[1];
rz(-2.0696938) q[2];
sx q[2];
rz(-1.1985156) q[2];
sx q[2];
rz(1.4494214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8755175) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(0.43550272) q[1];
x q[2];
rz(-3.1413583) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(-1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.7447727) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(-0.74950779) q[2];
sx q[2];
rz(-1.3386249) q[2];
sx q[2];
rz(-2.7809536) q[2];
rz(2.1288539) q[3];
sx q[3];
rz(-1.2663208) q[3];
sx q[3];
rz(-2.3940621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];