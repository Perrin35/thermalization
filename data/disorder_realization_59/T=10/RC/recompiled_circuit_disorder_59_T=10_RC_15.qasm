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
rz(1.6127365) q[2];
sx q[2];
rz(-1.1239927) q[2];
sx q[2];
rz(-0.97181335) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4692987) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(1.6507571) q[1];
rz(1.1997585) q[3];
sx q[3];
rz(-1.8098117) q[3];
sx q[3];
rz(2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(-1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(2.326139) q[1];
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
rz(1.2977807) q[2];
sx q[2];
rz(-2.2815939) q[2];
sx q[2];
rz(-3.0836011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38824575) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(-1.6080329) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4957982) q[3];
sx q[3];
rz(-1.3919953) q[3];
sx q[3];
rz(-1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(2.2667623) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-2.1420746) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8529352) q[0];
sx q[0];
rz(-1.4212199) q[0];
sx q[0];
rz(1.268671) q[0];
x q[1];
rz(-0.91018422) q[2];
sx q[2];
rz(-2.0816457) q[2];
sx q[2];
rz(-0.78188932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0083025) q[1];
sx q[1];
rz(-1.8255207) q[1];
sx q[1];
rz(2.6363274) q[1];
rz(-1.4286832) q[3];
sx q[3];
rz(-2.6283773) q[3];
sx q[3];
rz(1.3877447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6040566) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(2.5615454) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(2.8667563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3954454) q[0];
sx q[0];
rz(-1.4822042) q[0];
sx q[0];
rz(-0.079061411) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0312112) q[2];
sx q[2];
rz(-1.5693671) q[2];
sx q[2];
rz(1.187385) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.91314544) q[1];
rz(-pi) q[2];
rz(1.4196017) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.5083195) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(2.143798) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.6246187) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18407962) q[0];
sx q[0];
rz(-2.3947869) q[0];
sx q[0];
rz(1.0190796) q[0];
x q[1];
rz(-0.60913182) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(-2.6513211) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4167037) q[1];
sx q[1];
rz(-2.5973338) q[1];
sx q[1];
rz(-0.7087899) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3715641) q[3];
sx q[3];
rz(-2.3838245) q[3];
sx q[3];
rz(-0.023035223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8020442) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.5312622) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(0.34067672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37939385) q[0];
sx q[0];
rz(-1.5389991) q[0];
sx q[0];
rz(-1.5048774) q[0];
x q[1];
rz(1.2773877) q[2];
sx q[2];
rz(-2.3050606) q[2];
sx q[2];
rz(0.66213995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50960474) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(1.8297086) q[1];
rz(2.7192781) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(2.5715668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(-3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(2.6766434) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5349605) q[0];
sx q[0];
rz(-2.6323942) q[0];
sx q[0];
rz(-1.9726994) q[0];
rz(-2.2642235) q[2];
sx q[2];
rz(-2.4420028) q[2];
sx q[2];
rz(2.1811461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9629434) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(-2.9220198) q[1];
rz(-pi) q[2];
rz(-3.1157324) q[3];
sx q[3];
rz(-1.6169294) q[3];
sx q[3];
rz(-2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(-0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0067622234) q[0];
sx q[0];
rz(-1.5015331) q[0];
sx q[0];
rz(0.80471054) q[0];
rz(2.9044754) q[2];
sx q[2];
rz(-1.8097005) q[2];
sx q[2];
rz(-0.6616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6137177) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(1.2916958) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7275229) q[3];
sx q[3];
rz(-1.5369475) q[3];
sx q[3];
rz(2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-2.890214) q[2];
rz(0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(-0.87402469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066814518) q[0];
sx q[0];
rz(-1.5621119) q[0];
sx q[0];
rz(2.3735223) q[0];
rz(2.7810077) q[2];
sx q[2];
rz(-1.9546095) q[2];
sx q[2];
rz(-1.6894481) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2487138) q[1];
sx q[1];
rz(-1.4993748) q[1];
sx q[1];
rz(2.0715269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6690528) q[3];
sx q[3];
rz(-2.4774385) q[3];
sx q[3];
rz(-0.00045517552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.1788517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9261242) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(-1.8490851) q[0];
rz(1.4775425) q[2];
sx q[2];
rz(-1.6972731) q[2];
sx q[2];
rz(1.7699514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19503838) q[1];
sx q[1];
rz(-0.95609162) q[1];
sx q[1];
rz(0.22919319) q[1];
x q[2];
rz(-1.7970656) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62109229) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(2.3836366) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-2.010871) q[2];
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
