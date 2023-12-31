OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(1.7083038) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(-1.9967611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.618453) q[0];
sx q[0];
rz(-1.3292399) q[0];
sx q[0];
rz(-3.0735344) q[0];
rz(-0.98923367) q[2];
sx q[2];
rz(-0.49058149) q[2];
sx q[2];
rz(-1.4197592) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45935985) q[1];
sx q[1];
rz(-2.0233005) q[1];
sx q[1];
rz(0.020999055) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54361312) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.9887259) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69764187) q[0];
sx q[0];
rz(-1.9662494) q[0];
sx q[0];
rz(-0.4801428) q[0];
x q[1];
rz(1.2752091) q[2];
sx q[2];
rz(-2.327773) q[2];
sx q[2];
rz(0.57397599) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3165247) q[1];
sx q[1];
rz(-2.1148588) q[1];
sx q[1];
rz(-1.4909418) q[1];
rz(2.2291525) q[3];
sx q[3];
rz(-0.53386253) q[3];
sx q[3];
rz(-0.57953366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(-2.4222597) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-0.57139325) q[0];
rz(-0.96673036) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29887154) q[0];
sx q[0];
rz(-1.5776331) q[0];
sx q[0];
rz(0.30928916) q[0];
rz(-pi) q[1];
rz(1.2450561) q[2];
sx q[2];
rz(-2.2176761) q[2];
sx q[2];
rz(1.5145472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5273683) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(0.36348344) q[1];
rz(-pi) q[2];
rz(-2.7819488) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(1.5639203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30806914) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(-1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546394) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(0.8262659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47047939) q[0];
sx q[0];
rz(-2.2320036) q[0];
sx q[0];
rz(0.83755042) q[0];
rz(-0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(-0.61566478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5845881) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(1.9083379) q[1];
x q[2];
rz(0.51059254) q[3];
sx q[3];
rz(-2.4435352) q[3];
sx q[3];
rz(1.5654246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-1.1337093) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78695801) q[0];
sx q[0];
rz(-1.5411396) q[0];
sx q[0];
rz(1.5251072) q[0];
rz(-pi) q[1];
rz(1.1281625) q[2];
sx q[2];
rz(-2.1652522) q[2];
sx q[2];
rz(2.9885459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6286271) q[1];
sx q[1];
rz(-1.291853) q[1];
sx q[1];
rz(2.1668424) q[1];
x q[2];
rz(-1.5836309) q[3];
sx q[3];
rz(-1.395263) q[3];
sx q[3];
rz(-0.83735355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0952547) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(-2.6020223) q[0];
rz(-pi) q[1];
rz(0.96713709) q[2];
sx q[2];
rz(-1.7269616) q[2];
sx q[2];
rz(-0.16317633) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99299586) q[1];
sx q[1];
rz(-0.33278782) q[1];
sx q[1];
rz(-0.51698835) q[1];
rz(3.1230934) q[3];
sx q[3];
rz(-1.0780932) q[3];
sx q[3];
rz(-2.8871356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-2.6426962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38567625) q[0];
sx q[0];
rz(-1.3323116) q[0];
sx q[0];
rz(3.0576865) q[0];
rz(1.1646541) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(-3.090812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(2.9700301) q[1];
rz(0.76544806) q[3];
sx q[3];
rz(-2.7141889) q[3];
sx q[3];
rz(-3.1020853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.6361902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7751986) q[0];
sx q[0];
rz(-1.1530071) q[0];
sx q[0];
rz(1.9154857) q[0];
x q[1];
rz(0.82377388) q[2];
sx q[2];
rz(-1.082755) q[2];
sx q[2];
rz(-1.4866231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38547388) q[1];
sx q[1];
rz(-1.6033152) q[1];
sx q[1];
rz(-2.4454152) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32780111) q[3];
sx q[3];
rz(-2.1468688) q[3];
sx q[3];
rz(-2.2097261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(0.15052477) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(2.495893) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(-1.8766778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0255233) q[0];
sx q[0];
rz(-0.86206305) q[0];
sx q[0];
rz(0.52225964) q[0];
x q[1];
rz(-1.132498) q[2];
sx q[2];
rz(-2.0813137) q[2];
sx q[2];
rz(2.8709708) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9727271) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(1.6070073) q[1];
rz(1.7725138) q[3];
sx q[3];
rz(-1.1551876) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-2.6169422) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.6773178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20477644) q[0];
sx q[0];
rz(-1.4602772) q[0];
sx q[0];
rz(2.9226581) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2400869) q[2];
sx q[2];
rz(-2.5291574) q[2];
sx q[2];
rz(-1.9233179) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97800868) q[1];
sx q[1];
rz(-2.2469236) q[1];
sx q[1];
rz(-0.34983695) q[1];
rz(-pi) q[2];
rz(1.3947992) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(-1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-0.71184288) q[2];
sx q[2];
rz(-1.1887475) q[2];
sx q[2];
rz(2.9954994) q[2];
rz(2.5682156) q[3];
sx q[3];
rz(-1.2497414) q[3];
sx q[3];
rz(-2.9336815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
