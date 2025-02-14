OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2744098) q[0];
sx q[0];
rz(-0.46719587) q[0];
sx q[0];
rz(2.245477) q[0];
rz(5.4652228) q[1];
sx q[1];
rz(4.735534) q[1];
sx q[1];
rz(8.4374333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28361904) q[0];
sx q[0];
rz(-1.19157) q[0];
sx q[0];
rz(0.88421706) q[0];
x q[1];
rz(-0.82168401) q[2];
sx q[2];
rz(-1.1237306) q[2];
sx q[2];
rz(0.38082235) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7973229) q[1];
sx q[1];
rz(-1.8472278) q[1];
sx q[1];
rz(1.7452507) q[1];
rz(-pi) q[2];
rz(-0.36611661) q[3];
sx q[3];
rz(-1.7185655) q[3];
sx q[3];
rz(-0.63945668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92153111) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(0.49911487) q[2];
rz(-2.3257997) q[3];
sx q[3];
rz(-0.59320265) q[3];
sx q[3];
rz(-0.35713404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4895184) q[0];
sx q[0];
rz(-0.3638142) q[0];
sx q[0];
rz(-2.9079085) q[0];
rz(-3.0417327) q[1];
sx q[1];
rz(-1.654511) q[1];
sx q[1];
rz(1.608009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8312992) q[0];
sx q[0];
rz(-1.6673894) q[0];
sx q[0];
rz(-1.6929168) q[0];
x q[1];
rz(-2.7614715) q[2];
sx q[2];
rz(-0.46762662) q[2];
sx q[2];
rz(1.1884226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69791693) q[1];
sx q[1];
rz(-1.6622826) q[1];
sx q[1];
rz(0.96443022) q[1];
rz(-pi) q[2];
rz(2.8321974) q[3];
sx q[3];
rz(-2.0338879) q[3];
sx q[3];
rz(0.17406305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.097215501) q[2];
sx q[2];
rz(-0.88572398) q[2];
sx q[2];
rz(3.0866747) q[2];
rz(-2.4954259) q[3];
sx q[3];
rz(-1.5271527) q[3];
sx q[3];
rz(-3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570479) q[0];
sx q[0];
rz(-0.2453198) q[0];
sx q[0];
rz(0.74485892) q[0];
rz(-1.2767731) q[1];
sx q[1];
rz(-2.1121912) q[1];
sx q[1];
rz(2.2778817) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.540032) q[0];
sx q[0];
rz(-1.1313725) q[0];
sx q[0];
rz(-1.5063398) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50541877) q[2];
sx q[2];
rz(-0.8552455) q[2];
sx q[2];
rz(-0.16929132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7924713) q[1];
sx q[1];
rz(-1.0949666) q[1];
sx q[1];
rz(-0.70648593) q[1];
rz(-pi) q[2];
rz(1.0705804) q[3];
sx q[3];
rz(-2.1693834) q[3];
sx q[3];
rz(0.14414302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.471571) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(2.7460597) q[2];
rz(2.1192571) q[3];
sx q[3];
rz(-0.46415713) q[3];
sx q[3];
rz(-0.49183229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763497) q[0];
sx q[0];
rz(-1.9597766) q[0];
sx q[0];
rz(-0.055971948) q[0];
rz(-2.5296027) q[1];
sx q[1];
rz(-1.5490218) q[1];
sx q[1];
rz(-0.8459808) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9214981) q[0];
sx q[0];
rz(-1.8290231) q[0];
sx q[0];
rz(-2.1058215) q[0];
rz(-0.34807713) q[2];
sx q[2];
rz(-0.98053369) q[2];
sx q[2];
rz(-1.8463617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31797925) q[1];
sx q[1];
rz(-0.3682963) q[1];
sx q[1];
rz(-0.55693926) q[1];
rz(-pi) q[2];
rz(-0.9153644) q[3];
sx q[3];
rz(-0.5753606) q[3];
sx q[3];
rz(1.2498472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1299639) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(-2.4212627) q[2];
rz(0.35308853) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(2.8928355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40097749) q[0];
sx q[0];
rz(-1.688513) q[0];
sx q[0];
rz(-2.154696) q[0];
rz(-1.1445716) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(2.1867337) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3026149) q[0];
sx q[0];
rz(-0.90388008) q[0];
sx q[0];
rz(0.033292183) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0097505) q[2];
sx q[2];
rz(-1.8510783) q[2];
sx q[2];
rz(0.55299475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92063078) q[1];
sx q[1];
rz(-0.60505962) q[1];
sx q[1];
rz(-2.4431538) q[1];
rz(-pi) q[2];
rz(1.3172888) q[3];
sx q[3];
rz(-1.714332) q[3];
sx q[3];
rz(1.9726799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5130875) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(1.0465735) q[2];
rz(0.70932499) q[3];
sx q[3];
rz(-1.9966634) q[3];
sx q[3];
rz(2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76383048) q[0];
sx q[0];
rz(-1.7338294) q[0];
sx q[0];
rz(2.8705226) q[0];
rz(-0.079004869) q[1];
sx q[1];
rz(-2.7075691) q[1];
sx q[1];
rz(-1.431042) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9535429) q[0];
sx q[0];
rz(-2.8859647) q[0];
sx q[0];
rz(0.60582692) q[0];
rz(-pi) q[1];
rz(-2.5779294) q[2];
sx q[2];
rz(-2.0437634) q[2];
sx q[2];
rz(1.2207292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.62785) q[1];
sx q[1];
rz(-1.3249036) q[1];
sx q[1];
rz(2.4469083) q[1];
rz(-pi) q[2];
rz(-1.1523139) q[3];
sx q[3];
rz(-1.7177204) q[3];
sx q[3];
rz(2.4148586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0926823) q[2];
sx q[2];
rz(-1.1834669) q[2];
sx q[2];
rz(3.1177706) q[2];
rz(1.228099) q[3];
sx q[3];
rz(-0.3796328) q[3];
sx q[3];
rz(0.14321271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.1517076) q[0];
sx q[0];
rz(-0.30549529) q[0];
sx q[0];
rz(0.09224961) q[0];
rz(2.8406738) q[1];
sx q[1];
rz(-1.9124799) q[1];
sx q[1];
rz(-1.0801962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0475104) q[0];
sx q[0];
rz(-1.2258397) q[0];
sx q[0];
rz(3.0408919) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2713234) q[2];
sx q[2];
rz(-2.610321) q[2];
sx q[2];
rz(2.9722555) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7799973) q[1];
sx q[1];
rz(-1.4228586) q[1];
sx q[1];
rz(-1.5725488) q[1];
rz(-1.2545414) q[3];
sx q[3];
rz(-2.2407534) q[3];
sx q[3];
rz(3.0681821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13176189) q[2];
sx q[2];
rz(-0.35085446) q[2];
sx q[2];
rz(-0.62823137) q[2];
rz(-0.96588165) q[3];
sx q[3];
rz(-1.5647669) q[3];
sx q[3];
rz(-2.8204744) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1107165) q[0];
sx q[0];
rz(-0.70699152) q[0];
sx q[0];
rz(1.734717) q[0];
rz(-0.12786099) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(2.0157287) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8517673) q[0];
sx q[0];
rz(-0.85677108) q[0];
sx q[0];
rz(1.8147574) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1229805) q[2];
sx q[2];
rz(-1.1191223) q[2];
sx q[2];
rz(-1.9960595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5347157) q[1];
sx q[1];
rz(-1.770853) q[1];
sx q[1];
rz(-1.6354531) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5223658) q[3];
sx q[3];
rz(-1.2017177) q[3];
sx q[3];
rz(-2.823373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0109978) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(-3.0478743) q[2];
rz(-1.4334009) q[3];
sx q[3];
rz(-0.283537) q[3];
sx q[3];
rz(-1.3933498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.4505287) q[0];
sx q[0];
rz(-1.0858303) q[0];
sx q[0];
rz(2.1375256) q[0];
rz(-1.1035236) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(1.7335266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0226949) q[0];
sx q[0];
rz(-1.0556882) q[0];
sx q[0];
rz(1.5638419) q[0];
rz(-pi) q[1];
rz(-2.8084178) q[2];
sx q[2];
rz(-2.8123724) q[2];
sx q[2];
rz(1.4179937) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4346659) q[1];
sx q[1];
rz(-1.3916755) q[1];
sx q[1];
rz(2.0113025) q[1];
rz(2.6226165) q[3];
sx q[3];
rz(-1.6138612) q[3];
sx q[3];
rz(1.8881292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7496926) q[2];
sx q[2];
rz(-1.5211279) q[2];
sx q[2];
rz(-2.5938972) q[2];
rz(-0.31217602) q[3];
sx q[3];
rz(-1.7269937) q[3];
sx q[3];
rz(1.7267905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80426973) q[0];
sx q[0];
rz(-1.4708568) q[0];
sx q[0];
rz(2.3626589) q[0];
rz(2.7670822) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(-0.83190727) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9269126) q[0];
sx q[0];
rz(-0.88145602) q[0];
sx q[0];
rz(1.9124407) q[0];
x q[1];
rz(0.7829297) q[2];
sx q[2];
rz(-2.1179869) q[2];
sx q[2];
rz(-1.7683355) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0147881) q[1];
sx q[1];
rz(-2.2378057) q[1];
sx q[1];
rz(2.2875525) q[1];
x q[2];
rz(2.34487) q[3];
sx q[3];
rz(-1.5490388) q[3];
sx q[3];
rz(0.12783229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-0.35934862) q[2];
sx q[2];
rz(3.1331151) q[2];
rz(-0.40555412) q[3];
sx q[3];
rz(-1.4529934) q[3];
sx q[3];
rz(-1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(0.99210284) q[0];
sx q[0];
rz(-1.6076417) q[0];
sx q[0];
rz(-2.3544307) q[0];
rz(1.0943195) q[1];
sx q[1];
rz(-2.7631187) q[1];
sx q[1];
rz(-0.43597058) q[1];
rz(-1.7827674) q[2];
sx q[2];
rz(-1.4478085) q[2];
sx q[2];
rz(0.87926452) q[2];
rz(-0.094219128) q[3];
sx q[3];
rz(-1.5757271) q[3];
sx q[3];
rz(-1.4477391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
