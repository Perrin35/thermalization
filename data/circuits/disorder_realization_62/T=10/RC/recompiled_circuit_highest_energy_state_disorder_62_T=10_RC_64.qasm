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
rz(-0.57617968) q[0];
sx q[0];
rz(4.8186995) q[0];
sx q[0];
rz(8.7719593) q[0];
rz(2.6989812) q[1];
sx q[1];
rz(-1.1224597) q[1];
sx q[1];
rz(-2.519156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25678062) q[0];
sx q[0];
rz(-2.8076165) q[0];
sx q[0];
rz(-1.9335174) q[0];
rz(-pi) q[1];
rz(-0.4680674) q[2];
sx q[2];
rz(-0.65535802) q[2];
sx q[2];
rz(-1.6201902) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.039923819) q[1];
sx q[1];
rz(-1.0743272) q[1];
sx q[1];
rz(2.9324004) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7951598) q[3];
sx q[3];
rz(-0.62073675) q[3];
sx q[3];
rz(-0.94514314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5114674) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(-0.54568616) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(-2.1849476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086455258) q[0];
sx q[0];
rz(-2.4392023) q[0];
sx q[0];
rz(2.0461244) q[0];
rz(-2.144004) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(-1.9445317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67162469) q[0];
sx q[0];
rz(-1.8238245) q[0];
sx q[0];
rz(-0.96970153) q[0];
rz(-pi) q[1];
rz(2.3648974) q[2];
sx q[2];
rz(-2.0414845) q[2];
sx q[2];
rz(-1.7913417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5250796) q[1];
sx q[1];
rz(-0.54821205) q[1];
sx q[1];
rz(-1.0031149) q[1];
rz(1.3838816) q[3];
sx q[3];
rz(-1.8946365) q[3];
sx q[3];
rz(1.327654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7257488) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(-1.7572629) q[2];
rz(-1.7978801) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(-1.869092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5196853) q[0];
sx q[0];
rz(-2.3540731) q[0];
sx q[0];
rz(-0.4271048) q[0];
rz(1.2318132) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(2.5274091) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012497525) q[0];
sx q[0];
rz(-1.3007727) q[0];
sx q[0];
rz(-0.34059033) q[0];
x q[1];
rz(-2.3947507) q[2];
sx q[2];
rz(-1.7201506) q[2];
sx q[2];
rz(-0.17490444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28084785) q[1];
sx q[1];
rz(-0.13361803) q[1];
sx q[1];
rz(1.0722776) q[1];
rz(-2.6355686) q[3];
sx q[3];
rz(-1.3576686) q[3];
sx q[3];
rz(-1.0552849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4517639) q[2];
sx q[2];
rz(-0.58388766) q[2];
sx q[2];
rz(-0.27099398) q[2];
rz(2.6594035) q[3];
sx q[3];
rz(-2.2995583) q[3];
sx q[3];
rz(-1.2838001) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2778306) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(0.41896391) q[0];
rz(0.22467443) q[1];
sx q[1];
rz(-1.2089665) q[1];
sx q[1];
rz(3.0640501) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60643911) q[0];
sx q[0];
rz(-0.61302671) q[0];
sx q[0];
rz(-2.0962578) q[0];
x q[1];
rz(0.77247844) q[2];
sx q[2];
rz(-1.315552) q[2];
sx q[2];
rz(2.8359063) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3540707) q[1];
sx q[1];
rz(-1.7131299) q[1];
sx q[1];
rz(-2.700129) q[1];
rz(-pi) q[2];
rz(0.29234286) q[3];
sx q[3];
rz(-1.3707531) q[3];
sx q[3];
rz(0.26250473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.055723995) q[2];
sx q[2];
rz(-2.8102977) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(0.27082768) q[3];
sx q[3];
rz(-1.2374249) q[3];
sx q[3];
rz(2.6868611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5772575) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(-1.7876392) q[0];
rz(-2.1063781) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(1.6114906) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28883067) q[0];
sx q[0];
rz(-2.0977328) q[0];
sx q[0];
rz(-0.88165347) q[0];
x q[1];
rz(-0.3283073) q[2];
sx q[2];
rz(-1.7264778) q[2];
sx q[2];
rz(2.1552483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4799696) q[1];
sx q[1];
rz(-1.9723867) q[1];
sx q[1];
rz(1.5127403) q[1];
rz(0.41095419) q[3];
sx q[3];
rz(-1.6235213) q[3];
sx q[3];
rz(-0.26782521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80663854) q[2];
sx q[2];
rz(-2.1636476) q[2];
sx q[2];
rz(-0.063892603) q[2];
rz(-0.70819267) q[3];
sx q[3];
rz(-1.4539366) q[3];
sx q[3];
rz(1.8074869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.719249) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(-1.0759906) q[0];
rz(0.05323449) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(2.7630973) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66010338) q[0];
sx q[0];
rz(-1.3162587) q[0];
sx q[0];
rz(1.4362364) q[0];
rz(-pi) q[1];
rz(-1.4827864) q[2];
sx q[2];
rz(-1.6122136) q[2];
sx q[2];
rz(2.8633022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15684756) q[1];
sx q[1];
rz(-1.7498921) q[1];
sx q[1];
rz(2.8791134) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1355576) q[3];
sx q[3];
rz(-2.2544686) q[3];
sx q[3];
rz(2.5292252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9516912) q[2];
sx q[2];
rz(-1.6804164) q[2];
sx q[2];
rz(-1.0339197) q[2];
rz(2.2377491) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(-0.044053642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164301) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(-0.12475363) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(0.4932901) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.54087) q[0];
sx q[0];
rz(-2.713795) q[0];
sx q[0];
rz(2.2912301) q[0];
x q[1];
rz(2.9134995) q[2];
sx q[2];
rz(-1.5465294) q[2];
sx q[2];
rz(2.9920141) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9700765) q[1];
sx q[1];
rz(-1.0708717) q[1];
sx q[1];
rz(-1.0830256) q[1];
rz(-pi) q[2];
rz(-1.4629355) q[3];
sx q[3];
rz(-1.2965805) q[3];
sx q[3];
rz(1.3317684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7571681) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(-1.1473131) q[2];
rz(2.3186963) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230474) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(0.4185032) q[0];
rz(-3.0283527) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(1.07771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54544577) q[0];
sx q[0];
rz(-1.8592863) q[0];
sx q[0];
rz(-3.122807) q[0];
rz(-0.32025614) q[2];
sx q[2];
rz(-0.82112193) q[2];
sx q[2];
rz(-2.1964354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5716963) q[1];
sx q[1];
rz(-1.2087617) q[1];
sx q[1];
rz(-0.052459474) q[1];
rz(-pi) q[2];
rz(-1.5111132) q[3];
sx q[3];
rz(-2.2264997) q[3];
sx q[3];
rz(-0.34464097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30617994) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(1.6861247) q[2];
rz(2.9742187) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(-1.0719871) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48285943) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(2.7400548) q[0];
rz(-2.7032779) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(-1.6848791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6922214) q[0];
sx q[0];
rz(-1.5367994) q[0];
sx q[0];
rz(0.0095937455) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4859547) q[2];
sx q[2];
rz(-1.3961424) q[2];
sx q[2];
rz(-1.2740335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46692013) q[1];
sx q[1];
rz(-0.90306696) q[1];
sx q[1];
rz(2.2064462) q[1];
rz(-pi) q[2];
rz(2.0450122) q[3];
sx q[3];
rz(-1.1062396) q[3];
sx q[3];
rz(-1.43917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7598286) q[2];
sx q[2];
rz(-2.9324053) q[2];
sx q[2];
rz(2.4299202) q[2];
rz(0.034218637) q[3];
sx q[3];
rz(-0.70845571) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3592247) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(-0.79886287) q[0];
rz(-2.0955739) q[1];
sx q[1];
rz(-1.9452399) q[1];
sx q[1];
rz(-0.23652133) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4170161) q[0];
sx q[0];
rz(-1.5484705) q[0];
sx q[0];
rz(1.5100368) q[0];
x q[1];
rz(-0.92774074) q[2];
sx q[2];
rz(-1.0272214) q[2];
sx q[2];
rz(-1.3613995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7812209) q[1];
sx q[1];
rz(-2.0016333) q[1];
sx q[1];
rz(1.9442417) q[1];
x q[2];
rz(1.4528689) q[3];
sx q[3];
rz(-0.8632568) q[3];
sx q[3];
rz(0.70048571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.11655) q[2];
sx q[2];
rz(-0.5566842) q[2];
sx q[2];
rz(0.64209783) q[2];
rz(-2.5721278) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(1.1915462) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(0.12267648) q[2];
sx q[2];
rz(-0.79570607) q[2];
sx q[2];
rz(-0.88655587) q[2];
rz(-0.33470086) q[3];
sx q[3];
rz(-1.9751493) q[3];
sx q[3];
rz(-2.5720664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
