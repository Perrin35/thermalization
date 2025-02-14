OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(-2.3655565) q[0];
rz(0.70932055) q[1];
sx q[1];
rz(-1.6242937) q[1];
sx q[1];
rz(-0.59901839) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020928362) q[0];
sx q[0];
rz(-2.0428223) q[0];
sx q[0];
rz(3.0042404) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9089977) q[2];
sx q[2];
rz(-0.71824348) q[2];
sx q[2];
rz(2.3423549) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30556074) q[1];
sx q[1];
rz(-1.2183883) q[1];
sx q[1];
rz(2.6120017) q[1];
rz(-pi) q[2];
rz(0.43384775) q[3];
sx q[3];
rz(-2.3126855) q[3];
sx q[3];
rz(2.4494107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1385931) q[2];
sx q[2];
rz(-1.1370167) q[2];
sx q[2];
rz(-2.9453759) q[2];
rz(1.044322) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(0.31061068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0394548) q[0];
sx q[0];
rz(-0.65334833) q[0];
sx q[0];
rz(0.85773221) q[0];
rz(-2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(-0.78261715) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1038541) q[0];
sx q[0];
rz(-1.4840992) q[0];
sx q[0];
rz(-0.95392761) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8939607) q[2];
sx q[2];
rz(-2.4944381) q[2];
sx q[2];
rz(2.0314856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.23797135) q[1];
sx q[1];
rz(-0.97717932) q[1];
sx q[1];
rz(-1.3376544) q[1];
rz(0.55763839) q[3];
sx q[3];
rz(-2.1774315) q[3];
sx q[3];
rz(2.5733657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9767849) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(-2.5505193) q[2];
rz(2.1137386) q[3];
sx q[3];
rz(-0.63575345) q[3];
sx q[3];
rz(-3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117821) q[0];
sx q[0];
rz(-0.52017838) q[0];
sx q[0];
rz(-2.0853364) q[0];
rz(2.1504869) q[1];
sx q[1];
rz(-1.3737563) q[1];
sx q[1];
rz(0.80449218) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4933518) q[0];
sx q[0];
rz(-1.8364779) q[0];
sx q[0];
rz(2.4881287) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7011004) q[2];
sx q[2];
rz(-0.63281239) q[2];
sx q[2];
rz(1.9595343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3884405) q[1];
sx q[1];
rz(-2.8049251) q[1];
sx q[1];
rz(-0.80516385) q[1];
rz(-pi) q[2];
rz(-0.25193416) q[3];
sx q[3];
rz(-2.478699) q[3];
sx q[3];
rz(-2.5980662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6911917) q[2];
sx q[2];
rz(-2.2914026) q[2];
sx q[2];
rz(2.5737393) q[2];
rz(1.19207) q[3];
sx q[3];
rz(-2.6046533) q[3];
sx q[3];
rz(0.15538628) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464226) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(-1.1871185) q[0];
rz(2.8861956) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(-0.38937169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1215243) q[0];
sx q[0];
rz(-0.25857718) q[0];
sx q[0];
rz(0.54988396) q[0];
x q[1];
rz(3.1115628) q[2];
sx q[2];
rz(-2.2976934) q[2];
sx q[2];
rz(-2.6384357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7070605) q[1];
sx q[1];
rz(-1.3816621) q[1];
sx q[1];
rz(0.906859) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3741174) q[3];
sx q[3];
rz(-1.4405319) q[3];
sx q[3];
rz(1.2127339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7369467) q[2];
sx q[2];
rz(-2.8179822) q[2];
sx q[2];
rz(-1.3673937) q[2];
rz(-0.5395475) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(1.8168943) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668199) q[0];
sx q[0];
rz(-0.19817752) q[0];
sx q[0];
rz(-0.5603801) q[0];
rz(-2.9226774) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(2.970649) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4881011) q[0];
sx q[0];
rz(-0.91918901) q[0];
sx q[0];
rz(0.96820684) q[0];
rz(-pi) q[1];
rz(2.937491) q[2];
sx q[2];
rz(-1.7095672) q[2];
sx q[2];
rz(-2.6134174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.51851398) q[1];
sx q[1];
rz(-0.85757035) q[1];
sx q[1];
rz(0.70116373) q[1];
x q[2];
rz(0.95052743) q[3];
sx q[3];
rz(-1.5135362) q[3];
sx q[3];
rz(-1.4365774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35974744) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(2.183059) q[2];
rz(-2.9790699) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(-1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.573134) q[0];
sx q[0];
rz(-3.077226) q[0];
sx q[0];
rz(-0.74813133) q[0];
rz(1.1185147) q[1];
sx q[1];
rz(-0.41901127) q[1];
sx q[1];
rz(-2.8245139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40561179) q[0];
sx q[0];
rz(-1.6389585) q[0];
sx q[0];
rz(2.8951485) q[0];
rz(2.9002764) q[2];
sx q[2];
rz(-0.20212999) q[2];
sx q[2];
rz(1.8692819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71419628) q[1];
sx q[1];
rz(-0.80679379) q[1];
sx q[1];
rz(0.098735672) q[1];
rz(-pi) q[2];
rz(-0.77692399) q[3];
sx q[3];
rz(-1.5424929) q[3];
sx q[3];
rz(1.5212718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3743484) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(-1.7283776) q[2];
rz(3.0610541) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(0.22182375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0291075) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(1.8674194) q[0];
rz(1.796465) q[1];
sx q[1];
rz(-2.3990264) q[1];
sx q[1];
rz(1.8003731) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2224436) q[0];
sx q[0];
rz(-2.0608927) q[0];
sx q[0];
rz(0.63700139) q[0];
x q[1];
rz(2.6244782) q[2];
sx q[2];
rz(-0.54143751) q[2];
sx q[2];
rz(1.1351897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6051424) q[1];
sx q[1];
rz(-1.7243694) q[1];
sx q[1];
rz(1.7642412) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0242434) q[3];
sx q[3];
rz(-0.92399358) q[3];
sx q[3];
rz(-2.5600159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93366569) q[2];
sx q[2];
rz(-1.1873446) q[2];
sx q[2];
rz(-0.46736091) q[2];
rz(-2.3815239) q[3];
sx q[3];
rz(-1.6642539) q[3];
sx q[3];
rz(-1.2490341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46045983) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(1.6946174) q[0];
rz(0.32577062) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(-1.6206585) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6353519) q[0];
sx q[0];
rz(-1.8367447) q[0];
sx q[0];
rz(-0.5269993) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0041107508) q[2];
sx q[2];
rz(-3.0499524) q[2];
sx q[2];
rz(0.81697538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7829166) q[1];
sx q[1];
rz(-2.7442928) q[1];
sx q[1];
rz(-2.0372169) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86674848) q[3];
sx q[3];
rz(-1.4102077) q[3];
sx q[3];
rz(0.15204568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(0.61526862) q[2];
rz(-1.2728914) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(-0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.0203005) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(2.2913388) q[0];
rz(-2.0203159) q[1];
sx q[1];
rz(-1.3669776) q[1];
sx q[1];
rz(-0.77176315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2126425) q[0];
sx q[0];
rz(-1.54689) q[0];
sx q[0];
rz(0.8769518) q[0];
rz(-pi) q[1];
rz(2.9608126) q[2];
sx q[2];
rz(-1.639699) q[2];
sx q[2];
rz(-3.1071752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62006751) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-3.0848178) q[1];
x q[2];
rz(-1.2466627) q[3];
sx q[3];
rz(-1.3132902) q[3];
sx q[3];
rz(0.83453629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0261592) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(2.3010632) q[2];
rz(-0.080951512) q[3];
sx q[3];
rz(-2.5637124) q[3];
sx q[3];
rz(-0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7131272) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(-1.1567098) q[0];
rz(1.0276065) q[1];
sx q[1];
rz(-1.3778069) q[1];
sx q[1];
rz(2.3311232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58662232) q[0];
sx q[0];
rz(-1.3548618) q[0];
sx q[0];
rz(0.76540375) q[0];
x q[1];
rz(-1.9058305) q[2];
sx q[2];
rz(-1.050569) q[2];
sx q[2];
rz(-2.4851785) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1406882) q[1];
sx q[1];
rz(-1.2869121) q[1];
sx q[1];
rz(-1.8769916) q[1];
rz(-1.719172) q[3];
sx q[3];
rz(-1.0142438) q[3];
sx q[3];
rz(-0.30239964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.345574) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(-1.5239117) q[2];
rz(-1.1003305) q[3];
sx q[3];
rz(-1.9610145) q[3];
sx q[3];
rz(2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1179467) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(-0.042451518) q[1];
sx q[1];
rz(-2.0274542) q[1];
sx q[1];
rz(-1.7995119) q[1];
rz(-1.1645198) q[2];
sx q[2];
rz(-1.6222519) q[2];
sx q[2];
rz(-1.9380515) q[2];
rz(0.45741882) q[3];
sx q[3];
rz(-1.6030967) q[3];
sx q[3];
rz(0.23489192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
