OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.175617) q[0];
sx q[0];
rz(4.7369851) q[0];
sx q[0];
rz(10.935258) q[0];
rz(-2.053396) q[1];
sx q[1];
rz(-1.8367986) q[1];
sx q[1];
rz(0.78707492) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2505456) q[0];
sx q[0];
rz(-2.4536687) q[0];
sx q[0];
rz(2.964818) q[0];
x q[1];
rz(3.0120968) q[2];
sx q[2];
rz(-0.39948764) q[2];
sx q[2];
rz(-2.0129888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2752645) q[1];
sx q[1];
rz(-1.5342848) q[1];
sx q[1];
rz(2.9300772) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6179469) q[3];
sx q[3];
rz(-2.0701346) q[3];
sx q[3];
rz(-1.7955347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8481019) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(0.21609406) q[2];
rz(-2.6381093) q[3];
sx q[3];
rz(-1.2516021) q[3];
sx q[3];
rz(-3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9965432) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-2.766372) q[0];
rz(0.61620617) q[1];
sx q[1];
rz(-2.1245978) q[1];
sx q[1];
rz(-0.18888758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2611885) q[0];
sx q[0];
rz(-1.9828026) q[0];
sx q[0];
rz(-0.82033495) q[0];
rz(-pi) q[1];
rz(0.99449956) q[2];
sx q[2];
rz(-2.3220571) q[2];
sx q[2];
rz(-2.799965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1950043) q[1];
sx q[1];
rz(-0.71041162) q[1];
sx q[1];
rz(-1.2571774) q[1];
rz(-pi) q[2];
rz(0.62633212) q[3];
sx q[3];
rz(-2.5181558) q[3];
sx q[3];
rz(1.3960736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6836493) q[2];
sx q[2];
rz(-1.318149) q[2];
sx q[2];
rz(-1.9981492) q[2];
rz(0.6360561) q[3];
sx q[3];
rz(-0.053074107) q[3];
sx q[3];
rz(1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627593) q[0];
sx q[0];
rz(-0.23985671) q[0];
sx q[0];
rz(1.973961) q[0];
rz(3.0150343) q[1];
sx q[1];
rz(-1.5153706) q[1];
sx q[1];
rz(2.5680241) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.180998) q[0];
sx q[0];
rz(-0.74686909) q[0];
sx q[0];
rz(-1.7043714) q[0];
x q[1];
rz(-1.838348) q[2];
sx q[2];
rz(-1.9256136) q[2];
sx q[2];
rz(-1.6521542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3237649) q[1];
sx q[1];
rz(-2.7080302) q[1];
sx q[1];
rz(-1.4563515) q[1];
rz(-0.6595936) q[3];
sx q[3];
rz(-1.2183739) q[3];
sx q[3];
rz(3.0291758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6557189) q[2];
sx q[2];
rz(-1.3081009) q[2];
sx q[2];
rz(1.6386848) q[2];
rz(1.8466628) q[3];
sx q[3];
rz(-2.86627) q[3];
sx q[3];
rz(-2.3143229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5824222) q[0];
sx q[0];
rz(-0.14635135) q[0];
sx q[0];
rz(-0.91648066) q[0];
rz(-0.16042635) q[1];
sx q[1];
rz(-1.8564936) q[1];
sx q[1];
rz(2.4442576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19577414) q[0];
sx q[0];
rz(-1.8911288) q[0];
sx q[0];
rz(0.26024241) q[0];
x q[1];
rz(0.38390145) q[2];
sx q[2];
rz(-0.6989726) q[2];
sx q[2];
rz(-1.1363824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8749344) q[1];
sx q[1];
rz(-1.2447847) q[1];
sx q[1];
rz(-0.74649109) q[1];
rz(0.13015161) q[3];
sx q[3];
rz(-1.5569382) q[3];
sx q[3];
rz(-2.7133301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15630284) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(-1.1227597) q[2];
rz(1.5924234) q[3];
sx q[3];
rz(-1.4607818) q[3];
sx q[3];
rz(-0.39337015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8321946) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(2.0124117) q[0];
rz(-2.26217) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(0.62854356) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.578772) q[0];
sx q[0];
rz(-1.6878753) q[0];
sx q[0];
rz(1.4463918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1206003) q[2];
sx q[2];
rz(-2.1537184) q[2];
sx q[2];
rz(-1.9457119) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5694738) q[1];
sx q[1];
rz(-0.55357305) q[1];
sx q[1];
rz(1.4941503) q[1];
x q[2];
rz(-1.699963) q[3];
sx q[3];
rz(-2.3979275) q[3];
sx q[3];
rz(0.63275933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2200372) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(1.3735695) q[2];
rz(-0.31275648) q[3];
sx q[3];
rz(-2.0650568) q[3];
sx q[3];
rz(-0.086418644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81109989) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(0.16264859) q[0];
rz(1.9074408) q[1];
sx q[1];
rz(-1.1470497) q[1];
sx q[1];
rz(-2.2208234) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4402078) q[0];
sx q[0];
rz(-0.21344859) q[0];
sx q[0];
rz(-2.0646854) q[0];
x q[1];
rz(-0.92451325) q[2];
sx q[2];
rz(-2.0343896) q[2];
sx q[2];
rz(1.6789503) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.858135) q[1];
sx q[1];
rz(-2.3871536) q[1];
sx q[1];
rz(-3.0467176) q[1];
x q[2];
rz(-2.9716039) q[3];
sx q[3];
rz(-2.3521479) q[3];
sx q[3];
rz(1.9758177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2938701) q[2];
sx q[2];
rz(-0.72124481) q[2];
sx q[2];
rz(-3.0976307) q[2];
rz(-1.4219159) q[3];
sx q[3];
rz(-2.7094262) q[3];
sx q[3];
rz(0.034984263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111572) q[0];
sx q[0];
rz(-0.80428094) q[0];
sx q[0];
rz(-0.093611896) q[0];
rz(1.0253819) q[1];
sx q[1];
rz(-1.5461812) q[1];
sx q[1];
rz(-2.3536918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2219395) q[0];
sx q[0];
rz(-1.0353695) q[0];
sx q[0];
rz(-0.8485283) q[0];
rz(-pi) q[1];
rz(-2.5727243) q[2];
sx q[2];
rz(-1.427257) q[2];
sx q[2];
rz(1.6678866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7514536) q[1];
sx q[1];
rz(-0.97242224) q[1];
sx q[1];
rz(1.5214439) q[1];
x q[2];
rz(-3.0209846) q[3];
sx q[3];
rz(-1.4292307) q[3];
sx q[3];
rz(-2.0742576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0755783) q[2];
sx q[2];
rz(-0.54328537) q[2];
sx q[2];
rz(0.96599609) q[2];
rz(-0.14723369) q[3];
sx q[3];
rz(-1.4598264) q[3];
sx q[3];
rz(-0.60432965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052893355) q[0];
sx q[0];
rz(-1.7558782) q[0];
sx q[0];
rz(-1.8889486) q[0];
rz(3.0839651) q[1];
sx q[1];
rz(-1.7689972) q[1];
sx q[1];
rz(-2.7431814) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3301901) q[0];
sx q[0];
rz(-0.63985642) q[0];
sx q[0];
rz(-0.29307239) q[0];
x q[1];
rz(1.4008825) q[2];
sx q[2];
rz(-2.0644098) q[2];
sx q[2];
rz(2.5371671) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0115464) q[1];
sx q[1];
rz(-0.1999487) q[1];
sx q[1];
rz(-1.5910831) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1755821) q[3];
sx q[3];
rz(-2.5915603) q[3];
sx q[3];
rz(-1.8635141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5742699) q[2];
sx q[2];
rz(-1.8696573) q[2];
sx q[2];
rz(-2.2197913) q[2];
rz(-0.70720339) q[3];
sx q[3];
rz(-1.0450398) q[3];
sx q[3];
rz(0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1929753) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(2.2344672) q[0];
rz(-0.47997296) q[1];
sx q[1];
rz(-1.0555121) q[1];
sx q[1];
rz(-2.5352246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84068882) q[0];
sx q[0];
rz(-0.32561526) q[0];
sx q[0];
rz(-2.8226844) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1290413) q[2];
sx q[2];
rz(-2.0336359) q[2];
sx q[2];
rz(1.3430581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1923101) q[1];
sx q[1];
rz(-2.3282166) q[1];
sx q[1];
rz(1.9793012) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0010292) q[3];
sx q[3];
rz(-0.66359164) q[3];
sx q[3];
rz(-2.666245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85764641) q[2];
sx q[2];
rz(-1.593677) q[2];
sx q[2];
rz(-0.66961163) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-2.0958869) q[3];
sx q[3];
rz(1.7414352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6794716) q[0];
sx q[0];
rz(-1.3422048) q[0];
sx q[0];
rz(0.8959499) q[0];
rz(-0.040291928) q[1];
sx q[1];
rz(-2.1998684) q[1];
sx q[1];
rz(-1.5641854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58916726) q[0];
sx q[0];
rz(-1.5660813) q[0];
sx q[0];
rz(1.5472359) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1383923) q[2];
sx q[2];
rz(-0.93195019) q[2];
sx q[2];
rz(2.2496719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5385754) q[1];
sx q[1];
rz(-1.9780206) q[1];
sx q[1];
rz(2.1291705) q[1];
rz(-2.1765616) q[3];
sx q[3];
rz(-1.26767) q[3];
sx q[3];
rz(2.3245963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8335874) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(1.3099111) q[2];
rz(-1.1854019) q[3];
sx q[3];
rz(-2.2535321) q[3];
sx q[3];
rz(-1.5230702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401405) q[0];
sx q[0];
rz(-1.3012713) q[0];
sx q[0];
rz(-0.95300994) q[0];
rz(-0.078527191) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-0.25987249) q[2];
sx q[2];
rz(-0.91283471) q[2];
sx q[2];
rz(1.2531493) q[2];
rz(2.4918741) q[3];
sx q[3];
rz(-0.25593723) q[3];
sx q[3];
rz(-1.2091936) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
