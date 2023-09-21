OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(-0.43963471) q[0];
sx q[0];
rz(3.0602732) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047526) q[0];
sx q[0];
rz(-1.820571) q[0];
sx q[0];
rz(-0.063637861) q[0];
rz(-pi) q[1];
rz(-0.15687234) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(-3.0770609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1374955) q[1];
sx q[1];
rz(-1.6988519) q[1];
sx q[1];
rz(-1.1346243) q[1];
x q[2];
rz(2.8257915) q[3];
sx q[3];
rz(-1.485802) q[3];
sx q[3];
rz(-2.1080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(-0.23981747) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0175184) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(2.571884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1755978) q[2];
sx q[2];
rz(-0.36075337) q[2];
sx q[2];
rz(-0.15168562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51057306) q[1];
sx q[1];
rz(-2.3595855) q[1];
sx q[1];
rz(1.6454979) q[1];
rz(0.014157045) q[3];
sx q[3];
rz(-1.5870968) q[3];
sx q[3];
rz(1.0316767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5043162) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(-1.3298539) q[2];
rz(-1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-0.66816107) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6324368) q[0];
sx q[0];
rz(-1.8189438) q[0];
sx q[0];
rz(-3.0993673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5553603) q[2];
sx q[2];
rz(-1.4484324) q[2];
sx q[2];
rz(-2.8107779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8733858) q[1];
sx q[1];
rz(-1.8516487) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-1.4170253) q[3];
sx q[3];
rz(-2.4745998) q[3];
sx q[3];
rz(0.94661498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-0.032547396) q[2];
rz(-0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(-1.9354405) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.6548086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2457755) q[0];
sx q[0];
rz(-0.44101199) q[0];
sx q[0];
rz(-2.4367024) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47866486) q[2];
sx q[2];
rz(-2.475111) q[2];
sx q[2];
rz(2.037231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.023991) q[1];
sx q[1];
rz(-1.4211434) q[1];
sx q[1];
rz(-0.28848044) q[1];
rz(-pi) q[2];
rz(-0.82641043) q[3];
sx q[3];
rz(-1.5658169) q[3];
sx q[3];
rz(1.5116364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.460357) q[0];
rz(1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(3.1255186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.21503) q[0];
sx q[0];
rz(-0.83850551) q[0];
sx q[0];
rz(-2.3123884) q[0];
rz(-pi) q[1];
rz(-3.1332544) q[2];
sx q[2];
rz(-2.001488) q[2];
sx q[2];
rz(2.8991933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9863723) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(2.1834454) q[1];
rz(-pi) q[2];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.9330977) q[3];
sx q[3];
rz(3.0294861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1092704) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(0.62189046) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46366102) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-0.93925516) q[0];
x q[1];
rz(2.9333026) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(1.8722033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3983706) q[1];
sx q[1];
rz(-0.9658769) q[1];
sx q[1];
rz(0.066992316) q[1];
x q[2];
rz(2.8858658) q[3];
sx q[3];
rz(-0.41318196) q[3];
sx q[3];
rz(-0.43858389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796824) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-2.4408128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38802559) q[0];
sx q[0];
rz(-1.1209079) q[0];
sx q[0];
rz(3.0168424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3406616) q[2];
sx q[2];
rz(-0.33764687) q[2];
sx q[2];
rz(1.9549191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.339387) q[1];
sx q[1];
rz(-1.6610258) q[1];
sx q[1];
rz(1.0046563) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0573781) q[3];
sx q[3];
rz(-2.299752) q[3];
sx q[3];
rz(-0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7769527) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(3.11943) q[2];
rz(2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448924) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-0.21960396) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58015577) q[2];
sx q[2];
rz(-1.9556502) q[2];
sx q[2];
rz(-0.58154026) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2878694) q[1];
sx q[1];
rz(-0.81765491) q[1];
sx q[1];
rz(-2.5649648) q[1];
x q[2];
rz(1.1710839) q[3];
sx q[3];
rz(-0.69972316) q[3];
sx q[3];
rz(1.8613929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(2.8083535) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.7766215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9997172) q[0];
sx q[0];
rz(-1.5293855) q[0];
sx q[0];
rz(-2.3977604) q[0];
rz(-1.1599837) q[2];
sx q[2];
rz(-1.64752) q[2];
sx q[2];
rz(-3.0551747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2862128) q[1];
sx q[1];
rz(-2.3533822) q[1];
sx q[1];
rz(0.49459313) q[1];
rz(1.2737208) q[3];
sx q[3];
rz(-1.7617412) q[3];
sx q[3];
rz(-1.2862658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68010846) q[0];
sx q[0];
rz(-0.52163863) q[0];
sx q[0];
rz(-2.3011544) q[0];
rz(-pi) q[1];
rz(-2.2553026) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(0.71000368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26607516) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(-2.7060899) q[1];
rz(-0.00023437436) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(-2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0633462) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(1.1408268) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-0.36623476) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(1.2583854) q[2];
sx q[2];
rz(-2.2956143) q[2];
sx q[2];
rz(2.14239) q[2];
rz(-1.0352186) q[3];
sx q[3];
rz(-0.62789161) q[3];
sx q[3];
rz(-0.37554489) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
