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
rz(5.8435506) q[0];
sx q[0];
rz(6.2018659) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(-2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5188132) q[0];
sx q[0];
rz(-0.25758994) q[0];
sx q[0];
rz(1.326484) q[0];
x q[1];
rz(0.15687234) q[2];
sx q[2];
rz(-1.9777159) q[2];
sx q[2];
rz(0.064531782) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1374955) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(2.0069684) q[1];
rz(-pi) q[2];
rz(2.8257915) q[3];
sx q[3];
rz(-1.485802) q[3];
sx q[3];
rz(1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2471182) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(2.0882864) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(2.9017752) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0175184) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(-0.56970861) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96599483) q[2];
sx q[2];
rz(-0.36075337) q[2];
sx q[2];
rz(-2.989907) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6310196) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(1.6454979) q[1];
rz(-pi) q[2];
rz(0.014157045) q[3];
sx q[3];
rz(-1.5544958) q[3];
sx q[3];
rz(2.1099159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(-1.8117388) q[2];
rz(-1.7999533) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(-1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(-0.66816107) q[0];
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
rz(-2.0695755) q[0];
sx q[0];
rz(-1.5298651) q[0];
sx q[0];
rz(-1.3224365) q[0];
rz(0.5553603) q[2];
sx q[2];
rz(-1.6931603) q[2];
sx q[2];
rz(2.8107779) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26820688) q[1];
sx q[1];
rz(-1.8516487) q[1];
sx q[1];
rz(-2.1594949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12001868) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(2.0002055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.187591) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(0.032547396) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86768326) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.6548086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2457755) q[0];
sx q[0];
rz(-0.44101199) q[0];
sx q[0];
rz(-2.4367024) q[0];
x q[1];
rz(-0.47866486) q[2];
sx q[2];
rz(-0.66648167) q[2];
sx q[2];
rz(2.037231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1292124) q[1];
sx q[1];
rz(-0.32402363) q[1];
sx q[1];
rz(-0.48735168) q[1];
rz(-2.3151822) q[3];
sx q[3];
rz(-1.5658169) q[3];
sx q[3];
rz(1.6299562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-3.0116459) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1059882) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(-1.6812356) q[0];
rz(1.5530855) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-3.1255186) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.21503) q[0];
sx q[0];
rz(-2.3030871) q[0];
sx q[0];
rz(2.3123884) q[0];
x q[1];
rz(-2.0015012) q[2];
sx q[2];
rz(-1.5632196) q[2];
sx q[2];
rz(1.8097144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4453455) q[1];
sx q[1];
rz(-1.0187341) q[1];
sx q[1];
rz(-0.50260431) q[1];
x q[2];
rz(-1.5662976) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(-0.11210657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1092704) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.2348256) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779316) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-2.2023375) q[0];
x q[1];
rz(2.0133063) q[2];
sx q[2];
rz(-1.3820717) q[2];
sx q[2];
rz(-0.21274266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3983706) q[1];
sx q[1];
rz(-2.1757158) q[1];
sx q[1];
rz(-3.0746003) q[1];
rz(-pi) q[2];
rz(-2.740432) q[3];
sx q[3];
rz(-1.4690555) q[3];
sx q[3];
rz(-1.3672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
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
rz(2.5612223) q[0];
rz(1.0549818) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(-2.4408128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472815) q[0];
sx q[0];
rz(-2.6758709) q[0];
sx q[0];
rz(1.3186243) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3406616) q[2];
sx q[2];
rz(-2.8039458) q[2];
sx q[2];
rz(1.1866736) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80220561) q[1];
sx q[1];
rz(-1.6610258) q[1];
sx q[1];
rz(2.1369364) q[1];
x q[2];
rz(0.084214597) q[3];
sx q[3];
rz(-2.299752) q[3];
sx q[3];
rz(-0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(-0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78684029) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.7250852) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211593) q[0];
sx q[0];
rz(-1.7283068) q[0];
sx q[0];
rz(2.349855) q[0];
rz(-pi) q[1];
rz(2.0217998) q[2];
sx q[2];
rz(-2.1037256) q[2];
sx q[2];
rz(1.9110796) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8430431) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(0.72960735) q[1];
x q[2];
rz(0.31657747) q[3];
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
rz(1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.7766215) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1418755) q[0];
sx q[0];
rz(-1.5293855) q[0];
sx q[0];
rz(-2.3977604) q[0];
rz(1.3806254) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(1.6585569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0627928) q[1];
sx q[1];
rz(-1.9140869) q[1];
sx q[1];
rz(-0.72453665) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98718231) q[3];
sx q[3];
rz(-0.35161388) q[3];
sx q[3];
rz(-0.83948638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.4779133) q[0];
rz(-1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(0.96819425) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010846) q[0];
sx q[0];
rz(-0.52163863) q[0];
sx q[0];
rz(-0.84043829) q[0];
rz(-pi) q[1];
rz(2.7231611) q[2];
sx q[2];
rz(-2.0327339) q[2];
sx q[2];
rz(0.07428169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.26607516) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(0.43550272) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2486542) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(0.9957046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(-3.1402804) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(-1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(-0.33404074) q[2];
sx q[2];
rz(-0.77790778) q[2];
sx q[2];
rz(1.6890656) q[2];
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