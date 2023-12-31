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
rz(-2.4913139) q[1];
sx q[1];
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37047526) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(-0.063637861) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9847203) q[2];
sx q[2];
rz(-1.9777159) q[2];
sx q[2];
rz(-3.0770609) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0040972) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(-2.0069684) q[1];
x q[2];
rz(2.8257915) q[3];
sx q[3];
rz(-1.6557906) q[3];
sx q[3];
rz(-1.0335361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2471182) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.2540934) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-0.23981747) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.996802) q[0];
sx q[0];
rz(-1.4142493) q[0];
sx q[0];
rz(-2.8926204) q[0];
rz(1.2698783) q[2];
sx q[2];
rz(-1.3687203) q[2];
sx q[2];
rz(-0.84504715) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6310196) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(-1.6454979) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5544942) q[3];
sx q[3];
rz(-1.5849515) q[3];
sx q[3];
rz(-2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6372765) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(-1.3298539) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(-1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-1.0659165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.802793) q[0];
sx q[0];
rz(-0.25164139) q[0];
sx q[0];
rz(1.7358857) q[0];
x q[1];
rz(1.714528) q[2];
sx q[2];
rz(-2.1215237) q[2];
sx q[2];
rz(1.8260173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4855811) q[1];
sx q[1];
rz(-1.0080358) q[1];
sx q[1];
rz(2.807711) q[1];
rz(-pi) q[2];
rz(1.7245674) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(3.1090453) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(-2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86768326) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(-1.4867841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8958172) q[0];
sx q[0];
rz(-0.44101199) q[0];
sx q[0];
rz(2.4367024) q[0];
x q[1];
rz(2.6629278) q[2];
sx q[2];
rz(-0.66648167) q[2];
sx q[2];
rz(2.037231) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1292124) q[1];
sx q[1];
rz(-2.817569) q[1];
sx q[1];
rz(-2.654241) q[1];
x q[2];
rz(-3.1348226) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(3.0778459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(2.13307) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0356045) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-3.1255186) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.21503) q[0];
sx q[0];
rz(-2.3030871) q[0];
sx q[0];
rz(0.82920427) q[0];
x q[1];
rz(1.5889421) q[2];
sx q[2];
rz(-0.43076736) q[2];
sx q[2];
rz(-0.22242966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6962471) q[1];
sx q[1];
rz(-1.0187341) q[1];
sx q[1];
rz(-2.6389883) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(0.11210657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.68836132) q[1];
sx q[1];
rz(-1.2247359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8794494) q[0];
sx q[0];
rz(-2.1855542) q[0];
sx q[0];
rz(0.26457796) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1282863) q[2];
sx q[2];
rz(-1.759521) q[2];
sx q[2];
rz(-0.21274266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9310301) q[1];
sx q[1];
rz(-1.6258874) q[1];
sx q[1];
rz(-0.964826) q[1];
rz(-2.740432) q[3];
sx q[3];
rz(-1.4690555) q[3];
sx q[3];
rz(-1.3672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(-2.2018946) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(-0.70077983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535671) q[0];
sx q[0];
rz(-2.0206847) q[0];
sx q[0];
rz(3.0168424) q[0];
x q[1];
rz(-0.3406616) q[2];
sx q[2];
rz(-2.8039458) q[2];
sx q[2];
rz(-1.1866736) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5138445) q[1];
sx q[1];
rz(-0.57250896) q[1];
sx q[1];
rz(-1.7379012) q[1];
rz(-pi) q[2];
rz(0.084214597) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78684029) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(-1.4165075) q[0];
rz(1.7658866) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4448924) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-2.9219887) q[0];
x q[1];
rz(0.6363836) q[2];
sx q[2];
rz(-0.6837662) q[2];
sx q[2];
rz(2.6725339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2878694) q[1];
sx q[1];
rz(-0.81765491) q[1];
sx q[1];
rz(-0.57662782) q[1];
x q[2];
rz(-1.9705087) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(1.2801998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(2.419557) q[0];
rz(2.8083535) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.7766215) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38395912) q[0];
sx q[0];
rz(-2.3968292) q[0];
sx q[0];
rz(3.0804759) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1599837) q[2];
sx q[2];
rz(-1.64752) q[2];
sx q[2];
rz(3.0551747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8553798) q[1];
sx q[1];
rz(-2.3533822) q[1];
sx q[1];
rz(2.6469995) q[1];
rz(-pi) q[2];
rz(2.9421259) q[3];
sx q[3];
rz(-1.27928) q[3];
sx q[3];
rz(0.22649543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(-0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-2.506822) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(2.1733984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68010846) q[0];
sx q[0];
rz(-2.619954) q[0];
sx q[0];
rz(-2.3011544) q[0];
rz(-pi) q[1];
rz(-2.2553026) q[2];
sx q[2];
rz(-0.61294014) q[2];
sx q[2];
rz(-0.71000368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72659513) q[1];
sx q[1];
rz(-1.2508878) q[1];
sx q[1];
rz(1.4160316) q[1];
x q[2];
rz(0.00023437436) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(-0.57516592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(-1.1408268) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(-1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(2.8075519) q[2];
sx q[2];
rz(-0.77790778) q[2];
sx q[2];
rz(1.6890656) q[2];
rz(1.0352186) q[3];
sx q[3];
rz(-2.513701) q[3];
sx q[3];
rz(2.7660478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
