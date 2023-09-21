OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(-2.3614376) q[1];
sx q[1];
rz(-1.0649788) q[1];
sx q[1];
rz(0.87632626) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(-2.602306) q[0];
x q[1];
rz(2.1562188) q[2];
sx q[2];
rz(-1.7314163) q[2];
sx q[2];
rz(-2.0656245) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7180011) q[1];
sx q[1];
rz(-2.8112486) q[1];
sx q[1];
rz(1.524339) q[1];
rz(-pi) q[2];
rz(-0.34378864) q[3];
sx q[3];
rz(-0.4292092) q[3];
sx q[3];
rz(-0.49376282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9709388) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(2.412964) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(-0.87444011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9782198) q[0];
sx q[0];
rz(-1.1436497) q[0];
sx q[0];
rz(0.33421974) q[0];
rz(2.9875056) q[2];
sx q[2];
rz(-1.0006957) q[2];
sx q[2];
rz(-0.64683435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9718711) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(0.68193087) q[1];
rz(-pi) q[2];
x q[2];
rz(3.087895) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(-0.43593105) q[2];
rz(2.46051) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(-2.0667734) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(2.7456465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.031035) q[0];
sx q[0];
rz(-1.9651411) q[0];
sx q[0];
rz(2.1179384) q[0];
rz(-pi) q[1];
rz(1.5244353) q[2];
sx q[2];
rz(-2.8472387) q[2];
sx q[2];
rz(2.4183395) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46519854) q[1];
sx q[1];
rz(-0.68323831) q[1];
sx q[1];
rz(3.1264683) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8860502) q[3];
sx q[3];
rz(-1.5934172) q[3];
sx q[3];
rz(0.25657755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6039156) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(0.23920693) q[2];
rz(-3.0662597) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.72162119) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(0.81992942) q[0];
rz(0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-0.23342361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963835) q[0];
sx q[0];
rz(-1.5102981) q[0];
sx q[0];
rz(-1.4738826) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-0.13194612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(0.48308785) q[1];
x q[2];
rz(-2.3573973) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0507811) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(-0.24211611) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500047) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-2.3805526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8262186) q[0];
sx q[0];
rz(-1.6008908) q[0];
sx q[0];
rz(1.8316359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76743482) q[2];
sx q[2];
rz(-0.27786294) q[2];
sx q[2];
rz(-2.1674736) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.099247301) q[1];
sx q[1];
rz(-1.5464916) q[1];
sx q[1];
rz(2.821032) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8499591) q[3];
sx q[3];
rz(-1.8431292) q[3];
sx q[3];
rz(-2.2945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(3.0495194) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(1.7592643) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(-2.81566) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1074166) q[0];
sx q[0];
rz(-1.0506442) q[0];
sx q[0];
rz(-2.0053894) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5599401) q[2];
sx q[2];
rz(-0.51917167) q[2];
sx q[2];
rz(-0.66228629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50800486) q[1];
sx q[1];
rz(-1.272401) q[1];
sx q[1];
rz(-1.7119346) q[1];
x q[2];
rz(-0.16485729) q[3];
sx q[3];
rz(-3.0391209) q[3];
sx q[3];
rz(-2.7285679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59763336) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.4298965) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(-2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(-3.022335) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28173791) q[0];
sx q[0];
rz(-0.87346948) q[0];
sx q[0];
rz(1.9182693) q[0];
rz(-pi) q[1];
rz(1.4085521) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(-0.42019368) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6357928) q[1];
sx q[1];
rz(-1.074082) q[1];
sx q[1];
rz(-0.15091166) q[1];
rz(-pi) q[2];
rz(-1.9414385) q[3];
sx q[3];
rz(-2.2439085) q[3];
sx q[3];
rz(-2.9692269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44935903) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-0.53708491) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83157241) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(-0.88561052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6016156) q[0];
sx q[0];
rz(-2.5936539) q[0];
sx q[0];
rz(0.45666306) q[0];
rz(-0.9912668) q[2];
sx q[2];
rz(-1.5577003) q[2];
sx q[2];
rz(0.0038557204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1427666) q[1];
sx q[1];
rz(-2.4803659) q[1];
sx q[1];
rz(-2.2426474) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86730154) q[3];
sx q[3];
rz(-0.9356381) q[3];
sx q[3];
rz(-2.9150073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4132335) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(2.0098861) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383012) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(-1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-2.0057604) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95490676) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(-2.9493939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7068107) q[2];
sx q[2];
rz(-1.5506622) q[2];
sx q[2];
rz(0.15205631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3867934) q[1];
sx q[1];
rz(-0.42897412) q[1];
sx q[1];
rz(0.11744833) q[1];
x q[2];
rz(0.42580749) q[3];
sx q[3];
rz(-2.1771181) q[3];
sx q[3];
rz(-3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0231126) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(2.2686968) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.89235598) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.409317) q[0];
sx q[0];
rz(-1.1336375) q[0];
sx q[0];
rz(-2.2073295) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1080152) q[2];
sx q[2];
rz(-0.26460755) q[2];
sx q[2];
rz(-2.3679784) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69233209) q[1];
sx q[1];
rz(-2.3093866) q[1];
sx q[1];
rz(-0.87528054) q[1];
rz(1.9479806) q[3];
sx q[3];
rz(-1.0918655) q[3];
sx q[3];
rz(2.9951028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(-2.3847053) q[2];
sx q[2];
rz(-2.6464528) q[2];
sx q[2];
rz(1.1337627) q[2];
rz(-0.79896169) q[3];
sx q[3];
rz(-0.81294717) q[3];
sx q[3];
rz(-3.0019928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];