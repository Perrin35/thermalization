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
rz(-2.4565826) q[0];
sx q[0];
rz(-0.69184875) q[0];
sx q[0];
rz(2.703171) q[0];
rz(-0.11998478) q[1];
sx q[1];
rz(-0.24155231) q[1];
sx q[1];
rz(0.99689364) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50580286) q[0];
sx q[0];
rz(-1.6576901) q[0];
sx q[0];
rz(1.723423) q[0];
rz(-pi) q[1];
rz(-0.14340372) q[2];
sx q[2];
rz(-0.273168) q[2];
sx q[2];
rz(0.96410018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.62694395) q[1];
sx q[1];
rz(-2.6730826) q[1];
sx q[1];
rz(-2.1283988) q[1];
x q[2];
rz(-1.3658872) q[3];
sx q[3];
rz(-1.8031604) q[3];
sx q[3];
rz(0.13162498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74242705) q[2];
sx q[2];
rz(-0.70968598) q[2];
sx q[2];
rz(0.7134552) q[2];
rz(0.82792884) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(-0.92638612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4543318) q[0];
sx q[0];
rz(-2.5255272) q[0];
sx q[0];
rz(1.2607505) q[0];
rz(-0.24185355) q[1];
sx q[1];
rz(-1.2956023) q[1];
sx q[1];
rz(2.5557925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5457348) q[0];
sx q[0];
rz(-0.35484609) q[0];
sx q[0];
rz(1.097358) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5620147) q[2];
sx q[2];
rz(-1.9987371) q[2];
sx q[2];
rz(-1.2195171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9487544) q[1];
sx q[1];
rz(-0.84872972) q[1];
sx q[1];
rz(2.7428736) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5942138) q[3];
sx q[3];
rz(-3.1298198) q[3];
sx q[3];
rz(-2.443832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5213617) q[2];
sx q[2];
rz(-0.62180454) q[2];
sx q[2];
rz(0.28548959) q[2];
rz(-2.1776958) q[3];
sx q[3];
rz(-2.4745092) q[3];
sx q[3];
rz(-2.9662761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46850878) q[0];
sx q[0];
rz(-0.54247576) q[0];
sx q[0];
rz(2.8520404) q[0];
rz(0.30329224) q[1];
sx q[1];
rz(-1.2723609) q[1];
sx q[1];
rz(-1.2264651) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4604386) q[0];
sx q[0];
rz(-1.1484206) q[0];
sx q[0];
rz(0.85737164) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7711839) q[2];
sx q[2];
rz(-2.5287712) q[2];
sx q[2];
rz(-2.1343663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7715622) q[1];
sx q[1];
rz(-0.65415934) q[1];
sx q[1];
rz(-2.0473718) q[1];
x q[2];
rz(-2.3555894) q[3];
sx q[3];
rz(-1.6542098) q[3];
sx q[3];
rz(1.0939897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44846416) q[2];
sx q[2];
rz(-2.144564) q[2];
sx q[2];
rz(2.2779951) q[2];
rz(3.038285) q[3];
sx q[3];
rz(-1.7938675) q[3];
sx q[3];
rz(-0.13532713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0057664) q[0];
sx q[0];
rz(-0.62060452) q[0];
sx q[0];
rz(-0.045106877) q[0];
rz(-0.82334423) q[1];
sx q[1];
rz(-2.7594559) q[1];
sx q[1];
rz(-2.8270922) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5161834) q[0];
sx q[0];
rz(-2.4780031) q[0];
sx q[0];
rz(0.13701464) q[0];
rz(-pi) q[1];
rz(2.1067736) q[2];
sx q[2];
rz(-1.2064486) q[2];
sx q[2];
rz(-1.8672158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.040995251) q[1];
sx q[1];
rz(-1.081859) q[1];
sx q[1];
rz(-1.7799499) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8174386) q[3];
sx q[3];
rz(-1.7739033) q[3];
sx q[3];
rz(0.58989159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0663674) q[2];
sx q[2];
rz(-1.6547357) q[2];
sx q[2];
rz(2.5328947) q[2];
rz(1.6913951) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(0.47877413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4920014) q[0];
sx q[0];
rz(-2.2719125) q[0];
sx q[0];
rz(0.087652303) q[0];
rz(-1.2610669) q[1];
sx q[1];
rz(-1.7365716) q[1];
sx q[1];
rz(-2.2671949) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1058348) q[0];
sx q[0];
rz(-0.89983655) q[0];
sx q[0];
rz(0.041759848) q[0];
x q[1];
rz(-3.1071859) q[2];
sx q[2];
rz(-0.98034562) q[2];
sx q[2];
rz(-1.4077983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4276909) q[1];
sx q[1];
rz(-0.9966313) q[1];
sx q[1];
rz(-0.68820895) q[1];
rz(-pi) q[2];
rz(-1.9256853) q[3];
sx q[3];
rz(-0.57587934) q[3];
sx q[3];
rz(1.3368541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7586691) q[2];
sx q[2];
rz(-2.1047968) q[2];
sx q[2];
rz(0.60302889) q[2];
rz(-2.1488819) q[3];
sx q[3];
rz(-0.38645667) q[3];
sx q[3];
rz(-0.64811289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.8822766) q[0];
sx q[0];
rz(-2.3969605) q[0];
sx q[0];
rz(2.4477006) q[0];
rz(-2.4248185) q[1];
sx q[1];
rz(-2.0718772) q[1];
sx q[1];
rz(2.1778291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0345273) q[0];
sx q[0];
rz(-1.1363159) q[0];
sx q[0];
rz(1.4407506) q[0];
rz(-pi) q[1];
rz(1.7441673) q[2];
sx q[2];
rz(-1.7503305) q[2];
sx q[2];
rz(0.6703568) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55685654) q[1];
sx q[1];
rz(-1.1968166) q[1];
sx q[1];
rz(2.1316865) q[1];
rz(-pi) q[2];
rz(3.0984466) q[3];
sx q[3];
rz(-0.57278297) q[3];
sx q[3];
rz(1.5348032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11513772) q[2];
sx q[2];
rz(-0.28885767) q[2];
sx q[2];
rz(-3.0333983) q[2];
rz(-0.9555971) q[3];
sx q[3];
rz(-0.038766131) q[3];
sx q[3];
rz(2.5819216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3947802) q[0];
sx q[0];
rz(-2.96947) q[0];
sx q[0];
rz(-0.10699233) q[0];
rz(-2.6394898) q[1];
sx q[1];
rz(-1.5109477) q[1];
sx q[1];
rz(-0.14920251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049292663) q[0];
sx q[0];
rz(-2.9098616) q[0];
sx q[0];
rz(0.85706237) q[0];
rz(2.0890498) q[2];
sx q[2];
rz(-2.0796806) q[2];
sx q[2];
rz(-0.38478068) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1625357) q[1];
sx q[1];
rz(-1.2675036) q[1];
sx q[1];
rz(2.3065662) q[1];
rz(-1.427569) q[3];
sx q[3];
rz(-1.1455854) q[3];
sx q[3];
rz(2.3452206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2627829) q[2];
sx q[2];
rz(-2.171319) q[2];
sx q[2];
rz(-0.61857569) q[2];
rz(-0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(0.98606199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2672284) q[0];
sx q[0];
rz(-1.769913) q[0];
sx q[0];
rz(-0.57269639) q[0];
rz(2.122208) q[1];
sx q[1];
rz(-0.23713325) q[1];
sx q[1];
rz(0.076489732) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.50216) q[0];
sx q[0];
rz(-1.0442088) q[0];
sx q[0];
rz(2.9798085) q[0];
rz(-1.3504845) q[2];
sx q[2];
rz(-1.2800084) q[2];
sx q[2];
rz(-0.55401582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44750252) q[1];
sx q[1];
rz(-1.4610485) q[1];
sx q[1];
rz(-0.90417273) q[1];
rz(-pi) q[2];
rz(-2.7143257) q[3];
sx q[3];
rz(-0.92822853) q[3];
sx q[3];
rz(1.515732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1028334) q[2];
sx q[2];
rz(-2.2566654) q[2];
sx q[2];
rz(0.49212512) q[2];
rz(-0.37880185) q[3];
sx q[3];
rz(-2.7087961) q[3];
sx q[3];
rz(2.2545599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73961306) q[0];
sx q[0];
rz(-2.6519863) q[0];
sx q[0];
rz(-0.42994764) q[0];
rz(0.49067378) q[1];
sx q[1];
rz(-0.48777598) q[1];
sx q[1];
rz(-2.1388163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96129629) q[0];
sx q[0];
rz(-1.3081828) q[0];
sx q[0];
rz(-2.0735334) q[0];
x q[1];
rz(0.90012365) q[2];
sx q[2];
rz(-1.5104093) q[2];
sx q[2];
rz(0.65298572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2247315) q[1];
sx q[1];
rz(-1.6186348) q[1];
sx q[1];
rz(0.46617723) q[1];
rz(-pi) q[2];
rz(0.53193502) q[3];
sx q[3];
rz(-0.77456512) q[3];
sx q[3];
rz(2.3331785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3873202) q[2];
sx q[2];
rz(-0.60101271) q[2];
sx q[2];
rz(-2.4499272) q[2];
rz(-3.0129041) q[3];
sx q[3];
rz(-1.5920937) q[3];
sx q[3];
rz(2.5531829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16828951) q[0];
sx q[0];
rz(-0.099270865) q[0];
sx q[0];
rz(-0.6231935) q[0];
rz(-2.9469931) q[1];
sx q[1];
rz(-2.01229) q[1];
sx q[1];
rz(-0.87619877) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19080432) q[0];
sx q[0];
rz(-1.5271458) q[0];
sx q[0];
rz(1.4556637) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4191462) q[2];
sx q[2];
rz(-2.306621) q[2];
sx q[2];
rz(-0.73260546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33223486) q[1];
sx q[1];
rz(-1.4351265) q[1];
sx q[1];
rz(-0.43121454) q[1];
rz(-2.629452) q[3];
sx q[3];
rz(-0.63203963) q[3];
sx q[3];
rz(1.8919945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96560043) q[2];
sx q[2];
rz(-0.24933641) q[2];
sx q[2];
rz(-0.636379) q[2];
rz(-1.6217568) q[3];
sx q[3];
rz(-2.2607925) q[3];
sx q[3];
rz(2.9145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31430055) q[0];
sx q[0];
rz(-1.5999595) q[0];
sx q[0];
rz(-0.40524361) q[0];
rz(1.7397407) q[1];
sx q[1];
rz(-1.6442465) q[1];
sx q[1];
rz(0.13784611) q[1];
rz(-2.9707303) q[2];
sx q[2];
rz(-2.5124585) q[2];
sx q[2];
rz(-0.72825904) q[2];
rz(-0.53608175) q[3];
sx q[3];
rz(-1.5403219) q[3];
sx q[3];
rz(2.2684682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
