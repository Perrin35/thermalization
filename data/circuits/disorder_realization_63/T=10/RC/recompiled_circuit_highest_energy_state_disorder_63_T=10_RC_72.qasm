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
rz(1.4242564) q[0];
sx q[0];
rz(-1.2196701) q[0];
sx q[0];
rz(2.5665459) q[0];
rz(-2.9991034) q[1];
sx q[1];
rz(-1.2187076) q[1];
sx q[1];
rz(-2.277318) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605715) q[0];
sx q[0];
rz(-0.29566524) q[0];
sx q[0];
rz(-0.85938248) q[0];
x q[1];
rz(1.7792542) q[2];
sx q[2];
rz(-2.499806) q[2];
sx q[2];
rz(-1.2520977) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6576523) q[1];
sx q[1];
rz(-2.6115656) q[1];
sx q[1];
rz(3.1375259) q[1];
rz(0.56700403) q[3];
sx q[3];
rz(-2.0006913) q[3];
sx q[3];
rz(0.56741949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.076866604) q[2];
sx q[2];
rz(-1.4996424) q[2];
sx q[2];
rz(-2.1323252) q[2];
rz(-2.3276954) q[3];
sx q[3];
rz(-0.32865694) q[3];
sx q[3];
rz(-2.8024659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0547884) q[0];
sx q[0];
rz(-1.4084933) q[0];
sx q[0];
rz(0.24222294) q[0];
rz(1.9107266) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(-0.79808527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4373598) q[0];
sx q[0];
rz(-0.95600545) q[0];
sx q[0];
rz(-1.4324709) q[0];
rz(-pi) q[1];
rz(2.9954916) q[2];
sx q[2];
rz(-1.2309936) q[2];
sx q[2];
rz(-0.81865849) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.07814503) q[1];
sx q[1];
rz(-1.1157827) q[1];
sx q[1];
rz(1.9719506) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7834218) q[3];
sx q[3];
rz(-1.4205975) q[3];
sx q[3];
rz(1.1800486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4855087) q[2];
sx q[2];
rz(-2.0723497) q[2];
sx q[2];
rz(0.83267027) q[2];
rz(-2.5101856) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(-0.00028636534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1239419) q[0];
sx q[0];
rz(-0.90921062) q[0];
sx q[0];
rz(2.9538474) q[0];
rz(0.84367696) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(0.63293308) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5143941) q[0];
sx q[0];
rz(-0.36982515) q[0];
sx q[0];
rz(-1.7677977) q[0];
rz(2.454921) q[2];
sx q[2];
rz(-1.8455659) q[2];
sx q[2];
rz(-2.6426154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5714076) q[1];
sx q[1];
rz(-1.8708036) q[1];
sx q[1];
rz(-2.6064998) q[1];
rz(-pi) q[2];
rz(-1.2237596) q[3];
sx q[3];
rz(-1.4489902) q[3];
sx q[3];
rz(2.1852139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8741499) q[2];
sx q[2];
rz(-2.1953857) q[2];
sx q[2];
rz(-1.8702742) q[2];
rz(-1.6455796) q[3];
sx q[3];
rz(-1.4964024) q[3];
sx q[3];
rz(0.97150826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520168) q[0];
sx q[0];
rz(-2.3887964) q[0];
sx q[0];
rz(-0.65943199) q[0];
rz(-1.3176428) q[1];
sx q[1];
rz(-1.3392071) q[1];
sx q[1];
rz(-0.79016322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20445828) q[0];
sx q[0];
rz(-0.35872981) q[0];
sx q[0];
rz(2.782194) q[0];
x q[1];
rz(1.0264364) q[2];
sx q[2];
rz(-3.1033278) q[2];
sx q[2];
rz(-2.1452034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8276538) q[1];
sx q[1];
rz(-2.0472333) q[1];
sx q[1];
rz(1.6446132) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2983367) q[3];
sx q[3];
rz(-2.108223) q[3];
sx q[3];
rz(0.3934653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.059375199) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(2.8455632) q[2];
rz(0.56240231) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(3.0276022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89428467) q[0];
sx q[0];
rz(-1.9534651) q[0];
sx q[0];
rz(2.9344015) q[0];
rz(1.0317135) q[1];
sx q[1];
rz(-1.1528015) q[1];
sx q[1];
rz(-1.4854887) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15940672) q[0];
sx q[0];
rz(-0.80538087) q[0];
sx q[0];
rz(-0.55434395) q[0];
rz(-0.29003044) q[2];
sx q[2];
rz(-2.1129015) q[2];
sx q[2];
rz(-2.5190767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1904162) q[1];
sx q[1];
rz(-0.37280478) q[1];
sx q[1];
rz(0.35721161) q[1];
x q[2];
rz(1.6909825) q[3];
sx q[3];
rz(-0.73420364) q[3];
sx q[3];
rz(2.7199573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1032054) q[2];
sx q[2];
rz(-2.3296671) q[2];
sx q[2];
rz(-2.0337598) q[2];
rz(0.92249089) q[3];
sx q[3];
rz(-1.731571) q[3];
sx q[3];
rz(-0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.51782) q[0];
sx q[0];
rz(-2.6277442) q[0];
sx q[0];
rz(-0.22860953) q[0];
rz(2.8672583) q[1];
sx q[1];
rz(-0.81293303) q[1];
sx q[1];
rz(-2.6913604) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4086356) q[0];
sx q[0];
rz(-1.3191603) q[0];
sx q[0];
rz(0.77045124) q[0];
x q[1];
rz(-1.3415496) q[2];
sx q[2];
rz(-1.563213) q[2];
sx q[2];
rz(0.46310616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6321559) q[1];
sx q[1];
rz(-0.56427279) q[1];
sx q[1];
rz(-0.48308259) q[1];
rz(-2.130533) q[3];
sx q[3];
rz(-0.478906) q[3];
sx q[3];
rz(0.68747179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7776514) q[2];
sx q[2];
rz(-0.55299091) q[2];
sx q[2];
rz(0.027776329) q[2];
rz(0.76644301) q[3];
sx q[3];
rz(-1.1602297) q[3];
sx q[3];
rz(-0.69507039) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22208333) q[0];
sx q[0];
rz(-1.5064025) q[0];
sx q[0];
rz(2.5191504) q[0];
rz(2.4355603) q[1];
sx q[1];
rz(-2.249554) q[1];
sx q[1];
rz(-2.5546254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76510274) q[0];
sx q[0];
rz(-2.4060632) q[0];
sx q[0];
rz(0.88659783) q[0];
rz(-pi) q[1];
rz(-2.6640011) q[2];
sx q[2];
rz(-0.53672635) q[2];
sx q[2];
rz(2.0411185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7365954) q[1];
sx q[1];
rz(-0.51144236) q[1];
sx q[1];
rz(0.14528017) q[1];
rz(-pi) q[2];
rz(-1.5587646) q[3];
sx q[3];
rz(-0.71206743) q[3];
sx q[3];
rz(-2.4648417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.075835) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(2.5471121) q[2];
rz(0.25932702) q[3];
sx q[3];
rz(-1.8766873) q[3];
sx q[3];
rz(1.9243141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.1031621) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(2.5627947) q[0];
rz(1.3102866) q[1];
sx q[1];
rz(-1.879004) q[1];
sx q[1];
rz(2.328918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8978764) q[0];
sx q[0];
rz(-1.0930632) q[0];
sx q[0];
rz(-1.7327139) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2590821) q[2];
sx q[2];
rz(-2.1171085) q[2];
sx q[2];
rz(-1.3457042) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26749565) q[1];
sx q[1];
rz(-2.2363642) q[1];
sx q[1];
rz(-2.0366621) q[1];
rz(-0.854354) q[3];
sx q[3];
rz(-0.70374792) q[3];
sx q[3];
rz(3.0643413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7909214) q[2];
sx q[2];
rz(-1.8427883) q[2];
sx q[2];
rz(2.217963) q[2];
rz(-2.1212497) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(-3.0237107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0608805) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(-1.6545779) q[0];
rz(-3.0912073) q[1];
sx q[1];
rz(-0.81704187) q[1];
sx q[1];
rz(0.64687669) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6481224) q[0];
sx q[0];
rz(-2.1645344) q[0];
sx q[0];
rz(2.6208682) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66604659) q[2];
sx q[2];
rz(-1.9589309) q[2];
sx q[2];
rz(1.1626409) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9641831) q[1];
sx q[1];
rz(-0.46176592) q[1];
sx q[1];
rz(-0.48514556) q[1];
rz(-0.96109747) q[3];
sx q[3];
rz(-0.25843474) q[3];
sx q[3];
rz(1.0040545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50463027) q[2];
sx q[2];
rz(-1.5865822) q[2];
sx q[2];
rz(0.75616765) q[2];
rz(1.470083) q[3];
sx q[3];
rz(-0.78592891) q[3];
sx q[3];
rz(0.80529958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.037755448) q[0];
sx q[0];
rz(-1.2498195) q[0];
sx q[0];
rz(-1.1081498) q[0];
rz(0.26421079) q[1];
sx q[1];
rz(-1.6725531) q[1];
sx q[1];
rz(-0.60233751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3964993) q[0];
sx q[0];
rz(-1.4317997) q[0];
sx q[0];
rz(1.2953912) q[0];
rz(-2.9649023) q[2];
sx q[2];
rz(-1.4193221) q[2];
sx q[2];
rz(0.02515153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.332676) q[1];
sx q[1];
rz(-2.4530468) q[1];
sx q[1];
rz(-1.4946412) q[1];
rz(-0.29903166) q[3];
sx q[3];
rz(-2.4686738) q[3];
sx q[3];
rz(0.80489075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5084874) q[2];
sx q[2];
rz(-2.5999531) q[2];
sx q[2];
rz(2.2701021) q[2];
rz(-1.2095215) q[3];
sx q[3];
rz(-2.8612374) q[3];
sx q[3];
rz(0.59396321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7021983) q[0];
sx q[0];
rz(-1.7560503) q[0];
sx q[0];
rz(-0.51881292) q[0];
rz(0.58204542) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(1.9459361) q[2];
sx q[2];
rz(-2.2715501) q[2];
sx q[2];
rz(3.0987433) q[2];
rz(0.20584917) q[3];
sx q[3];
rz(-0.69910739) q[3];
sx q[3];
rz(-1.7144937) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
