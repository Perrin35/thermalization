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
rz(1.9873729) q[0];
sx q[0];
rz(-1.6183102) q[0];
sx q[0];
rz(-0.14259882) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(3.7945336) q[1];
sx q[1];
rz(8.3795587) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3616989) q[0];
sx q[0];
rz(-2.4068659) q[0];
sx q[0];
rz(2.421342) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93306834) q[2];
sx q[2];
rz(-1.5412113) q[2];
sx q[2];
rz(-0.78630182) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75946158) q[1];
sx q[1];
rz(-0.18759094) q[1];
sx q[1];
rz(-1.0008903) q[1];
x q[2];
rz(1.0020489) q[3];
sx q[3];
rz(-1.4706597) q[3];
sx q[3];
rz(1.1796234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1656701) q[2];
sx q[2];
rz(-1.0079577) q[2];
sx q[2];
rz(-0.21318501) q[2];
rz(-0.91607696) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857472) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(-0.47072738) q[0];
rz(2.0384906) q[1];
sx q[1];
rz(-2.6328937) q[1];
sx q[1];
rz(0.57977605) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3382028) q[0];
sx q[0];
rz(-1.0472327) q[0];
sx q[0];
rz(-1.6399327) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2389555) q[2];
sx q[2];
rz(-2.08925) q[2];
sx q[2];
rz(-1.7331074) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.55706876) q[1];
sx q[1];
rz(-1.1864206) q[1];
sx q[1];
rz(-2.1904883) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9409431) q[3];
sx q[3];
rz(-1.3251628) q[3];
sx q[3];
rz(0.64526886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(0.095005438) q[2];
rz(-2.5659918) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(-1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59548241) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(-0.27221671) q[0];
rz(-0.1046003) q[1];
sx q[1];
rz(-1.0403386) q[1];
sx q[1];
rz(1.6786172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9538077) q[0];
sx q[0];
rz(-1.6208795) q[0];
sx q[0];
rz(0.84651504) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3562695) q[2];
sx q[2];
rz(-1.5253272) q[2];
sx q[2];
rz(-0.77982219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8045878) q[1];
sx q[1];
rz(-2.2175601) q[1];
sx q[1];
rz(1.2156562) q[1];
rz(-1.5103064) q[3];
sx q[3];
rz(-1.20245) q[3];
sx q[3];
rz(-2.9591536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0649123) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(-0.40985516) q[2];
rz(3.0900433) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(2.4079017) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9728397) q[0];
sx q[0];
rz(-2.3641455) q[0];
sx q[0];
rz(-0.69946104) q[0];
rz(-2.2109924) q[1];
sx q[1];
rz(-1.7761296) q[1];
sx q[1];
rz(-2.0420989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.612815) q[0];
sx q[0];
rz(-1.8474425) q[0];
sx q[0];
rz(2.1601281) q[0];
rz(-2.900057) q[2];
sx q[2];
rz(-2.1061828) q[2];
sx q[2];
rz(-2.8098742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9605254) q[1];
sx q[1];
rz(-1.8904422) q[1];
sx q[1];
rz(-1.8962217) q[1];
rz(-1.9774262) q[3];
sx q[3];
rz(-0.96897954) q[3];
sx q[3];
rz(-0.072657771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76116556) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(2.579651) q[3];
sx q[3];
rz(-2.1886531) q[3];
sx q[3];
rz(-0.57653069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1171653) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(1.0166136) q[0];
rz(0.92574614) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(3.064916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0177855) q[0];
sx q[0];
rz(-1.9636256) q[0];
sx q[0];
rz(0.17946243) q[0];
rz(-0.20598866) q[2];
sx q[2];
rz(-1.811337) q[2];
sx q[2];
rz(2.5702916) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5591084) q[1];
sx q[1];
rz(-0.6000207) q[1];
sx q[1];
rz(-1.1901072) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2604976) q[3];
sx q[3];
rz(-2.3633133) q[3];
sx q[3];
rz(0.33184127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3735247) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(-0.12316556) q[2];
rz(0.67433107) q[3];
sx q[3];
rz(-2.6682523) q[3];
sx q[3];
rz(2.3136852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8318361) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-2.6605666) q[0];
rz(-2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(-0.13490881) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19770007) q[0];
sx q[0];
rz(-2.0257844) q[0];
sx q[0];
rz(-2.3030192) q[0];
rz(2.9685855) q[2];
sx q[2];
rz(-2.5998625) q[2];
sx q[2];
rz(-2.2587905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1926769) q[1];
sx q[1];
rz(-0.8444311) q[1];
sx q[1];
rz(2.7706258) q[1];
rz(2.6916299) q[3];
sx q[3];
rz(-1.6358161) q[3];
sx q[3];
rz(1.3715594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2065108) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(-1.3235462) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-2.2836253) q[3];
sx q[3];
rz(2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(-2.4830699) q[0];
rz(-1.5497442) q[1];
sx q[1];
rz(-1.0127944) q[1];
sx q[1];
rz(2.6351567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69909912) q[0];
sx q[0];
rz(-0.37537071) q[0];
sx q[0];
rz(0.18113776) q[0];
rz(1.6295048) q[2];
sx q[2];
rz(-2.65911) q[2];
sx q[2];
rz(2.5168632) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68489198) q[1];
sx q[1];
rz(-1.646767) q[1];
sx q[1];
rz(-0.09346813) q[1];
rz(0.063252216) q[3];
sx q[3];
rz(-0.58025415) q[3];
sx q[3];
rz(1.751251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8062313) q[2];
sx q[2];
rz(-2.3956617) q[2];
sx q[2];
rz(1.2696421) q[2];
rz(-1.7757724) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(-2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7954471) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(2.7224097) q[0];
rz(3.0508793) q[1];
sx q[1];
rz(-0.9181298) q[1];
sx q[1];
rz(2.1144287) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61495435) q[0];
sx q[0];
rz(-2.0169584) q[0];
sx q[0];
rz(-0.63132186) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1187657) q[2];
sx q[2];
rz(-1.5628624) q[2];
sx q[2];
rz(1.4102907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60909437) q[1];
sx q[1];
rz(-0.23819345) q[1];
sx q[1];
rz(-1.1042117) q[1];
rz(-pi) q[2];
rz(-2.2563491) q[3];
sx q[3];
rz(-1.5285057) q[3];
sx q[3];
rz(1.919618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1586228) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(-0.99009222) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7541499) q[0];
sx q[0];
rz(-0.70873547) q[0];
sx q[0];
rz(-2.5842174) q[0];
rz(2.7793461) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-0.16709669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8734735) q[0];
sx q[0];
rz(-0.41129204) q[0];
sx q[0];
rz(2.4231572) q[0];
rz(-pi) q[1];
rz(-1.3873151) q[2];
sx q[2];
rz(-2.263732) q[2];
sx q[2];
rz(0.82833457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36896781) q[1];
sx q[1];
rz(-2.4030622) q[1];
sx q[1];
rz(-0.46486295) q[1];
rz(-0.88730335) q[3];
sx q[3];
rz(-1.512882) q[3];
sx q[3];
rz(-1.37552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(-2.7131405) q[2];
rz(0.7306478) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44883248) q[0];
sx q[0];
rz(-1.6535783) q[0];
sx q[0];
rz(-1.1195419) q[0];
rz(0.66850942) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(0.25984919) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4762806) q[0];
sx q[0];
rz(-0.57054936) q[0];
sx q[0];
rz(-1.0744988) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.75801) q[2];
sx q[2];
rz(-2.0077939) q[2];
sx q[2];
rz(1.1326157) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4230792) q[1];
sx q[1];
rz(-1.2260409) q[1];
sx q[1];
rz(0.43889795) q[1];
x q[2];
rz(-2.2606528) q[3];
sx q[3];
rz(-1.4487113) q[3];
sx q[3];
rz(-1.087904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6707637) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(1.0864351) q[2];
rz(2.6492665) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(2.4194748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(-1.8863574) q[1];
sx q[1];
rz(-1.8240758) q[1];
sx q[1];
rz(0.12771894) q[1];
rz(0.095070953) q[2];
sx q[2];
rz(-1.9683206) q[2];
sx q[2];
rz(-1.5299601) q[2];
rz(-2.2333748) q[3];
sx q[3];
rz(-1.5208018) q[3];
sx q[3];
rz(3.0326299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
