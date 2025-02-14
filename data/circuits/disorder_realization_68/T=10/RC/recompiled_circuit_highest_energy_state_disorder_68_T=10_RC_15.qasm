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
rz(0.23213586) q[0];
sx q[0];
rz(-0.27096662) q[0];
sx q[0];
rz(-1.1722857) q[0];
rz(-1.6074033) q[1];
sx q[1];
rz(-1.489403) q[1];
sx q[1];
rz(-0.54816562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70943862) q[0];
sx q[0];
rz(-2.84542) q[0];
sx q[0];
rz(-1.1838566) q[0];
rz(-pi) q[1];
rz(-0.51082533) q[2];
sx q[2];
rz(-2.0693972) q[2];
sx q[2];
rz(-0.076534903) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2506574) q[1];
sx q[1];
rz(-1.0708628) q[1];
sx q[1];
rz(1.2229162) q[1];
x q[2];
rz(1.4413766) q[3];
sx q[3];
rz(-2.6886422) q[3];
sx q[3];
rz(-0.85080409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24070172) q[2];
sx q[2];
rz(-1.9343932) q[2];
sx q[2];
rz(1.5473676) q[2];
rz(2.2089925) q[3];
sx q[3];
rz(-0.91351944) q[3];
sx q[3];
rz(-2.3368321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43657434) q[0];
sx q[0];
rz(-2.5387634) q[0];
sx q[0];
rz(0.21827179) q[0];
rz(-1.6328579) q[1];
sx q[1];
rz(-0.873133) q[1];
sx q[1];
rz(1.1223209) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22161417) q[0];
sx q[0];
rz(-1.99619) q[0];
sx q[0];
rz(-0.90817389) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.400945) q[2];
sx q[2];
rz(-1.4221995) q[2];
sx q[2];
rz(2.1987178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3914267) q[1];
sx q[1];
rz(-1.1216333) q[1];
sx q[1];
rz(2.5264747) q[1];
rz(-1.8863986) q[3];
sx q[3];
rz(-1.8132993) q[3];
sx q[3];
rz(0.24569233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0108769) q[2];
sx q[2];
rz(-1.9515832) q[2];
sx q[2];
rz(-2.7575764) q[2];
rz(-2.5859517) q[3];
sx q[3];
rz(-0.58755392) q[3];
sx q[3];
rz(0.89239341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4246849) q[0];
sx q[0];
rz(-0.45775828) q[0];
sx q[0];
rz(2.6440788) q[0];
rz(2.6566907) q[1];
sx q[1];
rz(-1.5496016) q[1];
sx q[1];
rz(-2.7000694) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0245666) q[0];
sx q[0];
rz(-0.68156201) q[0];
sx q[0];
rz(1.3016537) q[0];
x q[1];
rz(1.9633358) q[2];
sx q[2];
rz(-0.52670331) q[2];
sx q[2];
rz(1.3046427) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42558266) q[1];
sx q[1];
rz(-2.0507318) q[1];
sx q[1];
rz(2.369057) q[1];
rz(0.27350815) q[3];
sx q[3];
rz(-1.8778749) q[3];
sx q[3];
rz(2.8321137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32450822) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(0.95477742) q[2];
rz(0.020141007) q[3];
sx q[3];
rz(-2.3432799) q[3];
sx q[3];
rz(2.8153343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666053) q[0];
sx q[0];
rz(-2.6191481) q[0];
sx q[0];
rz(-3.1368384) q[0];
rz(0.44113723) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(1.7316679) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49214464) q[0];
sx q[0];
rz(-2.7810367) q[0];
sx q[0];
rz(1.9515522) q[0];
x q[1];
rz(-0.84216046) q[2];
sx q[2];
rz(-1.3628863) q[2];
sx q[2];
rz(0.72170695) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.04206229) q[1];
sx q[1];
rz(-2.5914609) q[1];
sx q[1];
rz(-2.0127556) q[1];
x q[2];
rz(-3.1143777) q[3];
sx q[3];
rz(-1.3821954) q[3];
sx q[3];
rz(-1.0116497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8498174) q[2];
sx q[2];
rz(-1.8287649) q[2];
sx q[2];
rz(0.69668359) q[2];
rz(-1.5289395) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(-2.4618885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.9686389) q[0];
sx q[0];
rz(-2.8346859) q[0];
sx q[0];
rz(-0.90721834) q[0];
rz(0.13547678) q[1];
sx q[1];
rz(-1.8699346) q[1];
sx q[1];
rz(2.4438593) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0402148) q[0];
sx q[0];
rz(-1.7003978) q[0];
sx q[0];
rz(1.1762397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7656828) q[2];
sx q[2];
rz(-1.5729674) q[2];
sx q[2];
rz(-1.2819829) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2574324) q[1];
sx q[1];
rz(-1.2425649) q[1];
sx q[1];
rz(-1.753367) q[1];
rz(-1.5248564) q[3];
sx q[3];
rz(-1.5410863) q[3];
sx q[3];
rz(1.3256095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9542784) q[2];
sx q[2];
rz(-0.66437393) q[2];
sx q[2];
rz(0.3453671) q[2];
rz(0.49515381) q[3];
sx q[3];
rz(-2.5457355) q[3];
sx q[3];
rz(-0.12721795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0536163) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(1.1548868) q[0];
rz(0.79479533) q[1];
sx q[1];
rz(-1.844901) q[1];
sx q[1];
rz(0.48387873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34209868) q[0];
sx q[0];
rz(-2.3364107) q[0];
sx q[0];
rz(2.1495887) q[0];
rz(1.2829928) q[2];
sx q[2];
rz(-0.51454267) q[2];
sx q[2];
rz(-0.3798011) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72493785) q[1];
sx q[1];
rz(-2.0339801) q[1];
sx q[1];
rz(2.928726) q[1];
rz(-pi) q[2];
rz(-1.9779786) q[3];
sx q[3];
rz(-1.0308427) q[3];
sx q[3];
rz(-2.0288717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7530219) q[2];
sx q[2];
rz(-1.2349671) q[2];
sx q[2];
rz(-1.1080326) q[2];
rz(-1.5293416) q[3];
sx q[3];
rz(-3.022091) q[3];
sx q[3];
rz(-0.39335462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18688467) q[0];
sx q[0];
rz(-1.0857546) q[0];
sx q[0];
rz(2.4672274) q[0];
rz(-2.2652594) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(-3.1092627) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040213765) q[0];
sx q[0];
rz(-1.4353936) q[0];
sx q[0];
rz(0.12286326) q[0];
x q[1];
rz(1.7252847) q[2];
sx q[2];
rz(-1.3488028) q[2];
sx q[2];
rz(-2.7374008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8994042) q[1];
sx q[1];
rz(-2.5061878) q[1];
sx q[1];
rz(-0.59235488) q[1];
rz(0.97801925) q[3];
sx q[3];
rz(-1.4109352) q[3];
sx q[3];
rz(-0.39409742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.18675599) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(2.1046861) q[2];
rz(-1.9131276) q[3];
sx q[3];
rz(-2.2322673) q[3];
sx q[3];
rz(-0.19074805) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73795885) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(-0.24895689) q[0];
rz(2.9320993) q[1];
sx q[1];
rz(-1.8003576) q[1];
sx q[1];
rz(0.6627717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4512206) q[0];
sx q[0];
rz(-1.5826384) q[0];
sx q[0];
rz(1.5607587) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0515819) q[2];
sx q[2];
rz(-2.2741246) q[2];
sx q[2];
rz(2.472773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6697093) q[1];
sx q[1];
rz(-0.85195573) q[1];
sx q[1];
rz(2.813176) q[1];
rz(-pi) q[2];
rz(1.2609188) q[3];
sx q[3];
rz(-2.204748) q[3];
sx q[3];
rz(0.11989633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2304307) q[2];
sx q[2];
rz(-1.5713567) q[2];
sx q[2];
rz(2.7479808) q[2];
rz(2.7502934) q[3];
sx q[3];
rz(-2.6063882) q[3];
sx q[3];
rz(2.3894501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0126208) q[0];
sx q[0];
rz(-0.18652815) q[0];
sx q[0];
rz(0.57998002) q[0];
rz(2.773556) q[1];
sx q[1];
rz(-0.8152222) q[1];
sx q[1];
rz(3.0267402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2942685) q[0];
sx q[0];
rz(-2.2912187) q[0];
sx q[0];
rz(-0.17981932) q[0];
x q[1];
rz(0.1505974) q[2];
sx q[2];
rz(-2.0080749) q[2];
sx q[2];
rz(-2.8018746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65832135) q[1];
sx q[1];
rz(-1.8206925) q[1];
sx q[1];
rz(2.8037594) q[1];
rz(-pi) q[2];
rz(-0.40267085) q[3];
sx q[3];
rz(-0.91710233) q[3];
sx q[3];
rz(-1.634906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7241235) q[2];
sx q[2];
rz(-2.8026411) q[2];
sx q[2];
rz(0.69619703) q[2];
rz(-1.1160858) q[3];
sx q[3];
rz(-0.79056549) q[3];
sx q[3];
rz(0.74151403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29257947) q[0];
sx q[0];
rz(-2.1722023) q[0];
sx q[0];
rz(2.868929) q[0];
rz(0.69372454) q[1];
sx q[1];
rz(-2.0246181) q[1];
sx q[1];
rz(-2.5781217) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7359223) q[0];
sx q[0];
rz(-1.5338681) q[0];
sx q[0];
rz(-1.4775299) q[0];
x q[1];
rz(0.23691249) q[2];
sx q[2];
rz(-1.713551) q[2];
sx q[2];
rz(2.3430489) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8565791) q[1];
sx q[1];
rz(-1.6520629) q[1];
sx q[1];
rz(1.8829499) q[1];
rz(-pi) q[2];
rz(-0.93592398) q[3];
sx q[3];
rz(-2.1800123) q[3];
sx q[3];
rz(3.0425231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1571265) q[2];
sx q[2];
rz(-2.0968585) q[2];
sx q[2];
rz(-2.5926479) q[2];
rz(0.99556154) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(-0.63001776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7231049) q[0];
sx q[0];
rz(-2.1514819) q[0];
sx q[0];
rz(2.6958382) q[0];
rz(0.86722974) q[1];
sx q[1];
rz(-2.0384616) q[1];
sx q[1];
rz(-1.9834317) q[1];
rz(-0.42648496) q[2];
sx q[2];
rz(-0.51649649) q[2];
sx q[2];
rz(-2.6994097) q[2];
rz(-0.23602939) q[3];
sx q[3];
rz(-0.56047816) q[3];
sx q[3];
rz(-2.6666835) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
