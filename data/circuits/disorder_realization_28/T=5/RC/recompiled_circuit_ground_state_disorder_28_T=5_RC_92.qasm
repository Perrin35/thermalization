OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(-1.5225141) q[0];
sx q[0];
rz(-0.063152753) q[0];
rz(1.8237279) q[1];
sx q[1];
rz(-0.39931077) q[1];
sx q[1];
rz(2.3656486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.522811) q[0];
sx q[0];
rz(-1.543985) q[0];
sx q[0];
rz(-0.045700886) q[0];
rz(-2.2065483) q[2];
sx q[2];
rz(-1.9823977) q[2];
sx q[2];
rz(-1.7807478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9358238) q[1];
sx q[1];
rz(-0.73038126) q[1];
sx q[1];
rz(1.303506) q[1];
rz(-1.3315866) q[3];
sx q[3];
rz(-2.0612217) q[3];
sx q[3];
rz(1.8338721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6905288) q[2];
sx q[2];
rz(-1.7454001) q[2];
sx q[2];
rz(-0.53831354) q[2];
rz(2.8484143) q[3];
sx q[3];
rz(-1.5436951) q[3];
sx q[3];
rz(1.258491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20741589) q[0];
sx q[0];
rz(-0.84686142) q[0];
sx q[0];
rz(2.9127981) q[0];
rz(-3.0711807) q[1];
sx q[1];
rz(-1.2733302) q[1];
sx q[1];
rz(1.061903) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.966514) q[0];
sx q[0];
rz(-2.0850777) q[0];
sx q[0];
rz(1.8755336) q[0];
x q[1];
rz(-2.0309762) q[2];
sx q[2];
rz(-1.1374344) q[2];
sx q[2];
rz(3.1040807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50256601) q[1];
sx q[1];
rz(-0.4430534) q[1];
sx q[1];
rz(2.4094617) q[1];
rz(-pi) q[2];
rz(-2.6999989) q[3];
sx q[3];
rz(-0.46735763) q[3];
sx q[3];
rz(0.87477126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.215302) q[2];
sx q[2];
rz(-2.6965202) q[2];
sx q[2];
rz(3.0253547) q[2];
rz(0.91533533) q[3];
sx q[3];
rz(-1.6719619) q[3];
sx q[3];
rz(-2.5168929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4794469) q[0];
sx q[0];
rz(-1.5818469) q[0];
sx q[0];
rz(0.33552718) q[0];
rz(-1.7299995) q[1];
sx q[1];
rz(-0.84452191) q[1];
sx q[1];
rz(-1.9427049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480961) q[0];
sx q[0];
rz(-2.5272547) q[0];
sx q[0];
rz(1.0636281) q[0];
rz(-pi) q[1];
rz(1.185965) q[2];
sx q[2];
rz(-1.0116546) q[2];
sx q[2];
rz(-2.2988479) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2553067) q[1];
sx q[1];
rz(-0.77091187) q[1];
sx q[1];
rz(-0.16718276) q[1];
rz(2.8923404) q[3];
sx q[3];
rz(-2.0468759) q[3];
sx q[3];
rz(0.92633807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18232839) q[2];
sx q[2];
rz(-0.79572833) q[2];
sx q[2];
rz(1.1129334) q[2];
rz(0.5111323) q[3];
sx q[3];
rz(-1.7685578) q[3];
sx q[3];
rz(-2.6326877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48753259) q[0];
sx q[0];
rz(-3.0656116) q[0];
sx q[0];
rz(2.7677166) q[0];
rz(-1.182425) q[1];
sx q[1];
rz(-1.1610718) q[1];
sx q[1];
rz(0.16071308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492225) q[0];
sx q[0];
rz(-2.9806031) q[0];
sx q[0];
rz(0.54246728) q[0];
x q[1];
rz(-0.79419996) q[2];
sx q[2];
rz(-2.080285) q[2];
sx q[2];
rz(1.2603261) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.070878167) q[1];
sx q[1];
rz(-1.1859815) q[1];
sx q[1];
rz(-1.0139731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75467189) q[3];
sx q[3];
rz(-0.42413482) q[3];
sx q[3];
rz(-1.6744259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.60761991) q[2];
sx q[2];
rz(-1.1090358) q[2];
sx q[2];
rz(2.4268761) q[2];
rz(2.886582) q[3];
sx q[3];
rz(-1.3304354) q[3];
sx q[3];
rz(2.3431006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0890927) q[0];
sx q[0];
rz(-2.8369501) q[0];
sx q[0];
rz(-0.76328817) q[0];
rz(0.25453645) q[1];
sx q[1];
rz(-1.3724047) q[1];
sx q[1];
rz(-1.6932142) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096222046) q[0];
sx q[0];
rz(-1.9217446) q[0];
sx q[0];
rz(-0.36349067) q[0];
x q[1];
rz(2.9853742) q[2];
sx q[2];
rz(-0.82619452) q[2];
sx q[2];
rz(2.790839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9824286) q[1];
sx q[1];
rz(-0.64627534) q[1];
sx q[1];
rz(-0.29781945) q[1];
x q[2];
rz(-0.45381399) q[3];
sx q[3];
rz(-0.76440135) q[3];
sx q[3];
rz(1.3168471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.280507) q[2];
sx q[2];
rz(-1.0887159) q[2];
sx q[2];
rz(-2.1121934) q[2];
rz(-1.5329125) q[3];
sx q[3];
rz(-2.2423988) q[3];
sx q[3];
rz(-2.5077584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3935881) q[0];
sx q[0];
rz(-0.91359502) q[0];
sx q[0];
rz(-2.9314801) q[0];
rz(-2.4402319) q[1];
sx q[1];
rz(-1.3664093) q[1];
sx q[1];
rz(2.0393541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0118679) q[0];
sx q[0];
rz(-1.1107011) q[0];
sx q[0];
rz(0.97614093) q[0];
x q[1];
rz(-0.66221018) q[2];
sx q[2];
rz(-1.5726316) q[2];
sx q[2];
rz(0.83787943) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21747281) q[1];
sx q[1];
rz(-2.365781) q[1];
sx q[1];
rz(0.45507832) q[1];
rz(-0.26138283) q[3];
sx q[3];
rz(-2.5845575) q[3];
sx q[3];
rz(-1.6121395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8474951) q[2];
sx q[2];
rz(-1.558344) q[2];
sx q[2];
rz(-1.6451277) q[2];
rz(-1.3930813) q[3];
sx q[3];
rz(-2.2025509) q[3];
sx q[3];
rz(0.67162544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74712718) q[0];
sx q[0];
rz(-1.8165996) q[0];
sx q[0];
rz(-2.6649244) q[0];
rz(-1.3564302) q[1];
sx q[1];
rz(-1.8442267) q[1];
sx q[1];
rz(2.2043601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8097157) q[0];
sx q[0];
rz(-2.2599155) q[0];
sx q[0];
rz(-0.046781311) q[0];
rz(0.71800443) q[2];
sx q[2];
rz(-1.5239626) q[2];
sx q[2];
rz(0.61344922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46980935) q[1];
sx q[1];
rz(-1.9783164) q[1];
sx q[1];
rz(1.9923366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9415575) q[3];
sx q[3];
rz(-2.6851401) q[3];
sx q[3];
rz(-0.70213041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53200191) q[2];
sx q[2];
rz(-1.8014182) q[2];
sx q[2];
rz(1.0756005) q[2];
rz(0.39014751) q[3];
sx q[3];
rz(-1.8307999) q[3];
sx q[3];
rz(-2.3384317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6049062) q[0];
sx q[0];
rz(-2.0669879) q[0];
sx q[0];
rz(-0.54501504) q[0];
rz(0.99264985) q[1];
sx q[1];
rz(-2.1558709) q[1];
sx q[1];
rz(1.4535905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6777991) q[0];
sx q[0];
rz(-2.2132769) q[0];
sx q[0];
rz(0.55046659) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7033061) q[2];
sx q[2];
rz(-1.6269267) q[2];
sx q[2];
rz(2.0126337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29991789) q[1];
sx q[1];
rz(-2.3280488) q[1];
sx q[1];
rz(2.9857789) q[1];
rz(1.6176295) q[3];
sx q[3];
rz(-1.6308271) q[3];
sx q[3];
rz(1.5200652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7327026) q[2];
sx q[2];
rz(-1.7444892) q[2];
sx q[2];
rz(0.2962386) q[2];
rz(-2.564751) q[3];
sx q[3];
rz(-0.39671612) q[3];
sx q[3];
rz(-1.860994) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883564) q[0];
sx q[0];
rz(-2.3529973) q[0];
sx q[0];
rz(-2.9175135) q[0];
rz(-0.60399404) q[1];
sx q[1];
rz(-2.150034) q[1];
sx q[1];
rz(1.4858861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1224269) q[0];
sx q[0];
rz(-0.44394025) q[0];
sx q[0];
rz(0.60444085) q[0];
x q[1];
rz(-1.3732128) q[2];
sx q[2];
rz(-1.7181529) q[2];
sx q[2];
rz(1.4214398) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2961041) q[1];
sx q[1];
rz(-1.8252686) q[1];
sx q[1];
rz(0.62800447) q[1];
rz(0.4918488) q[3];
sx q[3];
rz(-1.3054732) q[3];
sx q[3];
rz(2.6387174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6622582) q[2];
sx q[2];
rz(-2.0749638) q[2];
sx q[2];
rz(0.38702854) q[2];
rz(-0.28338638) q[3];
sx q[3];
rz(-2.3895388) q[3];
sx q[3];
rz(-0.40248218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1166444) q[0];
sx q[0];
rz(-0.270917) q[0];
sx q[0];
rz(-2.0181632) q[0];
rz(1.4943538) q[1];
sx q[1];
rz(-1.6554183) q[1];
sx q[1];
rz(-2.2656238) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1946404) q[0];
sx q[0];
rz(-0.84492296) q[0];
sx q[0];
rz(1.915917) q[0];
x q[1];
rz(0.45166679) q[2];
sx q[2];
rz(-2.0658138) q[2];
sx q[2];
rz(-1.5908013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.068066464) q[1];
sx q[1];
rz(-1.3519533) q[1];
sx q[1];
rz(2.2576451) q[1];
rz(1.3920636) q[3];
sx q[3];
rz(-1.1471841) q[3];
sx q[3];
rz(2.7152747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7312077) q[2];
sx q[2];
rz(-1.7816252) q[2];
sx q[2];
rz(-3.0044921) q[2];
rz(-1.7588663) q[3];
sx q[3];
rz(-2.2189249) q[3];
sx q[3];
rz(-1.9130798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.726534) q[0];
sx q[0];
rz(-2.6714323) q[0];
sx q[0];
rz(2.5191125) q[0];
rz(-0.38901916) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(-0.71826886) q[2];
sx q[2];
rz(-2.5296938) q[2];
sx q[2];
rz(-1.1863248) q[2];
rz(2.9403654) q[3];
sx q[3];
rz(-1.22898) q[3];
sx q[3];
rz(-2.1635319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
