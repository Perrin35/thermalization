OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.03117938) q[0];
sx q[0];
rz(-0.92209446) q[0];
sx q[0];
rz(2.3715012) q[0];
rz(0.71647477) q[1];
sx q[1];
rz(4.0840277) q[1];
sx q[1];
rz(10.066636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72015136) q[0];
sx q[0];
rz(-0.23166616) q[0];
sx q[0];
rz(2.6257444) q[0];
x q[1];
rz(-0.0014512295) q[2];
sx q[2];
rz(-1.5696758) q[2];
sx q[2];
rz(1.4933153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0986276) q[1];
sx q[1];
rz(-2.7471099) q[1];
sx q[1];
rz(-2.7261655) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3083616) q[3];
sx q[3];
rz(-1.1166443) q[3];
sx q[3];
rz(1.3564701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98051071) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(0.80992997) q[2];
rz(3.1071281) q[3];
sx q[3];
rz(-2.4813215) q[3];
sx q[3];
rz(3.0380429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.506839) q[0];
sx q[0];
rz(-0.50951183) q[0];
sx q[0];
rz(-0.65938812) q[0];
rz(-1.4393282) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(0.48092458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.07671) q[0];
sx q[0];
rz(-2.0079029) q[0];
sx q[0];
rz(-2.9124898) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4056817) q[2];
sx q[2];
rz(-2.1440268) q[2];
sx q[2];
rz(0.57511759) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3673076) q[1];
sx q[1];
rz(-1.3626422) q[1];
sx q[1];
rz(-3.1104196) q[1];
x q[2];
rz(-1.102785) q[3];
sx q[3];
rz(-1.0429405) q[3];
sx q[3];
rz(-0.46609391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7646358) q[2];
sx q[2];
rz(-2.3048293) q[2];
sx q[2];
rz(-0.62977201) q[2];
rz(-1.1791641) q[3];
sx q[3];
rz(-2.4383014) q[3];
sx q[3];
rz(1.111697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-0.010882219) q[0];
sx q[0];
rz(-0.96845281) q[0];
rz(0.14006607) q[1];
sx q[1];
rz(-1.7882971) q[1];
sx q[1];
rz(2.5586939) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2067277) q[0];
sx q[0];
rz(-1.4040274) q[0];
sx q[0];
rz(-1.5251446) q[0];
rz(-2.9045461) q[2];
sx q[2];
rz(-2.2893527) q[2];
sx q[2];
rz(-1.7463373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0576833) q[1];
sx q[1];
rz(-1.0094125) q[1];
sx q[1];
rz(2.3777005) q[1];
x q[2];
rz(1.2517334) q[3];
sx q[3];
rz(-1.8382065) q[3];
sx q[3];
rz(-1.2280994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6110903) q[2];
sx q[2];
rz(-3.0641596) q[2];
sx q[2];
rz(-0.21406847) q[2];
rz(0.30238447) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(-0.42588699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9532303) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(-1.0444214) q[0];
rz(-2.2638679) q[1];
sx q[1];
rz(-1.5889771) q[1];
sx q[1];
rz(-2.9728319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2917738) q[0];
sx q[0];
rz(-0.55233228) q[0];
sx q[0];
rz(2.0122347) q[0];
rz(-1.020476) q[2];
sx q[2];
rz(-2.762971) q[2];
sx q[2];
rz(-2.052321) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9738439) q[1];
sx q[1];
rz(-0.97626057) q[1];
sx q[1];
rz(-1.3092625) q[1];
rz(0.041958001) q[3];
sx q[3];
rz(-2.2760512) q[3];
sx q[3];
rz(0.77159568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7793444) q[2];
sx q[2];
rz(-1.2835953) q[2];
sx q[2];
rz(-1.0603504) q[2];
rz(2.849071) q[3];
sx q[3];
rz(-2.4157603) q[3];
sx q[3];
rz(3.0574851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58920687) q[0];
sx q[0];
rz(-2.4942106) q[0];
sx q[0];
rz(0.86828434) q[0];
rz(0.6262511) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.3583604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3247899) q[0];
sx q[0];
rz(-0.30203095) q[0];
sx q[0];
rz(0.30257757) q[0];
x q[1];
rz(0.38154885) q[2];
sx q[2];
rz(-1.7722436) q[2];
sx q[2];
rz(-0.82009456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.080278668) q[1];
sx q[1];
rz(-1.3484203) q[1];
sx q[1];
rz(0.23134065) q[1];
x q[2];
rz(1.129398) q[3];
sx q[3];
rz(-2.1611161) q[3];
sx q[3];
rz(-2.2143827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9354349) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(-0.52331501) q[2];
rz(-2.7715136) q[3];
sx q[3];
rz(-0.78032929) q[3];
sx q[3];
rz(-1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
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
rz(3.0321781) q[0];
sx q[0];
rz(-0.21571708) q[0];
sx q[0];
rz(-2.7874462) q[0];
rz(2.1996563) q[1];
sx q[1];
rz(-1.4553921) q[1];
sx q[1];
rz(-1.8249493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8982074) q[0];
sx q[0];
rz(-1.6750209) q[0];
sx q[0];
rz(2.6263155) q[0];
rz(-pi) q[1];
rz(-0.022158547) q[2];
sx q[2];
rz(-2.5922734) q[2];
sx q[2];
rz(-1.6946799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.020654708) q[1];
sx q[1];
rz(-1.5096501) q[1];
sx q[1];
rz(-0.17423363) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7735071) q[3];
sx q[3];
rz(-0.67782611) q[3];
sx q[3];
rz(0.2673291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8159304) q[2];
sx q[2];
rz(-0.38620913) q[2];
sx q[2];
rz(-2.8032934) q[2];
rz(-2.6650688) q[3];
sx q[3];
rz(-0.75348133) q[3];
sx q[3];
rz(-0.34059718) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4395831) q[0];
sx q[0];
rz(-1.5832573) q[0];
sx q[0];
rz(2.6038792) q[0];
rz(-1.733755) q[1];
sx q[1];
rz(-0.45999637) q[1];
sx q[1];
rz(-0.62526155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475554) q[0];
sx q[0];
rz(-0.68935822) q[0];
sx q[0];
rz(-0.58527845) q[0];
x q[1];
rz(-0.53886063) q[2];
sx q[2];
rz(-1.6340573) q[2];
sx q[2];
rz(2.6841595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.96237732) q[1];
sx q[1];
rz(-1.7853702) q[1];
sx q[1];
rz(-2.5291165) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37276778) q[3];
sx q[3];
rz(-0.81171747) q[3];
sx q[3];
rz(-1.852664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9739146) q[2];
sx q[2];
rz(-0.46678552) q[2];
sx q[2];
rz(1.7612877) q[2];
rz(0.4365094) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(-2.4280587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9647144) q[0];
sx q[0];
rz(-2.5202993) q[0];
sx q[0];
rz(-0.51625133) q[0];
rz(-2.3618354) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(-2.0947184) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16653331) q[0];
sx q[0];
rz(-1.2724845) q[0];
sx q[0];
rz(-1.7214283) q[0];
rz(-pi) q[1];
rz(0.70430906) q[2];
sx q[2];
rz(-2.4302135) q[2];
sx q[2];
rz(-0.34397438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2566764) q[1];
sx q[1];
rz(-3.0312928) q[1];
sx q[1];
rz(1.0156471) q[1];
x q[2];
rz(-0.47886301) q[3];
sx q[3];
rz(-1.9267462) q[3];
sx q[3];
rz(-2.1321116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9376935) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(-2.8474478) q[3];
sx q[3];
rz(-2.011994) q[3];
sx q[3];
rz(-0.2778151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99943632) q[0];
sx q[0];
rz(-2.9111828) q[0];
sx q[0];
rz(-2.9545422) q[0];
rz(2.710178) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(1.6492708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6607912) q[0];
sx q[0];
rz(-1.4993164) q[0];
sx q[0];
rz(-3.0707703) q[0];
x q[1];
rz(1.236294) q[2];
sx q[2];
rz(-1.3906456) q[2];
sx q[2];
rz(-2.1366773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3826661) q[1];
sx q[1];
rz(-1.9793643) q[1];
sx q[1];
rz(-2.9620785) q[1];
x q[2];
rz(1.0338551) q[3];
sx q[3];
rz(-2.9599422) q[3];
sx q[3];
rz(2.8665115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46818587) q[2];
sx q[2];
rz(-0.91172051) q[2];
sx q[2];
rz(2.1569596) q[2];
rz(2.6268688) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(1.908186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.89727) q[0];
sx q[0];
rz(-1.4800625) q[0];
sx q[0];
rz(0.75862128) q[0];
rz(1.1933391) q[1];
sx q[1];
rz(-1.9494282) q[1];
sx q[1];
rz(-1.6903445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1275996) q[0];
sx q[0];
rz(-1.3422478) q[0];
sx q[0];
rz(2.9319113) q[0];
x q[1];
rz(-0.042165857) q[2];
sx q[2];
rz(-1.3683967) q[2];
sx q[2];
rz(0.037994904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32827863) q[1];
sx q[1];
rz(-0.61610389) q[1];
sx q[1];
rz(-1.6170943) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91345738) q[3];
sx q[3];
rz(-11/(7*pi)) q[3];
sx q[3];
rz(-0.72266912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54185581) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(-2.3403781) q[2];
rz(0.93252212) q[3];
sx q[3];
rz(-1.2152117) q[3];
sx q[3];
rz(2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5781317) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(0.32147944) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(0.16531113) q[2];
sx q[2];
rz(-1.4480235) q[2];
sx q[2];
rz(1.5245246) q[2];
rz(-0.88820171) q[3];
sx q[3];
rz(-2.785688) q[3];
sx q[3];
rz(1.9733081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
