OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57103676) q[0];
sx q[0];
rz(-2.4693627) q[0];
sx q[0];
rz(1.6663405) q[0];
rz(0.9530468) q[1];
sx q[1];
rz(3.3068125) q[1];
sx q[1];
rz(10.693476) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8707298) q[0];
sx q[0];
rz(-2.6755973) q[0];
sx q[0];
rz(2.6129338) q[0];
rz(-0.4920636) q[2];
sx q[2];
rz(-0.97717077) q[2];
sx q[2];
rz(-0.77897859) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7172473) q[1];
sx q[1];
rz(-1.0043036) q[1];
sx q[1];
rz(2.9866108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8016863) q[3];
sx q[3];
rz(-0.82726631) q[3];
sx q[3];
rz(0.32318599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98422009) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(-0.79322195) q[2];
rz(-2.4685517) q[3];
sx q[3];
rz(-2.0561736) q[3];
sx q[3];
rz(1.2696772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9840045) q[0];
sx q[0];
rz(-0.80261153) q[0];
sx q[0];
rz(2.5322835) q[0];
rz(-2.9318103) q[1];
sx q[1];
rz(-0.79133004) q[1];
sx q[1];
rz(0.9598859) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5426838) q[0];
sx q[0];
rz(-3.050878) q[0];
sx q[0];
rz(0.51341052) q[0];
x q[1];
rz(3.0716607) q[2];
sx q[2];
rz(-1.0834104) q[2];
sx q[2];
rz(-2.5068381) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9896796) q[1];
sx q[1];
rz(-1.508058) q[1];
sx q[1];
rz(0.83338085) q[1];
rz(-pi) q[2];
rz(2.0783362) q[3];
sx q[3];
rz(-1.9572658) q[3];
sx q[3];
rz(-1.8563351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7858872) q[2];
sx q[2];
rz(-2.3023119) q[2];
sx q[2];
rz(-0.4915702) q[2];
rz(0.0056886557) q[3];
sx q[3];
rz(-1.0892884) q[3];
sx q[3];
rz(0.51386851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095057644) q[0];
sx q[0];
rz(-0.53545606) q[0];
sx q[0];
rz(0.74263483) q[0];
rz(-0.62478089) q[1];
sx q[1];
rz(-2.2014328) q[1];
sx q[1];
rz(2.0754441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0381706) q[0];
sx q[0];
rz(-2.8763412) q[0];
sx q[0];
rz(-2.1384504) q[0];
rz(2.9511388) q[2];
sx q[2];
rz(-1.8862269) q[2];
sx q[2];
rz(-2.4078232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1101493) q[1];
sx q[1];
rz(-2.717579) q[1];
sx q[1];
rz(1.8808238) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0793909) q[3];
sx q[3];
rz(-1.3771476) q[3];
sx q[3];
rz(-2.1313063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2523969) q[2];
sx q[2];
rz(-0.19813457) q[2];
sx q[2];
rz(-0.55257094) q[2];
rz(-0.50734723) q[3];
sx q[3];
rz(-2.0523968) q[3];
sx q[3];
rz(2.5721917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4813389) q[0];
sx q[0];
rz(-0.57732552) q[0];
sx q[0];
rz(-1.8950155) q[0];
rz(-1.735911) q[1];
sx q[1];
rz(-0.77189267) q[1];
sx q[1];
rz(0.71800047) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3858741) q[0];
sx q[0];
rz(-1.0654952) q[0];
sx q[0];
rz(2.2875209) q[0];
rz(-pi) q[1];
rz(1.0077072) q[2];
sx q[2];
rz(-1.3674842) q[2];
sx q[2];
rz(0.85693923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93054616) q[1];
sx q[1];
rz(-0.77096838) q[1];
sx q[1];
rz(-2.9322185) q[1];
rz(-pi) q[2];
rz(-0.15851373) q[3];
sx q[3];
rz(-1.6974546) q[3];
sx q[3];
rz(1.4617006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2672853) q[2];
sx q[2];
rz(-0.75338805) q[2];
sx q[2];
rz(0.70449746) q[2];
rz(-2.7590175) q[3];
sx q[3];
rz(-1.7035328) q[3];
sx q[3];
rz(2.2883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48069561) q[0];
sx q[0];
rz(-3.1013885) q[0];
sx q[0];
rz(0.30886343) q[0];
rz(2.1924696) q[1];
sx q[1];
rz(-2.821065) q[1];
sx q[1];
rz(-0.55108756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26924414) q[0];
sx q[0];
rz(-1.1820894) q[0];
sx q[0];
rz(1.7354911) q[0];
rz(-pi) q[1];
x q[1];
rz(0.068821235) q[2];
sx q[2];
rz(-0.91409475) q[2];
sx q[2];
rz(0.35730413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6092564) q[1];
sx q[1];
rz(-0.21363959) q[1];
sx q[1];
rz(2.3268301) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9948009) q[3];
sx q[3];
rz(-2.3596968) q[3];
sx q[3];
rz(-1.1679389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71517313) q[2];
sx q[2];
rz(-0.74843633) q[2];
sx q[2];
rz(2.5227762) q[2];
rz(-2.5185781) q[3];
sx q[3];
rz(-2.2996733) q[3];
sx q[3];
rz(0.4272517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1022559) q[0];
sx q[0];
rz(-0.66348851) q[0];
sx q[0];
rz(-0.46184194) q[0];
rz(0.49679187) q[1];
sx q[1];
rz(-1.8180314) q[1];
sx q[1];
rz(-1.7740446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0858618) q[0];
sx q[0];
rz(-1.6501969) q[0];
sx q[0];
rz(-0.88831304) q[0];
rz(-pi) q[1];
x q[1];
rz(1.073764) q[2];
sx q[2];
rz(-1.4364527) q[2];
sx q[2];
rz(1.6430935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9744579) q[1];
sx q[1];
rz(-2.8394545) q[1];
sx q[1];
rz(-1.4767417) q[1];
rz(-pi) q[2];
rz(0.9487834) q[3];
sx q[3];
rz(-1.7025196) q[3];
sx q[3];
rz(1.6831116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10609145) q[2];
sx q[2];
rz(-1.1320628) q[2];
sx q[2];
rz(2.3512225) q[2];
rz(-2.5370989) q[3];
sx q[3];
rz(-2.1166182) q[3];
sx q[3];
rz(0.42009556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6919493) q[0];
sx q[0];
rz(-2.255891) q[0];
sx q[0];
rz(-2.7272136) q[0];
rz(0.62037933) q[1];
sx q[1];
rz(-1.2216156) q[1];
sx q[1];
rz(1.5827804) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5511502) q[0];
sx q[0];
rz(-1.6550468) q[0];
sx q[0];
rz(0.1131375) q[0];
rz(-0.16838603) q[2];
sx q[2];
rz(-0.10659519) q[2];
sx q[2];
rz(-1.3063947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1852946) q[1];
sx q[1];
rz(-1.5130338) q[1];
sx q[1];
rz(0.0094780427) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9314194) q[3];
sx q[3];
rz(-0.40935959) q[3];
sx q[3];
rz(1.5081852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15535007) q[2];
sx q[2];
rz(-2.4537931) q[2];
sx q[2];
rz(-2.9637994) q[2];
rz(-2.6627461) q[3];
sx q[3];
rz(-2.0279453) q[3];
sx q[3];
rz(3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.14684045) q[0];
sx q[0];
rz(-2.1629592) q[0];
sx q[0];
rz(-0.11257182) q[0];
rz(-0.62537891) q[1];
sx q[1];
rz(-2.0140779) q[1];
sx q[1];
rz(0.88923997) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13876943) q[0];
sx q[0];
rz(-0.089931503) q[0];
sx q[0];
rz(1.0732717) q[0];
rz(-pi) q[1];
rz(2.8577789) q[2];
sx q[2];
rz(-2.8655911) q[2];
sx q[2];
rz(2.2021879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4496838) q[1];
sx q[1];
rz(-2.0912716) q[1];
sx q[1];
rz(-2.8694105) q[1];
rz(-pi) q[2];
rz(1.8440644) q[3];
sx q[3];
rz(-0.24431657) q[3];
sx q[3];
rz(2.6586146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2832977) q[2];
sx q[2];
rz(-3.0136643) q[2];
sx q[2];
rz(2.8683786) q[2];
rz(-2.9122399) q[3];
sx q[3];
rz(-1.5867686) q[3];
sx q[3];
rz(1.8730414) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6943618) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(0.91621512) q[0];
rz(0.18237309) q[1];
sx q[1];
rz(-0.20621754) q[1];
sx q[1];
rz(1.1590385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3834476) q[0];
sx q[0];
rz(-2.8383857) q[0];
sx q[0];
rz(-0.68443735) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5196475) q[2];
sx q[2];
rz(-1.6661834) q[2];
sx q[2];
rz(-2.9254828) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7024544) q[1];
sx q[1];
rz(-1.7299486) q[1];
sx q[1];
rz(-1.8887146) q[1];
x q[2];
rz(0.93709903) q[3];
sx q[3];
rz(-2.1917412) q[3];
sx q[3];
rz(1.4439195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3850022) q[2];
sx q[2];
rz(-2.341541) q[2];
sx q[2];
rz(2.8263367) q[2];
rz(2.5473525) q[3];
sx q[3];
rz(-0.88163328) q[3];
sx q[3];
rz(-2.5074904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2600128) q[0];
sx q[0];
rz(-0.6686815) q[0];
sx q[0];
rz(2.9823533) q[0];
rz(-1.0881933) q[1];
sx q[1];
rz(-2.2290778) q[1];
sx q[1];
rz(2.870627) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9342459) q[0];
sx q[0];
rz(-1.1109612) q[0];
sx q[0];
rz(2.8608972) q[0];
rz(-pi) q[1];
rz(3.0130544) q[2];
sx q[2];
rz(-0.67691278) q[2];
sx q[2];
rz(0.9062137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3828363) q[1];
sx q[1];
rz(-0.74302948) q[1];
sx q[1];
rz(-0.50937517) q[1];
rz(-pi) q[2];
rz(1.9064398) q[3];
sx q[3];
rz(-2.0291532) q[3];
sx q[3];
rz(-0.021307746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8012041) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(2.8188952) q[2];
rz(0.05803756) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(-2.8431852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651481) q[0];
sx q[0];
rz(-1.6315176) q[0];
sx q[0];
rz(2.3254707) q[0];
rz(0.25794087) q[1];
sx q[1];
rz(-2.0057269) q[1];
sx q[1];
rz(-1.534091) q[1];
rz(-0.77946812) q[2];
sx q[2];
rz(-2.7457954) q[2];
sx q[2];
rz(2.376062) q[2];
rz(0.54775379) q[3];
sx q[3];
rz(-0.73368254) q[3];
sx q[3];
rz(2.0077326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
