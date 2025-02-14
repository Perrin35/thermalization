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
rz(2.9958148) q[0];
sx q[0];
rz(-1.79359) q[0];
sx q[0];
rz(-0.3682799) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(-2.2987125) q[1];
sx q[1];
rz(0.24922961) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1220142) q[0];
sx q[0];
rz(-1.6824357) q[0];
sx q[0];
rz(0.48099244) q[0];
x q[1];
rz(0.31166844) q[2];
sx q[2];
rz(-1.0578007) q[2];
sx q[2];
rz(-2.5530961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6492086) q[1];
sx q[1];
rz(-1.8064152) q[1];
sx q[1];
rz(-0.71028965) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29549142) q[3];
sx q[3];
rz(-2.7682332) q[3];
sx q[3];
rz(-0.96366027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92980802) q[2];
sx q[2];
rz(-0.70694184) q[2];
sx q[2];
rz(0.24016538) q[2];
rz(-2.5943878) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(-2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7996247) q[0];
sx q[0];
rz(-0.37536055) q[0];
sx q[0];
rz(0.64390916) q[0];
rz(2.0945235) q[1];
sx q[1];
rz(-1.1771076) q[1];
sx q[1];
rz(2.8783096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7264284) q[0];
sx q[0];
rz(-1.5405354) q[0];
sx q[0];
rz(1.5432158) q[0];
rz(0.77575923) q[2];
sx q[2];
rz(-3.0490766) q[2];
sx q[2];
rz(-0.55769071) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66052478) q[1];
sx q[1];
rz(-0.98036375) q[1];
sx q[1];
rz(0.51312889) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41384048) q[3];
sx q[3];
rz(-2.1713421) q[3];
sx q[3];
rz(-2.5418856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(-2.4158884) q[2];
rz(-2.8343685) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(-1.0544302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211408) q[0];
sx q[0];
rz(-1.6571925) q[0];
sx q[0];
rz(2.3725574) q[0];
rz(-0.10737315) q[1];
sx q[1];
rz(-1.702405) q[1];
sx q[1];
rz(2.4420787) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61460068) q[0];
sx q[0];
rz(-1.8934947) q[0];
sx q[0];
rz(-1.2100156) q[0];
rz(-pi) q[1];
rz(1.5753463) q[2];
sx q[2];
rz(-2.1855436) q[2];
sx q[2];
rz(-2.4640623) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4710083) q[1];
sx q[1];
rz(-2.0194897) q[1];
sx q[1];
rz(3.0506794) q[1];
x q[2];
rz(0.75306692) q[3];
sx q[3];
rz(-0.99374607) q[3];
sx q[3];
rz(-2.6073539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37381441) q[2];
sx q[2];
rz(-2.5641597) q[2];
sx q[2];
rz(2.2334297) q[2];
rz(-0.57458893) q[3];
sx q[3];
rz(-1.2535973) q[3];
sx q[3];
rz(0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3278367) q[0];
sx q[0];
rz(-1.3264553) q[0];
sx q[0];
rz(2.6326219) q[0];
rz(0.058188997) q[1];
sx q[1];
rz(-1.4570844) q[1];
sx q[1];
rz(-1.0675272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021085652) q[0];
sx q[0];
rz(-2.3160546) q[0];
sx q[0];
rz(2.6322281) q[0];
rz(-pi) q[1];
rz(3.0661461) q[2];
sx q[2];
rz(-1.1598931) q[2];
sx q[2];
rz(0.50830799) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0587973) q[1];
sx q[1];
rz(-1.4230295) q[1];
sx q[1];
rz(3.0331103) q[1];
x q[2];
rz(-1.3467784) q[3];
sx q[3];
rz(-0.30034143) q[3];
sx q[3];
rz(-3.0413922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1099781) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(-0.99008375) q[2];
rz(1.3958942) q[3];
sx q[3];
rz(-0.11002222) q[3];
sx q[3];
rz(-1.2466189) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2331959) q[0];
sx q[0];
rz(-1.3968503) q[0];
sx q[0];
rz(-2.0821849) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.474294) q[1];
sx q[1];
rz(1.7105191) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1170332) q[0];
sx q[0];
rz(-2.9142671) q[0];
sx q[0];
rz(2.2973799) q[0];
rz(-pi) q[1];
rz(2.0756948) q[2];
sx q[2];
rz(-1.3422728) q[2];
sx q[2];
rz(-1.8375719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0986745) q[1];
sx q[1];
rz(-2.201115) q[1];
sx q[1];
rz(-1.8712256) q[1];
rz(1.4214755) q[3];
sx q[3];
rz(-1.3589371) q[3];
sx q[3];
rz(2.6184591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.3502256) q[2];
sx q[2];
rz(-0.75608772) q[2];
sx q[2];
rz(1.4678601) q[2];
rz(-2.098846) q[3];
sx q[3];
rz(-1.0576893) q[3];
sx q[3];
rz(2.3500672) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4333711) q[0];
sx q[0];
rz(-1.8480166) q[0];
sx q[0];
rz(2.9300387) q[0];
rz(2.4987706) q[1];
sx q[1];
rz(-2.0170409) q[1];
sx q[1];
rz(1.0483673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12473561) q[0];
sx q[0];
rz(-0.87364679) q[0];
sx q[0];
rz(-1.7378896) q[0];
rz(-2.9625234) q[2];
sx q[2];
rz(-1.6229651) q[2];
sx q[2];
rz(-1.4526001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.051474007) q[1];
sx q[1];
rz(-1.5475382) q[1];
sx q[1];
rz(-0.6384797) q[1];
rz(1.5023175) q[3];
sx q[3];
rz(-2.5279495) q[3];
sx q[3];
rz(2.5023823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50531498) q[2];
sx q[2];
rz(-1.9376398) q[2];
sx q[2];
rz(-2.8875276) q[2];
rz(-2.9680179) q[3];
sx q[3];
rz(-0.55145276) q[3];
sx q[3];
rz(1.180163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.5467065) q[0];
sx q[0];
rz(-2.7017024) q[0];
sx q[0];
rz(0.79676262) q[0];
rz(-2.2548389) q[1];
sx q[1];
rz(-1.795307) q[1];
sx q[1];
rz(-0.60060445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89413801) q[0];
sx q[0];
rz(-2.8464212) q[0];
sx q[0];
rz(-0.18991332) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1372631) q[2];
sx q[2];
rz(-2.3676927) q[2];
sx q[2];
rz(1.3900666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4485885) q[1];
sx q[1];
rz(-1.4088165) q[1];
sx q[1];
rz(2.2605091) q[1];
rz(-pi) q[2];
rz(-2.7891385) q[3];
sx q[3];
rz(-2.0766933) q[3];
sx q[3];
rz(2.231489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1116921) q[2];
sx q[2];
rz(-1.4375968) q[2];
sx q[2];
rz(1.0924443) q[2];
rz(1.7732636) q[3];
sx q[3];
rz(-2.5320801) q[3];
sx q[3];
rz(0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4261037) q[0];
sx q[0];
rz(-2.0273209) q[0];
sx q[0];
rz(-2.6901167) q[0];
rz(0.077797912) q[1];
sx q[1];
rz(-1.0323689) q[1];
sx q[1];
rz(2.480004) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1790893) q[0];
sx q[0];
rz(-2.668619) q[0];
sx q[0];
rz(-1.0731368) q[0];
x q[1];
rz(1.0724474) q[2];
sx q[2];
rz(-1.5094286) q[2];
sx q[2];
rz(0.59624404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80967951) q[1];
sx q[1];
rz(-0.47985489) q[1];
sx q[1];
rz(1.859568) q[1];
rz(-pi) q[2];
rz(-1.8351848) q[3];
sx q[3];
rz(-1.9548023) q[3];
sx q[3];
rz(2.6155909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0852069) q[2];
sx q[2];
rz(-1.8737917) q[2];
sx q[2];
rz(-2.3351604) q[2];
rz(-1.8993529) q[3];
sx q[3];
rz(-2.9759088) q[3];
sx q[3];
rz(-3.0012259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8388782) q[0];
sx q[0];
rz(-2.0063722) q[0];
sx q[0];
rz(2.7022527) q[0];
rz(1.9413403) q[1];
sx q[1];
rz(-0.83963436) q[1];
sx q[1];
rz(2.1913948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59965682) q[0];
sx q[0];
rz(-2.1996351) q[0];
sx q[0];
rz(-2.6758183) q[0];
x q[1];
rz(-2.9885811) q[2];
sx q[2];
rz(-1.8540314) q[2];
sx q[2];
rz(0.52545122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5193351) q[1];
sx q[1];
rz(-0.37168113) q[1];
sx q[1];
rz(-1.2787914) q[1];
x q[2];
rz(-0.40117694) q[3];
sx q[3];
rz(-2.0755575) q[3];
sx q[3];
rz(-1.1324319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83303678) q[2];
sx q[2];
rz(-2.6335282) q[2];
sx q[2];
rz(-1.9527831) q[2];
rz(-0.10558072) q[3];
sx q[3];
rz(-2.2002386) q[3];
sx q[3];
rz(2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.062716) q[0];
sx q[0];
rz(-0.3083516) q[0];
sx q[0];
rz(-2.4421332) q[0];
rz(-1.4106916) q[1];
sx q[1];
rz(-0.78838333) q[1];
sx q[1];
rz(2.9404822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58696207) q[0];
sx q[0];
rz(-0.86269683) q[0];
sx q[0];
rz(-2.5434341) q[0];
rz(-2.1664326) q[2];
sx q[2];
rz(-1.9819519) q[2];
sx q[2];
rz(1.3871297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.048677532) q[1];
sx q[1];
rz(-2.4138422) q[1];
sx q[1];
rz(-1.2341094) q[1];
rz(-pi) q[2];
rz(-0.17989017) q[3];
sx q[3];
rz(-1.6597865) q[3];
sx q[3];
rz(-0.51489553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19645277) q[2];
sx q[2];
rz(-1.8818734) q[2];
sx q[2];
rz(-0.067616612) q[2];
rz(1.128528) q[3];
sx q[3];
rz(-1.0849378) q[3];
sx q[3];
rz(1.5170001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6832798) q[0];
sx q[0];
rz(-2.0742317) q[0];
sx q[0];
rz(-1.7924894) q[0];
rz(1.6364527) q[1];
sx q[1];
rz(-0.22947336) q[1];
sx q[1];
rz(-0.089692399) q[1];
rz(-2.8436974) q[2];
sx q[2];
rz(-2.5174601) q[2];
sx q[2];
rz(-2.8909825) q[2];
rz(-1.296386) q[3];
sx q[3];
rz(-1.2518223) q[3];
sx q[3];
rz(2.8444461) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
