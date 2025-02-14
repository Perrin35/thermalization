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
rz(-0.14577785) q[0];
sx q[0];
rz(-1.3480027) q[0];
sx q[0];
rz(0.3682799) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(0.84288016) q[1];
sx q[1];
rz(9.1755484) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6347353) q[0];
sx q[0];
rz(-1.0930499) q[0];
sx q[0];
rz(1.696582) q[0];
rz(-0.31166844) q[2];
sx q[2];
rz(-1.0578007) q[2];
sx q[2];
rz(-0.58849653) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.492384) q[1];
sx q[1];
rz(-1.8064152) q[1];
sx q[1];
rz(0.71028965) q[1];
rz(-1.6843819) q[3];
sx q[3];
rz(-1.2143597) q[3];
sx q[3];
rz(0.64768744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92980802) q[2];
sx q[2];
rz(-0.70694184) q[2];
sx q[2];
rz(-0.24016538) q[2];
rz(-0.54720488) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3419679) q[0];
sx q[0];
rz(-2.7662321) q[0];
sx q[0];
rz(0.64390916) q[0];
rz(1.0470692) q[1];
sx q[1];
rz(-1.1771076) q[1];
sx q[1];
rz(0.26328304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8174648) q[0];
sx q[0];
rz(-3.1006515) q[0];
sx q[0];
rz(0.73887478) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3658334) q[2];
sx q[2];
rz(-0.09251602) q[2];
sx q[2];
rz(-2.5839019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.689641) q[1];
sx q[1];
rz(-2.3799689) q[1];
sx q[1];
rz(-2.2030001) q[1];
rz(-2.2130737) q[3];
sx q[3];
rz(-1.9089724) q[3];
sx q[3];
rz(2.4137792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(-2.4158884) q[2];
rz(0.30722412) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(2.0871625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82045186) q[0];
sx q[0];
rz(-1.4844002) q[0];
sx q[0];
rz(2.3725574) q[0];
rz(-3.0342195) q[1];
sx q[1];
rz(-1.702405) q[1];
sx q[1];
rz(-2.4420787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25722593) q[0];
sx q[0];
rz(-0.47927893) q[0];
sx q[0];
rz(2.3291161) q[0];
rz(-pi) q[1];
rz(1.5662464) q[2];
sx q[2];
rz(-2.1855436) q[2];
sx q[2];
rz(-0.67753032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87772885) q[1];
sx q[1];
rz(-2.6843964) q[1];
sx q[1];
rz(-1.3844107) q[1];
rz(-pi) q[2];
rz(0.76074184) q[3];
sx q[3];
rz(-0.91289744) q[3];
sx q[3];
rz(2.6321312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7677782) q[2];
sx q[2];
rz(-2.5641597) q[2];
sx q[2];
rz(-0.90816298) q[2];
rz(2.5670037) q[3];
sx q[3];
rz(-1.8879954) q[3];
sx q[3];
rz(-0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.813756) q[0];
sx q[0];
rz(-1.8151374) q[0];
sx q[0];
rz(2.6326219) q[0];
rz(-3.0834037) q[1];
sx q[1];
rz(-1.4570844) q[1];
sx q[1];
rz(2.0740654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71004407) q[0];
sx q[0];
rz(-2.2674034) q[0];
sx q[0];
rz(2.0569274) q[0];
rz(-pi) q[1];
rz(1.1588489) q[2];
sx q[2];
rz(-1.6399523) q[2];
sx q[2];
rz(2.1092887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6696251) q[1];
sx q[1];
rz(-1.6780919) q[1];
sx q[1];
rz(-1.7194242) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8640609) q[3];
sx q[3];
rz(-1.505027) q[3];
sx q[3];
rz(1.6848967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1099781) q[2];
sx q[2];
rz(-1.5974533) q[2];
sx q[2];
rz(2.1515089) q[2];
rz(1.7456985) q[3];
sx q[3];
rz(-0.11002222) q[3];
sx q[3];
rz(-1.8949738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083967) q[0];
sx q[0];
rz(-1.3968503) q[0];
sx q[0];
rz(-1.0594077) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.474294) q[1];
sx q[1];
rz(-1.4310736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0245594) q[0];
sx q[0];
rz(-2.9142671) q[0];
sx q[0];
rz(0.84421279) q[0];
rz(-pi) q[1];
rz(1.0658979) q[2];
sx q[2];
rz(-1.7993198) q[2];
sx q[2];
rz(-1.8375719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0986745) q[1];
sx q[1];
rz(-0.94047767) q[1];
sx q[1];
rz(1.8712256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7201172) q[3];
sx q[3];
rz(-1.7826555) q[3];
sx q[3];
rz(2.6184591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7913671) q[2];
sx q[2];
rz(-2.3855049) q[2];
sx q[2];
rz(-1.6737326) q[2];
rz(-2.098846) q[3];
sx q[3];
rz(-1.0576893) q[3];
sx q[3];
rz(2.3500672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4333711) q[0];
sx q[0];
rz(-1.2935761) q[0];
sx q[0];
rz(-0.21155393) q[0];
rz(-0.64282203) q[1];
sx q[1];
rz(-1.1245518) q[1];
sx q[1];
rz(-1.0483673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38163588) q[0];
sx q[0];
rz(-0.71361976) q[0];
sx q[0];
rz(2.9455393) q[0];
rz(2.8564212) q[2];
sx q[2];
rz(-0.18643555) q[2];
sx q[2];
rz(2.9793019) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0901186) q[1];
sx q[1];
rz(-1.5475382) q[1];
sx q[1];
rz(0.6384797) q[1];
rz(0.95825742) q[3];
sx q[3];
rz(-1.5313834) q[3];
sx q[3];
rz(2.266021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50531498) q[2];
sx q[2];
rz(-1.9376398) q[2];
sx q[2];
rz(-0.25406507) q[2];
rz(2.9680179) q[3];
sx q[3];
rz(-2.5901399) q[3];
sx q[3];
rz(1.180163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5467065) q[0];
sx q[0];
rz(-0.43989023) q[0];
sx q[0];
rz(-2.34483) q[0];
rz(-2.2548389) q[1];
sx q[1];
rz(-1.795307) q[1];
sx q[1];
rz(-0.60060445) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89413801) q[0];
sx q[0];
rz(-2.8464212) q[0];
sx q[0];
rz(2.9516793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76917665) q[2];
sx q[2];
rz(-1.666579) q[2];
sx q[2];
rz(2.8623919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.010041865) q[1];
sx q[1];
rz(-2.2497592) q[1];
sx q[1];
rz(0.20874397) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7891385) q[3];
sx q[3];
rz(-2.0766933) q[3];
sx q[3];
rz(-2.231489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1116921) q[2];
sx q[2];
rz(-1.7039958) q[2];
sx q[2];
rz(-1.0924443) q[2];
rz(1.368329) q[3];
sx q[3];
rz(-2.5320801) q[3];
sx q[3];
rz(-0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4261037) q[0];
sx q[0];
rz(-1.1142718) q[0];
sx q[0];
rz(-0.45147595) q[0];
rz(0.077797912) q[1];
sx q[1];
rz(-1.0323689) q[1];
sx q[1];
rz(-0.66158867) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104354) q[0];
sx q[0];
rz(-1.9826188) q[0];
sx q[0];
rz(-2.9020082) q[0];
rz(0.069839283) q[2];
sx q[2];
rz(-2.0681212) q[2];
sx q[2];
rz(0.94119173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.132838) q[1];
sx q[1];
rz(-2.0292205) q[1];
sx q[1];
rz(0.14713344) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57406942) q[3];
sx q[3];
rz(-0.46246734) q[3];
sx q[3];
rz(-1.9898349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0563858) q[2];
sx q[2];
rz(-1.2678009) q[2];
sx q[2];
rz(0.80643225) q[2];
rz(1.8993529) q[3];
sx q[3];
rz(-0.16568383) q[3];
sx q[3];
rz(-3.0012259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8388782) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(2.7022527) q[0];
rz(1.9413403) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(-2.1913948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0341102) q[0];
sx q[0];
rz(-2.378298) q[0];
sx q[0];
rz(2.123968) q[0];
rz(0.15301159) q[2];
sx q[2];
rz(-1.2875612) q[2];
sx q[2];
rz(2.6161414) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31017329) q[1];
sx q[1];
rz(-1.2155639) q[1];
sx q[1];
rz(-3.0298477) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7404157) q[3];
sx q[3];
rz(-1.0660352) q[3];
sx q[3];
rz(-1.1324319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83303678) q[2];
sx q[2];
rz(-2.6335282) q[2];
sx q[2];
rz(-1.1888095) q[2];
rz(3.0360119) q[3];
sx q[3];
rz(-0.9413541) q[3];
sx q[3];
rz(-2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.062716) q[0];
sx q[0];
rz(-0.3083516) q[0];
sx q[0];
rz(0.69945949) q[0];
rz(-1.4106916) q[1];
sx q[1];
rz(-2.3532093) q[1];
sx q[1];
rz(-2.9404822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58696207) q[0];
sx q[0];
rz(-0.86269683) q[0];
sx q[0];
rz(2.5434341) q[0];
rz(-pi) q[1];
rz(-0.91014782) q[2];
sx q[2];
rz(-0.70932271) q[2];
sx q[2];
rz(-2.7924479) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2664908) q[1];
sx q[1];
rz(-1.7923585) q[1];
sx q[1];
rz(-2.2699577) q[1];
x q[2];
rz(-1.6612382) q[3];
sx q[3];
rz(-1.3916257) q[3];
sx q[3];
rz(-2.1018525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19645277) q[2];
sx q[2];
rz(-1.2597193) q[2];
sx q[2];
rz(-0.067616612) q[2];
rz(1.128528) q[3];
sx q[3];
rz(-2.0566548) q[3];
sx q[3];
rz(1.6245925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45831281) q[0];
sx q[0];
rz(-2.0742317) q[0];
sx q[0];
rz(-1.7924894) q[0];
rz(-1.6364527) q[1];
sx q[1];
rz(-2.9121193) q[1];
sx q[1];
rz(3.0519003) q[1];
rz(0.2978953) q[2];
sx q[2];
rz(-2.5174601) q[2];
sx q[2];
rz(-2.8909825) q[2];
rz(2.811089) q[3];
sx q[3];
rz(-1.3105583) q[3];
sx q[3];
rz(1.185598) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
