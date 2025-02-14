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
rz(-2.2987125) q[1];
sx q[1];
rz(-2.892363) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1220142) q[0];
sx q[0];
rz(-1.459157) q[0];
sx q[0];
rz(0.48099244) q[0];
rz(-pi) q[1];
rz(0.31166844) q[2];
sx q[2];
rz(-2.083792) q[2];
sx q[2];
rz(-0.58849653) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8650444) q[1];
sx q[1];
rz(-2.2575592) q[1];
sx q[1];
rz(-1.264132) q[1];
rz(1.4572108) q[3];
sx q[3];
rz(-1.2143597) q[3];
sx q[3];
rz(-2.4939052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2117846) q[2];
sx q[2];
rz(-0.70694184) q[2];
sx q[2];
rz(2.9014273) q[2];
rz(2.5943878) q[3];
sx q[3];
rz(-1.5144843) q[3];
sx q[3];
rz(0.98946324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419679) q[0];
sx q[0];
rz(-2.7662321) q[0];
sx q[0];
rz(-0.64390916) q[0];
rz(-2.0945235) q[1];
sx q[1];
rz(-1.1771076) q[1];
sx q[1];
rz(0.26328304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8174648) q[0];
sx q[0];
rz(-3.1006515) q[0];
sx q[0];
rz(-0.73887478) q[0];
x q[1];
rz(0.77575923) q[2];
sx q[2];
rz(-3.0490766) q[2];
sx q[2];
rz(-0.55769071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.689641) q[1];
sx q[1];
rz(-0.76162377) q[1];
sx q[1];
rz(-2.2030001) q[1];
rz(-pi) q[2];
rz(-0.41384048) q[3];
sx q[3];
rz(-2.1713421) q[3];
sx q[3];
rz(2.5418856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.457966) q[2];
sx q[2];
rz(-1.7650975) q[2];
sx q[2];
rz(2.4158884) q[2];
rz(2.8343685) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(1.0544302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82045186) q[0];
sx q[0];
rz(-1.6571925) q[0];
sx q[0];
rz(0.76903525) q[0];
rz(-0.10737315) q[1];
sx q[1];
rz(-1.702405) q[1];
sx q[1];
rz(-0.69951397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0663129) q[0];
sx q[0];
rz(-1.9121721) q[0];
sx q[0];
rz(-2.7983448) q[0];
rz(1.5753463) q[2];
sx q[2];
rz(-0.9560491) q[2];
sx q[2];
rz(2.4640623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86068639) q[1];
sx q[1];
rz(-1.4889034) q[1];
sx q[1];
rz(-1.1204835) q[1];
rz(-pi) q[2];
rz(0.84229509) q[3];
sx q[3];
rz(-2.1811322) q[3];
sx q[3];
rz(-1.6323324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7677782) q[2];
sx q[2];
rz(-0.57743293) q[2];
sx q[2];
rz(-0.90816298) q[2];
rz(0.57458893) q[3];
sx q[3];
rz(-1.2535973) q[3];
sx q[3];
rz(2.7558034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.813756) q[0];
sx q[0];
rz(-1.8151374) q[0];
sx q[0];
rz(2.6326219) q[0];
rz(3.0834037) q[1];
sx q[1];
rz(-1.4570844) q[1];
sx q[1];
rz(1.0675272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71004407) q[0];
sx q[0];
rz(-0.87418927) q[0];
sx q[0];
rz(-1.0846653) q[0];
x q[1];
rz(3.0661461) q[2];
sx q[2];
rz(-1.9816996) q[2];
sx q[2];
rz(2.6332847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6696251) q[1];
sx q[1];
rz(-1.6780919) q[1];
sx q[1];
rz(1.7194242) q[1];
rz(-pi) q[2];
rz(-1.8640609) q[3];
sx q[3];
rz(-1.6365657) q[3];
sx q[3];
rz(1.4566959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.031614583) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(-0.99008375) q[2];
rz(-1.7456985) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(1.2466189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083967) q[0];
sx q[0];
rz(-1.7447423) q[0];
sx q[0];
rz(2.0821849) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.6672986) q[1];
sx q[1];
rz(-1.7105191) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16742584) q[0];
sx q[0];
rz(-1.7210809) q[0];
sx q[0];
rz(1.7420064) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25974187) q[2];
sx q[2];
rz(-2.0613823) q[2];
sx q[2];
rz(-0.39133137) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7943253) q[1];
sx q[1];
rz(-1.3293943) q[1];
sx q[1];
rz(2.4893087) q[1];
rz(2.5364872) q[3];
sx q[3];
rz(-0.25854585) q[3];
sx q[3];
rz(-0.097878067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3502256) q[2];
sx q[2];
rz(-2.3855049) q[2];
sx q[2];
rz(-1.6737326) q[2];
rz(-1.0427467) q[3];
sx q[3];
rz(-2.0839033) q[3];
sx q[3];
rz(-0.79152542) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4333711) q[0];
sx q[0];
rz(-1.8480166) q[0];
sx q[0];
rz(0.21155393) q[0];
rz(-2.4987706) q[1];
sx q[1];
rz(-2.0170409) q[1];
sx q[1];
rz(2.0932253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38163588) q[0];
sx q[0];
rz(-2.4279729) q[0];
sx q[0];
rz(2.9455393) q[0];
rz(-2.9625234) q[2];
sx q[2];
rz(-1.6229651) q[2];
sx q[2];
rz(1.6889926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.590946) q[1];
sx q[1];
rz(-0.63884402) q[1];
sx q[1];
rz(0.039012564) q[1];
x q[2];
rz(-2.1833352) q[3];
sx q[3];
rz(-1.6102092) q[3];
sx q[3];
rz(0.87557169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6362777) q[2];
sx q[2];
rz(-1.2039528) q[2];
sx q[2];
rz(0.25406507) q[2];
rz(-2.9680179) q[3];
sx q[3];
rz(-2.5901399) q[3];
sx q[3];
rz(1.9614296) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59488615) q[0];
sx q[0];
rz(-2.7017024) q[0];
sx q[0];
rz(-2.34483) q[0];
rz(-2.2548389) q[1];
sx q[1];
rz(-1.795307) q[1];
sx q[1];
rz(-0.60060445) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85854218) q[0];
sx q[0];
rz(-1.5158537) q[0];
sx q[0];
rz(-2.8514329) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4378635) q[2];
sx q[2];
rz(-0.80604751) q[2];
sx q[2];
rz(-1.9423167) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6930042) q[1];
sx q[1];
rz(-1.4088165) q[1];
sx q[1];
rz(2.2605091) q[1];
rz(2.1280471) q[3];
sx q[3];
rz(-2.5339014) q[3];
sx q[3];
rz(0.26085873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1116921) q[2];
sx q[2];
rz(-1.4375968) q[2];
sx q[2];
rz(2.0491484) q[2];
rz(1.368329) q[3];
sx q[3];
rz(-2.5320801) q[3];
sx q[3];
rz(-0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4261037) q[0];
sx q[0];
rz(-1.1142718) q[0];
sx q[0];
rz(2.6901167) q[0];
rz(-0.077797912) q[1];
sx q[1];
rz(-1.0323689) q[1];
sx q[1];
rz(-2.480004) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104354) q[0];
sx q[0];
rz(-1.9826188) q[0];
sx q[0];
rz(2.9020082) q[0];
rz(-pi) q[1];
rz(-0.069839283) q[2];
sx q[2];
rz(-2.0681212) q[2];
sx q[2];
rz(2.2004009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.132838) q[1];
sx q[1];
rz(-2.0292205) q[1];
sx q[1];
rz(-2.9944592) q[1];
rz(-pi) q[2];
rz(2.745146) q[3];
sx q[3];
rz(-1.8155193) q[3];
sx q[3];
rz(-0.94371599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0852069) q[2];
sx q[2];
rz(-1.8737917) q[2];
sx q[2];
rz(-0.80643225) q[2];
rz(1.8993529) q[3];
sx q[3];
rz(-0.16568383) q[3];
sx q[3];
rz(0.14036673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8388782) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(0.43933991) q[0];
rz(-1.2002523) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(-2.1913948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10748246) q[0];
sx q[0];
rz(-2.378298) q[0];
sx q[0];
rz(1.0176246) q[0];
rz(-pi) q[1];
rz(0.15301159) q[2];
sx q[2];
rz(-1.8540314) q[2];
sx q[2];
rz(0.52545122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2216144) q[1];
sx q[1];
rz(-1.6755381) q[1];
sx q[1];
rz(1.2135189) q[1];
x q[2];
rz(-2.1113273) q[3];
sx q[3];
rz(-1.221962) q[3];
sx q[3];
rz(0.23603786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3085559) q[2];
sx q[2];
rz(-0.50806442) q[2];
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
rz(3.062716) q[0];
sx q[0];
rz(-2.833241) q[0];
sx q[0];
rz(0.69945949) q[0];
rz(-1.4106916) q[1];
sx q[1];
rz(-0.78838333) q[1];
sx q[1];
rz(2.9404822) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5546306) q[0];
sx q[0];
rz(-2.2788958) q[0];
sx q[0];
rz(-0.59815852) q[0];
rz(-0.91014782) q[2];
sx q[2];
rz(-2.4322699) q[2];
sx q[2];
rz(-0.34914474) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48702588) q[1];
sx q[1];
rz(-2.2495765) q[1];
sx q[1];
rz(0.28622932) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9617025) q[3];
sx q[3];
rz(-1.6597865) q[3];
sx q[3];
rz(2.6266971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9451399) q[2];
sx q[2];
rz(-1.2597193) q[2];
sx q[2];
rz(0.067616612) q[2];
rz(-2.0130646) q[3];
sx q[3];
rz(-1.0849378) q[3];
sx q[3];
rz(-1.6245925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
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
rz(1.5051399) q[1];
sx q[1];
rz(-2.9121193) q[1];
sx q[1];
rz(3.0519003) q[1];
rz(-2.8436974) q[2];
sx q[2];
rz(-2.5174601) q[2];
sx q[2];
rz(-2.8909825) q[2];
rz(-0.33050362) q[3];
sx q[3];
rz(-1.3105583) q[3];
sx q[3];
rz(1.185598) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
