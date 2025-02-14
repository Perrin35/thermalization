OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6177144) q[0];
sx q[0];
rz(-0.59433794) q[0];
sx q[0];
rz(0.30881) q[0];
rz(0.29769695) q[1];
sx q[1];
rz(4.3990064) q[1];
sx q[1];
rz(10.036751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64498211) q[0];
sx q[0];
rz(-0.57304806) q[0];
sx q[0];
rz(1.4047506) q[0];
rz(0.16024477) q[2];
sx q[2];
rz(-1.8547684) q[2];
sx q[2];
rz(2.430973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0667324) q[1];
sx q[1];
rz(-1.7111527) q[1];
sx q[1];
rz(1.235515) q[1];
x q[2];
rz(1.520492) q[3];
sx q[3];
rz(-2.2626749) q[3];
sx q[3];
rz(1.8722201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.97418) q[2];
sx q[2];
rz(-1.9638502) q[2];
sx q[2];
rz(2.5766032) q[2];
rz(2.6307093) q[3];
sx q[3];
rz(-0.23879819) q[3];
sx q[3];
rz(1.5272944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086455) q[0];
sx q[0];
rz(-0.2810418) q[0];
sx q[0];
rz(0.99579048) q[0];
rz(0.58798724) q[1];
sx q[1];
rz(-0.43991393) q[1];
sx q[1];
rz(-1.3194552) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73896995) q[0];
sx q[0];
rz(-1.8085294) q[0];
sx q[0];
rz(-2.0031702) q[0];
x q[1];
rz(0.075602268) q[2];
sx q[2];
rz(-1.8972862) q[2];
sx q[2];
rz(-1.9734427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.118757) q[1];
sx q[1];
rz(-2.3110227) q[1];
sx q[1];
rz(2.8603795) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6922582) q[3];
sx q[3];
rz(-2.0503042) q[3];
sx q[3];
rz(-3.1288655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1002645) q[2];
sx q[2];
rz(-2.1375956) q[2];
sx q[2];
rz(-0.022493258) q[2];
rz(3.0371173) q[3];
sx q[3];
rz(-1.6040809) q[3];
sx q[3];
rz(0.87048602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81258881) q[0];
sx q[0];
rz(-2.1109695) q[0];
sx q[0];
rz(1.6382244) q[0];
rz(-3.0606048) q[1];
sx q[1];
rz(-0.65892017) q[1];
sx q[1];
rz(1.9714877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2161897) q[0];
sx q[0];
rz(-1.4099551) q[0];
sx q[0];
rz(0.58165929) q[0];
rz(1.3721714) q[2];
sx q[2];
rz(-2.2057057) q[2];
sx q[2];
rz(0.18066192) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8997716) q[1];
sx q[1];
rz(-1.5994497) q[1];
sx q[1];
rz(1.0706399) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.891344) q[3];
sx q[3];
rz(-2.0648533) q[3];
sx q[3];
rz(1.9039465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66204232) q[2];
sx q[2];
rz(-0.68632555) q[2];
sx q[2];
rz(0.043206841) q[2];
rz(2.9442545) q[3];
sx q[3];
rz(-2.1102326) q[3];
sx q[3];
rz(0.432338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8869121) q[0];
sx q[0];
rz(-1.0858902) q[0];
sx q[0];
rz(-0.32549724) q[0];
rz(0.85572851) q[1];
sx q[1];
rz(-2.7266462) q[1];
sx q[1];
rz(1.0861446) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551712) q[0];
sx q[0];
rz(-2.060411) q[0];
sx q[0];
rz(-1.7167164) q[0];
rz(2.3065673) q[2];
sx q[2];
rz(-2.1673598) q[2];
sx q[2];
rz(1.8955829) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4068674) q[1];
sx q[1];
rz(-1.735209) q[1];
sx q[1];
rz(-2.0937992) q[1];
rz(0.086768199) q[3];
sx q[3];
rz(-1.4283071) q[3];
sx q[3];
rz(-0.87228197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2073652) q[2];
sx q[2];
rz(-0.97565979) q[2];
sx q[2];
rz(-2.9577241) q[2];
rz(-1.81987) q[3];
sx q[3];
rz(-0.70122856) q[3];
sx q[3];
rz(0.13300657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4399399) q[0];
sx q[0];
rz(-0.30783215) q[0];
sx q[0];
rz(-0.93691784) q[0];
rz(2.2266455) q[1];
sx q[1];
rz(-2.2888384) q[1];
sx q[1];
rz(1.6818887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2258218) q[0];
sx q[0];
rz(-2.2647595) q[0];
sx q[0];
rz(3.0208605) q[0];
rz(-pi) q[1];
rz(-2.5499623) q[2];
sx q[2];
rz(-1.0359633) q[2];
sx q[2];
rz(0.66885751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25681695) q[1];
sx q[1];
rz(-0.76013541) q[1];
sx q[1];
rz(1.7124618) q[1];
x q[2];
rz(0.57444467) q[3];
sx q[3];
rz(-0.46009053) q[3];
sx q[3];
rz(-2.121832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0437643) q[2];
sx q[2];
rz(-1.8057258) q[2];
sx q[2];
rz(-0.14399993) q[2];
rz(1.0675659) q[3];
sx q[3];
rz(-0.34393603) q[3];
sx q[3];
rz(2.4501154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095793) q[0];
sx q[0];
rz(-0.80736512) q[0];
sx q[0];
rz(2.373234) q[0];
rz(-0.8575303) q[1];
sx q[1];
rz(-1.0464959) q[1];
sx q[1];
rz(0.55364496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.377877) q[0];
sx q[0];
rz(-1.3175497) q[0];
sx q[0];
rz(-1.1035155) q[0];
rz(-pi) q[1];
rz(-2.3403999) q[2];
sx q[2];
rz(-1.7320219) q[2];
sx q[2];
rz(3.1163505) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.044477) q[1];
sx q[1];
rz(-0.71050853) q[1];
sx q[1];
rz(2.249975) q[1];
x q[2];
rz(2.4870212) q[3];
sx q[3];
rz(-2.304545) q[3];
sx q[3];
rz(-0.76920054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7553317) q[2];
sx q[2];
rz(-0.43060455) q[2];
sx q[2];
rz(-0.44580305) q[2];
rz(1.2465994) q[3];
sx q[3];
rz(-1.3695025) q[3];
sx q[3];
rz(2.2915452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069139473) q[0];
sx q[0];
rz(-0.9599762) q[0];
sx q[0];
rz(3.0015216) q[0];
rz(-2.2135997) q[1];
sx q[1];
rz(-1.3368139) q[1];
sx q[1];
rz(1.450052) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194013) q[0];
sx q[0];
rz(-2.9474576) q[0];
sx q[0];
rz(2.1737264) q[0];
x q[1];
rz(-1.9533402) q[2];
sx q[2];
rz(-1.3556644) q[2];
sx q[2];
rz(2.6184788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90277744) q[1];
sx q[1];
rz(-1.1401145) q[1];
sx q[1];
rz(-0.68044739) q[1];
rz(-0.18168707) q[3];
sx q[3];
rz(-1.4188926) q[3];
sx q[3];
rz(0.10703281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.037584) q[2];
sx q[2];
rz(-0.17825492) q[2];
sx q[2];
rz(-0.083871052) q[2];
rz(0.9907848) q[3];
sx q[3];
rz(-1.1889941) q[3];
sx q[3];
rz(1.2880464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.24638076) q[0];
sx q[0];
rz(-1.3061433) q[0];
sx q[0];
rz(0.35686785) q[0];
rz(-1.6995947) q[1];
sx q[1];
rz(-0.60209638) q[1];
sx q[1];
rz(0.009036202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4215338) q[0];
sx q[0];
rz(-1.3429317) q[0];
sx q[0];
rz(0.22703815) q[0];
x q[1];
rz(2.4414012) q[2];
sx q[2];
rz(-1.2692361) q[2];
sx q[2];
rz(-1.4789326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0865514) q[1];
sx q[1];
rz(-2.8541871) q[1];
sx q[1];
rz(1.2432008) q[1];
rz(2.6983016) q[3];
sx q[3];
rz(-1.6066243) q[3];
sx q[3];
rz(1.9323424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8845727) q[2];
sx q[2];
rz(-1.7947861) q[2];
sx q[2];
rz(-2.436077) q[2];
rz(0.26933119) q[3];
sx q[3];
rz(-1.0764542) q[3];
sx q[3];
rz(-2.1721325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89328289) q[0];
sx q[0];
rz(-0.80535424) q[0];
sx q[0];
rz(2.4321108) q[0];
rz(-0.64208883) q[1];
sx q[1];
rz(-2.700192) q[1];
sx q[1];
rz(-2.716224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4152483) q[0];
sx q[0];
rz(-1.9972902) q[0];
sx q[0];
rz(0.70856673) q[0];
rz(-pi) q[1];
rz(0.98871059) q[2];
sx q[2];
rz(-1.9217689) q[2];
sx q[2];
rz(-2.7711353) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44008419) q[1];
sx q[1];
rz(-0.45508859) q[1];
sx q[1];
rz(0.20642682) q[1];
rz(-1.2669417) q[3];
sx q[3];
rz(-1.0179449) q[3];
sx q[3];
rz(1.0859717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2702668) q[2];
sx q[2];
rz(-0.17921236) q[2];
sx q[2];
rz(0.47879177) q[2];
rz(-0.35017961) q[3];
sx q[3];
rz(-1.1778573) q[3];
sx q[3];
rz(0.42696264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66788524) q[0];
sx q[0];
rz(-0.29692867) q[0];
sx q[0];
rz(2.3042451) q[0];
rz(0.26652023) q[1];
sx q[1];
rz(-1.8217249) q[1];
sx q[1];
rz(-1.7841608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0251009) q[0];
sx q[0];
rz(-3.0550346) q[0];
sx q[0];
rz(1.5309912) q[0];
x q[1];
rz(-0.063284782) q[2];
sx q[2];
rz(-1.4522168) q[2];
sx q[2];
rz(2.0381387) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7559133) q[1];
sx q[1];
rz(-2.2230621) q[1];
sx q[1];
rz(1.5550809) q[1];
rz(-pi) q[2];
rz(-2.1641047) q[3];
sx q[3];
rz(-2.5989669) q[3];
sx q[3];
rz(-3.0491587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9799161) q[2];
sx q[2];
rz(-2.6075173) q[2];
sx q[2];
rz(-2.7034289) q[2];
rz(2.2257889) q[3];
sx q[3];
rz(-1.5235498) q[3];
sx q[3];
rz(2.3384136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61951471) q[0];
sx q[0];
rz(-1.2418455) q[0];
sx q[0];
rz(-1.0157304) q[0];
rz(-0.10454128) q[1];
sx q[1];
rz(-1.3019982) q[1];
sx q[1];
rz(-1.7560538) q[1];
rz(-1.9011433) q[2];
sx q[2];
rz(-1.6295682) q[2];
sx q[2];
rz(-0.30420797) q[2];
rz(-2.2661891) q[3];
sx q[3];
rz(-2.2973552) q[3];
sx q[3];
rz(-2.4612343) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
