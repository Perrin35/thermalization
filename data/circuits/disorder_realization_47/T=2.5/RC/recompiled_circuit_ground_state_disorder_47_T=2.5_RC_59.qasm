OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.2831777) q[0];
sx q[0];
rz(-1.5812961) q[0];
sx q[0];
rz(-0.66891447) q[0];
rz(-2.7954697) q[1];
sx q[1];
rz(-0.82155138) q[1];
sx q[1];
rz(-2.2023444) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4209983) q[0];
sx q[0];
rz(-1.2828553) q[0];
sx q[0];
rz(-1.9701411) q[0];
rz(-pi) q[1];
rz(-1.16667) q[2];
sx q[2];
rz(-1.6138645) q[2];
sx q[2];
rz(1.7801746) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2841442) q[1];
sx q[1];
rz(-1.836187) q[1];
sx q[1];
rz(-1.1737203) q[1];
x q[2];
rz(1.2143308) q[3];
sx q[3];
rz(-2.2495396) q[3];
sx q[3];
rz(-1.1671305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5758489) q[2];
sx q[2];
rz(-1.0172903) q[2];
sx q[2];
rz(3.1080833) q[2];
rz(-1.9401898) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(2.0638154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031438436) q[0];
sx q[0];
rz(-2.8903676) q[0];
sx q[0];
rz(2.187619) q[0];
rz(-1.6961478) q[1];
sx q[1];
rz(-2.0868389) q[1];
sx q[1];
rz(1.9704069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7273219) q[0];
sx q[0];
rz(-2.1571299) q[0];
sx q[0];
rz(0.81887736) q[0];
rz(-pi) q[1];
rz(-2.5043152) q[2];
sx q[2];
rz(-1.8101781) q[2];
sx q[2];
rz(2.0429037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5038472) q[1];
sx q[1];
rz(-0.23023573) q[1];
sx q[1];
rz(3.0158325) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0061139) q[3];
sx q[3];
rz(-1.2655228) q[3];
sx q[3];
rz(2.0729617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90616068) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(-1.99235) q[2];
rz(0.033128459) q[3];
sx q[3];
rz(-1.6413611) q[3];
sx q[3];
rz(-2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0572492) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(-1.0302011) q[0];
rz(1.239981) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(0.99303594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5169553) q[0];
sx q[0];
rz(-2.0258396) q[0];
sx q[0];
rz(-1.6242473) q[0];
x q[1];
rz(-0.66752357) q[2];
sx q[2];
rz(-1.4771531) q[2];
sx q[2];
rz(0.33622959) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5883623) q[1];
sx q[1];
rz(-0.85236406) q[1];
sx q[1];
rz(0.78939446) q[1];
rz(-pi) q[2];
rz(2.7168324) q[3];
sx q[3];
rz(-0.83213193) q[3];
sx q[3];
rz(0.55603851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.21732907) q[2];
sx q[2];
rz(-1.1161048) q[2];
sx q[2];
rz(2.7868311) q[2];
rz(2.1445856) q[3];
sx q[3];
rz(-1.6544147) q[3];
sx q[3];
rz(-2.2009489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724801) q[0];
sx q[0];
rz(-1.7153772) q[0];
sx q[0];
rz(2.0347563) q[0];
rz(-0.86889443) q[1];
sx q[1];
rz(-0.51742253) q[1];
sx q[1];
rz(0.69721627) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2223059) q[0];
sx q[0];
rz(-0.59658748) q[0];
sx q[0];
rz(0.35937341) q[0];
rz(-pi) q[1];
rz(-3.0142586) q[2];
sx q[2];
rz(-0.72605726) q[2];
sx q[2];
rz(-2.6301094) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2576506) q[1];
sx q[1];
rz(-1.202794) q[1];
sx q[1];
rz(-1.9872527) q[1];
x q[2];
rz(0.02454464) q[3];
sx q[3];
rz(-1.595579) q[3];
sx q[3];
rz(1.9056232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8640459) q[2];
sx q[2];
rz(-2.9329381) q[2];
sx q[2];
rz(2.8475658) q[2];
rz(-2.0781519) q[3];
sx q[3];
rz(-1.4592417) q[3];
sx q[3];
rz(-2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71642891) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(0.46318769) q[0];
rz(2.7376392) q[1];
sx q[1];
rz(-1.633176) q[1];
sx q[1];
rz(2.892866) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8578019) q[0];
sx q[0];
rz(-1.3732855) q[0];
sx q[0];
rz(3.1398612) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8291446) q[2];
sx q[2];
rz(-1.0444006) q[2];
sx q[2];
rz(2.0451982) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7686004) q[1];
sx q[1];
rz(-1.2443719) q[1];
sx q[1];
rz(1.1821163) q[1];
rz(-pi) q[2];
rz(-2.4234613) q[3];
sx q[3];
rz(-0.87801027) q[3];
sx q[3];
rz(-0.75961514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8467466) q[2];
sx q[2];
rz(-1.3174026) q[2];
sx q[2];
rz(1.741629) q[2];
rz(-1.8585662) q[3];
sx q[3];
rz(-1.7581519) q[3];
sx q[3];
rz(0.2259026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95752174) q[0];
sx q[0];
rz(-0.83860832) q[0];
sx q[0];
rz(-2.7958909) q[0];
rz(0.050994571) q[1];
sx q[1];
rz(-1.9146999) q[1];
sx q[1];
rz(0.7720224) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2564023) q[0];
sx q[0];
rz(-1.9213543) q[0];
sx q[0];
rz(1.8392483) q[0];
x q[1];
rz(0.096944158) q[2];
sx q[2];
rz(-2.2205995) q[2];
sx q[2];
rz(0.61496269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4881106) q[1];
sx q[1];
rz(-2.1381492) q[1];
sx q[1];
rz(-0.385872) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0431248) q[3];
sx q[3];
rz(-1.3195795) q[3];
sx q[3];
rz(2.4285897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2766075) q[2];
sx q[2];
rz(-1.7165246) q[2];
sx q[2];
rz(0.15743206) q[2];
rz(0.43846798) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(2.1854775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601783) q[0];
sx q[0];
rz(-2.0745451) q[0];
sx q[0];
rz(9/(11*pi)) q[0];
rz(-0.57834894) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-2.9885805) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48400797) q[0];
sx q[0];
rz(-1.0562684) q[0];
sx q[0];
rz(-2.5637676) q[0];
rz(-pi) q[1];
rz(-3.0566759) q[2];
sx q[2];
rz(-1.6667637) q[2];
sx q[2];
rz(1.3451479) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.93219705) q[1];
sx q[1];
rz(-1.160566) q[1];
sx q[1];
rz(2.5305629) q[1];
rz(-pi) q[2];
rz(2.1649394) q[3];
sx q[3];
rz(-2.6043731) q[3];
sx q[3];
rz(-1.9013167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.962062) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(2.8774101) q[2];
rz(-1.8386486) q[3];
sx q[3];
rz(-1.3774201) q[3];
sx q[3];
rz(-2.0343659) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1846979) q[0];
sx q[0];
rz(-0.95320025) q[0];
sx q[0];
rz(-1.6860929) q[0];
rz(-2.1799344) q[1];
sx q[1];
rz(-1.3788297) q[1];
sx q[1];
rz(0.89967322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6379999) q[0];
sx q[0];
rz(-1.3595306) q[0];
sx q[0];
rz(-0.46163033) q[0];
rz(0.2238621) q[2];
sx q[2];
rz(-2.2831585) q[2];
sx q[2];
rz(1.7162042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65708651) q[1];
sx q[1];
rz(-2.5528347) q[1];
sx q[1];
rz(1.4935054) q[1];
rz(-pi) q[2];
rz(1.7877696) q[3];
sx q[3];
rz(-1.5343175) q[3];
sx q[3];
rz(0.36408261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4453033) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(2.3882833) q[2];
rz(-2.8478012) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(-2.0297348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927354) q[0];
sx q[0];
rz(-0.98965544) q[0];
sx q[0];
rz(-0.85187546) q[0];
rz(1.4082255) q[1];
sx q[1];
rz(-1.4829166) q[1];
sx q[1];
rz(-0.3826938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9155884) q[0];
sx q[0];
rz(-1.8126285) q[0];
sx q[0];
rz(-2.7126724) q[0];
rz(1.3415975) q[2];
sx q[2];
rz(-2.3275073) q[2];
sx q[2];
rz(2.2687721) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1074658) q[1];
sx q[1];
rz(-1.4728947) q[1];
sx q[1];
rz(2.8393954) q[1];
rz(-0.15573536) q[3];
sx q[3];
rz(-2.6504575) q[3];
sx q[3];
rz(-2.8409426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7644299) q[2];
sx q[2];
rz(-0.96401507) q[2];
sx q[2];
rz(0.60297472) q[2];
rz(-1.0682586) q[3];
sx q[3];
rz(-1.4385185) q[3];
sx q[3];
rz(-0.8409797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0076440796) q[0];
sx q[0];
rz(-1.6772062) q[0];
sx q[0];
rz(0.35368791) q[0];
rz(-1.1324646) q[1];
sx q[1];
rz(-2.1734889) q[1];
sx q[1];
rz(-2.7521334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6521527) q[0];
sx q[0];
rz(-1.3394233) q[0];
sx q[0];
rz(-2.6315297) q[0];
rz(-pi) q[1];
rz(-1.9612938) q[2];
sx q[2];
rz(-0.26528851) q[2];
sx q[2];
rz(-0.92084322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.25915256) q[1];
sx q[1];
rz(-2.5969567) q[1];
sx q[1];
rz(-0.55014054) q[1];
rz(-pi) q[2];
rz(-0.52822379) q[3];
sx q[3];
rz(-1.6162795) q[3];
sx q[3];
rz(2.8128565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0103717) q[2];
sx q[2];
rz(-1.6348569) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(1.5366588) q[3];
sx q[3];
rz(-2.6585572) q[3];
sx q[3];
rz(-0.35227942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1491886) q[0];
sx q[0];
rz(-2.438899) q[0];
sx q[0];
rz(-2.6371523) q[0];
rz(-3.0814677) q[1];
sx q[1];
rz(-2.6838214) q[1];
sx q[1];
rz(2.6837742) q[1];
rz(3.0679476) q[2];
sx q[2];
rz(-1.4630058) q[2];
sx q[2];
rz(-1.502682) q[2];
rz(1.6144013) q[3];
sx q[3];
rz(-0.674896) q[3];
sx q[3];
rz(1.9622635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
