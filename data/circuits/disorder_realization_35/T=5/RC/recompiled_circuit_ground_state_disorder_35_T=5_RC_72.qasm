OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(-1.2195725) q[0];
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(0.9019444) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.886698) q[0];
sx q[0];
rz(-2.5097846) q[0];
sx q[0];
rz(-2.8796701) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7576394) q[2];
sx q[2];
rz(-1.9451491) q[2];
sx q[2];
rz(2.6461305) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8512501) q[1];
sx q[1];
rz(-2.800673) q[1];
sx q[1];
rz(-2.1147637) q[1];
x q[2];
rz(-0.73027421) q[3];
sx q[3];
rz(-1.0515107) q[3];
sx q[3];
rz(1.9744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1787662) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(2.315305) q[2];
rz(-1.4055584) q[3];
sx q[3];
rz(-0.30705753) q[3];
sx q[3];
rz(-0.74339408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58251441) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(1.7412809) q[0];
rz(-1.1236313) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(-2.0434911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207635) q[0];
sx q[0];
rz(-1.6801103) q[0];
sx q[0];
rz(2.9846322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2400572) q[2];
sx q[2];
rz(-1.8914701) q[2];
sx q[2];
rz(1.641524) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2929165) q[1];
sx q[1];
rz(-1.6458938) q[1];
sx q[1];
rz(-0.1212705) q[1];
rz(-pi) q[2];
rz(1.9722749) q[3];
sx q[3];
rz(-2.7352834) q[3];
sx q[3];
rz(0.10029785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0850247) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(-1.2497831) q[2];
rz(-1.362484) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(1.484163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054841) q[0];
sx q[0];
rz(-0.69378575) q[0];
sx q[0];
rz(-0.27534494) q[0];
rz(0.45922008) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(-3.088248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33063525) q[0];
sx q[0];
rz(-2.914408) q[0];
sx q[0];
rz(-0.35463984) q[0];
rz(1.533639) q[2];
sx q[2];
rz(-0.46122631) q[2];
sx q[2];
rz(-0.30468582) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1054105) q[1];
sx q[1];
rz(-2.1222669) q[1];
sx q[1];
rz(1.5230194) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2417913) q[3];
sx q[3];
rz(-1.7681554) q[3];
sx q[3];
rz(-2.1579735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3205545) q[2];
sx q[2];
rz(-2.7802763) q[2];
sx q[2];
rz(-1.6374755) q[2];
rz(1.641168) q[3];
sx q[3];
rz(-1.0502366) q[3];
sx q[3];
rz(2.8659081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7212873) q[0];
sx q[0];
rz(-2.5581701) q[0];
sx q[0];
rz(-0.48459184) q[0];
rz(-1.2424319) q[1];
sx q[1];
rz(-0.74598765) q[1];
sx q[1];
rz(0.34908435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1734245) q[0];
sx q[0];
rz(-1.5281046) q[0];
sx q[0];
rz(1.0220362) q[0];
rz(-pi) q[1];
rz(3.066272) q[2];
sx q[2];
rz(-1.2284797) q[2];
sx q[2];
rz(1.7750193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1368833) q[1];
sx q[1];
rz(-1.7994969) q[1];
sx q[1];
rz(-1.7654224) q[1];
rz(-pi) q[2];
rz(-2.146017) q[3];
sx q[3];
rz(-0.68505008) q[3];
sx q[3];
rz(-0.12888651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35859534) q[2];
sx q[2];
rz(-1.2025669) q[2];
sx q[2];
rz(-0.45004582) q[2];
rz(-2.3027244) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(-0.1990327) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49804509) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(0.099844649) q[0];
rz(1.8210583) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(2.5340396) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927101) q[0];
sx q[0];
rz(-1.6537871) q[0];
sx q[0];
rz(-1.6458428) q[0];
x q[1];
rz(-1.4133006) q[2];
sx q[2];
rz(-2.11065) q[2];
sx q[2];
rz(-1.6423722) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6441356) q[1];
sx q[1];
rz(-1.0007528) q[1];
sx q[1];
rz(1.8291147) q[1];
rz(-pi) q[2];
rz(-1.2149548) q[3];
sx q[3];
rz(-0.27493335) q[3];
sx q[3];
rz(2.7235018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1943835) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-2.6109931) q[2];
rz(0.44140205) q[3];
sx q[3];
rz(-0.97359052) q[3];
sx q[3];
rz(0.94943625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7049103) q[0];
sx q[0];
rz(-1.3446151) q[0];
sx q[0];
rz(-2.5675024) q[0];
rz(-2.2615945) q[1];
sx q[1];
rz(-2.0206385) q[1];
sx q[1];
rz(1.3425739) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3898252) q[0];
sx q[0];
rz(-0.51503599) q[0];
sx q[0];
rz(-1.9866605) q[0];
rz(-pi) q[1];
rz(1.1583352) q[2];
sx q[2];
rz(-1.165373) q[2];
sx q[2];
rz(-2.6087922) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3664535) q[1];
sx q[1];
rz(-0.43700686) q[1];
sx q[1];
rz(-0.98458146) q[1];
x q[2];
rz(-0.66304147) q[3];
sx q[3];
rz(-1.6812857) q[3];
sx q[3];
rz(2.7821531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67419702) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(-0.44710844) q[2];
rz(2.1111264) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7717188) q[0];
sx q[0];
rz(-1.8444703) q[0];
sx q[0];
rz(0.41279992) q[0];
rz(-1.075047) q[1];
sx q[1];
rz(-1.1378891) q[1];
sx q[1];
rz(2.3379751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5301054) q[0];
sx q[0];
rz(-0.75615935) q[0];
sx q[0];
rz(0.79571457) q[0];
rz(-pi) q[1];
rz(-2.2718524) q[2];
sx q[2];
rz(-0.77574965) q[2];
sx q[2];
rz(2.9091331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0405443) q[1];
sx q[1];
rz(-0.78436995) q[1];
sx q[1];
rz(3.1215828) q[1];
x q[2];
rz(-0.17937029) q[3];
sx q[3];
rz(-1.495788) q[3];
sx q[3];
rz(-1.4574432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86218086) q[2];
sx q[2];
rz(-1.923424) q[2];
sx q[2];
rz(-1.6406406) q[2];
rz(0.62266478) q[3];
sx q[3];
rz(-2.3707135) q[3];
sx q[3];
rz(-1.6525432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.415446) q[0];
sx q[0];
rz(-1.5203238) q[0];
sx q[0];
rz(-2.5836482) q[0];
rz(-0.56124148) q[1];
sx q[1];
rz(-1.7702425) q[1];
sx q[1];
rz(-1.4161313) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5197882) q[0];
sx q[0];
rz(-1.1237696) q[0];
sx q[0];
rz(2.8267953) q[0];
rz(-pi) q[1];
rz(-0.92500706) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(-3.1104308) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0603283) q[1];
sx q[1];
rz(-0.34354478) q[1];
sx q[1];
rz(1.8650123) q[1];
rz(-pi) q[2];
rz(-0.75545488) q[3];
sx q[3];
rz(-1.7208365) q[3];
sx q[3];
rz(-2.9767075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2531565) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-2.1412795) q[2];
rz(3.0194164) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(2.9848671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4633453) q[0];
sx q[0];
rz(-1.954701) q[0];
sx q[0];
rz(-0.004322411) q[0];
rz(-0.16383544) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(-1.77553) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28986606) q[0];
sx q[0];
rz(-2.4938739) q[0];
sx q[0];
rz(-1.6169548) q[0];
rz(-pi) q[1];
rz(2.4539656) q[2];
sx q[2];
rz(-0.4677217) q[2];
sx q[2];
rz(-0.91910884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8722265) q[1];
sx q[1];
rz(-2.8878643) q[1];
sx q[1];
rz(-0.14464186) q[1];
x q[2];
rz(2.4716464) q[3];
sx q[3];
rz(-0.94196999) q[3];
sx q[3];
rz(2.3184702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(-2.5223993) q[2];
rz(-2.5891417) q[3];
sx q[3];
rz(-1.3055472) q[3];
sx q[3];
rz(-0.14510554) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904952) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(2.6628394) q[0];
rz(-1.9888196) q[1];
sx q[1];
rz(-2.8515127) q[1];
sx q[1];
rz(-2.1943888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21175948) q[0];
sx q[0];
rz(-0.18778983) q[0];
sx q[0];
rz(1.7540356) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98299111) q[2];
sx q[2];
rz(-1.3653367) q[2];
sx q[2];
rz(-3.0794155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45631623) q[1];
sx q[1];
rz(-2.3446977) q[1];
sx q[1];
rz(2.0048098) q[1];
rz(-pi) q[2];
rz(-2.5371505) q[3];
sx q[3];
rz(-0.42057188) q[3];
sx q[3];
rz(-2.4673691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3843711) q[2];
sx q[2];
rz(-2.9360076) q[2];
sx q[2];
rz(-2.8749386) q[2];
rz(-1.0673149) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(-2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08854475) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(-2.0883941) q[1];
sx q[1];
rz(-1.9269301) q[1];
sx q[1];
rz(-2.3763837) q[1];
rz(0.24198738) q[2];
sx q[2];
rz(-1.568071) q[2];
sx q[2];
rz(-0.65721401) q[2];
rz(-2.8070634) q[3];
sx q[3];
rz(-2.5361037) q[3];
sx q[3];
rz(-1.1246214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
