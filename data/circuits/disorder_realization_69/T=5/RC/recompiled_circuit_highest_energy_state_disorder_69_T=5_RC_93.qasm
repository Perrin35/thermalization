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
rz(-2.7954571) q[0];
sx q[0];
rz(-1.6383642) q[0];
sx q[0];
rz(-0.7315973) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(2.7076758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.159891) q[0];
sx q[0];
rz(-0.74112391) q[0];
sx q[0];
rz(1.6228049) q[0];
rz(-pi) q[1];
rz(0.97831867) q[2];
sx q[2];
rz(-1.0772395) q[2];
sx q[2];
rz(1.6101642) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9508236) q[1];
sx q[1];
rz(-0.94591037) q[1];
sx q[1];
rz(-0.48887078) q[1];
rz(-pi) q[2];
rz(1.7619261) q[3];
sx q[3];
rz(-1.6014631) q[3];
sx q[3];
rz(0.14740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4291541) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(-1.0038556) q[2];
rz(2.6058274) q[3];
sx q[3];
rz(-0.86447132) q[3];
sx q[3];
rz(-0.0011477688) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67343229) q[0];
sx q[0];
rz(-2.4591481) q[0];
sx q[0];
rz(-2.4144507) q[0];
rz(3.0739821) q[1];
sx q[1];
rz(-1.7865684) q[1];
sx q[1];
rz(1.1236069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0824995) q[0];
sx q[0];
rz(-1.8525665) q[0];
sx q[0];
rz(-0.49651329) q[0];
rz(-pi) q[1];
rz(1.4482165) q[2];
sx q[2];
rz(-1.5059587) q[2];
sx q[2];
rz(-1.4139869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.935168) q[1];
sx q[1];
rz(-1.670253) q[1];
sx q[1];
rz(1.2897435) q[1];
rz(-pi) q[2];
rz(1.0960078) q[3];
sx q[3];
rz(-1.1752306) q[3];
sx q[3];
rz(2.6717348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2850538) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(0.48173586) q[2];
rz(2.5887865) q[3];
sx q[3];
rz(-2.063664) q[3];
sx q[3];
rz(3.0520181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8168617) q[0];
sx q[0];
rz(-2.6664) q[0];
sx q[0];
rz(-0.09356308) q[0];
rz(1.7612673) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(-0.63724744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68604159) q[0];
sx q[0];
rz(-2.0860057) q[0];
sx q[0];
rz(-2.0705219) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31934505) q[2];
sx q[2];
rz(-0.63268662) q[2];
sx q[2];
rz(-1.3038532) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.6031076) q[1];
sx q[1];
rz(-2.7432495) q[1];
sx q[1];
rz(-1.6377691) q[1];
rz(-pi) q[2];
rz(-2.9985688) q[3];
sx q[3];
rz(-0.29353729) q[3];
sx q[3];
rz(2.0484201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20643413) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(-2.4578102) q[2];
rz(0.23009662) q[3];
sx q[3];
rz(-1.4068539) q[3];
sx q[3];
rz(-0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0919331) q[0];
sx q[0];
rz(-0.50685087) q[0];
sx q[0];
rz(1.7012713) q[0];
rz(2.6662042) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(2.5604274) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064539486) q[0];
sx q[0];
rz(-0.7270027) q[0];
sx q[0];
rz(-1.4420271) q[0];
rz(-pi) q[1];
rz(2.0330795) q[2];
sx q[2];
rz(-1.7812742) q[2];
sx q[2];
rz(2.0049948) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5161203) q[1];
sx q[1];
rz(-0.56159084) q[1];
sx q[1];
rz(-2.549704) q[1];
x q[2];
rz(1.4122333) q[3];
sx q[3];
rz(-2.057029) q[3];
sx q[3];
rz(2.9610046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0541957) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(2.6206214) q[2];
rz(1.6446796) q[3];
sx q[3];
rz(-1.729634) q[3];
sx q[3];
rz(1.6269256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5882551) q[0];
sx q[0];
rz(-2.8659358) q[0];
sx q[0];
rz(2.0113373) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.1993473) q[1];
sx q[1];
rz(-1.1111396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6040627) q[0];
sx q[0];
rz(-2.0274693) q[0];
sx q[0];
rz(-2.9567493) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5264185) q[2];
sx q[2];
rz(-1.0414972) q[2];
sx q[2];
rz(0.30454301) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6592404) q[1];
sx q[1];
rz(-2.2063046) q[1];
sx q[1];
rz(-2.9457432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5938528) q[3];
sx q[3];
rz(-0.99036694) q[3];
sx q[3];
rz(-1.2810027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7378716) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(-2.1007288) q[2];
rz(-2.3216085) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(-2.5177054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89989221) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(-2.4851121) q[0];
rz(-1.7680602) q[1];
sx q[1];
rz(-2.3479925) q[1];
sx q[1];
rz(2.0097282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8478125) q[0];
sx q[0];
rz(-0.72350693) q[0];
sx q[0];
rz(0.30059867) q[0];
rz(-pi) q[1];
rz(-2.4651338) q[2];
sx q[2];
rz(-2.0148811) q[2];
sx q[2];
rz(-1.9537587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51323503) q[1];
sx q[1];
rz(-1.5984922) q[1];
sx q[1];
rz(0.85978924) q[1];
rz(-pi) q[2];
rz(-1.5751714) q[3];
sx q[3];
rz(-1.7212015) q[3];
sx q[3];
rz(0.37585092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6953096) q[2];
sx q[2];
rz(-0.60574564) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(-0.22932912) q[3];
sx q[3];
rz(-2.1181483) q[3];
sx q[3];
rz(-0.91013175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0685773) q[0];
sx q[0];
rz(-1.1443161) q[0];
sx q[0];
rz(-1.2314433) q[0];
rz(-3.0629509) q[1];
sx q[1];
rz(-1.4631203) q[1];
sx q[1];
rz(2.9873649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336727) q[0];
sx q[0];
rz(-1.3929954) q[0];
sx q[0];
rz(2.7975797) q[0];
rz(-pi) q[1];
rz(-2.8381589) q[2];
sx q[2];
rz(-1.9643485) q[2];
sx q[2];
rz(0.21912609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5104766) q[1];
sx q[1];
rz(-1.0923903) q[1];
sx q[1];
rz(1.1523184) q[1];
x q[2];
rz(-1.3580136) q[3];
sx q[3];
rz(-1.1682421) q[3];
sx q[3];
rz(-0.29197955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8057033) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(1.4914782) q[2];
rz(1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(1.570805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0989646) q[0];
sx q[0];
rz(-2.0476116) q[0];
sx q[0];
rz(1.6866823) q[0];
rz(0.97575724) q[1];
sx q[1];
rz(-2.0819596) q[1];
sx q[1];
rz(1.3386493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8737301) q[0];
sx q[0];
rz(-0.81276816) q[0];
sx q[0];
rz(1.0229179) q[0];
rz(2.1326124) q[2];
sx q[2];
rz(-0.57719165) q[2];
sx q[2];
rz(-1.4711589) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4563096) q[1];
sx q[1];
rz(-1.7830308) q[1];
sx q[1];
rz(0.93385277) q[1];
x q[2];
rz(-1.9429132) q[3];
sx q[3];
rz(-0.94148472) q[3];
sx q[3];
rz(-0.019364186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51155382) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(-2.3021728) q[2];
rz(-1.2342341) q[3];
sx q[3];
rz(-2.0525565) q[3];
sx q[3];
rz(-3.1101826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3473174) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(0.41912249) q[0];
rz(-1.4979111) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(0.43325123) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33746359) q[0];
sx q[0];
rz(-0.47352284) q[0];
sx q[0];
rz(1.3086523) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1243043) q[2];
sx q[2];
rz(-2.4845036) q[2];
sx q[2];
rz(-1.8985871) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.309714) q[1];
sx q[1];
rz(-2.0818748) q[1];
sx q[1];
rz(-0.32541497) q[1];
x q[2];
rz(2.8178704) q[3];
sx q[3];
rz(-0.60327134) q[3];
sx q[3];
rz(0.7654741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7822632) q[2];
sx q[2];
rz(-1.1507582) q[2];
sx q[2];
rz(-0.19598728) q[2];
rz(-1.1228784) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(1.4518729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48542431) q[0];
sx q[0];
rz(-0.66101414) q[0];
sx q[0];
rz(-2.0458903) q[0];
rz(0.51086673) q[1];
sx q[1];
rz(-0.26900649) q[1];
sx q[1];
rz(-2.9313472) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73186648) q[0];
sx q[0];
rz(-1.7667701) q[0];
sx q[0];
rz(0.98742296) q[0];
rz(-3.0265507) q[2];
sx q[2];
rz(-0.91918531) q[2];
sx q[2];
rz(2.5082805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0970891) q[1];
sx q[1];
rz(-1.0877481) q[1];
sx q[1];
rz(1.0853069) q[1];
x q[2];
rz(-2.5737567) q[3];
sx q[3];
rz(-2.46008) q[3];
sx q[3];
rz(2.5208254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11253396) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(1.0851592) q[2];
rz(0.30760136) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-0.65348452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6899684) q[0];
sx q[0];
rz(-1.1545447) q[0];
sx q[0];
rz(1.92323) q[0];
rz(0.65144173) q[1];
sx q[1];
rz(-2.1919498) q[1];
sx q[1];
rz(0.776074) q[1];
rz(-3.1262911) q[2];
sx q[2];
rz(-2.2616784) q[2];
sx q[2];
rz(-0.26605284) q[2];
rz(-2.2163556) q[3];
sx q[3];
rz(-2.400077) q[3];
sx q[3];
rz(-2.526068) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
