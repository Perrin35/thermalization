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
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6990841) q[0];
sx q[0];
rz(-2.6538349) q[0];
sx q[0];
rz(0.91983025) q[0];
rz(0.046836179) q[2];
sx q[2];
rz(-1.9745262) q[2];
sx q[2];
rz(2.9138034) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85744845) q[1];
sx q[1];
rz(-1.3054056) q[1];
sx q[1];
rz(1.1737203) q[1];
rz(-pi) q[2];
rz(1.2143308) q[3];
sx q[3];
rz(-0.89205304) q[3];
sx q[3];
rz(-1.9744622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5758489) q[2];
sx q[2];
rz(-2.1243024) q[2];
sx q[2];
rz(-0.033509342) q[2];
rz(1.9401898) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(-2.0638154) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1101542) q[0];
sx q[0];
rz(-2.8903676) q[0];
sx q[0];
rz(0.95397368) q[0];
rz(1.4454449) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(-1.9704069) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508043) q[0];
sx q[0];
rz(-2.224824) q[0];
sx q[0];
rz(0.79933856) q[0];
rz(-pi) q[1];
rz(-1.2759802) q[2];
sx q[2];
rz(-0.95446842) q[2];
sx q[2];
rz(-2.4957239) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3747159) q[1];
sx q[1];
rz(-1.3424113) q[1];
sx q[1];
rz(1.5414052) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33447845) q[3];
sx q[3];
rz(-1.1568767) q[3];
sx q[3];
rz(2.778307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90616068) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(1.99235) q[2];
rz(-0.033128459) q[3];
sx q[3];
rz(-1.6413611) q[3];
sx q[3];
rz(-0.59605014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0572492) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(1.0302011) q[0];
rz(1.239981) q[1];
sx q[1];
rz(-1.7844618) q[1];
sx q[1];
rz(2.1485567) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1112615) q[0];
sx q[0];
rz(-1.5227888) q[0];
sx q[0];
rz(-0.45560776) q[0];
x q[1];
rz(0.15056653) q[2];
sx q[2];
rz(-0.67306256) q[2];
sx q[2];
rz(1.3526431) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60336411) q[1];
sx q[1];
rz(-2.1346655) q[1];
sx q[1];
rz(-0.67810525) q[1];
rz(1.1458323) q[3];
sx q[3];
rz(-0.83163762) q[3];
sx q[3];
rz(-1.1475565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9242636) q[2];
sx q[2];
rz(-2.0254878) q[2];
sx q[2];
rz(-2.7868311) q[2];
rz(-2.1445856) q[3];
sx q[3];
rz(-1.487178) q[3];
sx q[3];
rz(-2.2009489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9724801) q[0];
sx q[0];
rz(-1.4262154) q[0];
sx q[0];
rz(-2.0347563) q[0];
rz(2.2726982) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(-0.69721627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91928673) q[0];
sx q[0];
rz(-2.5450052) q[0];
sx q[0];
rz(-0.35937341) q[0];
x q[1];
rz(-0.72202335) q[2];
sx q[2];
rz(-1.486384) q[2];
sx q[2];
rz(-1.9868324) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7719318) q[1];
sx q[1];
rz(-0.54851739) q[1];
sx q[1];
rz(2.3322075) q[1];
rz(-2.3512164) q[3];
sx q[3];
rz(-0.034878313) q[3];
sx q[3];
rz(0.45524516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8640459) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-0.29402688) q[2];
rz(2.0781519) q[3];
sx q[3];
rz(-1.682351) q[3];
sx q[3];
rz(-2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251637) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(0.46318769) q[0];
rz(-0.40395346) q[1];
sx q[1];
rz(-1.5084167) q[1];
sx q[1];
rz(-2.892866) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489784) q[0];
sx q[0];
rz(-0.19751829) q[0];
sx q[0];
rz(-1.5794483) q[0];
x q[1];
rz(-1.8291446) q[2];
sx q[2];
rz(-1.0444006) q[2];
sx q[2];
rz(1.0963944) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46649964) q[1];
sx q[1];
rz(-2.6394301) q[1];
sx q[1];
rz(2.2999022) q[1];
rz(-pi) q[2];
rz(-2.4234613) q[3];
sx q[3];
rz(-0.87801027) q[3];
sx q[3];
rz(-0.75961514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2948461) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(-1.3999636) q[2];
rz(-1.8585662) q[3];
sx q[3];
rz(-1.3834407) q[3];
sx q[3];
rz(-0.2259026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1840709) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(-2.7958909) q[0];
rz(-0.050994571) q[1];
sx q[1];
rz(-1.9146999) q[1];
sx q[1];
rz(2.3695703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22020082) q[0];
sx q[0];
rz(-1.3190375) q[0];
sx q[0];
rz(2.779106) q[0];
rz(-pi) q[1];
rz(-1.4441024) q[2];
sx q[2];
rz(-2.4856353) q[2];
sx q[2];
rz(-0.7743338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1322546) q[1];
sx q[1];
rz(-1.8937832) q[1];
sx q[1];
rz(-2.1732974) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28874044) q[3];
sx q[3];
rz(-2.0802698) q[3];
sx q[3];
rz(-2.4276707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(-2.9841606) q[2];
rz(0.43846798) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(2.1854775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0601783) q[0];
sx q[0];
rz(-2.0745451) q[0];
sx q[0];
rz(2.881158) q[0];
rz(-0.57834894) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-2.9885805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3653202) q[0];
sx q[0];
rz(-2.0663102) q[0];
sx q[0];
rz(0.97712626) q[0];
rz(3.0566759) q[2];
sx q[2];
rz(-1.474829) q[2];
sx q[2];
rz(1.3451479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2305829) q[1];
sx q[1];
rz(-1.0167767) q[1];
sx q[1];
rz(2.0589214) q[1];
rz(-pi) q[2];
rz(2.1649394) q[3];
sx q[3];
rz(-2.6043731) q[3];
sx q[3];
rz(1.2402759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.962062) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(2.8774101) q[2];
rz(-1.8386486) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(-1.1072268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1846979) q[0];
sx q[0];
rz(-2.1883924) q[0];
sx q[0];
rz(-1.6860929) q[0];
rz(-0.9616583) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(-2.2419194) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1048516) q[0];
sx q[0];
rz(-1.1201968) q[0];
sx q[0];
rz(-1.3356871) q[0];
x q[1];
rz(-1.3191965) q[2];
sx q[2];
rz(-0.74079047) q[2];
sx q[2];
rz(2.0513926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3916495) q[1];
sx q[1];
rz(-2.1575621) q[1];
sx q[1];
rz(3.090078) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4028655) q[3];
sx q[3];
rz(-2.9216218) q[3];
sx q[3];
rz(-2.0988362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6962894) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(2.3882833) q[2];
rz(2.8478012) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(-1.1118579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927354) q[0];
sx q[0];
rz(-2.1519372) q[0];
sx q[0];
rz(0.85187546) q[0];
rz(1.4082255) q[1];
sx q[1];
rz(-1.6586761) q[1];
sx q[1];
rz(-2.7588989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2792042) q[0];
sx q[0];
rz(-0.48868256) q[0];
sx q[0];
rz(-0.53532289) q[0];
rz(-pi) q[1];
rz(1.7999951) q[2];
sx q[2];
rz(-2.3275073) q[2];
sx q[2];
rz(0.87282055) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.574461) q[1];
sx q[1];
rz(-1.2700915) q[1];
sx q[1];
rz(-1.4682795) q[1];
rz(-pi) q[2];
rz(-2.6555041) q[3];
sx q[3];
rz(-1.6440142) q[3];
sx q[3];
rz(-1.7338727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7644299) q[2];
sx q[2];
rz(-2.1775776) q[2];
sx q[2];
rz(-0.60297472) q[2];
rz(-2.073334) q[3];
sx q[3];
rz(-1.4385185) q[3];
sx q[3];
rz(-2.300613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0076440796) q[0];
sx q[0];
rz(-1.4643865) q[0];
sx q[0];
rz(2.7879047) q[0];
rz(-1.1324646) q[1];
sx q[1];
rz(-2.1734889) q[1];
sx q[1];
rz(-2.7521334) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48943993) q[0];
sx q[0];
rz(-1.8021694) q[0];
sx q[0];
rz(0.51006298) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0385397) q[2];
sx q[2];
rz(-1.8156869) q[2];
sx q[2];
rz(-2.6239397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25915256) q[1];
sx q[1];
rz(-0.54463592) q[1];
sx q[1];
rz(0.55014054) q[1];
rz(-0.52822379) q[3];
sx q[3];
rz(-1.6162795) q[3];
sx q[3];
rz(-0.32873617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13122095) q[2];
sx q[2];
rz(-1.6348569) q[2];
sx q[2];
rz(1.1466675) q[2];
rz(1.5366588) q[3];
sx q[3];
rz(-0.48303548) q[3];
sx q[3];
rz(-2.7893132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99240408) q[0];
sx q[0];
rz(-0.70269361) q[0];
sx q[0];
rz(0.50444034) q[0];
rz(-3.0814677) q[1];
sx q[1];
rz(-2.6838214) q[1];
sx q[1];
rz(2.6837742) q[1];
rz(1.4627152) q[2];
sx q[2];
rz(-1.6440132) q[2];
sx q[2];
rz(-3.0655412) q[2];
rz(-1.6144013) q[3];
sx q[3];
rz(-2.4666967) q[3];
sx q[3];
rz(-1.1793292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
