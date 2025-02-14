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
rz(0.54230827) q[0];
sx q[0];
rz(3.0071654) q[0];
sx q[0];
rz(10.472049) q[0];
rz(-0.32416999) q[1];
sx q[1];
rz(-0.14961641) q[1];
sx q[1];
rz(0.9447929) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0243264) q[0];
sx q[0];
rz(-0.42792441) q[0];
sx q[0];
rz(-0.86366349) q[0];
x q[1];
rz(2.9080897) q[2];
sx q[2];
rz(-0.93793195) q[2];
sx q[2];
rz(1.3695182) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8103247) q[1];
sx q[1];
rz(-2.4400418) q[1];
sx q[1];
rz(-0.74358209) q[1];
x q[2];
rz(2.0534001) q[3];
sx q[3];
rz(-1.4632097) q[3];
sx q[3];
rz(-1.3914598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2904539) q[2];
sx q[2];
rz(-1.6893427) q[2];
sx q[2];
rz(-1.5770844) q[2];
rz(-2.5751298) q[3];
sx q[3];
rz(-2.5201859) q[3];
sx q[3];
rz(1.8433146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8241149) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(2.4216006) q[0];
rz(0.27436817) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(2.5879587) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.276727) q[0];
sx q[0];
rz(-1.454852) q[0];
sx q[0];
rz(1.9254294) q[0];
rz(-2.738094) q[2];
sx q[2];
rz(-1.8581693) q[2];
sx q[2];
rz(2.0593875) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1914042) q[1];
sx q[1];
rz(-2.2850415) q[1];
sx q[1];
rz(-2.6493303) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6194044) q[3];
sx q[3];
rz(-2.0117674) q[3];
sx q[3];
rz(0.14312927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1729892) q[2];
sx q[2];
rz(-1.2022688) q[2];
sx q[2];
rz(0.36390641) q[2];
rz(0.11161741) q[3];
sx q[3];
rz(-0.93897396) q[3];
sx q[3];
rz(-2.3811471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97742057) q[0];
sx q[0];
rz(-2.3168679) q[0];
sx q[0];
rz(-3.136193) q[0];
rz(0.81575704) q[1];
sx q[1];
rz(-1.1382256) q[1];
sx q[1];
rz(-2.7556748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2337795) q[0];
sx q[0];
rz(-2.695437) q[0];
sx q[0];
rz(2.4491007) q[0];
rz(-pi) q[1];
rz(2.7806769) q[2];
sx q[2];
rz(-1.0385761) q[2];
sx q[2];
rz(0.82718933) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5623916) q[1];
sx q[1];
rz(-1.1913956) q[1];
sx q[1];
rz(-3.0048278) q[1];
x q[2];
rz(3.0073062) q[3];
sx q[3];
rz(-1.1597753) q[3];
sx q[3];
rz(0.81601303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58891121) q[2];
sx q[2];
rz(-2.3083355) q[2];
sx q[2];
rz(-1.3820648) q[2];
rz(2.9803993) q[3];
sx q[3];
rz(-1.7101945) q[3];
sx q[3];
rz(0.62895697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2564119) q[0];
sx q[0];
rz(-1.2762524) q[0];
sx q[0];
rz(1.8026344) q[0];
rz(-1.8244052) q[1];
sx q[1];
rz(-0.97511292) q[1];
sx q[1];
rz(-0.29979527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9343963) q[0];
sx q[0];
rz(-2.3863433) q[0];
sx q[0];
rz(2.62225) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1607397) q[2];
sx q[2];
rz(-1.9741657) q[2];
sx q[2];
rz(1.174675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.055017415) q[1];
sx q[1];
rz(-2.7624532) q[1];
sx q[1];
rz(1.192548) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13238975) q[3];
sx q[3];
rz(-0.93919884) q[3];
sx q[3];
rz(-1.6211444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51912159) q[2];
sx q[2];
rz(-2.2032479) q[2];
sx q[2];
rz(2.8224714) q[2];
rz(-0.090911344) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(-1.7849785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153445) q[0];
sx q[0];
rz(-0.81979668) q[0];
sx q[0];
rz(-0.0023284624) q[0];
rz(1.5709411) q[1];
sx q[1];
rz(-0.80051533) q[1];
sx q[1];
rz(-1.194582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1462458) q[0];
sx q[0];
rz(-1.4527391) q[0];
sx q[0];
rz(0.15583584) q[0];
rz(-1.6171239) q[2];
sx q[2];
rz(-2.4232658) q[2];
sx q[2];
rz(-2.3026932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1135322) q[1];
sx q[1];
rz(-0.36021566) q[1];
sx q[1];
rz(-1.8647782) q[1];
rz(-pi) q[2];
rz(0.43554775) q[3];
sx q[3];
rz(-0.60127944) q[3];
sx q[3];
rz(-1.5454436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2583367) q[2];
sx q[2];
rz(-0.77815762) q[2];
sx q[2];
rz(3.0641595) q[2];
rz(0.94610131) q[3];
sx q[3];
rz(-0.48282048) q[3];
sx q[3];
rz(2.6879123) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2133863) q[0];
sx q[0];
rz(-0.086406924) q[0];
sx q[0];
rz(1.7267831) q[0];
rz(-2.4746223) q[1];
sx q[1];
rz(-1.357115) q[1];
sx q[1];
rz(-0.40474969) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0721829) q[0];
sx q[0];
rz(-1.8526005) q[0];
sx q[0];
rz(2.826716) q[0];
rz(-pi) q[1];
rz(-2.8799492) q[2];
sx q[2];
rz(-0.930013) q[2];
sx q[2];
rz(-1.4705758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.09913832) q[1];
sx q[1];
rz(-2.6746866) q[1];
sx q[1];
rz(-1.7068638) q[1];
rz(-pi) q[2];
rz(2.1920106) q[3];
sx q[3];
rz(-0.57692617) q[3];
sx q[3];
rz(2.2872137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82379597) q[2];
sx q[2];
rz(-2.3522289) q[2];
sx q[2];
rz(0.44343597) q[2];
rz(-0.41380841) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(0.88576353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7498748) q[0];
sx q[0];
rz(-1.7789919) q[0];
sx q[0];
rz(-0.66746563) q[0];
rz(-0.04774566) q[1];
sx q[1];
rz(-1.1404488) q[1];
sx q[1];
rz(-2.023229) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76343173) q[0];
sx q[0];
rz(-1.973828) q[0];
sx q[0];
rz(2.5600938) q[0];
x q[1];
rz(2.1492858) q[2];
sx q[2];
rz(-1.6725146) q[2];
sx q[2];
rz(-2.5755239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1350007) q[1];
sx q[1];
rz(-1.9891095) q[1];
sx q[1];
rz(-0.63726421) q[1];
rz(-pi) q[2];
rz(0.067906109) q[3];
sx q[3];
rz(-2.0693681) q[3];
sx q[3];
rz(2.8907421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.045791693) q[2];
sx q[2];
rz(-1.8818776) q[2];
sx q[2];
rz(-0.0087139159) q[2];
rz(1.9376612) q[3];
sx q[3];
rz(-2.7976076) q[3];
sx q[3];
rz(0.022139942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1190204) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(0.88985306) q[0];
rz(1.285137) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(0.13592517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724919) q[0];
sx q[0];
rz(-1.5234424) q[0];
sx q[0];
rz(0.22914178) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1077584) q[2];
sx q[2];
rz(-1.0091678) q[2];
sx q[2];
rz(-1.461535) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4299851) q[1];
sx q[1];
rz(-1.9031591) q[1];
sx q[1];
rz(-1.3620939) q[1];
rz(-pi) q[2];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.7673878) q[3];
sx q[3];
rz(-2.1458613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70302427) q[2];
sx q[2];
rz(-0.88006222) q[2];
sx q[2];
rz(-1.0620037) q[2];
rz(0.92987531) q[3];
sx q[3];
rz(-0.78571856) q[3];
sx q[3];
rz(2.3532531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7787665) q[0];
sx q[0];
rz(-0.39149785) q[0];
sx q[0];
rz(2.7407001) q[0];
rz(-2.3800384) q[1];
sx q[1];
rz(-1.6490033) q[1];
sx q[1];
rz(0.78899312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3991206) q[0];
sx q[0];
rz(-1.6982887) q[0];
sx q[0];
rz(-1.5869058) q[0];
x q[1];
rz(-0.15536552) q[2];
sx q[2];
rz(-0.90165686) q[2];
sx q[2];
rz(0.63442996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9826451) q[1];
sx q[1];
rz(-0.98474307) q[1];
sx q[1];
rz(1.7141377) q[1];
x q[2];
rz(-1.5212359) q[3];
sx q[3];
rz(-1.0465099) q[3];
sx q[3];
rz(1.6848642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57364982) q[2];
sx q[2];
rz(-1.2929792) q[2];
sx q[2];
rz(-2.4207777) q[2];
rz(-1.6503664) q[3];
sx q[3];
rz(-2.3590915) q[3];
sx q[3];
rz(1.3097552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023130527) q[0];
sx q[0];
rz(-3.0266302) q[0];
sx q[0];
rz(-2.3638828) q[0];
rz(-2.481781) q[1];
sx q[1];
rz(-2.2751364) q[1];
sx q[1];
rz(0.39001098) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84108739) q[0];
sx q[0];
rz(-1.705637) q[0];
sx q[0];
rz(3.10411) q[0];
rz(-pi) q[1];
rz(-0.89945515) q[2];
sx q[2];
rz(-0.94591606) q[2];
sx q[2];
rz(2.867183) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6811878) q[1];
sx q[1];
rz(-1.8680267) q[1];
sx q[1];
rz(1.1497208) q[1];
x q[2];
rz(-2.635119) q[3];
sx q[3];
rz(-2.1586543) q[3];
sx q[3];
rz(-1.7755058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68221349) q[2];
sx q[2];
rz(-0.43872681) q[2];
sx q[2];
rz(0.48975804) q[2];
rz(-0.30719906) q[3];
sx q[3];
rz(-2.2178853) q[3];
sx q[3];
rz(-1.7507621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3601892) q[0];
sx q[0];
rz(-1.4908981) q[0];
sx q[0];
rz(1.4733462) q[0];
rz(-1.0745984) q[1];
sx q[1];
rz(-0.85292024) q[1];
sx q[1];
rz(-1.9180752) q[1];
rz(0.75847738) q[2];
sx q[2];
rz(-0.52634326) q[2];
sx q[2];
rz(1.1545622) q[2];
rz(-0.73978872) q[3];
sx q[3];
rz(-2.3935912) q[3];
sx q[3];
rz(-1.145515) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
