OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0924858) q[0];
sx q[0];
rz(-1.2993113) q[0];
sx q[0];
rz(-0.045825034) q[0];
rz(4.4412208) q[1];
sx q[1];
rz(4.9443359) q[1];
sx q[1];
rz(8.1716552) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9535256) q[0];
sx q[0];
rz(-1.3261516) q[0];
sx q[0];
rz(1.1172386) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9401993) q[2];
sx q[2];
rz(-1.5279576) q[2];
sx q[2];
rz(-0.34319718) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33838446) q[1];
sx q[1];
rz(-0.82370355) q[1];
sx q[1];
rz(1.4080774) q[1];
x q[2];
rz(-1.6982444) q[3];
sx q[3];
rz(-1.5747254) q[3];
sx q[3];
rz(1.1636101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8200298) q[2];
sx q[2];
rz(-1.4048615) q[2];
sx q[2];
rz(-2.826214) q[2];
rz(3.0553715) q[3];
sx q[3];
rz(-0.092970522) q[3];
sx q[3];
rz(0.84403795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4013937) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(-0.33541086) q[0];
rz(2.7566578) q[1];
sx q[1];
rz(-0.79610577) q[1];
sx q[1];
rz(-1.8523432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1803363) q[0];
sx q[0];
rz(-1.5085876) q[0];
sx q[0];
rz(0.38952413) q[0];
rz(-2.2122967) q[2];
sx q[2];
rz(-1.6469064) q[2];
sx q[2];
rz(-2.2141505) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7988879) q[1];
sx q[1];
rz(-1.2884166) q[1];
sx q[1];
rz(-2.7434231) q[1];
rz(-1.0689911) q[3];
sx q[3];
rz(-2.8659389) q[3];
sx q[3];
rz(0.75934809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5424767) q[2];
sx q[2];
rz(-1.7581538) q[2];
sx q[2];
rz(1.4625589) q[2];
rz(-0.8820495) q[3];
sx q[3];
rz(-2.4817011) q[3];
sx q[3];
rz(-1.5006458) q[3];
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
rz(-1.3682692) q[0];
sx q[0];
rz(-0.48650807) q[0];
sx q[0];
rz(-0.64273709) q[0];
rz(-2.2679988) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(-0.98683039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5499134) q[0];
sx q[0];
rz(-2.4351623) q[0];
sx q[0];
rz(1.3384992) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0796649) q[2];
sx q[2];
rz(-1.0871097) q[2];
sx q[2];
rz(0.40555039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9171608) q[1];
sx q[1];
rz(-0.81714918) q[1];
sx q[1];
rz(-0.78950531) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1961206) q[3];
sx q[3];
rz(-1.5244532) q[3];
sx q[3];
rz(0.99907334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4387536) q[2];
sx q[2];
rz(-0.025391014) q[2];
sx q[2];
rz(2.0849483) q[2];
rz(-1.7993401) q[3];
sx q[3];
rz(-1.4547576) q[3];
sx q[3];
rz(-2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037647) q[0];
sx q[0];
rz(-0.28452888) q[0];
sx q[0];
rz(1.9985265) q[0];
rz(-0.68483886) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(0.24982223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1764602) q[0];
sx q[0];
rz(-2.0270545) q[0];
sx q[0];
rz(-0.81911032) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1353829) q[2];
sx q[2];
rz(-0.57423019) q[2];
sx q[2];
rz(-1.8956313) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4433684) q[1];
sx q[1];
rz(-2.0495546) q[1];
sx q[1];
rz(2.4633292) q[1];
x q[2];
rz(-1.5644844) q[3];
sx q[3];
rz(-2.8996116) q[3];
sx q[3];
rz(0.93130563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41986856) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(1.8458337) q[2];
rz(-1.7660247) q[3];
sx q[3];
rz(-1.4294759) q[3];
sx q[3];
rz(-0.099055722) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7527723) q[0];
sx q[0];
rz(-1.7701912) q[0];
sx q[0];
rz(0.8465299) q[0];
rz(1.1835774) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(1.2394946) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13279937) q[0];
sx q[0];
rz(-1.6911048) q[0];
sx q[0];
rz(-1.8254721) q[0];
x q[1];
rz(-2.7991512) q[2];
sx q[2];
rz(-1.3991809) q[2];
sx q[2];
rz(1.1521074) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.073055) q[1];
sx q[1];
rz(-3.0843711) q[1];
sx q[1];
rz(1.7806609) q[1];
x q[2];
rz(-1.4828478) q[3];
sx q[3];
rz(-2.5426358) q[3];
sx q[3];
rz(-2.9943337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6401297) q[2];
sx q[2];
rz(-1.5093466) q[2];
sx q[2];
rz(2.787369) q[2];
rz(2.4223478) q[3];
sx q[3];
rz(-2.5651599) q[3];
sx q[3];
rz(2.332212) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3767913) q[0];
sx q[0];
rz(-0.23623315) q[0];
sx q[0];
rz(2.7948622) q[0];
rz(-2.4193343) q[1];
sx q[1];
rz(-1.6112695) q[1];
sx q[1];
rz(2.9647656) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2472025) q[0];
sx q[0];
rz(-2.6567334) q[0];
sx q[0];
rz(1.8762183) q[0];
rz(-pi) q[1];
rz(-3.0484285) q[2];
sx q[2];
rz(-1.0736862) q[2];
sx q[2];
rz(2.1456631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36784259) q[1];
sx q[1];
rz(-2.9108139) q[1];
sx q[1];
rz(2.6683067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7309639) q[3];
sx q[3];
rz(-2.5729542) q[3];
sx q[3];
rz(1.7344432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56883812) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(-1.0432165) q[2];
rz(0.15626945) q[3];
sx q[3];
rz(-2.0782491) q[3];
sx q[3];
rz(-3.0472896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909376) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(-1.2283196) q[0];
rz(-0.43371513) q[1];
sx q[1];
rz(-0.80077306) q[1];
sx q[1];
rz(1.0450276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.341808) q[0];
sx q[0];
rz(-1.5130454) q[0];
sx q[0];
rz(1.7736797) q[0];
x q[1];
rz(-2.3190193) q[2];
sx q[2];
rz(-1.9299091) q[2];
sx q[2];
rz(-0.99701478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1365388) q[1];
sx q[1];
rz(-2.1834259) q[1];
sx q[1];
rz(0.1392675) q[1];
rz(-pi) q[2];
rz(-0.77680334) q[3];
sx q[3];
rz(-1.9597199) q[3];
sx q[3];
rz(0.17059205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9016483) q[2];
sx q[2];
rz(-2.1176691) q[2];
sx q[2];
rz(-0.18590064) q[2];
rz(-1.9013885) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(-2.127229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0094725322) q[0];
sx q[0];
rz(-1.8931696) q[0];
sx q[0];
rz(-0.15952071) q[0];
rz(-0.80687833) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(-0.0085011403) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5497564) q[0];
sx q[0];
rz(-1.7248132) q[0];
sx q[0];
rz(1.1598253) q[0];
rz(-0.10430704) q[2];
sx q[2];
rz(-2.3416714) q[2];
sx q[2];
rz(1.4664354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2841683) q[1];
sx q[1];
rz(-1.4393974) q[1];
sx q[1];
rz(-0.22716503) q[1];
rz(-pi) q[2];
rz(1.0866585) q[3];
sx q[3];
rz(-1.7268506) q[3];
sx q[3];
rz(0.54359667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8622417) q[2];
sx q[2];
rz(-1.3631577) q[2];
sx q[2];
rz(0.19374338) q[2];
rz(-1.5396317) q[3];
sx q[3];
rz(-1.964147) q[3];
sx q[3];
rz(1.1295454) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047456) q[0];
sx q[0];
rz(-1.756825) q[0];
sx q[0];
rz(-2.4054476) q[0];
rz(2.2537117) q[1];
sx q[1];
rz(-2.6887951) q[1];
sx q[1];
rz(-0.050051659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5745478) q[0];
sx q[0];
rz(-2.7549681) q[0];
sx q[0];
rz(0.48166122) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6696981) q[2];
sx q[2];
rz(-1.8840232) q[2];
sx q[2];
rz(3.1196496) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1117275) q[1];
sx q[1];
rz(-1.7741412) q[1];
sx q[1];
rz(2.8257269) q[1];
rz(2.5727651) q[3];
sx q[3];
rz(-1.0320559) q[3];
sx q[3];
rz(-2.8728268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2729518) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(0.0090553332) q[2];
rz(-2.791259) q[3];
sx q[3];
rz(-0.30600268) q[3];
sx q[3];
rz(-2.8370324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7250799) q[0];
sx q[0];
rz(-2.0820936) q[0];
sx q[0];
rz(-2.1685725) q[0];
rz(1.579772) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(-0.80950338) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79671395) q[0];
sx q[0];
rz(-2.2911706) q[0];
sx q[0];
rz(-1.0864837) q[0];
rz(-pi) q[1];
rz(0.19301069) q[2];
sx q[2];
rz(-1.3456315) q[2];
sx q[2];
rz(-1.4373154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7743446) q[1];
sx q[1];
rz(-1.4092236) q[1];
sx q[1];
rz(2.0366686) q[1];
rz(-pi) q[2];
rz(0.46790917) q[3];
sx q[3];
rz(-2.10026) q[3];
sx q[3];
rz(0.61509815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96106625) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(2.6194438) q[2];
rz(0.3565878) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(-3.0780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4840354) q[0];
sx q[0];
rz(-2.7900896) q[0];
sx q[0];
rz(-1.8650613) q[0];
rz(-0.36318489) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(0.54375081) q[2];
sx q[2];
rz(-2.2866572) q[2];
sx q[2];
rz(0.86845749) q[2];
rz(1.2467066) q[3];
sx q[3];
rz(-2.1264599) q[3];
sx q[3];
rz(2.651941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
