OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.858415) q[0];
sx q[0];
rz(-1.5602966) q[0];
sx q[0];
rz(-2.4726782) q[0];
rz(-2.7954697) q[1];
sx q[1];
rz(-0.82155138) q[1];
sx q[1];
rz(0.93924826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4209983) q[0];
sx q[0];
rz(-1.2828553) q[0];
sx q[0];
rz(-1.1714515) q[0];
rz(3.0947565) q[2];
sx q[2];
rz(-1.9745262) q[2];
sx q[2];
rz(0.22778928) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8693162) q[1];
sx q[1];
rz(-0.47364911) q[1];
sx q[1];
rz(0.95817153) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40832728) q[3];
sx q[3];
rz(-0.75330594) q[3];
sx q[3];
rz(1.7024794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5758489) q[2];
sx q[2];
rz(-1.0172903) q[2];
sx q[2];
rz(-3.1080833) q[2];
rz(-1.2014028) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(1.0777773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031438436) q[0];
sx q[0];
rz(-0.25122508) q[0];
sx q[0];
rz(0.95397368) q[0];
rz(-1.4454449) q[1];
sx q[1];
rz(-2.0868389) q[1];
sx q[1];
rz(-1.9704069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69078835) q[0];
sx q[0];
rz(-2.224824) q[0];
sx q[0];
rz(0.79933856) q[0];
rz(-pi) q[1];
rz(-0.38925217) q[2];
sx q[2];
rz(-0.67485038) q[2];
sx q[2];
rz(-2.9794326) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.18942434) q[1];
sx q[1];
rz(-1.599424) q[1];
sx q[1];
rz(-0.22848033) q[1];
rz(-pi) q[2];
rz(-0.33447845) q[3];
sx q[3];
rz(-1.984716) q[3];
sx q[3];
rz(2.778307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90616068) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(1.1492427) q[2];
rz(3.1084642) q[3];
sx q[3];
rz(-1.5002316) q[3];
sx q[3];
rz(0.59605014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0843435) q[0];
sx q[0];
rz(-2.3488022) q[0];
sx q[0];
rz(1.0302011) q[0];
rz(-1.239981) q[1];
sx q[1];
rz(-1.7844618) q[1];
sx q[1];
rz(0.99303594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5035) q[0];
sx q[0];
rz(-0.4579535) q[0];
sx q[0];
rz(-3.0328337) q[0];
rz(-pi) q[1];
rz(-1.4517752) q[2];
sx q[2];
rz(-0.90672158) q[2];
sx q[2];
rz(-1.1609951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5532304) q[1];
sx q[1];
rz(-0.85236406) q[1];
sx q[1];
rz(-0.78939446) q[1];
rz(-pi) q[2];
rz(1.9957603) q[3];
sx q[3];
rz(-0.83163762) q[3];
sx q[3];
rz(1.1475565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9242636) q[2];
sx q[2];
rz(-1.1161048) q[2];
sx q[2];
rz(-0.35476157) q[2];
rz(2.1445856) q[3];
sx q[3];
rz(-1.6544147) q[3];
sx q[3];
rz(0.94064373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16911258) q[0];
sx q[0];
rz(-1.4262154) q[0];
sx q[0];
rz(1.1068363) q[0];
rz(-0.86889443) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(2.4443764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2223059) q[0];
sx q[0];
rz(-2.5450052) q[0];
sx q[0];
rz(-0.35937341) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72202335) q[2];
sx q[2];
rz(-1.6552087) q[2];
sx q[2];
rz(1.9868324) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47095151) q[1];
sx q[1];
rz(-1.1837675) q[1];
sx q[1];
rz(0.3989889) q[1];
x q[2];
rz(-1.5460062) q[3];
sx q[3];
rz(-1.5953334) q[3];
sx q[3];
rz(-0.33543521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2775468) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(2.8475658) q[2];
rz(2.0781519) q[3];
sx q[3];
rz(-1.4592417) q[3];
sx q[3];
rz(2.3991876) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71642891) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(2.678405) q[0];
rz(0.40395346) q[1];
sx q[1];
rz(-1.633176) q[1];
sx q[1];
rz(-2.892866) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2866658) q[0];
sx q[0];
rz(-1.5690985) q[0];
sx q[0];
rz(-1.3732852) q[0];
rz(-0.41422768) q[2];
sx q[2];
rz(-2.5606206) q[2];
sx q[2];
rz(1.5805949) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.675093) q[1];
sx q[1];
rz(-0.50216253) q[1];
sx q[1];
rz(2.2999022) q[1];
rz(-pi) q[2];
rz(-2.4234613) q[3];
sx q[3];
rz(-2.2635824) q[3];
sx q[3];
rz(0.75961514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8467466) q[2];
sx q[2];
rz(-1.3174026) q[2];
sx q[2];
rz(-1.3999636) q[2];
rz(1.8585662) q[3];
sx q[3];
rz(-1.3834407) q[3];
sx q[3];
rz(0.2259026) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1840709) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(0.34570178) q[0];
rz(-0.050994571) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(0.7720224) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22020082) q[0];
sx q[0];
rz(-1.3190375) q[0];
sx q[0];
rz(0.36248668) q[0];
x q[1];
rz(0.096944158) q[2];
sx q[2];
rz(-2.2205995) q[2];
sx q[2];
rz(-2.52663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4881106) q[1];
sx q[1];
rz(-1.0034434) q[1];
sx q[1];
rz(0.385872) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0431248) q[3];
sx q[3];
rz(-1.3195795) q[3];
sx q[3];
rz(0.71300292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(0.15743206) q[2];
rz(2.7031247) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(0.95611519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601783) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(2.881158) q[0];
rz(0.57834894) q[1];
sx q[1];
rz(-1.7702421) q[1];
sx q[1];
rz(-2.9885805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333014) q[0];
sx q[0];
rz(-0.75364764) q[0];
sx q[0];
rz(2.338999) q[0];
rz(-pi) q[1];
rz(3.0566759) q[2];
sx q[2];
rz(-1.474829) q[2];
sx q[2];
rz(-1.7964448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91100973) q[1];
sx q[1];
rz(-1.0167767) q[1];
sx q[1];
rz(-1.0826712) q[1];
rz(-pi) q[2];
rz(2.1649394) q[3];
sx q[3];
rz(-0.53721957) q[3];
sx q[3];
rz(-1.2402759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1795307) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(0.2641826) q[2];
rz(-1.3029441) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(1.1072268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1846979) q[0];
sx q[0];
rz(-2.1883924) q[0];
sx q[0];
rz(-1.4554998) q[0];
rz(-0.9616583) q[1];
sx q[1];
rz(-1.3788297) q[1];
sx q[1];
rz(2.2419194) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5035928) q[0];
sx q[0];
rz(-1.3595306) q[0];
sx q[0];
rz(2.6799623) q[0];
x q[1];
rz(-1.3191965) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(1.0902001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.051514678) q[1];
x q[2];
rz(3.1042388) q[3];
sx q[3];
rz(-1.787623) q[3];
sx q[3];
rz(1.2147533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4453033) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(-2.3882833) q[2];
rz(2.8478012) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(-1.1118579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488572) q[0];
sx q[0];
rz(-2.1519372) q[0];
sx q[0];
rz(-0.85187546) q[0];
rz(1.7333671) q[1];
sx q[1];
rz(-1.6586761) q[1];
sx q[1];
rz(-0.3826938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2792042) q[0];
sx q[0];
rz(-2.6529101) q[0];
sx q[0];
rz(0.53532289) q[0];
rz(-2.3716454) q[2];
sx q[2];
rz(-1.4048409) q[2];
sx q[2];
rz(-0.5391575) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5671317) q[1];
sx q[1];
rz(-1.2700915) q[1];
sx q[1];
rz(1.4682795) q[1];
x q[2];
rz(2.9858573) q[3];
sx q[3];
rz(-0.49113516) q[3];
sx q[3];
rz(-0.30065003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7644299) q[2];
sx q[2];
rz(-2.1775776) q[2];
sx q[2];
rz(-2.5386179) q[2];
rz(1.0682586) q[3];
sx q[3];
rz(-1.4385185) q[3];
sx q[3];
rz(-2.300613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339486) q[0];
sx q[0];
rz(-1.6772062) q[0];
sx q[0];
rz(-0.35368791) q[0];
rz(-2.0091281) q[1];
sx q[1];
rz(-2.1734889) q[1];
sx q[1];
rz(2.7521334) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4492466) q[0];
sx q[0];
rz(-0.55583891) q[0];
sx q[0];
rz(-2.6920094) q[0];
x q[1];
rz(1.816941) q[2];
sx q[2];
rz(-1.6707641) q[2];
sx q[2];
rz(1.0280746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3130541) q[1];
sx q[1];
rz(-1.2964998) q[1];
sx q[1];
rz(-2.6649339) q[1];
rz(-pi) q[2];
rz(-3.0515303) q[3];
sx q[3];
rz(-0.52999338) q[3];
sx q[3];
rz(1.3198205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0103717) q[2];
sx q[2];
rz(-1.6348569) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(1.5366588) q[3];
sx q[3];
rz(-2.6585572) q[3];
sx q[3];
rz(2.7893132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1491886) q[0];
sx q[0];
rz(-0.70269361) q[0];
sx q[0];
rz(0.50444034) q[0];
rz(-0.06012499) q[1];
sx q[1];
rz(-0.45777121) q[1];
sx q[1];
rz(-0.4578185) q[1];
rz(2.1679466) q[2];
sx q[2];
rz(-0.13046593) q[2];
sx q[2];
rz(2.2400357) q[2];
rz(-1.5271913) q[3];
sx q[3];
rz(-0.674896) q[3];
sx q[3];
rz(1.9622635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
