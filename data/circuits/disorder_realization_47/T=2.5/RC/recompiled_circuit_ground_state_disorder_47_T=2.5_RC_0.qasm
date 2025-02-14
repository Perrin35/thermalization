OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.4247704) q[0];
sx q[0];
rz(4.7228887) q[0];
sx q[0];
rz(6.9520998) q[0];
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72059435) q[0];
sx q[0];
rz(-1.8587374) q[0];
sx q[0];
rz(-1.1714515) q[0];
rz(-pi) q[1];
rz(1.9749227) q[2];
sx q[2];
rz(-1.6138645) q[2];
sx q[2];
rz(-1.3614181) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8693162) q[1];
sx q[1];
rz(-2.6679435) q[1];
sx q[1];
rz(-2.1834211) q[1];
rz(-1.2143308) q[3];
sx q[3];
rz(-0.89205304) q[3];
sx q[3];
rz(-1.1671305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5758489) q[2];
sx q[2];
rz(-2.1243024) q[2];
sx q[2];
rz(3.1080833) q[2];
rz(-1.9401898) q[3];
sx q[3];
rz(-1.7824495) q[3];
sx q[3];
rz(1.0777773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031438436) q[0];
sx q[0];
rz(-2.8903676) q[0];
sx q[0];
rz(0.95397368) q[0];
rz(1.4454449) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(1.1711858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508043) q[0];
sx q[0];
rz(-2.224824) q[0];
sx q[0];
rz(0.79933856) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38925217) q[2];
sx q[2];
rz(-2.4667423) q[2];
sx q[2];
rz(2.9794326) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7668768) q[1];
sx q[1];
rz(-1.3424113) q[1];
sx q[1];
rz(1.5414052) q[1];
rz(2.0061139) q[3];
sx q[3];
rz(-1.8760699) q[3];
sx q[3];
rz(1.0686309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.235432) q[2];
sx q[2];
rz(-1.7110598) q[2];
sx q[2];
rz(-1.1492427) q[2];
rz(3.1084642) q[3];
sx q[3];
rz(-1.5002316) q[3];
sx q[3];
rz(-2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0843435) q[0];
sx q[0];
rz(-0.79279041) q[0];
sx q[0];
rz(-1.0302011) q[0];
rz(-1.9016117) q[1];
sx q[1];
rz(-1.7844618) q[1];
sx q[1];
rz(-0.99303594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6246374) q[0];
sx q[0];
rz(-2.0258396) q[0];
sx q[0];
rz(-1.5173453) q[0];
x q[1];
rz(1.6898174) q[2];
sx q[2];
rz(-2.2348711) q[2];
sx q[2];
rz(1.1609951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5883623) q[1];
sx q[1];
rz(-0.85236406) q[1];
sx q[1];
rz(-2.3521982) q[1];
x q[2];
rz(0.42476023) q[3];
sx q[3];
rz(-2.3094607) q[3];
sx q[3];
rz(0.55603851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9242636) q[2];
sx q[2];
rz(-2.0254878) q[2];
sx q[2];
rz(0.35476157) q[2];
rz(2.1445856) q[3];
sx q[3];
rz(-1.487178) q[3];
sx q[3];
rz(2.2009489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16911258) q[0];
sx q[0];
rz(-1.7153772) q[0];
sx q[0];
rz(1.1068363) q[0];
rz(0.86889443) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(0.69721627) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95283629) q[0];
sx q[0];
rz(-1.3719014) q[0];
sx q[0];
rz(2.5752978) q[0];
rz(-pi) q[1];
rz(3.0142586) q[2];
sx q[2];
rz(-0.72605726) q[2];
sx q[2];
rz(-0.51148326) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6706411) q[1];
sx q[1];
rz(-1.9578251) q[1];
sx q[1];
rz(-2.7426038) q[1];
rz(-pi) q[2];
rz(3.117048) q[3];
sx q[3];
rz(-1.5460137) q[3];
sx q[3];
rz(-1.2359695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2775468) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-2.8475658) q[2];
rz(-2.0781519) q[3];
sx q[3];
rz(-1.4592417) q[3];
sx q[3];
rz(0.74240509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.71642891) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(-0.46318769) q[0];
rz(-0.40395346) q[1];
sx q[1];
rz(-1.633176) q[1];
sx q[1];
rz(-0.24872669) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8578019) q[0];
sx q[0];
rz(-1.3732855) q[0];
sx q[0];
rz(0.0017314712) q[0];
rz(-pi) q[1];
rz(2.727365) q[2];
sx q[2];
rz(-0.58097208) q[2];
sx q[2];
rz(-1.5805949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.675093) q[1];
sx q[1];
rz(-0.50216253) q[1];
sx q[1];
rz(-0.84169047) q[1];
rz(-pi) q[2];
rz(-2.4234613) q[3];
sx q[3];
rz(-2.2635824) q[3];
sx q[3];
rz(-2.3819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8467466) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(1.741629) q[2];
rz(1.8585662) q[3];
sx q[3];
rz(-1.7581519) q[3];
sx q[3];
rz(2.9156901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(0.95752174) q[0];
sx q[0];
rz(-2.3029843) q[0];
sx q[0];
rz(-0.34570178) q[0];
rz(0.050994571) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(2.3695703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9213918) q[0];
sx q[0];
rz(-1.3190375) q[0];
sx q[0];
rz(0.36248668) q[0];
x q[1];
rz(-2.2228681) q[2];
sx q[2];
rz(-1.6479392) q[2];
sx q[2];
rz(-0.89706286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.009338) q[1];
sx q[1];
rz(-1.2478095) q[1];
sx q[1];
rz(2.1732974) q[1];
x q[2];
rz(-1.0994301) q[3];
sx q[3];
rz(-2.5623218) q[3];
sx q[3];
rz(-1.260965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8649851) q[2];
sx q[2];
rz(-1.7165246) q[2];
sx q[2];
rz(2.9841606) q[2];
rz(2.7031247) q[3];
sx q[3];
rz(-0.74123588) q[3];
sx q[3];
rz(-0.95611519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081414374) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(9/(11*pi)) q[0];
rz(-0.57834894) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(0.15301212) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48400797) q[0];
sx q[0];
rz(-2.0853242) q[0];
sx q[0];
rz(-0.57782509) q[0];
rz(-0.084916755) q[2];
sx q[2];
rz(-1.6667637) q[2];
sx q[2];
rz(-1.3451479) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0205903) q[1];
sx q[1];
rz(-0.7210702) q[1];
sx q[1];
rz(-2.4929558) q[1];
rz(-2.1649394) q[3];
sx q[3];
rz(-0.53721957) q[3];
sx q[3];
rz(1.2402759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.962062) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(0.2641826) q[2];
rz(-1.8386486) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(-1.1072268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95689479) q[0];
sx q[0];
rz(-2.1883924) q[0];
sx q[0];
rz(1.6860929) q[0];
rz(0.9616583) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(2.2419194) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6379999) q[0];
sx q[0];
rz(-1.3595306) q[0];
sx q[0];
rz(2.6799623) q[0];
x q[1];
rz(2.9177306) q[2];
sx q[2];
rz(-0.85843411) q[2];
sx q[2];
rz(1.7162042) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2921999) q[1];
sx q[1];
rz(-1.527904) q[1];
sx q[1];
rz(2.1581743) q[1];
rz(1.353823) q[3];
sx q[3];
rz(-1.6072751) q[3];
sx q[3];
rz(0.36408261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6962894) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(2.3882833) q[2];
rz(-2.8478012) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(1.1118579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.4488572) q[0];
sx q[0];
rz(-2.1519372) q[0];
sx q[0];
rz(2.2897172) q[0];
rz(1.7333671) q[1];
sx q[1];
rz(-1.4829166) q[1];
sx q[1];
rz(0.3826938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4538759) q[0];
sx q[0];
rz(-1.9864489) q[0];
sx q[0];
rz(1.835653) q[0];
rz(-1.7999951) q[2];
sx q[2];
rz(-0.8140854) q[2];
sx q[2];
rz(-2.2687721) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.908757) q[1];
sx q[1];
rz(-0.31719734) q[1];
sx q[1];
rz(-0.31875191) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6535651) q[3];
sx q[3];
rz(-1.0861229) q[3];
sx q[3];
rz(3.017149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3771628) q[2];
sx q[2];
rz(-0.96401507) q[2];
sx q[2];
rz(2.5386179) q[2];
rz(1.0682586) q[3];
sx q[3];
rz(-1.4385185) q[3];
sx q[3];
rz(0.8409797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0076440796) q[0];
sx q[0];
rz(-1.6772062) q[0];
sx q[0];
rz(-0.35368791) q[0];
rz(-2.0091281) q[1];
sx q[1];
rz(-0.9681038) q[1];
sx q[1];
rz(-2.7521334) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2089473) q[0];
sx q[0];
rz(-2.0660127) q[0];
sx q[0];
rz(-1.8344648) q[0];
rz(-0.10305291) q[2];
sx q[2];
rz(-1.3259058) q[2];
sx q[2];
rz(-2.6239397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82853854) q[1];
sx q[1];
rz(-1.2964998) q[1];
sx q[1];
rz(2.6649339) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6234446) q[3];
sx q[3];
rz(-2.0984167) q[3];
sx q[3];
rz(-1.9260581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13122095) q[2];
sx q[2];
rz(-1.5067357) q[2];
sx q[2];
rz(1.1466675) q[2];
rz(-1.6049339) q[3];
sx q[3];
rz(-2.6585572) q[3];
sx q[3];
rz(2.7893132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.1491886) q[0];
sx q[0];
rz(-0.70269361) q[0];
sx q[0];
rz(0.50444034) q[0];
rz(0.06012499) q[1];
sx q[1];
rz(-2.6838214) q[1];
sx q[1];
rz(2.6837742) q[1];
rz(0.9736461) q[2];
sx q[2];
rz(-3.0111267) q[2];
sx q[2];
rz(-0.901557) q[2];
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
