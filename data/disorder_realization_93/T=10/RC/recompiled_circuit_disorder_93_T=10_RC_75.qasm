OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38102725) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(0.3737803) q[0];
rz(-pi) q[1];
rz(1.107723) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(-2.3067834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1957789) q[1];
sx q[1];
rz(-1.073277) q[1];
sx q[1];
rz(1.0832018) q[1];
rz(-pi) q[2];
rz(1.9507017) q[3];
sx q[3];
rz(-1.1067179) q[3];
sx q[3];
rz(2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.1967999) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.686036) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3128132) q[0];
sx q[0];
rz(-1.4038741) q[0];
sx q[0];
rz(0.98915999) q[0];
rz(2.935264) q[2];
sx q[2];
rz(-2.7596139) q[2];
sx q[2];
rz(1.0155592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93745366) q[1];
sx q[1];
rz(-1.0480282) q[1];
sx q[1];
rz(-0.015603113) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0807642) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.3519752) q[2];
rz(2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(-0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.9690537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5867509) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(2.5462333) q[0];
rz(-2.3861888) q[2];
sx q[2];
rz(-1.3284151) q[2];
sx q[2];
rz(2.541045) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2492003) q[1];
sx q[1];
rz(-1.8141659) q[1];
sx q[1];
rz(-1.0358441) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3200687) q[3];
sx q[3];
rz(-1.513895) q[3];
sx q[3];
rz(-2.1123321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-0.52350837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2098171) q[0];
sx q[0];
rz(-2.1486001) q[0];
sx q[0];
rz(1.125976) q[0];
rz(-pi) q[1];
rz(0.24387118) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(1.3432168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88672968) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(2.8121594) q[1];
x q[2];
rz(-0.13682271) q[3];
sx q[3];
rz(-2.1577583) q[3];
sx q[3];
rz(2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.4915285) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(0.99194828) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46734738) q[0];
sx q[0];
rz(-1.6773946) q[0];
sx q[0];
rz(2.9664413) q[0];
rz(-pi) q[1];
rz(-1.2679342) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(2.0451343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93047749) q[1];
sx q[1];
rz(-2.8301297) q[1];
sx q[1];
rz(-3.009797) q[1];
rz(-1.7080073) q[3];
sx q[3];
rz(-0.76223323) q[3];
sx q[3];
rz(-1.2587794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0118959) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(-2.8707855) q[2];
rz(-2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(1.7927992) q[0];
rz(0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7871008) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(1.8255193) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8259949) q[2];
sx q[2];
rz(-1.4824502) q[2];
sx q[2];
rz(1.7905854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52275601) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(1.418581) q[1];
rz(-pi) q[2];
rz(-2.2819727) q[3];
sx q[3];
rz(-0.73577995) q[3];
sx q[3];
rz(1.7914634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069698378) q[0];
sx q[0];
rz(-2.0928203) q[0];
sx q[0];
rz(-2.4126023) q[0];
rz(0.49634883) q[2];
sx q[2];
rz(-0.90663547) q[2];
sx q[2];
rz(1.6830483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2591178) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(0.74101733) q[1];
rz(-pi) q[2];
rz(1.20007) q[3];
sx q[3];
rz(-2.126174) q[3];
sx q[3];
rz(0.45645082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-2.3186671) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(-1.3051422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179498) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(-2.6108517) q[0];
rz(-0.70456409) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(1.2003843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0481938) q[1];
sx q[1];
rz(-1.1953925) q[1];
sx q[1];
rz(2.4673389) q[1];
x q[2];
rz(1.9100902) q[3];
sx q[3];
rz(-1.5929856) q[3];
sx q[3];
rz(2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1398853) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2625933) q[0];
sx q[0];
rz(-1.9290553) q[0];
sx q[0];
rz(0.4155638) q[0];
rz(-3.1019194) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(-0.85658011) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.326509) q[1];
sx q[1];
rz(-1.0352967) q[1];
sx q[1];
rz(-2.9498847) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95548198) q[3];
sx q[3];
rz(-1.2947047) q[3];
sx q[3];
rz(-2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(1.3396324) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.8113177) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9194473) q[0];
sx q[0];
rz(-1.3627909) q[0];
sx q[0];
rz(-0.39214765) q[0];
rz(2.8154545) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(-1.067576) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2965282) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(-1.8370085) q[1];
rz(-pi) q[2];
rz(0.43379421) q[3];
sx q[3];
rz(-1.1489831) q[3];
sx q[3];
rz(-2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(2.4516104) q[2];
sx q[2];
rz(-2.1524515) q[2];
sx q[2];
rz(-0.099302789) q[2];
rz(0.91924304) q[3];
sx q[3];
rz(-1.6332492) q[3];
sx q[3];
rz(1.9132683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
