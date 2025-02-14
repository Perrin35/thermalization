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
rz(1.4752969) q[0];
sx q[0];
rz(-1.2694321) q[0];
sx q[0];
rz(-2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(0.7575922) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.282895) q[0];
sx q[0];
rz(-3.0083249) q[0];
sx q[0];
rz(2.1758276) q[0];
x q[1];
rz(-1.937061) q[2];
sx q[2];
rz(-2.0431113) q[2];
sx q[2];
rz(-2.3031395) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54851531) q[1];
sx q[1];
rz(-0.47359772) q[1];
sx q[1];
rz(-2.6294336) q[1];
rz(-2.4663062) q[3];
sx q[3];
rz(-0.8199586) q[3];
sx q[3];
rz(0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0932833) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(1.3422356) q[2];
rz(-1.3679158) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(-1.5526519) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20464483) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(0.94888765) q[0];
rz(-0.10143796) q[1];
sx q[1];
rz(-2.0815492) q[1];
sx q[1];
rz(0.92996517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0895871) q[0];
sx q[0];
rz(-1.6723335) q[0];
sx q[0];
rz(0.2137645) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1279966) q[2];
sx q[2];
rz(-2.6363723) q[2];
sx q[2];
rz(-1.1824974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21446642) q[1];
sx q[1];
rz(-7*pi/15) q[1];
sx q[1];
rz(0.16204496) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5245776) q[3];
sx q[3];
rz(-0.91147826) q[3];
sx q[3];
rz(2.3544745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0019504) q[2];
sx q[2];
rz(-1.6763687) q[2];
sx q[2];
rz(0.742221) q[2];
rz(2.6774075) q[3];
sx q[3];
rz(-1.3451385) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6716229) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(1.6999014) q[0];
rz(-0.16547671) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(0.86732078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7169859) q[0];
sx q[0];
rz(-0.0005670003) q[0];
sx q[0];
rz(1.6419771) q[0];
rz(-pi) q[1];
rz(2.8702552) q[2];
sx q[2];
rz(-1.0808766) q[2];
sx q[2];
rz(0.64583954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9359303) q[1];
sx q[1];
rz(-1.2453658) q[1];
sx q[1];
rz(1.8673351) q[1];
rz(-pi) q[2];
rz(-0.8115143) q[3];
sx q[3];
rz(-1.9841188) q[3];
sx q[3];
rz(2.1275525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9637588) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(2.3826694) q[2];
rz(2.3376076) q[3];
sx q[3];
rz(-2.1920429) q[3];
sx q[3];
rz(-2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3114965) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(-2.6575644) q[0];
rz(-1.5361891) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(1.7162292) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.037402) q[0];
sx q[0];
rz(-1.0604825) q[0];
sx q[0];
rz(-1.181382) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3431802) q[2];
sx q[2];
rz(-1.8755442) q[2];
sx q[2];
rz(2.9272542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4933136) q[1];
sx q[1];
rz(-1.3505304) q[1];
sx q[1];
rz(0.92695285) q[1];
rz(-1.9462162) q[3];
sx q[3];
rz(-1.9945126) q[3];
sx q[3];
rz(-2.8590607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75590762) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(2.8483086) q[2];
rz(-3.0934603) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(-0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3047979) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(-1.1037214) q[0];
rz(2.2301105) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(0.10890659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2809362) q[0];
sx q[0];
rz(-2.0514279) q[0];
sx q[0];
rz(-2.3970766) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6066066) q[2];
sx q[2];
rz(-1.4145383) q[2];
sx q[2];
rz(-1.5076758) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8629294) q[1];
sx q[1];
rz(-2.4554774) q[1];
sx q[1];
rz(1.8300232) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6960228) q[3];
sx q[3];
rz(-2.516149) q[3];
sx q[3];
rz(-1.7991964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(-0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(-0.55735731) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3968762) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(1.9816403) q[0];
rz(0.57811919) q[1];
sx q[1];
rz(-1.6696397) q[1];
sx q[1];
rz(0.32346183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9834546) q[0];
sx q[0];
rz(-2.0120512) q[0];
sx q[0];
rz(2.6652314) q[0];
x q[1];
rz(-1.3176962) q[2];
sx q[2];
rz(-2.0539114) q[2];
sx q[2];
rz(-0.3699257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0929133) q[1];
sx q[1];
rz(-1.3447176) q[1];
sx q[1];
rz(2.5552555) q[1];
rz(-1.8927073) q[3];
sx q[3];
rz(-1.1309147) q[3];
sx q[3];
rz(1.3706576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9385927) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(0.8482376) q[2];
rz(1.8741459) q[3];
sx q[3];
rz(-1.4013314) q[3];
sx q[3];
rz(0.69376865) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271725) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(2.2681336) q[0];
rz(1.9050441) q[1];
sx q[1];
rz(-0.97573391) q[1];
sx q[1];
rz(0.083560856) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3665584) q[0];
sx q[0];
rz(-1.2511484) q[0];
sx q[0];
rz(1.699422) q[0];
rz(-pi) q[1];
rz(-1.0413076) q[2];
sx q[2];
rz(-1.363596) q[2];
sx q[2];
rz(0.71646128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2703184) q[1];
sx q[1];
rz(-1.9650998) q[1];
sx q[1];
rz(2.32347) q[1];
rz(1.8734345) q[3];
sx q[3];
rz(-1.916496) q[3];
sx q[3];
rz(1.8283591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7053335) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(-1.3667038) q[2];
rz(-2.8691835) q[3];
sx q[3];
rz(-1.0206157) q[3];
sx q[3];
rz(2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90683872) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(0.93836623) q[0];
rz(-0.80288184) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(0.82069194) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4107738) q[0];
sx q[0];
rz(-1.9389429) q[0];
sx q[0];
rz(2.9981722) q[0];
rz(-pi) q[1];
rz(2.9740779) q[2];
sx q[2];
rz(-1.2021087) q[2];
sx q[2];
rz(0.14376727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.332224) q[1];
sx q[1];
rz(-1.6767772) q[1];
sx q[1];
rz(-2.1952704) q[1];
rz(0.16444178) q[3];
sx q[3];
rz(-1.1749845) q[3];
sx q[3];
rz(2.7705517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20251033) q[2];
sx q[2];
rz(-0.15190092) q[2];
sx q[2];
rz(-1.067777) q[2];
rz(-2.0104525) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.985567) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(-0.40618968) q[0];
rz(-1.73229) q[1];
sx q[1];
rz(-2.1235762) q[1];
sx q[1];
rz(-1.0850151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4403518) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(1.7847654) q[0];
rz(-pi) q[1];
rz(2.5312838) q[2];
sx q[2];
rz(-1.5670735) q[2];
sx q[2];
rz(3.0368254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.93714) q[1];
sx q[1];
rz(-1.8984183) q[1];
sx q[1];
rz(0.83468584) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5625251) q[3];
sx q[3];
rz(-1.0985615) q[3];
sx q[3];
rz(-0.079426182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5857508) q[2];
sx q[2];
rz(-0.61183524) q[2];
sx q[2];
rz(2.2108868) q[2];
rz(-1.3961004) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(0.11955587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68468204) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(1.8344301) q[0];
rz(0.3859418) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-2.1070259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2725814) q[0];
sx q[0];
rz(-0.20244652) q[0];
sx q[0];
rz(-1.0342717) q[0];
x q[1];
rz(-0.46073353) q[2];
sx q[2];
rz(-1.7560395) q[2];
sx q[2];
rz(-0.68250193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9683796) q[1];
sx q[1];
rz(-0.82260859) q[1];
sx q[1];
rz(0.93507018) q[1];
rz(-pi) q[2];
rz(-1.7004556) q[3];
sx q[3];
rz(-0.76414062) q[3];
sx q[3];
rz(-1.0636927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6941541) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(-2.7745957) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-0.39803353) q[3];
sx q[3];
rz(-1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.3601396) q[0];
sx q[0];
rz(-1.1893138) q[0];
sx q[0];
rz(-1.0839533) q[0];
rz(3.0837334) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(2.9686676) q[2];
sx q[2];
rz(-1.4739707) q[2];
sx q[2];
rz(-2.9116918) q[2];
rz(-2.7230311) q[3];
sx q[3];
rz(-2.3971315) q[3];
sx q[3];
rz(-2.286398) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
