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
rz(0.67212927) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(-2.3840005) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286228) q[0];
sx q[0];
rz(-1.4951473) q[0];
sx q[0];
rz(1.4609758) q[0];
rz(-pi) q[1];
rz(1.2045317) q[2];
sx q[2];
rz(-1.0984813) q[2];
sx q[2];
rz(2.3031395) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1266844) q[1];
sx q[1];
rz(-1.9796625) q[1];
sx q[1];
rz(-1.3247299) q[1];
rz(-pi) q[2];
rz(0.67528649) q[3];
sx q[3];
rz(-2.3216341) q[3];
sx q[3];
rz(-0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0483094) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(1.7993571) q[2];
rz(-1.7736769) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(-1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9369478) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(2.192705) q[0];
rz(-0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(2.2116275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0596493) q[0];
sx q[0];
rz(-0.23632061) q[0];
sx q[0];
rz(-0.44775072) q[0];
rz(-1.5632767) q[2];
sx q[2];
rz(-2.0759656) q[2];
sx q[2];
rz(-1.9746321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8023493) q[1];
sx q[1];
rz(-1.7319458) q[1];
sx q[1];
rz(-1.6768964) q[1];
x q[2];
rz(-0.61701507) q[3];
sx q[3];
rz(-0.91147826) q[3];
sx q[3];
rz(-2.3544745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0019504) q[2];
sx q[2];
rz(-1.6763687) q[2];
sx q[2];
rz(-2.3993717) q[2];
rz(0.46418515) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6716229) q[0];
sx q[0];
rz(-2.9556584) q[0];
sx q[0];
rz(1.6999014) q[0];
rz(-0.16547671) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(0.86732078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.353426) q[0];
sx q[0];
rz(-1.5713619) q[0];
sx q[0];
rz(-3.1415523) q[0];
rz(-pi) q[1];
rz(2.8702552) q[2];
sx q[2];
rz(-2.0607161) q[2];
sx q[2];
rz(-0.64583954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8738385) q[1];
sx q[1];
rz(-1.8513362) q[1];
sx q[1];
rz(2.8024017) q[1];
x q[2];
rz(-1.0035361) q[3];
sx q[3];
rz(-2.2972882) q[3];
sx q[3];
rz(0.15641016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9637588) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(2.3826694) q[2];
rz(-0.80398503) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(-0.72833958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3114965) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-0.48402825) q[0];
rz(-1.6054035) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(-1.4253634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1041906) q[0];
sx q[0];
rz(-2.0811102) q[0];
sx q[0];
rz(1.9602106) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7984125) q[2];
sx q[2];
rz(-1.8755442) q[2];
sx q[2];
rz(-2.9272542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5023574) q[1];
sx q[1];
rz(-2.4662152) q[1];
sx q[1];
rz(1.2137854) q[1];
x q[2];
rz(-1.1953765) q[3];
sx q[3];
rz(-1.1470801) q[3];
sx q[3];
rz(0.28253192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.385685) q[2];
sx q[2];
rz(-1.0681095) q[2];
sx q[2];
rz(0.29328406) q[2];
rz(0.048132345) q[3];
sx q[3];
rz(-1.8354514) q[3];
sx q[3];
rz(-2.8369246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83679477) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(-1.1037214) q[0];
rz(0.91148218) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(-3.0326861) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1128588) q[0];
sx q[0];
rz(-2.2153531) q[0];
sx q[0];
rz(0.95405202) q[0];
rz(-pi) q[1];
rz(-2.8418903) q[2];
sx q[2];
rz(-2.5863918) q[2];
sx q[2];
rz(-0.19367684) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8629294) q[1];
sx q[1];
rz(-2.4554774) q[1];
sx q[1];
rz(1.3115694) q[1];
x q[2];
rz(2.5640574) q[3];
sx q[3];
rz(-1.3157237) q[3];
sx q[3];
rz(0.59781848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4917422) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(0.25807992) q[2];
rz(-0.38170013) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3968762) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(-1.9816403) q[0];
rz(0.57811919) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(-0.32346183) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37044034) q[0];
sx q[0];
rz(-1.9983294) q[0];
sx q[0];
rz(-1.0822791) q[0];
rz(-pi) q[1];
rz(-0.44536369) q[2];
sx q[2];
rz(-2.6008743) q[2];
sx q[2];
rz(-3.0034844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0929133) q[1];
sx q[1];
rz(-1.3447176) q[1];
sx q[1];
rz(-0.58633713) q[1];
rz(0.59185054) q[3];
sx q[3];
rz(-0.53880063) q[3];
sx q[3];
rz(2.0350698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20299992) q[2];
sx q[2];
rz(-2.3679569) q[2];
sx q[2];
rz(2.2933551) q[2];
rz(1.8741459) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271725) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(2.2681336) q[0];
rz(-1.9050441) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(0.083560856) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.386451) q[0];
sx q[0];
rz(-1.6928732) q[0];
sx q[0];
rz(0.32213078) q[0];
rz(1.1764063) q[2];
sx q[2];
rz(-2.5766234) q[2];
sx q[2];
rz(-0.51630563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.04490964) q[1];
sx q[1];
rz(-0.88777855) q[1];
sx q[1];
rz(0.51814305) q[1];
rz(1.2681581) q[3];
sx q[3];
rz(-1.916496) q[3];
sx q[3];
rz(1.3132335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7053335) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(1.7748888) q[2];
rz(2.8691835) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2347539) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(0.93836623) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(-0.82069194) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7308189) q[0];
sx q[0];
rz(-1.9389429) q[0];
sx q[0];
rz(0.14342043) q[0];
x q[1];
rz(0.16751473) q[2];
sx q[2];
rz(-1.2021087) q[2];
sx q[2];
rz(2.9978254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31466111) q[1];
sx q[1];
rz(-0.95035997) q[1];
sx q[1];
rz(-3.0112096) q[1];
rz(1.1701297) q[3];
sx q[3];
rz(-1.7224215) q[3];
sx q[3];
rz(-1.135863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9390823) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(-2.0738156) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(-0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15602569) q[0];
sx q[0];
rz(-1.1556867) q[0];
sx q[0];
rz(2.735403) q[0];
rz(1.73229) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(-1.0850151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.800348) q[0];
sx q[0];
rz(-2.4946176) q[0];
sx q[0];
rz(0.2917618) q[0];
rz(1.5753393) q[2];
sx q[2];
rz(-2.1811003) q[2];
sx q[2];
rz(-1.4686327) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1165365) q[1];
sx q[1];
rz(-2.3485314) q[1];
sx q[1];
rz(-1.1021986) q[1];
x q[2];
rz(-0.57906753) q[3];
sx q[3];
rz(-2.0430312) q[3];
sx q[3];
rz(0.079426182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5558418) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(0.93070585) q[2];
rz(-1.7454923) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(-0.11955587) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4569106) q[0];
sx q[0];
rz(-2.1645808) q[0];
sx q[0];
rz(1.8344301) q[0];
rz(0.3859418) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-2.1070259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32337727) q[0];
sx q[0];
rz(-1.744483) q[0];
sx q[0];
rz(-3.0370569) q[0];
rz(-2.7427086) q[2];
sx q[2];
rz(-2.6475057) q[2];
sx q[2];
rz(1.8979771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99920995) q[1];
sx q[1];
rz(-2.2015328) q[1];
sx q[1];
rz(-0.56908619) q[1];
rz(-pi) q[2];
rz(0.8108658) q[3];
sx q[3];
rz(-1.4812143) q[3];
sx q[3];
rz(0.60096622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4474386) q[2];
sx q[2];
rz(-2.4093781) q[2];
sx q[2];
rz(0.36699692) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(-1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.3601396) q[0];
sx q[0];
rz(-1.1893138) q[0];
sx q[0];
rz(-1.0839533) q[0];
rz(-3.0837334) q[1];
sx q[1];
rz(-1.2670988) q[1];
sx q[1];
rz(1.7115464) q[1];
rz(0.17292508) q[2];
sx q[2];
rz(-1.6676219) q[2];
sx q[2];
rz(0.22990083) q[2];
rz(1.9290942) q[3];
sx q[3];
rz(-0.90322106) q[3];
sx q[3];
rz(-1.7424104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
