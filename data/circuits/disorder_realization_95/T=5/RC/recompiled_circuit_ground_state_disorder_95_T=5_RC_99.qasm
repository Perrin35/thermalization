OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4841782) q[0];
sx q[0];
rz(-1.2793469) q[0];
sx q[0];
rz(-2.4980463) q[0];
rz(-2.9930288) q[1];
sx q[1];
rz(-0.64639503) q[1];
sx q[1];
rz(0.57599154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2781485) q[0];
sx q[0];
rz(-1.8810026) q[0];
sx q[0];
rz(2.1351867) q[0];
rz(-1.2744489) q[2];
sx q[2];
rz(-1.4055894) q[2];
sx q[2];
rz(-0.295754) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6580087) q[1];
sx q[1];
rz(-1.600126) q[1];
sx q[1];
rz(1.4714249) q[1];
x q[2];
rz(-2.9408119) q[3];
sx q[3];
rz(-0.83856486) q[3];
sx q[3];
rz(1.689336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5129471) q[2];
sx q[2];
rz(-0.60428667) q[2];
sx q[2];
rz(0.71913546) q[2];
rz(-2.5288845) q[3];
sx q[3];
rz(-0.49379525) q[3];
sx q[3];
rz(1.4601716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4145819) q[0];
sx q[0];
rz(-2.5717323) q[0];
sx q[0];
rz(-1.3563096) q[0];
rz(2.5510229) q[1];
sx q[1];
rz(-1.5868203) q[1];
sx q[1];
rz(2.2812567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3911678) q[0];
sx q[0];
rz(-0.92022395) q[0];
sx q[0];
rz(-2.9123996) q[0];
rz(0.27766547) q[2];
sx q[2];
rz(-2.5685453) q[2];
sx q[2];
rz(-2.3615357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7041143) q[1];
sx q[1];
rz(-1.8859004) q[1];
sx q[1];
rz(1.9060288) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11252304) q[3];
sx q[3];
rz(-0.72970245) q[3];
sx q[3];
rz(-2.2000748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4071953) q[2];
sx q[2];
rz(-2.037214) q[2];
sx q[2];
rz(0.26367903) q[2];
rz(0.10144357) q[3];
sx q[3];
rz(-2.0439309) q[3];
sx q[3];
rz(-1.6775848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2270671) q[0];
sx q[0];
rz(-2.3093746) q[0];
sx q[0];
rz(-1.9050003) q[0];
rz(-2.6446758) q[1];
sx q[1];
rz(-2.2221815) q[1];
sx q[1];
rz(-2.9317754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38302177) q[0];
sx q[0];
rz(-1.5514217) q[0];
sx q[0];
rz(-0.13067496) q[0];
rz(-pi) q[1];
rz(1.2926854) q[2];
sx q[2];
rz(-0.46574083) q[2];
sx q[2];
rz(-1.9139707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4821149) q[1];
sx q[1];
rz(-2.6365197) q[1];
sx q[1];
rz(1.6467846) q[1];
rz(-pi) q[2];
rz(0.83720343) q[3];
sx q[3];
rz(-1.687369) q[3];
sx q[3];
rz(2.0904515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1257552) q[2];
sx q[2];
rz(-2.0070984) q[2];
sx q[2];
rz(0.26975676) q[2];
rz(2.6848327) q[3];
sx q[3];
rz(-2.7220107) q[3];
sx q[3];
rz(2.8381798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.1064827) q[0];
sx q[0];
rz(-0.39316097) q[0];
sx q[0];
rz(3.0495354) q[0];
rz(-0.58097845) q[1];
sx q[1];
rz(-0.61994225) q[1];
sx q[1];
rz(2.7041025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1745146) q[0];
sx q[0];
rz(-2.1319492) q[0];
sx q[0];
rz(0.18995096) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0863414) q[2];
sx q[2];
rz(-2.8979282) q[2];
sx q[2];
rz(-0.59020868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60595861) q[1];
sx q[1];
rz(-1.8211963) q[1];
sx q[1];
rz(-1.1634924) q[1];
rz(-pi) q[2];
rz(-0.32467036) q[3];
sx q[3];
rz(-0.89327795) q[3];
sx q[3];
rz(-0.087740104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1979394) q[2];
sx q[2];
rz(-1.7655756) q[2];
sx q[2];
rz(2.9735238) q[2];
rz(2.4456612) q[3];
sx q[3];
rz(-1.1974502) q[3];
sx q[3];
rz(1.7456938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.593852) q[0];
sx q[0];
rz(-0.20741367) q[0];
sx q[0];
rz(2.4399624) q[0];
rz(0.54948366) q[1];
sx q[1];
rz(-1.0797078) q[1];
sx q[1];
rz(0.93542498) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0129378) q[0];
sx q[0];
rz(-1.8194981) q[0];
sx q[0];
rz(-1.1292581) q[0];
rz(-pi) q[1];
rz(-0.40570183) q[2];
sx q[2];
rz(-3.0076111) q[2];
sx q[2];
rz(-1.2359499) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7364663) q[1];
sx q[1];
rz(-1.163981) q[1];
sx q[1];
rz(-2.1975425) q[1];
x q[2];
rz(1.5059371) q[3];
sx q[3];
rz(-1.4993414) q[3];
sx q[3];
rz(2.2301729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51297411) q[2];
sx q[2];
rz(-1.3209891) q[2];
sx q[2];
rz(-0.67250195) q[2];
rz(0.77338141) q[3];
sx q[3];
rz(-2.0956109) q[3];
sx q[3];
rz(-2.8092303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8124939) q[0];
sx q[0];
rz(-2.2342873) q[0];
sx q[0];
rz(-1.6727653) q[0];
rz(-2.2724197) q[1];
sx q[1];
rz(-0.69907993) q[1];
sx q[1];
rz(-0.3259784) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7237968) q[0];
sx q[0];
rz(-1.7053331) q[0];
sx q[0];
rz(2.8409182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5766854) q[2];
sx q[2];
rz(-1.1400901) q[2];
sx q[2];
rz(-0.6503833) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80854844) q[1];
sx q[1];
rz(-0.99892975) q[1];
sx q[1];
rz(2.1480335) q[1];
x q[2];
rz(0.066727344) q[3];
sx q[3];
rz(-1.9983833) q[3];
sx q[3];
rz(-0.59837435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0526696) q[2];
sx q[2];
rz(-1.576705) q[2];
sx q[2];
rz(0.79724533) q[2];
rz(-3.0439324) q[3];
sx q[3];
rz(-0.74334136) q[3];
sx q[3];
rz(-1.2231479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.866211) q[0];
sx q[0];
rz(-0.99358639) q[0];
sx q[0];
rz(2.1589101) q[0];
rz(1.6536225) q[1];
sx q[1];
rz(-0.5653615) q[1];
sx q[1];
rz(0.32378325) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8651915) q[0];
sx q[0];
rz(-0.29769024) q[0];
sx q[0];
rz(0.28930347) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9672438) q[2];
sx q[2];
rz(-0.65047036) q[2];
sx q[2];
rz(0.54325587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.946879) q[1];
sx q[1];
rz(-2.4733798) q[1];
sx q[1];
rz(-0.18385033) q[1];
rz(-pi) q[2];
rz(0.57347943) q[3];
sx q[3];
rz(-2.5635168) q[3];
sx q[3];
rz(-0.40595192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4299778) q[2];
sx q[2];
rz(-2.6275676) q[2];
sx q[2];
rz(0.035710486) q[2];
rz(2.4014373) q[3];
sx q[3];
rz(-1.6868846) q[3];
sx q[3];
rz(-1.1470225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78298727) q[0];
sx q[0];
rz(-2.7645223) q[0];
sx q[0];
rz(2.9407035) q[0];
rz(-2.5783077) q[1];
sx q[1];
rz(-1.7726243) q[1];
sx q[1];
rz(2.7422781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92420805) q[0];
sx q[0];
rz(-1.7151807) q[0];
sx q[0];
rz(-1.420279) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6496193) q[2];
sx q[2];
rz(-0.48187253) q[2];
sx q[2];
rz(-0.80480591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3914403) q[1];
sx q[1];
rz(-1.925577) q[1];
sx q[1];
rz(-1.9998235) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7186728) q[3];
sx q[3];
rz(-0.81133553) q[3];
sx q[3];
rz(-0.75046152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1839972) q[2];
sx q[2];
rz(-2.5338569) q[2];
sx q[2];
rz(3.0228534) q[2];
rz(-0.27058288) q[3];
sx q[3];
rz(-2.1515473) q[3];
sx q[3];
rz(-2.9873007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0305233) q[0];
sx q[0];
rz(-1.1708165) q[0];
sx q[0];
rz(0.4554553) q[0];
rz(-0.2332553) q[1];
sx q[1];
rz(-2.3979954) q[1];
sx q[1];
rz(-3.1366248) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7457897) q[0];
sx q[0];
rz(-3.0533049) q[0];
sx q[0];
rz(-0.39678617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.754141) q[2];
sx q[2];
rz(-1.8672647) q[2];
sx q[2];
rz(-2.7067892) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50434166) q[1];
sx q[1];
rz(-2.0310861) q[1];
sx q[1];
rz(1.00072) q[1];
rz(-pi) q[2];
rz(-0.016214215) q[3];
sx q[3];
rz(-1.1450279) q[3];
sx q[3];
rz(-2.2174284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7838955) q[2];
sx q[2];
rz(-1.8507439) q[2];
sx q[2];
rz(2.2515187) q[2];
rz(-2.2642073) q[3];
sx q[3];
rz(-0.93534094) q[3];
sx q[3];
rz(-0.64524159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18168618) q[0];
sx q[0];
rz(-3.063995) q[0];
sx q[0];
rz(1.0542057) q[0];
rz(-2.8554845) q[1];
sx q[1];
rz(-2.095463) q[1];
sx q[1];
rz(-1.8792763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12114502) q[0];
sx q[0];
rz(-0.95709267) q[0];
sx q[0];
rz(-1.2200939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91210552) q[2];
sx q[2];
rz(-0.3276796) q[2];
sx q[2];
rz(-0.35057783) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8827829) q[1];
sx q[1];
rz(-1.38816) q[1];
sx q[1];
rz(-0.44328494) q[1];
rz(-pi) q[2];
rz(-1.0159053) q[3];
sx q[3];
rz(-1.4994748) q[3];
sx q[3];
rz(1.2268905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.012764843) q[2];
sx q[2];
rz(-2.3598857) q[2];
sx q[2];
rz(0.26415602) q[2];
rz(2.9014897) q[3];
sx q[3];
rz(-2.0930347) q[3];
sx q[3];
rz(-0.31606328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8964597) q[0];
sx q[0];
rz(-2.3871853) q[0];
sx q[0];
rz(2.2375794) q[0];
rz(-0.28884197) q[1];
sx q[1];
rz(-1.387351) q[1];
sx q[1];
rz(3.0659061) q[1];
rz(-2.4905229) q[2];
sx q[2];
rz(-2.747703) q[2];
sx q[2];
rz(2.2892164) q[2];
rz(-3.0700923) q[3];
sx q[3];
rz(-1.6115474) q[3];
sx q[3];
rz(-1.1781884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
