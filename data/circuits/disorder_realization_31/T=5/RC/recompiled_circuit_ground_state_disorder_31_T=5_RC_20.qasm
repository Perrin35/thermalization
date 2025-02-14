OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6346729) q[0];
sx q[0];
rz(-1.7393751) q[0];
sx q[0];
rz(2.3911067) q[0];
rz(-3.0783202) q[1];
sx q[1];
rz(-0.81875357) q[1];
sx q[1];
rz(-1.0190581) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9813596) q[0];
sx q[0];
rz(-2.5447114) q[0];
sx q[0];
rz(0.33905466) q[0];
x q[1];
rz(0.86111492) q[2];
sx q[2];
rz(-1.7758596) q[2];
sx q[2];
rz(0.67180639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77324919) q[1];
sx q[1];
rz(-1.7071525) q[1];
sx q[1];
rz(-2.5565867) q[1];
rz(-0.54368162) q[3];
sx q[3];
rz(-1.8988109) q[3];
sx q[3];
rz(-1.3222493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84746209) q[2];
sx q[2];
rz(-2.1475466) q[2];
sx q[2];
rz(-1.6729986) q[2];
rz(-0.20644203) q[3];
sx q[3];
rz(-1.8136181) q[3];
sx q[3];
rz(-2.4895721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1714529) q[0];
sx q[0];
rz(-1.411922) q[0];
sx q[0];
rz(-0.10301244) q[0];
rz(-0.39762321) q[1];
sx q[1];
rz(-0.42639521) q[1];
sx q[1];
rz(0.85223371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6838339) q[0];
sx q[0];
rz(-1.6310098) q[0];
sx q[0];
rz(1.0369861) q[0];
rz(-1.909341) q[2];
sx q[2];
rz(-1.8391092) q[2];
sx q[2];
rz(1.878083) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.045984775) q[1];
sx q[1];
rz(-1.5846445) q[1];
sx q[1];
rz(1.0471837) q[1];
rz(-1.6144595) q[3];
sx q[3];
rz(-1.3779197) q[3];
sx q[3];
rz(-1.5571644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6323866) q[2];
sx q[2];
rz(-2.1337324) q[2];
sx q[2];
rz(-0.81727916) q[2];
rz(-0.20538524) q[3];
sx q[3];
rz(-0.21586625) q[3];
sx q[3];
rz(0.5425905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83194724) q[0];
sx q[0];
rz(-2.0537856) q[0];
sx q[0];
rz(-0.52016869) q[0];
rz(-2.4179516) q[1];
sx q[1];
rz(-1.8579204) q[1];
sx q[1];
rz(2.9626194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163612) q[0];
sx q[0];
rz(-2.0746914) q[0];
sx q[0];
rz(-1.0625307) q[0];
rz(-pi) q[1];
rz(2.6232324) q[2];
sx q[2];
rz(-2.7438965) q[2];
sx q[2];
rz(-2.8373371) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6295669) q[1];
sx q[1];
rz(-0.24228748) q[1];
sx q[1];
rz(-2.971388) q[1];
x q[2];
rz(-2.4839999) q[3];
sx q[3];
rz(-1.7879855) q[3];
sx q[3];
rz(1.515102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2554864) q[2];
sx q[2];
rz(-2.8612558) q[2];
sx q[2];
rz(-2.172016) q[2];
rz(-0.3420091) q[3];
sx q[3];
rz(-1.0229144) q[3];
sx q[3];
rz(1.3268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.2372357) q[0];
sx q[0];
rz(-2.3396753) q[0];
sx q[0];
rz(-3.0661769) q[0];
rz(-0.30741179) q[1];
sx q[1];
rz(-1.1630029) q[1];
sx q[1];
rz(-2.7154162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722011) q[0];
sx q[0];
rz(-1.0474023) q[0];
sx q[0];
rz(2.4358569) q[0];
rz(-pi) q[1];
rz(1.8521339) q[2];
sx q[2];
rz(-2.0232213) q[2];
sx q[2];
rz(1.4624925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7245009) q[1];
sx q[1];
rz(-2.1414653) q[1];
sx q[1];
rz(2.9076734) q[1];
rz(3.0846023) q[3];
sx q[3];
rz(-2.7354623) q[3];
sx q[3];
rz(-0.95088357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9093466) q[2];
sx q[2];
rz(-1.1541977) q[2];
sx q[2];
rz(0.92440355) q[2];
rz(-1.1030446) q[3];
sx q[3];
rz(-2.0274935) q[3];
sx q[3];
rz(-1.2801142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1396609) q[0];
sx q[0];
rz(-2.6502471) q[0];
sx q[0];
rz(-2.370148) q[0];
rz(-1.3123243) q[1];
sx q[1];
rz(-1.7701365) q[1];
sx q[1];
rz(1.2244474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3712922) q[0];
sx q[0];
rz(-1.5988358) q[0];
sx q[0];
rz(-1.7078778) q[0];
rz(-pi) q[1];
rz(-0.0026826492) q[2];
sx q[2];
rz(-0.61842881) q[2];
sx q[2];
rz(-1.2346877) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2001554) q[1];
sx q[1];
rz(-1.9552045) q[1];
sx q[1];
rz(3.0430493) q[1];
rz(-pi) q[2];
rz(2.6825887) q[3];
sx q[3];
rz(-2.3025666) q[3];
sx q[3];
rz(2.9162198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9839342) q[2];
sx q[2];
rz(-2.5658786) q[2];
sx q[2];
rz(-2.0470587) q[2];
rz(-3.0648699) q[3];
sx q[3];
rz(-0.50944296) q[3];
sx q[3];
rz(-2.4491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.43129608) q[0];
sx q[0];
rz(-2.4831979) q[0];
sx q[0];
rz(-0.21507138) q[0];
rz(2.1521425) q[1];
sx q[1];
rz(-0.96885252) q[1];
sx q[1];
rz(1.4873803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0519245) q[0];
sx q[0];
rz(-1.4082272) q[0];
sx q[0];
rz(-2.6211365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2770259) q[2];
sx q[2];
rz(-1.1907628) q[2];
sx q[2];
rz(-0.88949163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1035005) q[1];
sx q[1];
rz(-1.2088782) q[1];
sx q[1];
rz(-1.1236217) q[1];
rz(-pi) q[2];
rz(-0.59923895) q[3];
sx q[3];
rz(-0.75746292) q[3];
sx q[3];
rz(0.77764192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.4789077) q[2];
sx q[2];
rz(-1.7272793) q[2];
sx q[2];
rz(-1.5642081) q[2];
rz(-0.94791895) q[3];
sx q[3];
rz(-1.0736059) q[3];
sx q[3];
rz(-1.7387559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74744144) q[0];
sx q[0];
rz(-1.6524557) q[0];
sx q[0];
rz(2.5519651) q[0];
rz(1.6472752) q[1];
sx q[1];
rz(-1.7604156) q[1];
sx q[1];
rz(3.0628915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704961) q[0];
sx q[0];
rz(-1.9645623) q[0];
sx q[0];
rz(-1.257819) q[0];
x q[1];
rz(1.3316657) q[2];
sx q[2];
rz(-2.0239365) q[2];
sx q[2];
rz(-1.5117642) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5906395) q[1];
sx q[1];
rz(-1.7788243) q[1];
sx q[1];
rz(2.8880638) q[1];
x q[2];
rz(0.77459333) q[3];
sx q[3];
rz(-2.4586185) q[3];
sx q[3];
rz(-0.81688629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9239203) q[2];
sx q[2];
rz(-0.87409002) q[2];
sx q[2];
rz(0.11202845) q[2];
rz(0.40515408) q[3];
sx q[3];
rz(-1.3978037) q[3];
sx q[3];
rz(2.9356094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21796313) q[0];
sx q[0];
rz(-2.1681652) q[0];
sx q[0];
rz(1.4189036) q[0];
rz(-1.0796684) q[1];
sx q[1];
rz(-0.38366145) q[1];
sx q[1];
rz(1.4071646) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0100219) q[0];
sx q[0];
rz(-2.6162806) q[0];
sx q[0];
rz(-2.4729687) q[0];
rz(-pi) q[1];
rz(-0.17282829) q[2];
sx q[2];
rz(-2.7044074) q[2];
sx q[2];
rz(0.83807105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8720334) q[1];
sx q[1];
rz(-0.57231325) q[1];
sx q[1];
rz(2.3523) q[1];
rz(2.3039303) q[3];
sx q[3];
rz(-0.3627844) q[3];
sx q[3];
rz(-2.5183122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5088523) q[2];
sx q[2];
rz(-0.036157046) q[2];
sx q[2];
rz(0.73842326) q[2];
rz(-2.9449055) q[3];
sx q[3];
rz(-2.5236712) q[3];
sx q[3];
rz(0.64016199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7554373) q[0];
sx q[0];
rz(-0.25179395) q[0];
sx q[0];
rz(-3.1016896) q[0];
rz(-2.9207322) q[1];
sx q[1];
rz(-1.618229) q[1];
sx q[1];
rz(1.1562645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98469767) q[0];
sx q[0];
rz(-1.3067208) q[0];
sx q[0];
rz(-0.15021245) q[0];
rz(-pi) q[1];
rz(2.4100163) q[2];
sx q[2];
rz(-0.28645624) q[2];
sx q[2];
rz(-2.4713304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6917518) q[1];
sx q[1];
rz(-0.66206703) q[1];
sx q[1];
rz(-1.0074703) q[1];
x q[2];
rz(1.6764033) q[3];
sx q[3];
rz(-1.4560513) q[3];
sx q[3];
rz(0.70699837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7463189) q[2];
sx q[2];
rz(-2.7776182) q[2];
sx q[2];
rz(-0.054595646) q[2];
rz(0.78628457) q[3];
sx q[3];
rz(-1.1966642) q[3];
sx q[3];
rz(1.1597077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.38128024) q[0];
sx q[0];
rz(-2.3535643) q[0];
sx q[0];
rz(0.81083167) q[0];
rz(-2.6122818) q[1];
sx q[1];
rz(-1.4365173) q[1];
sx q[1];
rz(1.1108105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68303525) q[0];
sx q[0];
rz(-2.1330569) q[0];
sx q[0];
rz(-1.109615) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48804897) q[2];
sx q[2];
rz(-1.0637317) q[2];
sx q[2];
rz(1.6973073) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6880092) q[1];
sx q[1];
rz(-2.3613259) q[1];
sx q[1];
rz(-0.23201402) q[1];
rz(-1.9930028) q[3];
sx q[3];
rz(-0.65311382) q[3];
sx q[3];
rz(-2.4401668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3436053) q[2];
sx q[2];
rz(-0.56861773) q[2];
sx q[2];
rz(-0.24151754) q[2];
rz(-0.24164116) q[3];
sx q[3];
rz(-1.6499949) q[3];
sx q[3];
rz(2.9901166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3342313) q[0];
sx q[0];
rz(-1.3483028) q[0];
sx q[0];
rz(-2.469113) q[0];
rz(0.45202759) q[1];
sx q[1];
rz(-2.1990135) q[1];
sx q[1];
rz(2.6262851) q[1];
rz(-1.128059) q[2];
sx q[2];
rz(-1.5166225) q[2];
sx q[2];
rz(-0.33467334) q[2];
rz(0.47361278) q[3];
sx q[3];
rz(-2.0423741) q[3];
sx q[3];
rz(-3.1033677) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
