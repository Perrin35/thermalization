OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(-0.24554645) q[0];
sx q[0];
rz(3.088933) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(2.8007562) q[1];
sx q[1];
rz(10.513289) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3153305) q[0];
sx q[0];
rz(-1.1754987) q[0];
sx q[0];
rz(-0.068282337) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7126738) q[2];
sx q[2];
rz(-1.4617702) q[2];
sx q[2];
rz(3.1224868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46228483) q[1];
sx q[1];
rz(-2.4801697) q[1];
sx q[1];
rz(-2.8712621) q[1];
rz(-pi) q[2];
x q[2];
rz(1.151122) q[3];
sx q[3];
rz(-1.5670781) q[3];
sx q[3];
rz(0.17824717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6344305) q[2];
sx q[2];
rz(-1.7796702) q[2];
sx q[2];
rz(0.6081028) q[2];
rz(-1.1242584) q[3];
sx q[3];
rz(-1.9367155) q[3];
sx q[3];
rz(-0.89157909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41481498) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(2.538105) q[0];
rz(-0.51015774) q[1];
sx q[1];
rz(-2.3699103) q[1];
sx q[1];
rz(-2.1147494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054508) q[0];
sx q[0];
rz(-1.742995) q[0];
sx q[0];
rz(-1.4132176) q[0];
rz(2.449069) q[2];
sx q[2];
rz(-2.1998458) q[2];
sx q[2];
rz(2.4170827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4812117) q[1];
sx q[1];
rz(-0.91970316) q[1];
sx q[1];
rz(2.2269985) q[1];
rz(-pi) q[2];
rz(-2.233674) q[3];
sx q[3];
rz(-2.0072122) q[3];
sx q[3];
rz(-1.392257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5828731) q[2];
sx q[2];
rz(-1.1483973) q[2];
sx q[2];
rz(-0.64644512) q[2];
rz(1.2785814) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(1.6897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2090787) q[0];
sx q[0];
rz(-3.0113853) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(-0.34700829) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(-0.86403799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0542674) q[0];
sx q[0];
rz(-2.1313166) q[0];
sx q[0];
rz(1.818324) q[0];
rz(-1.8225602) q[2];
sx q[2];
rz(-1.1472817) q[2];
sx q[2];
rz(2.8734795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5123295) q[1];
sx q[1];
rz(-2.1844668) q[1];
sx q[1];
rz(-2.3113109) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7517463) q[3];
sx q[3];
rz(-2.6084508) q[3];
sx q[3];
rz(2.0651434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.049909264) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(-1.4570215) q[2];
rz(-1.0157887) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(-2.8252025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24093534) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(1.578791) q[0];
rz(2.0511625) q[1];
sx q[1];
rz(-0.54160392) q[1];
sx q[1];
rz(1.1515559) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59127677) q[0];
sx q[0];
rz(-1.341686) q[0];
sx q[0];
rz(0.59058779) q[0];
rz(-2.0456373) q[2];
sx q[2];
rz(-2.2730388) q[2];
sx q[2];
rz(2.4202731) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2673064) q[1];
sx q[1];
rz(-1.6055672) q[1];
sx q[1];
rz(1.0465996) q[1];
x q[2];
rz(-2.1776803) q[3];
sx q[3];
rz(-1.5480032) q[3];
sx q[3];
rz(-0.43302872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84732032) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(-1.6418183) q[2];
rz(-3.0070987) q[3];
sx q[3];
rz(-1.2765063) q[3];
sx q[3];
rz(-2.3597778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0791557) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(-0.4656747) q[0];
rz(1.8648719) q[1];
sx q[1];
rz(-2.1379505) q[1];
sx q[1];
rz(1.0328971) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55556923) q[0];
sx q[0];
rz(-0.15842552) q[0];
sx q[0];
rz(1.5044295) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14012248) q[2];
sx q[2];
rz(-1.4714141) q[2];
sx q[2];
rz(-1.6665719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3440282) q[1];
sx q[1];
rz(-0.92272341) q[1];
sx q[1];
rz(-1.2540884) q[1];
x q[2];
rz(0.9744076) q[3];
sx q[3];
rz(-2.5112409) q[3];
sx q[3];
rz(-2.3977365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5501962) q[2];
sx q[2];
rz(-1.952221) q[2];
sx q[2];
rz(2.2097394) q[2];
rz(-0.094680928) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(1.3100821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4513627) q[0];
sx q[0];
rz(-2.9424423) q[0];
sx q[0];
rz(3.0354101) q[0];
rz(-0.55446082) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(-1.5415446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0743177) q[0];
sx q[0];
rz(-1.5948926) q[0];
sx q[0];
rz(2.7510608) q[0];
rz(-pi) q[1];
rz(-1.0511398) q[2];
sx q[2];
rz(-2.6259632) q[2];
sx q[2];
rz(1.3517429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5145077) q[1];
sx q[1];
rz(-0.70620757) q[1];
sx q[1];
rz(0.18893623) q[1];
rz(-0.24644773) q[3];
sx q[3];
rz(-1.6250027) q[3];
sx q[3];
rz(-2.0230296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9083378) q[2];
sx q[2];
rz(-2.4007863) q[2];
sx q[2];
rz(-1.4976658) q[2];
rz(-2.6805367) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(-1.8259995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(2.0627956) q[0];
rz(2.5476593) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(2.914391) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543252) q[0];
sx q[0];
rz(-2.9552412) q[0];
sx q[0];
rz(0.92599286) q[0];
rz(-pi) q[1];
rz(-0.14102139) q[2];
sx q[2];
rz(-1.2407836) q[2];
sx q[2];
rz(-0.38244707) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0409483) q[1];
sx q[1];
rz(-0.64926636) q[1];
sx q[1];
rz(-1.0430286) q[1];
rz(-pi) q[2];
rz(1.9428857) q[3];
sx q[3];
rz(-1.6642906) q[3];
sx q[3];
rz(-1.0013415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6195153) q[2];
sx q[2];
rz(-0.41831133) q[2];
sx q[2];
rz(0.72959161) q[2];
rz(1.6884165) q[3];
sx q[3];
rz(-1.4540693) q[3];
sx q[3];
rz(2.5280473) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59298092) q[0];
sx q[0];
rz(-0.91658968) q[0];
sx q[0];
rz(-2.7899637) q[0];
rz(2.8278606) q[1];
sx q[1];
rz(-1.9351363) q[1];
sx q[1];
rz(-1.3765913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2257654) q[0];
sx q[0];
rz(-1.6179913) q[0];
sx q[0];
rz(2.7412358) q[0];
x q[1];
rz(-1.4046762) q[2];
sx q[2];
rz(-1.6848279) q[2];
sx q[2];
rz(1.5064552) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9326068) q[1];
sx q[1];
rz(-2.5521005) q[1];
sx q[1];
rz(-2.8578651) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3198765) q[3];
sx q[3];
rz(-2.7610215) q[3];
sx q[3];
rz(-1.6264834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54479638) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(-1.6205622) q[2];
rz(-0.14791402) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(-1.773905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.77383298) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.897478) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(2.921385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7871683) q[0];
sx q[0];
rz(-0.97499135) q[0];
sx q[0];
rz(0.45136772) q[0];
x q[1];
rz(-2.979109) q[2];
sx q[2];
rz(-2.0318609) q[2];
sx q[2];
rz(1.7562255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0988834) q[1];
sx q[1];
rz(-2.115887) q[1];
sx q[1];
rz(2.6429151) q[1];
rz(-pi) q[2];
rz(-2.3784786) q[3];
sx q[3];
rz(-3.0060351) q[3];
sx q[3];
rz(0.28597304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1844909) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(1.5000878) q[2];
rz(0.084065048) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(-0.070092289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641007) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(-0.41859928) q[0];
rz(0.94340008) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(0.30280608) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43234149) q[0];
sx q[0];
rz(-0.21638566) q[0];
sx q[0];
rz(-0.88309137) q[0];
rz(2.3934028) q[2];
sx q[2];
rz(-2.8837969) q[2];
sx q[2];
rz(2.2141544) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.23754072) q[1];
sx q[1];
rz(-1.6142323) q[1];
sx q[1];
rz(0.082708184) q[1];
x q[2];
rz(1.2262086) q[3];
sx q[3];
rz(-1.6804964) q[3];
sx q[3];
rz(-1.4668224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6454978) q[2];
sx q[2];
rz(-2.2690513) q[2];
sx q[2];
rz(2.3910451) q[2];
rz(1.4613072) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(2.6245978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208329) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(-1.9137406) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(1.3881793) q[2];
sx q[2];
rz(-2.1299405) q[2];
sx q[2];
rz(-2.401939) q[2];
rz(1.0629366) q[3];
sx q[3];
rz(-2.776317) q[3];
sx q[3];
rz(-1.1767514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
