OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7115241) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(0.67396069) q[0];
rz(-3.4586973) q[1];
sx q[1];
rz(4.7749333) q[1];
sx q[1];
rz(7.0778579) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7793517) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(-2.7428521) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48912666) q[2];
sx q[2];
rz(-1.4406275) q[2];
sx q[2];
rz(-0.39718539) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7197345) q[1];
sx q[1];
rz(-1.0654447) q[1];
sx q[1];
rz(-1.964633) q[1];
rz(2.0178954) q[3];
sx q[3];
rz(-1.6959794) q[3];
sx q[3];
rz(3.0931635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(-2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(2.5141292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9976881) q[0];
sx q[0];
rz(-1.5630957) q[0];
sx q[0];
rz(1.5639923) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10208315) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(-1.4174457) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97034772) q[1];
sx q[1];
rz(-1.8106318) q[1];
sx q[1];
rz(0.14648267) q[1];
rz(3.1291817) q[3];
sx q[3];
rz(-1.59613) q[3];
sx q[3];
rz(-0.94143553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.6905789) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(-1.6987945) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028570024) q[0];
sx q[0];
rz(-2.0355814) q[0];
sx q[0];
rz(1.3010498) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.548544) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(0.74860886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29979953) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(1.6714765) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.090431902) q[3];
sx q[3];
rz(-1.7548429) q[3];
sx q[3];
rz(2.4726766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(2.9262503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1911083) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(0.82975181) q[0];
rz(-1.8183069) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(1.9844696) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1064062) q[1];
sx q[1];
rz(-1.7654997) q[1];
sx q[1];
rz(3.0343642) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68656355) q[3];
sx q[3];
rz(-0.6696223) q[3];
sx q[3];
rz(1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5061491) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895806) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(2.1767298) q[0];
rz(-3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(2.4749277) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(1.9712649) q[0];
rz(-1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(-2.0999883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9469229) q[1];
sx q[1];
rz(-1.8538845) q[1];
sx q[1];
rz(-2.4592295) q[1];
x q[2];
rz(-1.8673709) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(0.34542686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91606402) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1154293) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.6113575) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513072) q[0];
sx q[0];
rz(-1.4978652) q[0];
sx q[0];
rz(0.42071995) q[0];
rz(-0.62263454) q[2];
sx q[2];
rz(-0.52402516) q[2];
sx q[2];
rz(-3.0571836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6758319) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(-1.0632221) q[1];
rz(-1.3176765) q[3];
sx q[3];
rz(-1.5217921) q[3];
sx q[3];
rz(-0.24310902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7851012) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-2.9587865) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(-2.6307154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683838) q[0];
sx q[0];
rz(-1.3084992) q[0];
sx q[0];
rz(-1.6845076) q[0];
rz(0.73166087) q[2];
sx q[2];
rz(-2.1795142) q[2];
sx q[2];
rz(-0.70639709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70799815) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(-0.57452332) q[1];
rz(-pi) q[2];
rz(-2.8227311) q[3];
sx q[3];
rz(-1.0705035) q[3];
sx q[3];
rz(1.6047516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-0.62057173) q[2];
rz(-0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-0.63527766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23112049) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(1.2088103) q[0];
x q[1];
rz(-2.9674171) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(0.27506405) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5105671) q[1];
sx q[1];
rz(-2.4706565) q[1];
sx q[1];
rz(2.7061694) q[1];
rz(-1.876272) q[3];
sx q[3];
rz(-0.59738041) q[3];
sx q[3];
rz(1.3437831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0044331) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(-2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(2.231853) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37659392) q[0];
sx q[0];
rz(-1.24876) q[0];
sx q[0];
rz(1.8565208) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54833834) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-0.48834947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4917131) q[1];
sx q[1];
rz(-0.24398206) q[1];
sx q[1];
rz(2.6073543) q[1];
rz(-pi) q[2];
rz(0.25760381) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(-1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(2.6249028) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-0.89458481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0588194) q[0];
sx q[0];
rz(-1.6130157) q[0];
sx q[0];
rz(-2.204338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9616227) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(0.48197907) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4641061) q[1];
sx q[1];
rz(-1.9312526) q[1];
sx q[1];
rz(3.0739215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15595943) q[3];
sx q[3];
rz(-2.2150063) q[3];
sx q[3];
rz(-0.19176602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(1.0220035) q[2];
rz(-1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5932896) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-0.54626089) q[2];
sx q[2];
rz(-2.1585474) q[2];
sx q[2];
rz(-1.3331158) q[2];
rz(0.048687497) q[3];
sx q[3];
rz(-1.2702474) q[3];
sx q[3];
rz(-2.4297759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
