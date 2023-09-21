OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(-0.67396069) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90097839) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(-0.76529495) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8698679) q[2];
sx q[2];
rz(-0.50479111) q[2];
sx q[2];
rz(-1.7287849) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1311878) q[1];
sx q[1];
rz(-2.511575) q[1];
sx q[1];
rz(2.5351934) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0029293) q[3];
sx q[3];
rz(-1.1274459) q[3];
sx q[3];
rz(-1.5594304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6661466) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-0.95970884) q[0];
rz(1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9976881) q[0];
sx q[0];
rz(-1.578497) q[0];
sx q[0];
rz(1.5639923) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5487899) q[2];
sx q[2];
rz(-1.7823879) q[2];
sx q[2];
rz(1.6197268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5261425) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(-1.0326833) q[1];
rz(-pi) q[2];
rz(-2.0262358) q[3];
sx q[3];
rz(-0.028209837) q[3];
sx q[3];
rz(-1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6015357) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(2.9648932) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(1.7938991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5603148) q[0];
sx q[0];
rz(-0.53239765) q[0];
sx q[0];
rz(2.6530932) q[0];
x q[1];
rz(-1.548544) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(0.74860886) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2587535) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-0.12211166) q[1];
rz(-pi) q[2];
rz(-1.755581) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(-2.256306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.1536095) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.233498) q[1];
sx q[1];
rz(0.21534236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0523473) q[0];
sx q[0];
rz(-2.2627292) q[0];
sx q[0];
rz(-2.7034876) q[0];
rz(-pi) q[1];
rz(-2.1583546) q[2];
sx q[2];
rz(-1.431627) q[2];
sx q[2];
rz(0.61901865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.035186471) q[1];
sx q[1];
rz(-1.7654997) q[1];
sx q[1];
rz(0.10722843) q[1];
rz(2.4550291) q[3];
sx q[3];
rz(-0.6696223) q[3];
sx q[3];
rz(1.7945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5061491) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(0.62600342) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(2.4749277) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0502888) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(-1.9712649) q[0];
rz(0.35258106) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(2.7378766) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1946698) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(0.6823632) q[1];
rz(-pi) q[2];
rz(-1.2742217) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(-0.34542686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.8886245) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1154293) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.6113575) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513072) q[0];
sx q[0];
rz(-1.6437274) q[0];
sx q[0];
rz(-0.42071995) q[0];
rz(-2.7026664) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(-0.93026464) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6758319) q[1];
sx q[1];
rz(-0.45468802) q[1];
sx q[1];
rz(-2.0783706) q[1];
rz(1.8239162) q[3];
sx q[3];
rz(-1.6198006) q[3];
sx q[3];
rz(-2.8984836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(-1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(-2.706066) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(2.6307154) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2832527) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(-0.39985379) q[0];
rz(-pi) q[1];
rz(2.3235544) q[2];
sx q[2];
rz(-2.15089) q[2];
sx q[2];
rz(1.8028508) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.43607298) q[1];
sx q[1];
rz(-1.1735859) q[1];
sx q[1];
rz(0.7049837) q[1];
rz(-2.0914518) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(-1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997416) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(1.0808806) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-0.63527766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23112049) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(1.2088103) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.073457) q[2];
sx q[2];
rz(-1.4173696) q[2];
sx q[2];
rz(1.9286326) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(2.5177588) q[1];
rz(1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(1.7978096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(-2.2837158) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(-2.7776921) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(0.90973967) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1014935) q[0];
sx q[0];
rz(-1.8414521) q[0];
sx q[0];
rz(-2.8069242) q[0];
rz(-2.5932543) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-0.48834947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(-2.9305305) q[1];
rz(-pi) q[2];
rz(2.8839888) q[3];
sx q[3];
rz(-1.8766878) q[3];
sx q[3];
rz(-2.5739939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(2.6249028) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(2.2470078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54297011) q[0];
sx q[0];
rz(-2.2036836) q[0];
sx q[0];
rz(0.052368725) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1012494) q[2];
sx q[2];
rz(-1.7551127) q[2];
sx q[2];
rz(-1.7061526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4641061) q[1];
sx q[1];
rz(-1.2103401) q[1];
sx q[1];
rz(3.0739215) q[1];
rz(2.22088) q[3];
sx q[3];
rz(-1.6953141) q[3];
sx q[3];
rz(1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(-1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548303) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(0.54626089) q[2];
sx q[2];
rz(-0.98304521) q[2];
sx q[2];
rz(1.8084768) q[2];
rz(-1.7265504) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
