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
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7793517) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(-2.7428521) q[0];
rz(1.7180213) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(-1.1046315) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1311878) q[1];
sx q[1];
rz(-2.511575) q[1];
sx q[1];
rz(-2.5351934) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0178954) q[3];
sx q[3];
rz(-1.4456133) q[3];
sx q[3];
rz(3.0931635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(0.62746343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146485) q[0];
sx q[0];
rz(-1.5639925) q[0];
sx q[0];
rz(3.1338918) q[0];
rz(-0.21164125) q[2];
sx q[2];
rz(-1.5492808) q[2];
sx q[2];
rz(3.0880398) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97034772) q[1];
sx q[1];
rz(-1.3309609) q[1];
sx q[1];
rz(-2.99511) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0262358) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(-1.7445604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6015357) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.7938991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5603148) q[0];
sx q[0];
rz(-0.53239765) q[0];
sx q[0];
rz(-2.6530932) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10920306) q[2];
sx q[2];
rz(-2.9401527) q[2];
sx q[2];
rz(-2.2815454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8417931) q[1];
sx q[1];
rz(-1.6922957) q[1];
sx q[1];
rz(-1.6714765) q[1];
rz(-1.1191145) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(-1.1289489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(-1.1536095) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(0.56458449) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(2.9262503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72258654) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(-1.0976085) q[0];
rz(-pi) q[1];
rz(-1.8183069) q[2];
sx q[2];
rz(-0.60192054) q[2];
sx q[2];
rz(1.9844696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54289651) q[1];
sx q[1];
rz(-2.9196432) q[1];
sx q[1];
rz(-1.0735682) q[1];
x q[2];
rz(-1.1057042) q[3];
sx q[3];
rz(-1.0700738) q[3];
sx q[3];
rz(-0.53946686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85201207) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(0.96486282) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-2.4749277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.6717523) q[0];
sx q[0];
rz(-1.9712649) q[0];
rz(-pi) q[1];
rz(-1.0235734) q[2];
sx q[2];
rz(-1.2663411) q[2];
sx q[2];
rz(-0.98642245) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95496817) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(-2.7093922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2090962) q[3];
sx q[3];
rz(-2.7780048) q[3];
sx q[3];
rz(-0.61908412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.2529681) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.7793659) q[0];
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
rz(-2.0936733) q[0];
sx q[0];
rz(-1.1512655) q[0];
sx q[0];
rz(-1.6506667) q[0];
rz(-pi) q[1];
rz(0.43892626) q[2];
sx q[2];
rz(-1.274684) q[2];
sx q[2];
rz(0.93026464) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6758319) q[1];
sx q[1];
rz(-0.45468802) q[1];
sx q[1];
rz(-1.0632221) q[1];
x q[2];
rz(-0.050614428) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(-1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(1.6507089) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(-2.6307154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85833997) q[0];
sx q[0];
rz(-0.28536277) q[0];
sx q[0];
rz(2.7417389) q[0];
x q[1];
rz(2.3235544) q[2];
sx q[2];
rz(-0.99070264) q[2];
sx q[2];
rz(1.3387418) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70799815) q[1];
sx q[1];
rz(-0.79213789) q[1];
sx q[1];
rz(2.5670693) q[1];
x q[2];
rz(2.8227311) q[3];
sx q[3];
rz(-1.0705035) q[3];
sx q[3];
rz(-1.6047516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(-0.62057173) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(-1.0808806) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(-2.506315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23112049) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(1.9327823) q[0];
x q[1];
rz(1.8842472) q[2];
sx q[2];
rz(-0.51856504) q[2];
sx q[2];
rz(0.6322565) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63102555) q[1];
sx q[1];
rz(-0.67093611) q[1];
sx q[1];
rz(2.7061694) q[1];
rz(-pi) q[2];
x q[2];
rz(1.876272) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(-1.7978096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(-2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(2.231853) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1247524) q[0];
sx q[0];
rz(-2.7144055) q[0];
sx q[0];
rz(-0.70144002) q[0];
rz(2.5932543) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-2.6532432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7418356) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(0.21106212) q[1];
rz(-pi) q[2];
rz(2.8839888) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(2.5739939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412823) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(2.6249028) q[0];
rz(1.5453045) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986225) q[0];
sx q[0];
rz(-2.2036836) q[0];
sx q[0];
rz(3.0892239) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0403433) q[2];
sx q[2];
rz(-1.7551127) q[2];
sx q[2];
rz(1.4354401) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91720944) q[1];
sx q[1];
rz(-1.6341126) q[1];
sx q[1];
rz(-1.2095832) q[1];
rz(-2.22088) q[3];
sx q[3];
rz(-1.4462785) q[3];
sx q[3];
rz(1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14780012) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548303) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-0.9085761) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
rz(-1.4150423) q[3];
sx q[3];
rz(-2.8372436) q[3];
sx q[3];
rz(-2.5929034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];