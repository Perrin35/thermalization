OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.5716612) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(10.098739) q[0];
rz(-3.4586973) q[1];
sx q[1];
rz(4.7749333) q[1];
sx q[1];
rz(7.0778579) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3622409) q[0];
sx q[0];
rz(-1.2138214) q[0];
sx q[0];
rz(0.39874052) q[0];
x q[1];
rz(1.4235714) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(1.1046315) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7197345) q[1];
sx q[1];
rz(-1.0654447) q[1];
sx q[1];
rz(-1.964633) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13866339) q[3];
sx q[3];
rz(-2.0141467) q[3];
sx q[3];
rz(-1.5594304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871724) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7213631) q[0];
sx q[0];
rz(-3.1313167) q[0];
sx q[0];
rz(2.4179439) q[0];
rz(-pi) q[1];
rz(1.5928028) q[2];
sx q[2];
rz(-1.7823879) q[2];
sx q[2];
rz(1.6197268) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5761766) q[1];
sx q[1];
rz(-1.7130573) q[1];
sx q[1];
rz(1.3284645) q[1];
rz(-pi) q[2];
rz(-1.1153568) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(1.7445604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6015357) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.6905789) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
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
rz(-1.9243762) q[1];
sx q[1];
rz(-1.7938991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5603148) q[0];
sx q[0];
rz(-2.609195) q[0];
sx q[0];
rz(-0.48849948) q[0];
rz(-pi) q[1];
rz(-1.548544) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(0.74860886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8828391) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-0.12211166) q[1];
x q[2];
rz(-1.755581) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(-2.256306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.1536095) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(0.56458449) q[0];
rz(-3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-0.21534236) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9504844) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(-2.3118408) q[0];
rz(-pi) q[1];
rz(0.16673659) q[2];
sx q[2];
rz(-0.98966375) q[2];
sx q[2];
rz(2.2819448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1064062) q[1];
sx q[1];
rz(-1.7654997) q[1];
sx q[1];
rz(3.0343642) q[1];
rz(-2.4550291) q[3];
sx q[3];
rz(-2.4719704) q[3];
sx q[3];
rz(-1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5061491) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(0.60194683) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(0.66666493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0502888) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(-1.9712649) q[0];
x q[1];
rz(1.0235734) q[2];
sx q[2];
rz(-1.2663411) q[2];
sx q[2];
rz(-2.1551702) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1866245) q[1];
sx q[1];
rz(-2.4116227) q[1];
sx q[1];
rz(2.7093922) q[1];
rz(-pi) q[2];
rz(-0.93249647) q[3];
sx q[3];
rz(-2.7780048) q[3];
sx q[3];
rz(0.61908412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(1.3145674) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(-1.6113575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0936733) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(1.6506667) q[0];
x q[1];
rz(-0.43892626) q[2];
sx q[2];
rz(-1.274684) q[2];
sx q[2];
rz(-0.93026464) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1215026) q[1];
sx q[1];
rz(-1.1768747) q[1];
sx q[1];
rz(-0.23328383) q[1];
rz(3.0909782) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(1.3150172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(-1.2221378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71762639) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(2.9587865) q[0];
rz(-2.706066) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85833997) q[0];
sx q[0];
rz(-0.28536277) q[0];
sx q[0];
rz(2.7417389) q[0];
x q[1];
rz(-0.81803825) q[2];
sx q[2];
rz(-0.99070264) q[2];
sx q[2];
rz(1.3387418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43607298) q[1];
sx q[1];
rz(-1.1735859) q[1];
sx q[1];
rz(2.436609) q[1];
x q[2];
rz(2.0914518) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(-2.1396162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(-0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.871199) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(2.0525232) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(0.63527766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3978183) q[0];
sx q[0];
rz(-1.213307) q[0];
sx q[0];
rz(-2.9767569) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17417553) q[2];
sx q[2];
rz(-1.0798228) q[2];
sx q[2];
rz(0.27506405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5105671) q[1];
sx q[1];
rz(-0.67093611) q[1];
sx q[1];
rz(-0.4354233) q[1];
x q[2];
rz(-1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(1.3437831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(-0.00017246406) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8453318) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(-2.231853) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7649987) q[0];
sx q[0];
rz(-1.8928327) q[0];
sx q[0];
rz(1.8565208) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54833834) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-2.6532432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.197387) q[1];
sx q[1];
rz(-1.7802317) q[1];
sx q[1];
rz(-1.6968813) q[1];
x q[2];
rz(1.886456) q[3];
sx q[3];
rz(-1.8161895) q[3];
sx q[3];
rz(-1.0823702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6870118) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(-1.6420341) q[0];
x q[1];
rz(-1.1012494) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(1.7061526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67748653) q[1];
sx q[1];
rz(-1.9312526) q[1];
sx q[1];
rz(-3.0739215) q[1];
x q[2];
rz(-0.92071269) q[3];
sx q[3];
rz(-1.6953141) q[3];
sx q[3];
rz(1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.2330166) q[2];
sx q[2];
rz(-0.77975811) q[2];
sx q[2];
rz(-2.1644885) q[2];
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
