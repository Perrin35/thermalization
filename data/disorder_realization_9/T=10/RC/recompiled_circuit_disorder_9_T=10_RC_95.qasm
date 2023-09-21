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
rz(2.467632) q[0];
rz(-3.4586973) q[1];
sx q[1];
rz(4.7749333) q[1];
sx q[1];
rz(7.0778579) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2406143) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(2.3762977) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27172471) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(-1.4128078) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0104048) q[1];
sx q[1];
rz(-2.511575) q[1];
sx q[1];
rz(-2.5351934) q[1];
rz(-pi) q[2];
rz(-1.1236973) q[3];
sx q[3];
rz(-1.4456133) q[3];
sx q[3];
rz(-3.0931635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-2.5141292) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146485) q[0];
sx q[0];
rz(-1.5776002) q[0];
sx q[0];
rz(0.0077008458) q[0];
x q[1];
rz(-0.21164125) q[2];
sx q[2];
rz(-1.5923119) q[2];
sx q[2];
rz(-3.0880398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1712449) q[1];
sx q[1];
rz(-1.3309609) q[1];
sx q[1];
rz(-2.99511) q[1];
x q[2];
rz(-1.1153568) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(-1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6015357) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(-1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5805805) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.7938991) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5603148) q[0];
sx q[0];
rz(-2.609195) q[0];
sx q[0];
rz(-0.48849948) q[0];
x q[1];
rz(1.5930487) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(-2.3929838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2587535) q[1];
sx q[1];
rz(-1.6707318) q[1];
sx q[1];
rz(-3.019481) q[1];
rz(1.1191145) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(-2.0126437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(1.9879831) q[2];
rz(2.7518318) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089245361) q[0];
sx q[0];
rz(-0.87886341) q[0];
sx q[0];
rz(-2.7034876) q[0];
rz(0.98323804) q[2];
sx q[2];
rz(-1.7099656) q[2];
sx q[2];
rz(2.522574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0680244) q[1];
x q[2];
rz(-0.68656355) q[3];
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
rz(-2.052867) q[2];
sx q[2];
rz(-0.60194683) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(2.1767298) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(2.4749277) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0502888) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(1.9712649) q[0];
rz(0.35258106) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(2.7378766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95496817) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(0.43220046) q[1];
x q[2];
rz(-0.22294754) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(1.851561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(-1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
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
rz(1.0479193) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(1.6506667) q[0];
x q[1];
rz(-2.5189581) q[2];
sx q[2];
rz(-0.52402516) q[2];
sx q[2];
rz(-0.084409075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46576071) q[1];
sx q[1];
rz(-0.45468802) q[1];
sx q[1];
rz(1.0632221) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7641894) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(1.5148439) q[3];
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
rz(3.1121758) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(-1.2221378) q[3];
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
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144007) q[0];
sx q[0];
rz(-1.6806024) q[0];
sx q[0];
rz(-0.26392428) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4099318) q[2];
sx q[2];
rz(-0.96207843) q[2];
sx q[2];
rz(-0.70639709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4335945) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(-2.5670693) q[1];
x q[2];
rz(0.31886153) q[3];
sx q[3];
rz(-1.0705035) q[3];
sx q[3];
rz(-1.5368411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9110979) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(-1.0808806) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3978183) q[0];
sx q[0];
rz(-1.213307) q[0];
sx q[0];
rz(-2.9767569) q[0];
x q[1];
rz(0.17417553) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(-2.8665286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0465225) q[1];
sx q[1];
rz(-2.1696643) q[1];
sx q[1];
rz(-1.2477161) q[1];
x q[2];
rz(1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(1.7978096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0044331) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37659392) q[0];
sx q[0];
rz(-1.24876) q[0];
sx q[0];
rz(-1.2850719) q[0];
rz(-pi) q[1];
rz(-0.34809525) q[2];
sx q[2];
rz(-1.3653792) q[2];
sx q[2];
rz(1.5945438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(2.9305305) q[1];
x q[2];
rz(-0.25760381) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(-0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.9177115) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-0.51668984) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(-0.89458481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45458083) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(1.6420341) q[0];
x q[1];
rz(-2.0403433) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(1.4354401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91720944) q[1];
sx q[1];
rz(-1.50748) q[1];
sx q[1];
rz(1.9320095) q[1];
rz(1.3668725) q[3];
sx q[3];
rz(-2.4813934) q[3];
sx q[3];
rz(-3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(-1.8036802) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(-0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-2.2330877) q[2];
sx q[2];
rz(-1.1237332) q[2];
sx q[2];
rz(-0.087469812) q[2];
rz(1.4150423) q[3];
sx q[3];
rz(-0.30434904) q[3];
sx q[3];
rz(0.54868922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];