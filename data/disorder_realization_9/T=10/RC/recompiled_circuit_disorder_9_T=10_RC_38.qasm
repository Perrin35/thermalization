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
rz(-0.79467264) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3622409) q[0];
sx q[0];
rz(-1.9277713) q[0];
sx q[0];
rz(-0.39874052) q[0];
x q[1];
rz(-2.652466) q[2];
sx q[2];
rz(-1.4406275) q[2];
sx q[2];
rz(2.7444073) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1311878) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(-0.6063993) q[1];
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
rz(pi/2) q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(-2.5141292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1439046) q[0];
sx q[0];
rz(-1.5630957) q[0];
sx q[0];
rz(1.5776004) q[0];
x q[1];
rz(-3.0395095) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(-1.724147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1712449) q[1];
sx q[1];
rz(-1.3309609) q[1];
sx q[1];
rz(2.99511) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1153568) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54005694) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(0.97989782) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.7938991) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130226) q[0];
sx q[0];
rz(-2.0355814) q[0];
sx q[0];
rz(-1.8405429) q[0];
rz(-pi) q[1];
rz(3.0323896) q[2];
sx q[2];
rz(-0.20143992) q[2];
sx q[2];
rz(-0.86004721) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2587535) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-0.12211166) q[1];
x q[2];
rz(2.0224781) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(-1.1289489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6692114) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4190061) q[0];
sx q[0];
rz(-0.79918396) q[0];
sx q[0];
rz(-1.0976085) q[0];
rz(-2.9748561) q[2];
sx q[2];
rz(-2.1519289) q[2];
sx q[2];
rz(0.8596479) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5986961) q[1];
sx q[1];
rz(-2.9196432) q[1];
sx q[1];
rz(2.0680244) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68656355) q[3];
sx q[3];
rz(-2.4719704) q[3];
sx q[3];
rz(-1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(0.60194683) q[2];
rz(-0.69058949) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-0.66666493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71320888) q[0];
sx q[0];
rz(-2.7292626) q[0];
sx q[0];
rz(-1.8250188) q[0];
x q[1];
rz(0.35258106) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(2.7378766) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95496817) q[1];
sx q[1];
rz(-0.72997) q[1];
sx q[1];
rz(-0.43220046) q[1];
rz(-0.22294754) q[3];
sx q[3];
rz(-1.8604391) q[3];
sx q[3];
rz(1.2900316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-1.130828) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.6113575) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0479193) q[0];
sx q[0];
rz(-1.9903272) q[0];
sx q[0];
rz(1.6506667) q[0];
x q[1];
rz(1.8958695) q[2];
sx q[2];
rz(-1.1522066) q[2];
sx q[2];
rz(-2.3649154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1215026) q[1];
sx q[1];
rz(-1.964718) q[1];
sx q[1];
rz(-0.23328383) q[1];
x q[2];
rz(1.3176765) q[3];
sx q[3];
rz(-1.6198006) q[3];
sx q[3];
rz(2.8984836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71762639) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683838) q[0];
sx q[0];
rz(-1.8330935) q[0];
sx q[0];
rz(-1.6845076) q[0];
rz(-pi) q[1];
rz(-2.3350231) q[2];
sx q[2];
rz(-0.91420805) q[2];
sx q[2];
rz(2.8442531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4527013) q[1];
sx q[1];
rz(-2.211262) q[1];
sx q[1];
rz(-2.074261) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0914518) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(-0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1875293) q[0];
sx q[0];
rz(-2.7494193) q[0];
sx q[0];
rz(-1.1568882) q[0];
rz(-1.8842472) q[2];
sx q[2];
rz(-0.51856504) q[2];
sx q[2];
rz(-0.6322565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2891846) q[1];
sx q[1];
rz(-1.305456) q[1];
sx q[1];
rz(2.5177588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99526309) q[3];
sx q[3];
rz(-1.4008153) q[3];
sx q[3];
rz(0.02804027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13715956) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(0.90973967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7649987) q[0];
sx q[0];
rz(-1.24876) q[0];
sx q[0];
rz(1.2850719) q[0];
rz(-pi) q[1];
rz(0.34809525) q[2];
sx q[2];
rz(-1.7762134) q[2];
sx q[2];
rz(1.5945438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(0.21106212) q[1];
rz(-pi) q[2];
rz(0.89197253) q[3];
sx q[3];
rz(-0.39728764) q[3];
sx q[3];
rz(2.9904423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.9177115) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0588194) q[0];
sx q[0];
rz(-1.528577) q[0];
sx q[0];
rz(-0.93725462) q[0];
rz(-pi) q[1];
rz(-2.9355029) q[2];
sx q[2];
rz(-1.1098252) q[2];
sx q[2];
rz(-0.042630171) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6539508) q[1];
sx q[1];
rz(-2.7751121) q[1];
sx q[1];
rz(1.3932863) q[1];
rz(-pi) q[2];
rz(-0.92071269) q[3];
sx q[3];
rz(-1.4462785) q[3];
sx q[3];
rz(1.4731821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-2.5125304) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(2.2330166) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
rz(1.2699119) q[3];
sx q[3];
rz(-1.5242929) q[3];
sx q[3];
rz(2.268189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
