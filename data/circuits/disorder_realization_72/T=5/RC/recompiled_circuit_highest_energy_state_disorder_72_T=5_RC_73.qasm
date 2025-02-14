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
rz(-2.8895145) q[0];
sx q[0];
rz(-0.57555389) q[0];
sx q[0];
rz(-3.0192896) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(-3.1032739) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0487028) q[0];
sx q[0];
rz(-0.99703353) q[0];
sx q[0];
rz(2.696373) q[0];
x q[1];
rz(-2.3549809) q[2];
sx q[2];
rz(-0.16737882) q[2];
sx q[2];
rz(-2.6712863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.029509228) q[1];
sx q[1];
rz(-1.2473628) q[1];
sx q[1];
rz(-2.1548163) q[1];
rz(-pi) q[2];
rz(0.081101426) q[3];
sx q[3];
rz(-2.8934294) q[3];
sx q[3];
rz(-2.1901772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4965839) q[2];
sx q[2];
rz(-1.714548) q[2];
sx q[2];
rz(1.4487779) q[2];
rz(2.951176) q[3];
sx q[3];
rz(-2.0864291) q[3];
sx q[3];
rz(0.14934389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377624) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(0.25800905) q[0];
rz(-1.5915271) q[1];
sx q[1];
rz(-1.240088) q[1];
sx q[1];
rz(2.9749427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44010559) q[0];
sx q[0];
rz(-2.1512554) q[0];
sx q[0];
rz(1.6850182) q[0];
x q[1];
rz(-0.3853674) q[2];
sx q[2];
rz(-1.8358942) q[2];
sx q[2];
rz(0.4601269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8434324) q[1];
sx q[1];
rz(-2.2800804) q[1];
sx q[1];
rz(-1.7357566) q[1];
x q[2];
rz(-2.3349668) q[3];
sx q[3];
rz(-0.72128937) q[3];
sx q[3];
rz(-1.862135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0698645) q[2];
sx q[2];
rz(-2.0191777) q[2];
sx q[2];
rz(1.2321164) q[2];
rz(-2.2169436) q[3];
sx q[3];
rz(-2.2053714) q[3];
sx q[3];
rz(-1.1054976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6919747) q[0];
sx q[0];
rz(-1.6716577) q[0];
sx q[0];
rz(0.0019419226) q[0];
rz(-0.034189668) q[1];
sx q[1];
rz(-1.2119774) q[1];
sx q[1];
rz(-1.5431822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7352075) q[0];
sx q[0];
rz(-1.8282561) q[0];
sx q[0];
rz(2.6763112) q[0];
rz(2.4056881) q[2];
sx q[2];
rz(-1.2703478) q[2];
sx q[2];
rz(-1.7395333) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1605121) q[1];
sx q[1];
rz(-2.7857669) q[1];
sx q[1];
rz(0.11248223) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9228306) q[3];
sx q[3];
rz(-1.3746972) q[3];
sx q[3];
rz(0.01250532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89692846) q[2];
sx q[2];
rz(-0.43411532) q[2];
sx q[2];
rz(1.0373235) q[2];
rz(1.1050998) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(2.0028116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26745519) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(1.4111891) q[0];
rz(-2.5022068) q[1];
sx q[1];
rz(-1.4566028) q[1];
sx q[1];
rz(-0.04714084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5306803) q[0];
sx q[0];
rz(-2.147104) q[0];
sx q[0];
rz(-1.2333825) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.669692) q[2];
sx q[2];
rz(-1.7444897) q[2];
sx q[2];
rz(1.398441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6005391) q[1];
sx q[1];
rz(-2.4324634) q[1];
sx q[1];
rz(-2.1313271) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.239226) q[3];
sx q[3];
rz(-2.3904388) q[3];
sx q[3];
rz(1.0855261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3186657) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(1.9226496) q[2];
rz(-2.8259891) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8197935) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(0.34859443) q[0];
rz(1.2527342) q[1];
sx q[1];
rz(-1.9786973) q[1];
sx q[1];
rz(1.8399651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0785261) q[0];
sx q[0];
rz(-1.491703) q[0];
sx q[0];
rz(0.80711676) q[0];
rz(-pi) q[1];
rz(-2.581004) q[2];
sx q[2];
rz(-1.212767) q[2];
sx q[2];
rz(2.5303417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0685454) q[1];
sx q[1];
rz(-1.3497735) q[1];
sx q[1];
rz(-0.18494341) q[1];
rz(2.1215277) q[3];
sx q[3];
rz(-1.2675084) q[3];
sx q[3];
rz(1.0098977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9584413) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(1.4754254) q[2];
rz(-2.4557377) q[3];
sx q[3];
rz(-0.97969046) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1118065) q[0];
sx q[0];
rz(-2.989558) q[0];
sx q[0];
rz(1.6492122) q[0];
rz(2.1624508) q[1];
sx q[1];
rz(-1.5359595) q[1];
sx q[1];
rz(1.917256) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92517744) q[0];
sx q[0];
rz(-1.6013711) q[0];
sx q[0];
rz(-2.9435817) q[0];
rz(-0.54927214) q[2];
sx q[2];
rz(-1.7584929) q[2];
sx q[2];
rz(0.21085462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8679132) q[1];
sx q[1];
rz(-2.1845594) q[1];
sx q[1];
rz(-0.015271827) q[1];
rz(-pi) q[2];
rz(2.200279) q[3];
sx q[3];
rz(-1.286003) q[3];
sx q[3];
rz(0.79508699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7225723) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(-0.9160308) q[2];
rz(-1.7570868) q[3];
sx q[3];
rz(-1.2576831) q[3];
sx q[3];
rz(-0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0093805669) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(0.95426553) q[0];
rz(-0.43680278) q[1];
sx q[1];
rz(-1.5780459) q[1];
sx q[1];
rz(3.0314235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6340652) q[0];
sx q[0];
rz(-1.8095836) q[0];
sx q[0];
rz(0.33383835) q[0];
rz(1.0344347) q[2];
sx q[2];
rz(-2.1948994) q[2];
sx q[2];
rz(-1.0076866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72790775) q[1];
sx q[1];
rz(-1.9164819) q[1];
sx q[1];
rz(3.137008) q[1];
x q[2];
rz(2.0131936) q[3];
sx q[3];
rz(-2.5964649) q[3];
sx q[3];
rz(-0.21827182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7097912) q[2];
sx q[2];
rz(-0.0063889901) q[2];
sx q[2];
rz(2.1978281) q[2];
rz(0.56378311) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(-2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875882) q[0];
sx q[0];
rz(-2.8546951) q[0];
sx q[0];
rz(2.7960844) q[0];
rz(-1.0954789) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(1.5798205) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51673543) q[0];
sx q[0];
rz(-2.3586732) q[0];
sx q[0];
rz(1.5262414) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.079633) q[2];
sx q[2];
rz(-2.2270348) q[2];
sx q[2];
rz(-1.7554754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4862389) q[1];
sx q[1];
rz(-1.2062688) q[1];
sx q[1];
rz(2.2097387) q[1];
rz(3.0828219) q[3];
sx q[3];
rz(-0.49859014) q[3];
sx q[3];
rz(2.9166247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3788508) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(1.4307107) q[2];
rz(-0.53449574) q[3];
sx q[3];
rz(-1.4662687) q[3];
sx q[3];
rz(-2.1889595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6119824) q[0];
sx q[0];
rz(-0.069267608) q[0];
sx q[0];
rz(-1.3105422) q[0];
rz(-2.2546841) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(-1.8720253) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16989141) q[0];
sx q[0];
rz(-1.9215688) q[0];
sx q[0];
rz(0.13668513) q[0];
x q[1];
rz(1.194095) q[2];
sx q[2];
rz(-1.4262876) q[2];
sx q[2];
rz(1.2928177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7901526) q[1];
sx q[1];
rz(-1.2987441) q[1];
sx q[1];
rz(-0.2403918) q[1];
rz(-pi) q[2];
rz(-2.9379775) q[3];
sx q[3];
rz(-1.609966) q[3];
sx q[3];
rz(0.43653442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.612192) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(0.23294918) q[2];
rz(-1.285078) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9832298) q[0];
sx q[0];
rz(-2.5322999) q[0];
sx q[0];
rz(0.097271517) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-0.709788) q[1];
sx q[1];
rz(0.21673094) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6410206) q[0];
sx q[0];
rz(-1.4002698) q[0];
sx q[0];
rz(3.0855623) q[0];
rz(-pi) q[1];
rz(0.31483908) q[2];
sx q[2];
rz(-2.9926569) q[2];
sx q[2];
rz(2.6268132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7864259) q[1];
sx q[1];
rz(-2.2192419) q[1];
sx q[1];
rz(0.72079682) q[1];
x q[2];
rz(-2.7085767) q[3];
sx q[3];
rz(-2.4320514) q[3];
sx q[3];
rz(0.39329986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8486166) q[2];
sx q[2];
rz(-0.28102195) q[2];
sx q[2];
rz(-2.6463553) q[2];
rz(-1.5589335) q[3];
sx q[3];
rz(-1.0844743) q[3];
sx q[3];
rz(-0.16286477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611298) q[0];
sx q[0];
rz(-1.6292138) q[0];
sx q[0];
rz(-0.52393352) q[0];
rz(1.0864661) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(-1.6025966) q[2];
sx q[2];
rz(-1.5780137) q[2];
sx q[2];
rz(0.93102166) q[2];
rz(-1.466352) q[3];
sx q[3];
rz(-0.71487311) q[3];
sx q[3];
rz(-2.5345595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
