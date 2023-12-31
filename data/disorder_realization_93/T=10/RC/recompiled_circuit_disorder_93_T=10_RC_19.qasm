OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(4.5259024) q[0];
sx q[0];
rz(10.685267) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(0.92372149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0246575) q[0];
sx q[0];
rz(-0.57514656) q[0];
sx q[0];
rz(-0.92098178) q[0];
x q[1];
rz(1.2878296) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(1.1793009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0337861) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(-2.4297907) q[1];
x q[2];
rz(1.9507017) q[3];
sx q[3];
rz(-1.1067179) q[3];
sx q[3];
rz(2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(-2.2062733) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.4555567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992422) q[0];
sx q[0];
rz(-2.143321) q[0];
sx q[0];
rz(2.9426129) q[0];
x q[1];
rz(0.37462072) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(2.3945216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93745366) q[1];
sx q[1];
rz(-1.0480282) q[1];
sx q[1];
rz(3.1259895) q[1];
rz(-0.060828408) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(-0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(-3.1128186) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5867509) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(-0.59535938) q[0];
x q[1];
rz(-0.34611361) q[2];
sx q[2];
rz(-2.3556404) q[2];
sx q[2];
rz(-1.9217938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2492003) q[1];
sx q[1];
rz(-1.8141659) q[1];
sx q[1];
rz(-1.0358441) q[1];
rz(1.8215239) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(2.1123321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-2.6180843) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1062766) q[0];
sx q[0];
rz(-1.9395394) q[0];
sx q[0];
rz(-0.62555255) q[0];
x q[1];
rz(1.4857616) q[2];
sx q[2];
rz(-1.8997314) q[2];
sx q[2];
rz(-1.0852244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88672968) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(0.32943326) q[1];
rz(-pi) q[2];
rz(1.773049) q[3];
sx q[3];
rz(-2.5407255) q[3];
sx q[3];
rz(2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52544242) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(-0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(2.1496444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46734738) q[0];
sx q[0];
rz(-1.464198) q[0];
sx q[0];
rz(0.17515134) q[0];
rz(-1.5422103) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(-2.6945393) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3494898) q[1];
sx q[1];
rz(-1.2621242) q[1];
sx q[1];
rz(1.5285138) q[1];
x q[2];
rz(-0.81327849) q[3];
sx q[3];
rz(-1.6653898) q[3];
sx q[3];
rz(0.41155848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7871008) q[0];
sx q[0];
rz(-2.0658501) q[0];
sx q[0];
rz(1.8255193) q[0];
rz(-pi) q[1];
rz(1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(-0.5459107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6188366) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(-1.418581) q[1];
rz(-pi) q[2];
rz(-0.96958843) q[3];
sx q[3];
rz(-1.117327) q[3];
sx q[3];
rz(-0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0075334) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-0.82180506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215295) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(-2.2277742) q[0];
rz(-pi) q[1];
rz(0.84341151) q[2];
sx q[2];
rz(-1.1864098) q[2];
sx q[2];
rz(2.9316528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0710443) q[1];
sx q[1];
rz(-0.76862915) q[1];
sx q[1];
rz(2.8119836) q[1];
rz(-pi) q[2];
rz(-2.5542198) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(1.3164933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(1.3051422) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236429) q[0];
sx q[0];
rz(-0.68730132) q[0];
sx q[0];
rz(-2.6108517) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9403946) q[2];
sx q[2];
rz(-2.2410789) q[2];
sx q[2];
rz(3.0073462) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95272428) q[1];
sx q[1];
rz(-0.75718588) q[1];
sx q[1];
rz(2.5785179) q[1];
x q[2];
rz(1.6373789) q[3];
sx q[3];
rz(-2.801602) q[3];
sx q[3];
rz(0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.6513599) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7548783) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(0.3219147) q[0];
rz(1.5362668) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(2.4386491) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2962869) q[0];
sx q[0];
rz(-1.9585113) q[0];
sx q[0];
rz(1.182343) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0614971) q[2];
sx q[2];
rz(-1.6054389) q[2];
sx q[2];
rz(2.408037) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1785537) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(-1.8815243) q[1];
x q[2];
rz(2.8076257) q[3];
sx q[3];
rz(-2.1595862) q[3];
sx q[3];
rz(-0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(-1.2257858) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(-2.6729565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2221453) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(-2.749445) q[0];
x q[1];
rz(-1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(2.7591443) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.66099) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(-1.3165228) q[1];
rz(-0.81810276) q[3];
sx q[3];
rz(-2.5460498) q[3];
sx q[3];
rz(0.18225741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(0.62190965) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-0.86482277) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(1.4680396) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
