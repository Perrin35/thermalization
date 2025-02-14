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
rz(2.8960462) q[0];
sx q[0];
rz(9.4774376) q[0];
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
rz(0.8262622) q[0];
sx q[0];
rz(-1.966094) q[0];
sx q[0];
rz(3.0733103) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4289189) q[2];
sx q[2];
rz(-1.6798225) q[2];
sx q[2];
rz(-3.1224868) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3237985) q[1];
sx q[1];
rz(-1.7355738) q[1];
sx q[1];
rz(2.4980157) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1375211) q[3];
sx q[3];
rz(-1.9904676) q[3];
sx q[3];
rz(1.3908902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6344305) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(-0.6081028) q[2];
rz(1.1242584) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(2.2500136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267777) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(0.60348764) q[0];
rz(0.51015774) q[1];
sx q[1];
rz(-0.77168232) q[1];
sx q[1];
rz(-2.1147494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.054508) q[0];
sx q[0];
rz(-1.3985976) q[0];
sx q[0];
rz(1.728375) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6925236) q[2];
sx q[2];
rz(-2.1998458) q[2];
sx q[2];
rz(2.4170827) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3852147) q[1];
sx q[1];
rz(-0.88895386) q[1];
sx q[1];
rz(-2.4663976) q[1];
rz(0.90791865) q[3];
sx q[3];
rz(-1.1343805) q[3];
sx q[3];
rz(-1.7493356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55871955) q[2];
sx q[2];
rz(-1.9931953) q[2];
sx q[2];
rz(2.4951475) q[2];
rz(1.8630113) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(-1.6897374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2090787) q[0];
sx q[0];
rz(-3.0113853) q[0];
sx q[0];
rz(-1.5694438) q[0];
rz(0.34700829) q[1];
sx q[1];
rz(-1.4748814) q[1];
sx q[1];
rz(2.2775547) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65007675) q[0];
sx q[0];
rz(-1.361761) q[0];
sx q[0];
rz(0.57455008) q[0];
rz(-pi) q[1];
rz(1.8225602) q[2];
sx q[2];
rz(-1.994311) q[2];
sx q[2];
rz(-0.26811312) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7622704) q[1];
sx q[1];
rz(-0.92298365) q[1];
sx q[1];
rz(-0.76384441) q[1];
rz(-pi) q[2];
rz(-3.0357843) q[3];
sx q[3];
rz(-2.0943301) q[3];
sx q[3];
rz(0.86712718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.049909264) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(1.6845711) q[2];
rz(-2.125804) q[3];
sx q[3];
rz(-0.53272811) q[3];
sx q[3];
rz(-2.8252025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006573) q[0];
sx q[0];
rz(-0.37739402) q[0];
sx q[0];
rz(1.578791) q[0];
rz(2.0511625) q[1];
sx q[1];
rz(-2.5999887) q[1];
sx q[1];
rz(-1.1515559) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59127677) q[0];
sx q[0];
rz(-1.341686) q[0];
sx q[0];
rz(2.5510049) q[0];
x q[1];
rz(0.49539613) q[2];
sx q[2];
rz(-2.317111) q[2];
sx q[2];
rz(-1.3935312) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71660605) q[1];
sx q[1];
rz(-2.0946436) q[1];
sx q[1];
rz(0.040158466) q[1];
rz(-pi) q[2];
rz(2.1776803) q[3];
sx q[3];
rz(-1.5480032) q[3];
sx q[3];
rz(-2.7085639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2942723) q[2];
sx q[2];
rz(-2.5966094) q[2];
sx q[2];
rz(1.6418183) q[2];
rz(-0.13449399) q[3];
sx q[3];
rz(-1.2765063) q[3];
sx q[3];
rz(-0.78181481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0791557) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(-2.675918) q[0];
rz(-1.8648719) q[1];
sx q[1];
rz(-2.1379505) q[1];
sx q[1];
rz(-1.0328971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55556923) q[0];
sx q[0];
rz(-0.15842552) q[0];
sx q[0];
rz(1.5044295) q[0];
rz(2.5215785) q[2];
sx q[2];
rz(-0.17159941) q[2];
sx q[2];
rz(2.624334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8414265) q[1];
sx q[1];
rz(-0.71117095) q[1];
sx q[1];
rz(0.39022846) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0276658) q[3];
sx q[3];
rz(-1.2333721) q[3];
sx q[3];
rz(1.3285445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5501962) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-2.2097394) q[2];
rz(0.094680928) q[3];
sx q[3];
rz(-0.90016142) q[3];
sx q[3];
rz(1.3100821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4513627) q[0];
sx q[0];
rz(-2.9424423) q[0];
sx q[0];
rz(-3.0354101) q[0];
rz(-0.55446082) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(1.6000481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743177) q[0];
sx q[0];
rz(-1.5948926) q[0];
sx q[0];
rz(0.39053183) q[0];
rz(-1.0511398) q[2];
sx q[2];
rz(-2.6259632) q[2];
sx q[2];
rz(1.3517429) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8733124) q[1];
sx q[1];
rz(-2.2619216) q[1];
sx q[1];
rz(-1.4119488) q[1];
rz(-pi) q[2];
rz(0.21884905) q[3];
sx q[3];
rz(-0.25222029) q[3];
sx q[3];
rz(2.9014933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23325486) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(-1.6439269) q[2];
rz(-2.6805367) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(-1.8259995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(-1.0787971) q[0];
rz(0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(-2.914391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086454) q[0];
sx q[0];
rz(-1.4221749) q[0];
sx q[0];
rz(0.11283837) q[0];
rz(-pi) q[1];
rz(1.237713) q[2];
sx q[2];
rz(-1.4374315) q[2];
sx q[2];
rz(1.9072717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.46887423) q[1];
sx q[1];
rz(-1.0212082) q[1];
sx q[1];
rz(2.7764715) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3183837) q[3];
sx q[3];
rz(-0.38312437) q[3];
sx q[3];
rz(-2.8069404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6195153) q[2];
sx q[2];
rz(-0.41831133) q[2];
sx q[2];
rz(0.72959161) q[2];
rz(1.6884165) q[3];
sx q[3];
rz(-1.6875234) q[3];
sx q[3];
rz(0.61354536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59298092) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(-2.7899637) q[0];
rz(2.8278606) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(-1.7650013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54395318) q[0];
sx q[0];
rz(-2.7386129) q[0];
sx q[0];
rz(0.12059327) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0259834) q[2];
sx q[2];
rz(-1.7358276) q[2];
sx q[2];
rz(0.045265667) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2089858) q[1];
sx q[1];
rz(-2.5521005) q[1];
sx q[1];
rz(-2.8578651) q[1];
rz(-pi) q[2];
x q[2];
rz(0.099011856) q[3];
sx q[3];
rz(-1.9388698) q[3];
sx q[3];
rz(1.2457445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54479638) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(1.6205622) q[2];
rz(-0.14791402) q[3];
sx q[3];
rz(-1.7971797) q[3];
sx q[3];
rz(1.773905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.3677597) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.2441147) q[0];
rz(-1.4338088) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(0.22020766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51802902) q[0];
sx q[0];
rz(-1.9401778) q[0];
sx q[0];
rz(2.2165038) q[0];
rz(-pi) q[1];
rz(-0.16248361) q[2];
sx q[2];
rz(-1.1097317) q[2];
sx q[2];
rz(-1.3853672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9088194) q[1];
sx q[1];
rz(-2.4202883) q[1];
sx q[1];
rz(-0.90296332) q[1];
rz(-pi) q[2];
x q[2];
rz(0.098252849) q[3];
sx q[3];
rz(-1.6643401) q[3];
sx q[3];
rz(-0.52631179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1844909) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(1.6415049) q[2];
rz(0.084065048) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(-0.070092289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641007) q[0];
sx q[0];
rz(-0.75834948) q[0];
sx q[0];
rz(0.41859928) q[0];
rz(-0.94340008) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(2.8387866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3269791) q[0];
sx q[0];
rz(-1.7075065) q[0];
sx q[0];
rz(1.4025406) q[0];
x q[1];
rz(2.9507016) q[2];
sx q[2];
rz(-1.396469) q[2];
sx q[2];
rz(-3.0534923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23754072) q[1];
sx q[1];
rz(-1.5273603) q[1];
sx q[1];
rz(0.082708184) q[1];
rz(3.0251011) q[3];
sx q[3];
rz(-1.9132275) q[3];
sx q[3];
rz(-0.064700944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4960949) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(-0.75054753) q[2];
rz(1.6802855) q[3];
sx q[3];
rz(-2.3716726) q[3];
sx q[3];
rz(2.6245978) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208329) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(1.2278521) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(-0.28245482) q[2];
sx q[2];
rz(-0.58517848) q[2];
sx q[2];
rz(-2.7369557) q[2];
rz(-1.8932992) q[3];
sx q[3];
rz(-1.7453946) q[3];
sx q[3];
rz(3.0562492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
