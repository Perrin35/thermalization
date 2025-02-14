OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(2.0164665) q[0];
rz(-5.2611051) q[1];
sx q[1];
rz(2.4740969) q[1];
sx q[1];
rz(11.752887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73972244) q[0];
sx q[0];
rz(-2.7474294) q[0];
sx q[0];
rz(1.1959082) q[0];
rz(-pi) q[1];
rz(0.44806077) q[2];
sx q[2];
rz(-1.2431113) q[2];
sx q[2];
rz(2.0304012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91650016) q[1];
sx q[1];
rz(-1.4218907) q[1];
sx q[1];
rz(2.6970106) q[1];
rz(-3.0510819) q[3];
sx q[3];
rz(-1.4078684) q[3];
sx q[3];
rz(-0.75608692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(-0.92864621) q[2];
rz(1.5422025) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(1.9699875) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20392513) q[0];
sx q[0];
rz(-1.3805905) q[0];
sx q[0];
rz(-0.12705886) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(0.7712706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5708904) q[0];
sx q[0];
rz(-1.2178019) q[0];
sx q[0];
rz(0.81309005) q[0];
rz(-pi) q[1];
rz(3.0808582) q[2];
sx q[2];
rz(-1.2337451) q[2];
sx q[2];
rz(-3.0063546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82133085) q[1];
sx q[1];
rz(-1.2148569) q[1];
sx q[1];
rz(-1.2004281) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3985474) q[3];
sx q[3];
rz(-2.4186169) q[3];
sx q[3];
rz(1.522775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1473006) q[2];
sx q[2];
rz(-1.9376126) q[2];
sx q[2];
rz(3.1331983) q[2];
rz(2.4781135) q[3];
sx q[3];
rz(-1.9050262) q[3];
sx q[3];
rz(2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0405149) q[0];
sx q[0];
rz(-0.82656693) q[0];
sx q[0];
rz(-2.7000632) q[0];
rz(2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(0.13557869) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39865935) q[0];
sx q[0];
rz(-0.017983111) q[0];
sx q[0];
rz(0.44264408) q[0];
rz(-pi) q[1];
rz(-0.71696059) q[2];
sx q[2];
rz(-1.8646984) q[2];
sx q[2];
rz(0.64168054) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8959055) q[1];
sx q[1];
rz(-2.062633) q[1];
sx q[1];
rz(-2.7373284) q[1];
rz(-pi) q[2];
rz(-2.3051585) q[3];
sx q[3];
rz(-1.7178917) q[3];
sx q[3];
rz(-0.35415748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93379891) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(0.03820339) q[2];
rz(-0.52538747) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(-2.4562522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66048375) q[0];
sx q[0];
rz(-2.3479192) q[0];
sx q[0];
rz(2.5307181) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.7673312) q[1];
sx q[1];
rz(0.23922051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4668381) q[0];
sx q[0];
rz(-1.309899) q[0];
sx q[0];
rz(-0.99715085) q[0];
rz(-pi) q[1];
rz(-0.47232136) q[2];
sx q[2];
rz(-1.2178253) q[2];
sx q[2];
rz(-1.2833088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4799616) q[1];
sx q[1];
rz(-1.1203655) q[1];
sx q[1];
rz(1.2110787) q[1];
x q[2];
rz(0.38762761) q[3];
sx q[3];
rz(-2.648733) q[3];
sx q[3];
rz(2.8475177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7188344) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(-3.139843) q[2];
rz(0.29378978) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(-2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2413498) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(-2.2985261) q[0];
rz(0.33755606) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(1.6036124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6809236) q[0];
sx q[0];
rz(-2.3516555) q[0];
sx q[0];
rz(-1.8415175) q[0];
x q[1];
rz(-1.6241239) q[2];
sx q[2];
rz(-1.7449208) q[2];
sx q[2];
rz(-1.8873896) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.027315779) q[1];
sx q[1];
rz(-1.3754579) q[1];
sx q[1];
rz(-2.8802425) q[1];
x q[2];
rz(-0.20974737) q[3];
sx q[3];
rz(-1.4919859) q[3];
sx q[3];
rz(-0.94890734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43188492) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(-1.0464) q[2];
rz(-0.31878582) q[3];
sx q[3];
rz(-2.4304978) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(-1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(-1.3290149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8662939) q[0];
sx q[0];
rz(-1.4614551) q[0];
sx q[0];
rz(-1.7924395) q[0];
x q[1];
rz(1.3565265) q[2];
sx q[2];
rz(-1.2870711) q[2];
sx q[2];
rz(-1.8510712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44101366) q[1];
sx q[1];
rz(-2.3668336) q[1];
sx q[1];
rz(-2.0225594) q[1];
rz(-pi) q[2];
rz(2.5446114) q[3];
sx q[3];
rz(-2.5425445) q[3];
sx q[3];
rz(2.0839276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38137388) q[2];
sx q[2];
rz(-1.2083283) q[2];
sx q[2];
rz(1.5488497) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(2.2729592) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(-2.1997531) q[0];
rz(-0.16009227) q[1];
sx q[1];
rz(-1.6611049) q[1];
sx q[1];
rz(-0.19217415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3775786) q[0];
sx q[0];
rz(-1.6068216) q[0];
sx q[0];
rz(-1.7948304) q[0];
rz(-pi) q[1];
rz(-1.6404387) q[2];
sx q[2];
rz(-2.6187839) q[2];
sx q[2];
rz(1.664134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77171626) q[1];
sx q[1];
rz(-1.2315005) q[1];
sx q[1];
rz(-2.3849065) q[1];
rz(-pi) q[2];
rz(-0.022054733) q[3];
sx q[3];
rz(-2.0682749) q[3];
sx q[3];
rz(-2.9201404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3840702) q[2];
sx q[2];
rz(-1.2382058) q[2];
sx q[2];
rz(-2.7080217) q[2];
rz(1.5844257) q[3];
sx q[3];
rz(-1.4945364) q[3];
sx q[3];
rz(2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.7973416) q[0];
rz(-1.2035707) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(-1.3628091) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8289611) q[0];
sx q[0];
rz(-1.5927218) q[0];
sx q[0];
rz(-1.5410822) q[0];
rz(-pi) q[1];
rz(1.7596471) q[2];
sx q[2];
rz(-2.3879693) q[2];
sx q[2];
rz(-0.058503956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2143933) q[1];
sx q[1];
rz(-1.7991369) q[1];
sx q[1];
rz(1.4424999) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1892494) q[3];
sx q[3];
rz(-0.7233215) q[3];
sx q[3];
rz(1.9017232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5726996) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(-0.9838689) q[2];
rz(0.24108663) q[3];
sx q[3];
rz(-2.2086996) q[3];
sx q[3];
rz(-2.4510395) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6702061) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(-2.8669226) q[0];
rz(-1.7999016) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(-2.0955657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90993308) q[0];
sx q[0];
rz(-0.69485661) q[0];
sx q[0];
rz(2.6207663) q[0];
rz(-pi) q[1];
rz(-3.0606424) q[2];
sx q[2];
rz(-1.7678723) q[2];
sx q[2];
rz(-2.112889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.14635) q[1];
sx q[1];
rz(-2.1533794) q[1];
sx q[1];
rz(-1.4928774) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6469429) q[3];
sx q[3];
rz(-2.6373495) q[3];
sx q[3];
rz(-1.7011736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(0.39946237) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-0.96650201) q[3];
sx q[3];
rz(0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852916) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(0.5823108) q[0];
rz(0.18320006) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(-2.3351672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9574653) q[0];
sx q[0];
rz(-0.13617198) q[0];
sx q[0];
rz(-0.86958142) q[0];
rz(-pi) q[1];
rz(-0.66566531) q[2];
sx q[2];
rz(-2.7241926) q[2];
sx q[2];
rz(1.1102499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90102531) q[1];
sx q[1];
rz(-1.7738924) q[1];
sx q[1];
rz(-0.74764772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7943198) q[3];
sx q[3];
rz(-1.6753622) q[3];
sx q[3];
rz(2.4109651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9091984) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(0.8052899) q[2];
rz(0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(-0.73827353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.88937) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(-1.5421142) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-3.0270544) q[2];
sx q[2];
rz(-2.2907612) q[2];
sx q[2];
rz(-1.9775122) q[2];
rz(-0.44381683) q[3];
sx q[3];
rz(-2.4469821) q[3];
sx q[3];
rz(-1.2367005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
