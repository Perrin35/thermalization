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
rz(1.4752969) q[0];
sx q[0];
rz(1.8721606) q[0];
sx q[0];
rz(8.7526487) q[0];
rz(-2.693306) q[1];
sx q[1];
rz(-1.6192133) q[1];
sx q[1];
rz(-0.7575922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8586977) q[0];
sx q[0];
rz(-3.0083249) q[0];
sx q[0];
rz(-0.96576502) q[0];
rz(-pi) q[1];
rz(2.640921) q[2];
sx q[2];
rz(-1.8953875) q[2];
sx q[2];
rz(-2.2364834) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1266844) q[1];
sx q[1];
rz(-1.9796625) q[1];
sx q[1];
rz(1.3247299) q[1];
rz(-pi) q[2];
rz(2.1610307) q[3];
sx q[3];
rz(-2.1781028) q[3];
sx q[3];
rz(0.052562873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0483094) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(-1.3422356) q[2];
rz(-1.3679158) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9369478) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(-0.94888765) q[0];
rz(0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(0.92996517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5407911) q[0];
sx q[0];
rz(-1.3581498) q[0];
sx q[0];
rz(-1.6746816) q[0];
x q[1];
rz(-0.013596046) q[2];
sx q[2];
rz(-2.6363723) q[2];
sx q[2];
rz(1.1824974) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9252214) q[1];
sx q[1];
rz(-0.19268806) q[1];
sx q[1];
rz(2.5641901) q[1];
rz(-pi) q[2];
rz(2.5245776) q[3];
sx q[3];
rz(-2.2301144) q[3];
sx q[3];
rz(-0.78711817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0019504) q[2];
sx q[2];
rz(-1.6763687) q[2];
sx q[2];
rz(0.742221) q[2];
rz(-0.46418515) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(-1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4699698) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(1.4416913) q[0];
rz(2.9761159) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(0.86732078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.353426) q[0];
sx q[0];
rz(-1.5713619) q[0];
sx q[0];
rz(3.1415523) q[0];
rz(-2.0365085) q[2];
sx q[2];
rz(-2.5869479) q[2];
sx q[2];
rz(1.1797734) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8738385) q[1];
sx q[1];
rz(-1.2902564) q[1];
sx q[1];
rz(-0.33919097) q[1];
rz(-pi) q[2];
rz(-2.1380566) q[3];
sx q[3];
rz(-0.84430443) q[3];
sx q[3];
rz(0.15641016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1778339) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(2.3826694) q[2];
rz(-0.80398503) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8300962) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(2.6575644) q[0];
rz(-1.6054035) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(-1.4253634) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4771627) q[0];
sx q[0];
rz(-1.2331404) q[0];
sx q[0];
rz(0.54414026) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3431802) q[2];
sx q[2];
rz(-1.8755442) q[2];
sx q[2];
rz(-2.9272542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64827907) q[1];
sx q[1];
rz(-1.7910622) q[1];
sx q[1];
rz(0.92695285) q[1];
x q[2];
rz(-0.68257733) q[3];
sx q[3];
rz(-2.5831476) q[3];
sx q[3];
rz(-2.6598833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75590762) q[2];
sx q[2];
rz(-1.0681095) q[2];
sx q[2];
rz(2.8483086) q[2];
rz(-0.048132345) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(0.3046681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3047979) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(-2.0378713) q[0];
rz(-0.91148218) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(-0.10890659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86065642) q[0];
sx q[0];
rz(-2.0514279) q[0];
sx q[0];
rz(-2.3970766) q[0];
rz(-pi) q[1];
rz(-2.6066066) q[2];
sx q[2];
rz(-1.7270544) q[2];
sx q[2];
rz(1.6339169) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6088691) q[1];
sx q[1];
rz(-0.91174284) q[1];
sx q[1];
rz(-2.9347035) q[1];
rz(-2.6960228) q[3];
sx q[3];
rz(-2.516149) q[3];
sx q[3];
rz(-1.3423962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(2.8835127) q[2];
rz(-0.38170013) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(-2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447164) q[0];
sx q[0];
rz(-2.1602186) q[0];
sx q[0];
rz(-1.9816403) q[0];
rz(0.57811919) q[1];
sx q[1];
rz(-1.4719529) q[1];
sx q[1];
rz(2.8181308) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1581381) q[0];
sx q[0];
rz(-1.1295415) q[0];
sx q[0];
rz(-2.6652314) q[0];
rz(1.3176962) q[2];
sx q[2];
rz(-2.0539114) q[2];
sx q[2];
rz(0.3699257) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9449294) q[1];
sx q[1];
rz(-2.5179836) q[1];
sx q[1];
rz(2.7476279) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59185054) q[3];
sx q[3];
rz(-0.53880063) q[3];
sx q[3];
rz(2.0350698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20299992) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(0.8482376) q[2];
rz(-1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11442014) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(2.2681336) q[0];
rz(-1.9050441) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(-3.0580318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3665584) q[0];
sx q[0];
rz(-1.2511484) q[0];
sx q[0];
rz(-1.4421706) q[0];
x q[1];
rz(-0.23891831) q[2];
sx q[2];
rz(-2.0878125) q[2];
sx q[2];
rz(0.97415249) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.04490964) q[1];
sx q[1];
rz(-2.2538141) q[1];
sx q[1];
rz(-0.51814305) q[1];
x q[2];
rz(2.450299) q[3];
sx q[3];
rz(-0.45540998) q[3];
sx q[3];
rz(0.56870715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7053335) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(-1.7748888) q[2];
rz(2.8691835) q[3];
sx q[3];
rz(-1.0206157) q[3];
sx q[3];
rz(0.75850707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347539) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(-0.93836623) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-0.96790853) q[1];
sx q[1];
rz(0.82069194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21194776) q[0];
sx q[0];
rz(-1.7045472) q[0];
sx q[0];
rz(1.1991713) q[0];
rz(-pi) q[1];
rz(-1.9442476) q[2];
sx q[2];
rz(-1.4146342) q[2];
sx q[2];
rz(1.3661623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31466111) q[1];
sx q[1];
rz(-0.95035997) q[1];
sx q[1];
rz(-3.0112096) q[1];
rz(-pi) q[2];
rz(-1.971463) q[3];
sx q[3];
rz(-1.7224215) q[3];
sx q[3];
rz(-1.135863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20251033) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(1.067777) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(-3.06156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.985567) q[0];
sx q[0];
rz(-1.9859059) q[0];
sx q[0];
rz(-0.40618968) q[0];
rz(1.73229) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(2.0565775) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4403518) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(-1.3568272) q[0];
rz(-pi) q[1];
rz(-3.135097) q[2];
sx q[2];
rz(-0.61031872) q[2];
sx q[2];
rz(1.6808866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.20445261) q[1];
sx q[1];
rz(-1.8984183) q[1];
sx q[1];
rz(-0.83468584) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3906443) q[3];
sx q[3];
rz(-2.4118858) q[3];
sx q[3];
rz(0.88353222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5857508) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(2.2108868) q[2];
rz(1.7454923) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4569106) q[0];
sx q[0];
rz(-2.1645808) q[0];
sx q[0];
rz(-1.3071625) q[0];
rz(2.7556509) q[1];
sx q[1];
rz(-1.7470876) q[1];
sx q[1];
rz(1.0345667) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912303) q[0];
sx q[0];
rz(-1.467839) q[0];
sx q[0];
rz(1.7454171) q[0];
x q[1];
rz(-1.7770281) q[2];
sx q[2];
rz(-1.1185371) q[2];
sx q[2];
rz(0.79712501) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1423827) q[1];
sx q[1];
rz(-2.2015328) q[1];
sx q[1];
rz(2.5725065) q[1];
rz(-0.12328445) q[3];
sx q[3];
rz(-2.3269231) q[3];
sx q[3];
rz(2.2565763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4474386) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(0.36699692) q[2];
rz(2.0224109) q[3];
sx q[3];
rz(-0.39803353) q[3];
sx q[3];
rz(1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.781453) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(-0.057859261) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(0.17292508) q[2];
sx q[2];
rz(-1.6676219) q[2];
sx q[2];
rz(0.22990083) q[2];
rz(-0.41856159) q[3];
sx q[3];
rz(-0.7444612) q[3];
sx q[3];
rz(0.85519467) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
