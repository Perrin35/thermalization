OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3988848) q[0];
sx q[0];
rz(-2.3595915) q[0];
sx q[0];
rz(-1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35933094) q[0];
sx q[0];
rz(-1.2027272) q[0];
sx q[0];
rz(-2.9666535) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6127365) q[2];
sx q[2];
rz(-2.0176) q[2];
sx q[2];
rz(-0.97181335) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4692987) q[1];
sx q[1];
rz(-1.0743595) q[1];
sx q[1];
rz(-1.4908355) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(-2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(2.326139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0464697) q[0];
sx q[0];
rz(-2.2182584) q[0];
sx q[0];
rz(-1.6460653) q[0];
rz(2.4120861) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(-1.6934998) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.193589) q[1];
sx q[1];
rz(-1.6063599) q[1];
sx q[1];
rz(2.8405632) q[1];
rz(-2.4957982) q[3];
sx q[3];
rz(-1.3919953) q[3];
sx q[3];
rz(1.4327232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(-0.30953428) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(0.99951807) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8130428) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(0.15655984) q[0];
rz(-pi) q[1];
rz(0.83061647) q[2];
sx q[2];
rz(-0.81095552) q[2];
sx q[2];
rz(-2.9142771) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(-1.2815777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7129094) q[3];
sx q[3];
rz(-0.51321533) q[3];
sx q[3];
rz(-1.7538479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6040566) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(2.3245658) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-0.27483637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4764458) q[0];
sx q[0];
rz(-3.0229212) q[0];
sx q[0];
rz(0.84400405) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1399973) q[2];
sx q[2];
rz(-2.0312107) q[2];
sx q[2];
rz(2.7574725) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53510016) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(-2.5938354) q[1];
x q[2];
rz(0.87423012) q[3];
sx q[3];
rz(-1.4733553) q[3];
sx q[3];
rz(1.3144573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.794902) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.5083195) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(-0.99779469) q[0];
rz(-2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.516974) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.957513) q[0];
sx q[0];
rz(-2.3947869) q[0];
sx q[0];
rz(-1.0190796) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.595572) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9285674) q[1];
sx q[1];
rz(-1.2270317) q[1];
sx q[1];
rz(-2.7108971) q[1];
rz(-pi) q[2];
rz(0.77002854) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.5312622) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76535392) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(-2.8009159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6402123) q[0];
sx q[0];
rz(-0.073177241) q[0];
sx q[0];
rz(1.120938) q[0];
x q[1];
rz(2.8315115) q[2];
sx q[2];
rz(-0.78044621) q[2];
sx q[2];
rz(0.23852894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6319879) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(1.8297086) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0814704) q[3];
sx q[3];
rz(-1.9483856) q[3];
sx q[3];
rz(-1.1946354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(-2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(-2.7959438) q[0];
rz(0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(2.6766434) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3195254) q[0];
sx q[0];
rz(-1.3789346) q[0];
sx q[0];
rz(-2.0454387) q[0];
rz(-pi) q[1];
rz(0.49352383) q[2];
sx q[2];
rz(-1.0527305) q[2];
sx q[2];
rz(1.7871737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3685776) q[1];
sx q[1];
rz(-1.3524719) q[1];
sx q[1];
rz(-1.6792084) q[1];
rz(-pi) q[2];
rz(-1.5246478) q[3];
sx q[3];
rz(-1.5966291) q[3];
sx q[3];
rz(-1.8125364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(-0.31759343) q[2];
rz(0.57146227) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(-0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(3.0632339) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6440455) q[0];
sx q[0];
rz(-0.80701485) q[0];
sx q[0];
rz(-3.0456196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23711726) q[2];
sx q[2];
rz(-1.3318921) q[2];
sx q[2];
rz(2.4799926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6137177) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(-1.2916958) q[1];
rz(-pi) q[2];
rz(-1.7844291) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(2.8182639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(2.890214) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-0.87402469) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286205) q[0];
sx q[0];
rz(-0.7681094) q[0];
sx q[0];
rz(-3.1290929) q[0];
x q[1];
rz(1.9782412) q[2];
sx q[2];
rz(-1.2375087) q[2];
sx q[2];
rz(-2.8826706) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8928788) q[1];
sx q[1];
rz(-1.4993748) q[1];
sx q[1];
rz(-1.0700657) q[1];
rz(-pi) q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.2864283) q[3];
sx q[3];
rz(-1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-0.61974636) q[2];
rz(-1.9571346) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26319474) q[0];
sx q[0];
rz(-2.7976755) q[0];
sx q[0];
rz(-2.2180936) q[0];
rz(-pi) q[1];
rz(1.6640501) q[2];
sx q[2];
rz(-1.4443195) q[2];
sx q[2];
rz(1.7699514) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.57943343) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(1.2594373) q[1];
rz(-pi) q[2];
rz(2.6188649) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(0.75795603) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(2.010871) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(-0.27579565) q[3];
sx q[3];
rz(-1.7134568) q[3];
sx q[3];
rz(-1.4207763) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];