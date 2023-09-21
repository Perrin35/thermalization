OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(-1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9278487) q[0];
sx q[0];
rz(-2.1670177) q[0];
sx q[0];
rz(2.5736789) q[0];
x q[1];
rz(2.6745559) q[2];
sx q[2];
rz(-0.39614284) q[2];
sx q[2];
rz(0.064183891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7429744) q[1];
sx q[1];
rz(-0.75098872) q[1];
sx q[1];
rz(-0.91575925) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7513566) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(-0.80871049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(1.4398549) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97025362) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.617584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47027662) q[0];
sx q[0];
rz(-0.19965262) q[0];
sx q[0];
rz(-1.5623564) q[0];
rz(-pi) q[1];
rz(1.6127869) q[2];
sx q[2];
rz(-0.4547555) q[2];
sx q[2];
rz(3.0595879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9873725) q[1];
sx q[1];
rz(-2.1402573) q[1];
sx q[1];
rz(3.0595368) q[1];
rz(-pi) q[2];
rz(1.1094692) q[3];
sx q[3];
rz(-0.20878775) q[3];
sx q[3];
rz(1.3267645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38561884) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.1478109) q[2];
rz(-1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(0.23708788) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-0.25207239) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35226563) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(0.69791039) q[0];
rz(-2.4263072) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(-0.31179024) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4430868) q[1];
sx q[1];
rz(-1.3830739) q[1];
sx q[1];
rz(2.1952573) q[1];
rz(0.73022233) q[3];
sx q[3];
rz(-2.632004) q[3];
sx q[3];
rz(0.93917055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(1.0954558) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(-2.0544255) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1336466) q[0];
sx q[0];
rz(-2.5292853) q[0];
sx q[0];
rz(0.72511073) q[0];
rz(-pi) q[1];
rz(-0.44666501) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(0.54892533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69869631) q[1];
sx q[1];
rz(-2.2203608) q[1];
sx q[1];
rz(-0.15028468) q[1];
x q[2];
rz(2.1726923) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(-0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(0.46009955) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(-1.3943577) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(-2.8869693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4655612) q[0];
sx q[0];
rz(-1.4916294) q[0];
sx q[0];
rz(-0.013750793) q[0];
rz(-pi) q[1];
rz(0.014572797) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(2.2968959) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4659781) q[1];
sx q[1];
rz(-1.8120159) q[1];
sx q[1];
rz(2.4720008) q[1];
x q[2];
rz(1.0765431) q[3];
sx q[3];
rz(-0.77270618) q[3];
sx q[3];
rz(-2.5172174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(-1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(2.9350231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20139192) q[0];
sx q[0];
rz(-2.3513146) q[0];
sx q[0];
rz(-1.9076365) q[0];
rz(-pi) q[1];
rz(1.0429522) q[2];
sx q[2];
rz(-1.2252508) q[2];
sx q[2];
rz(-2.8258459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5094604) q[1];
sx q[1];
rz(-1.3941947) q[1];
sx q[1];
rz(-2.280974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9566831) q[3];
sx q[3];
rz(-0.69291249) q[3];
sx q[3];
rz(-0.087156765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(-2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92641002) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(-1.8136224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0921811) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(1.1294424) q[0];
rz(-pi) q[1];
rz(2.1778657) q[2];
sx q[2];
rz(-2.4237195) q[2];
sx q[2];
rz(-1.3340064) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.59203397) q[1];
sx q[1];
rz(-2.4394819) q[1];
sx q[1];
rz(-1.3105427) q[1];
x q[2];
rz(1.1442723) q[3];
sx q[3];
rz(-2.0445619) q[3];
sx q[3];
rz(-1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61838377) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-2.1906733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5247105) q[0];
sx q[0];
rz(-1.2843772) q[0];
sx q[0];
rz(2.8911203) q[0];
rz(-pi) q[1];
rz(1.3046706) q[2];
sx q[2];
rz(-1.6309435) q[2];
sx q[2];
rz(-1.4820478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2086522) q[1];
sx q[1];
rz(-1.7793852) q[1];
sx q[1];
rz(1.5044466) q[1];
rz(2.7637134) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(0.32285238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(-0.65972796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10044554) q[0];
sx q[0];
rz(-2.1491096) q[0];
sx q[0];
rz(0.64097494) q[0];
x q[1];
rz(-1.3525891) q[2];
sx q[2];
rz(-2.3404684) q[2];
sx q[2];
rz(-1.4584691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5539726) q[1];
sx q[1];
rz(-2.4269322) q[1];
sx q[1];
rz(-2.0228533) q[1];
rz(-1.0321484) q[3];
sx q[3];
rz(-1.5218471) q[3];
sx q[3];
rz(-0.34070542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62347162) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(2.1155604) q[2];
rz(2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.3056668) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4982088) q[0];
sx q[0];
rz(-1.9783101) q[0];
sx q[0];
rz(1.2163175) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56634855) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(-1.7042421) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.021691572) q[1];
sx q[1];
rz(-2.7530115) q[1];
sx q[1];
rz(0.67967023) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3922937) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(-1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(-1.1516085) q[2];
rz(0.51268762) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.5079386) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(-0.097461854) q[2];
sx q[2];
rz(-0.41257358) q[2];
sx q[2];
rz(-1.67795) q[2];
rz(-0.012398331) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];