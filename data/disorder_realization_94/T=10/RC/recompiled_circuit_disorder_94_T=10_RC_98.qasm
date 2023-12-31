OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(-1.5885408) q[0];
sx q[0];
rz(1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1546254) q[0];
sx q[0];
rz(-2.0320503) q[0];
sx q[0];
rz(0.89303645) q[0];
x q[1];
rz(1.7569321) q[2];
sx q[2];
rz(-1.9225444) q[2];
sx q[2];
rz(2.577201) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2088036) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(2.6245481) q[1];
rz(-pi) q[2];
rz(-2.6184222) q[3];
sx q[3];
rz(-0.35029951) q[3];
sx q[3];
rz(-1.7821799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7011828) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(-2.8033076) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.679927) q[0];
sx q[0];
rz(-1.7704417) q[0];
sx q[0];
rz(-3.1398849) q[0];
rz(-pi) q[1];
x q[1];
rz(0.020521684) q[2];
sx q[2];
rz(-1.1164718) q[2];
sx q[2];
rz(0.12873912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9873725) q[1];
sx q[1];
rz(-2.1402573) q[1];
sx q[1];
rz(3.0595368) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3832983) q[3];
sx q[3];
rz(-1.6631931) q[3];
sx q[3];
rz(2.932991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(0.23708788) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-2.8258064) q[0];
rz(0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-2.8895203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35226563) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(-0.69791039) q[0];
rz(2.2601068) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(-1.2607247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0059800681) q[1];
sx q[1];
rz(-2.1826535) q[1];
sx q[1];
rz(2.9115885) q[1];
rz(-1.9275946) q[3];
sx q[3];
rz(-1.9427512) q[3];
sx q[3];
rz(1.7372204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0198274) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(-3.1075409) q[2];
rz(-3.1241336) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3084761) q[0];
sx q[0];
rz(-2.0154698) q[0];
sx q[0];
rz(-2.0067257) q[0];
x q[1];
rz(0.44666501) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(2.5926673) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45338079) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(1.7654256) q[1];
rz(-pi) q[2];
rz(-0.96890038) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(2.8440648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(-2.8989255) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3209155) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(0.25462338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.03527) q[0];
sx q[0];
rz(-1.584504) q[0];
sx q[0];
rz(-1.4916219) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1483634) q[2];
sx q[2];
rz(-1.5585871) q[2];
sx q[2];
rz(-0.73405594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9533206) q[1];
sx q[1];
rz(-2.4362262) q[1];
sx q[1];
rz(2.7642247) q[1];
rz(2.7084064) q[3];
sx q[3];
rz(-0.90900366) q[3];
sx q[3];
rz(1.8720686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(-1.3767892) q[2];
rz(1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(0.20656955) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6109989) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(0.80942746) q[0];
rz(-0.95025392) q[2];
sx q[2];
rz(-2.5197919) q[2];
sx q[2];
rz(-1.7813462) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.088591136) q[1];
sx q[1];
rz(-0.87391657) q[1];
sx q[1];
rz(-0.23115302) q[1];
x q[2];
rz(-2.8387186) q[3];
sx q[3];
rz(-2.2040963) q[3];
sx q[3];
rz(-0.57297046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(2.9747484) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92641002) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(0.38152951) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.3279703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8570003) q[0];
sx q[0];
rz(-0.45495957) q[0];
sx q[0];
rz(1.3083463) q[0];
x q[1];
rz(2.1778657) q[2];
sx q[2];
rz(-0.71787314) q[2];
sx q[2];
rz(1.3340064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59203397) q[1];
sx q[1];
rz(-2.4394819) q[1];
sx q[1];
rz(1.83105) q[1];
rz(1.9973203) q[3];
sx q[3];
rz(-1.0970308) q[3];
sx q[3];
rz(-1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(-0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-0.95091933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.788818) q[0];
sx q[0];
rz(-2.7633861) q[0];
sx q[0];
rz(2.2703855) q[0];
x q[1];
rz(1.3457001) q[2];
sx q[2];
rz(-0.27268073) q[2];
sx q[2];
rz(2.8358012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89815318) q[1];
sx q[1];
rz(-2.9228518) q[1];
sx q[1];
rz(2.8380413) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7637134) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(-2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(-0.075909464) q[2];
rz(-0.54801303) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(-2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(-2.7245522) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-0.65972796) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10044554) q[0];
sx q[0];
rz(-0.99248306) q[0];
sx q[0];
rz(-2.5006177) q[0];
x q[1];
rz(-1.3525891) q[2];
sx q[2];
rz(-0.80112427) q[2];
sx q[2];
rz(1.4584691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1253423) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(-2.7793105) q[1];
rz(0.057007313) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(1.9407335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.518121) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(-1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-1.0326577) q[3];
sx q[3];
rz(0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24348564) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(1.0356888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10782345) q[0];
sx q[0];
rz(-2.6080837) q[0];
sx q[0];
rz(0.67722042) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91949384) q[2];
sx q[2];
rz(-0.80923015) q[2];
sx q[2];
rz(2.2724255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1199011) q[1];
sx q[1];
rz(-0.38858116) q[1];
sx q[1];
rz(-2.4619224) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74929897) q[3];
sx q[3];
rz(-1.7710847) q[3];
sx q[3];
rz(-2.0405958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62853652) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(-0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(1.5079386) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-1.528231) q[2];
sx q[2];
rz(-1.1602989) q[2];
sx q[2];
rz(1.5699671) q[2];
rz(3.1291943) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
