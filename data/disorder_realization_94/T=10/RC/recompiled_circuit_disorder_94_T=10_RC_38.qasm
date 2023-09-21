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
rz(1.6341524) q[0];
rz(-1.545067) q[1];
sx q[1];
rz(-2.5453321) q[1];
sx q[1];
rz(2.526386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98696729) q[0];
sx q[0];
rz(-2.0320503) q[0];
sx q[0];
rz(0.89303645) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7841714) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-2.0703966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93278904) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(-2.6245481) q[1];
rz(-2.8350713) q[3];
sx q[3];
rz(-1.74311) q[3];
sx q[3];
rz(-2.4337208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(-2.8033076) q[2];
rz(1.4398549) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(0.030348226) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.5240086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.679927) q[0];
sx q[0];
rz(-1.7704417) q[0];
sx q[0];
rz(-0.0017077831) q[0];
x q[1];
rz(1.6127869) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(-3.0595879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46088947) q[1];
sx q[1];
rz(-1.5017121) q[1];
sx q[1];
rz(0.99980385) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1094692) q[3];
sx q[3];
rz(-2.9328049) q[3];
sx q[3];
rz(-1.3267645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7559738) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(2.8258064) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(2.8895203) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35226563) q[0];
sx q[0];
rz(-1.1061215) q[0];
sx q[0];
rz(-2.4436823) q[0];
rz(-pi) q[1];
rz(-0.75906934) q[2];
sx q[2];
rz(-2.1096863) q[2];
sx q[2];
rz(-2.3784504) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0059800681) q[1];
sx q[1];
rz(-0.95893919) q[1];
sx q[1];
rz(2.9115885) q[1];
x q[2];
rz(-1.2139981) q[3];
sx q[3];
rz(-1.9427512) q[3];
sx q[3];
rz(-1.7372204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(-3.1075409) q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(1.5267641) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0079460572) q[0];
sx q[0];
rz(-2.5292853) q[0];
sx q[0];
rz(0.72511073) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6949276) q[2];
sx q[2];
rz(-1.6840877) q[2];
sx q[2];
rz(-0.54892533) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.69869631) q[1];
sx q[1];
rz(-2.2203608) q[1];
sx q[1];
rz(-0.15028468) q[1];
x q[2];
rz(0.21541689) q[3];
sx q[3];
rz(-1.2691174) q[3];
sx q[3];
rz(0.33723436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.8206772) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(0.25462338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29339644) q[0];
sx q[0];
rz(-3.0612429) q[0];
sx q[0];
rz(1.7424165) q[0];
rz(3.1270199) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(-2.2968959) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9533206) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(-0.377368) q[1];
rz(0.43318627) q[3];
sx q[3];
rz(-2.232589) q[3];
sx q[3];
rz(-1.2695241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(-1.3767892) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(0.20656955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5305938) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(2.3321652) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1913387) q[2];
sx q[2];
rz(-2.5197919) q[2];
sx q[2];
rz(-1.3602464) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2628852) q[1];
sx q[1];
rz(-0.72808121) q[1];
sx q[1];
rz(1.8379777) q[1];
rz(-pi) q[2];
rz(-2.2264678) q[3];
sx q[3];
rz(-1.8135999) q[3];
sx q[3];
rz(1.1806928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(-0.28277961) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(0.16684428) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151826) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(2.5577257) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(1.8136224) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5751942) q[0];
sx q[0];
rz(-2.0090721) q[0];
sx q[0];
rz(3.0153494) q[0];
rz(-0.96372693) q[2];
sx q[2];
rz(-2.4237195) q[2];
sx q[2];
rz(1.8075862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77818645) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(2.2560675) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9973203) q[3];
sx q[3];
rz(-2.0445619) q[3];
sx q[3];
rz(1.1748479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(-1.7810812) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(-2.9597136) q[0];
rz(-0.47422844) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(0.95091933) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61688214) q[0];
sx q[0];
rz(-1.2843772) q[0];
sx q[0];
rz(-2.8911203) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0792564) q[2];
sx q[2];
rz(-1.3051635) q[2];
sx q[2];
rz(-3.0692284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2434395) q[1];
sx q[1];
rz(-0.21874084) q[1];
sx q[1];
rz(-2.8380413) q[1];
x q[2];
rz(-0.37787921) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(-2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(3.0656832) q[2];
rz(-0.54801303) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(-2.7245522) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(-0.65972796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83803672) q[0];
sx q[0];
rz(-0.8350026) q[0];
sx q[0];
rz(-2.3124218) q[0];
rz(-1.3525891) q[2];
sx q[2];
rz(-0.80112427) q[2];
sx q[2];
rz(-1.6831236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5539726) q[1];
sx q[1];
rz(-2.4269322) q[1];
sx q[1];
rz(-2.0228533) q[1];
x q[2];
rz(-2.1094443) q[3];
sx q[3];
rz(-1.6197455) q[3];
sx q[3];
rz(2.8008872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62347162) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(-2.1155604) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.3056668) q[0];
rz(-1.9650412) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730597) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(2.7102094) q[0];
rz(-pi) q[1];
rz(-2.5752441) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(-1.4373506) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.73951057) q[1];
sx q[1];
rz(-1.869919) q[1];
sx q[1];
rz(1.3189391) q[1];
rz(-pi) q[2];
rz(-1.8412748) q[3];
sx q[3];
rz(-0.83993739) q[3];
sx q[3];
rz(-0.65281103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(-1.9899842) q[2];
rz(-0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.37968996) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(1.6336541) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(0.097461854) q[2];
sx q[2];
rz(-2.7290191) q[2];
sx q[2];
rz(1.4636427) q[2];
rz(0.012398331) q[3];
sx q[3];
rz(-2.5231902) q[3];
sx q[3];
rz(-1.4004422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];