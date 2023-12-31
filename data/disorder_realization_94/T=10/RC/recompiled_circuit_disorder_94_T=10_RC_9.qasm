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
rz(1.5965257) q[1];
sx q[1];
rz(2.5453321) q[1];
sx q[1];
rz(8.8095713) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9278487) q[0];
sx q[0];
rz(-2.1670177) q[0];
sx q[0];
rz(0.56791373) q[0];
rz(-pi) q[1];
rz(-0.35742128) q[2];
sx q[2];
rz(-1.3961785) q[2];
sx q[2];
rz(-1.071196) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8021009) q[1];
sx q[1];
rz(-1.1420982) q[1];
sx q[1];
rz(2.2080253) q[1];
x q[2];
rz(2.6184222) q[3];
sx q[3];
rz(-2.7912931) q[3];
sx q[3];
rz(-1.7821799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7011828) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(2.8033076) q[2];
rz(-1.4398549) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(-1.617584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.679927) q[0];
sx q[0];
rz(-1.3711509) q[0];
sx q[0];
rz(-3.1398849) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.121071) q[2];
sx q[2];
rz(-1.1164718) q[2];
sx q[2];
rz(-3.0128535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6807032) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(0.99980385) q[1];
rz(0.094035427) q[3];
sx q[3];
rz(-1.7574851) q[3];
sx q[3];
rz(1.3446913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38561884) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(2.8258064) q[0];
rz(0.93859998) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(-0.25207239) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2825851) q[0];
sx q[0];
rz(-0.95882817) q[0];
sx q[0];
rz(0.99143272) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2601068) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(1.2607247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6985059) q[1];
sx q[1];
rz(-1.7585187) q[1];
sx q[1];
rz(0.94633533) q[1];
rz(1.2139981) q[3];
sx q[3];
rz(-1.1988415) q[3];
sx q[3];
rz(1.4043722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(-1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144192) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(1.5267641) q[0];
rz(1.0871672) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(0.70708752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93543816) q[0];
sx q[0];
rz(-1.9618789) q[0];
sx q[0];
rz(-0.48396707) q[0];
x q[1];
rz(-0.44666501) q[2];
sx q[2];
rz(-1.457505) q[2];
sx q[2];
rz(-0.54892533) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69869631) q[1];
sx q[1];
rz(-0.92123182) q[1];
sx q[1];
rz(2.991308) q[1];
rz(-pi) q[2];
rz(-0.96890038) q[3];
sx q[3];
rz(-2.7728191) q[3];
sx q[3];
rz(-0.2975279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84918555) q[2];
sx q[2];
rz(-1.9533998) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(-0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29339644) q[0];
sx q[0];
rz(-0.080349803) q[0];
sx q[0];
rz(-1.3991762) q[0];
rz(-pi) q[1];
rz(2.1483634) q[2];
sx q[2];
rz(-1.5830056) q[2];
sx q[2];
rz(-0.73405594) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.188272) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(0.377368) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0765431) q[3];
sx q[3];
rz(-2.3688865) q[3];
sx q[3];
rz(0.62437526) q[3];
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
rz(-1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(-0.29944637) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(-0.20656955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8812013) q[0];
sx q[0];
rz(-0.83575373) q[0];
sx q[0];
rz(-0.3221237) q[0];
rz(-pi) q[1];
rz(-0.39482306) q[2];
sx q[2];
rz(-1.0770814) q[2];
sx q[2];
rz(1.0600818) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8787074) q[1];
sx q[1];
rz(-0.72808121) q[1];
sx q[1];
rz(1.8379777) q[1];
rz(-pi) q[2];
rz(0.91512485) q[3];
sx q[3];
rz(-1.8135999) q[3];
sx q[3];
rz(-1.9608998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(-0.28277961) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(0.16684428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151826) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-2.5577257) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(-1.3279703) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751942) q[0];
sx q[0];
rz(-1.1325206) q[0];
sx q[0];
rz(-0.12624329) q[0];
x q[1];
rz(-2.6793924) q[2];
sx q[2];
rz(-2.1415347) q[2];
sx q[2];
rz(-2.5525023) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77818645) q[1];
sx q[1];
rz(-1.4038329) q[1];
sx q[1];
rz(-0.88552514) q[1];
x q[2];
rz(0.51289576) q[3];
sx q[3];
rz(-1.1937871) q[3];
sx q[3];
rz(0.60037724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-2.1906733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.788818) q[0];
sx q[0];
rz(-2.7633861) q[0];
sx q[0];
rz(-0.87120716) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3046706) q[2];
sx q[2];
rz(-1.6309435) q[2];
sx q[2];
rz(-1.6595449) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7656895) q[1];
sx q[1];
rz(-1.6357058) q[1];
sx q[1];
rz(-2.9325571) q[1];
x q[2];
rz(2.7637134) q[3];
sx q[3];
rz(-1.1933019) q[3];
sx q[3];
rz(-2.8187403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(0.075909464) q[2];
rz(-0.54801303) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(0.65972796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83803672) q[0];
sx q[0];
rz(-2.3065901) q[0];
sx q[0];
rz(-2.3124218) q[0];
x q[1];
rz(0.78166878) q[2];
sx q[2];
rz(-1.4147007) q[2];
sx q[2];
rz(-2.8761656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5539726) q[1];
sx q[1];
rz(-2.4269322) q[1];
sx q[1];
rz(2.0228533) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4755867) q[3];
sx q[3];
rz(-2.6009437) q[3];
sx q[3];
rz(1.3117865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62347162) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(2.1155604) q[2];
rz(-0.14885151) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.898107) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(1.0356888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.068533) q[0];
sx q[0];
rz(-1.8951299) q[0];
sx q[0];
rz(-2.7102094) q[0];
rz(-pi) q[1];
rz(-2.2220988) q[2];
sx q[2];
rz(-0.80923015) q[2];
sx q[2];
rz(0.86916718) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.90696883) q[1];
sx q[1];
rz(-1.3303489) q[1];
sx q[1];
rz(0.30827) q[1];
x q[2];
rz(-0.74929897) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(-1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62853652) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7619027) q[0];
sx q[0];
rz(-1.3544449) q[0];
sx q[0];
rz(0.83308573) q[0];
rz(1.5079386) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(-0.41082906) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(-1.5796173) q[3];
sx q[3];
rz(-2.1891441) q[3];
sx q[3];
rz(1.7563663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
