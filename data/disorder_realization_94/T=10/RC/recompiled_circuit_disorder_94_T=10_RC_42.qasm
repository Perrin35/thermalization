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
rz(-0.59626055) q[1];
sx q[1];
rz(0.61520666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9278487) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(-2.5736789) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46703672) q[2];
sx q[2];
rz(-0.39614284) q[2];
sx q[2];
rz(0.064183891) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33949172) q[1];
sx q[1];
rz(-1.9994945) q[1];
sx q[1];
rz(-2.2080253) q[1];
x q[2];
rz(-1.7513566) q[3];
sx q[3];
rz(-1.872634) q[3];
sx q[3];
rz(2.3328822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7011828) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(3.1112444) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.5240086) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47027662) q[0];
sx q[0];
rz(-2.94194) q[0];
sx q[0];
rz(-1.5792363) q[0];
x q[1];
rz(0.020521684) q[2];
sx q[2];
rz(-1.1164718) q[2];
sx q[2];
rz(-3.0128535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9873725) q[1];
sx q[1];
rz(-1.0013354) q[1];
sx q[1];
rz(-3.0595368) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0475572) q[3];
sx q[3];
rz(-1.7574851) q[3];
sx q[3];
rz(1.7969014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38561884) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(-1.8148445) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(0.25207239) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2825851) q[0];
sx q[0];
rz(-0.95882817) q[0];
sx q[0];
rz(-2.1501599) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2601068) q[2];
sx q[2];
rz(-2.2027317) q[2];
sx q[2];
rz(1.880868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38110106) q[1];
sx q[1];
rz(-2.4931506) q[1];
sx q[1];
rz(-1.8849424) q[1];
rz(-2.7470845) q[3];
sx q[3];
rz(-1.9022226) q[3];
sx q[3];
rz(-3.1098207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(-1.5267641) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2061545) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(-2.6576256) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44666501) q[2];
sx q[2];
rz(-1.457505) q[2];
sx q[2];
rz(-2.5926673) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1781588) q[1];
sx q[1];
rz(-1.6903094) q[1];
sx q[1];
rz(0.91576373) q[1];
rz(-pi) q[2];
rz(-1.8791734) q[3];
sx q[3];
rz(-1.3652507) q[3];
sx q[3];
rz(1.8431078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(-0.24266711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.8206772) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.7472349) q[0];
rz(-1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-2.8869693) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8481962) q[0];
sx q[0];
rz(-3.0612429) q[0];
sx q[0];
rz(-1.7424165) q[0];
rz(-pi) q[1];
rz(1.5931555) q[2];
sx q[2];
rz(-0.5776814) q[2];
sx q[2];
rz(2.3235842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8498807) q[1];
sx q[1];
rz(-0.92392081) q[1];
sx q[1];
rz(-1.8748267) q[1];
rz(2.0650495) q[3];
sx q[3];
rz(-2.3688865) q[3];
sx q[3];
rz(-2.5172174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1191117) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45143932) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(0.29944637) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(2.9350231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8812013) q[0];
sx q[0];
rz(-0.83575373) q[0];
sx q[0];
rz(-2.819469) q[0];
rz(-0.39482306) q[2];
sx q[2];
rz(-2.0645112) q[2];
sx q[2];
rz(2.0815108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8787074) q[1];
sx q[1];
rz(-2.4135114) q[1];
sx q[1];
rz(1.8379777) q[1];
rz(-pi) q[2];
rz(-2.2264678) q[3];
sx q[3];
rz(-1.8135999) q[3];
sx q[3];
rz(-1.9608998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(2.9747484) q[3];
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
x q[0];
rz(-pi) q[1];
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
rz(-0.54324141) q[1];
sx q[1];
rz(1.3279703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8570003) q[0];
sx q[0];
rz(-2.6866331) q[0];
sx q[0];
rz(-1.3083463) q[0];
x q[1];
rz(0.46220025) q[2];
sx q[2];
rz(-1.000058) q[2];
sx q[2];
rz(2.5525023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.92762891) q[1];
sx q[1];
rz(-0.8968401) q[1];
sx q[1];
rz(2.9272635) q[1];
rz(-pi) q[2];
rz(1.9973203) q[3];
sx q[3];
rz(-1.0970308) q[3];
sx q[3];
rz(-1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61838377) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(-1.3605114) q[2];
rz(1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1469864) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(2.1906733) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(2.8911203) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3046706) q[2];
sx q[2];
rz(-1.5106491) q[2];
sx q[2];
rz(1.6595449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2086522) q[1];
sx q[1];
rz(-1.7793852) q[1];
sx q[1];
rz(-1.637146) q[1];
rz(-pi) q[2];
rz(-0.82138737) q[3];
sx q[3];
rz(-0.52762023) q[3];
sx q[3];
rz(1.1452831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2237079) q[2];
sx q[2];
rz(-0.48030883) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(1.7653718) q[0];
rz(0.41704047) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(2.4818647) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0411471) q[0];
sx q[0];
rz(-0.99248306) q[0];
sx q[0];
rz(-0.64097494) q[0];
rz(-pi) q[1];
rz(-2.3599239) q[2];
sx q[2];
rz(-1.726892) q[2];
sx q[2];
rz(2.8761656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1253423) q[1];
sx q[1];
rz(-2.2012735) q[1];
sx q[1];
rz(0.36228212) q[1];
rz(-pi) q[2];
rz(1.0321484) q[3];
sx q[3];
rz(-1.5218471) q[3];
sx q[3];
rz(0.34070542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62347162) q[2];
sx q[2];
rz(-2.3770964) q[2];
sx q[2];
rz(2.1155604) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(-1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730597) q[0];
sx q[0];
rz(-1.2464628) q[0];
sx q[0];
rz(-2.7102094) q[0];
rz(2.5752441) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(1.4373506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2346238) q[1];
sx q[1];
rz(-1.3303489) q[1];
sx q[1];
rz(-2.8333227) q[1];
rz(-pi) q[2];
rz(2.3922937) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(-1.1009969) q[3];
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
rz(-1.1516085) q[2];
rz(0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
rz(1.6133616) q[2];
sx q[2];
rz(-1.1602989) q[2];
sx q[2];
rz(1.5699671) q[2];
rz(0.61836615) q[3];
sx q[3];
rz(-1.5779839) q[3];
sx q[3];
rz(0.18045651) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];