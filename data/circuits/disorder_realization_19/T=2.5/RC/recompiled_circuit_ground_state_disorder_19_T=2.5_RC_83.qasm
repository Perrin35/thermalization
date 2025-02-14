OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1524393) q[0];
sx q[0];
rz(-1.5751155) q[0];
sx q[0];
rz(-2.0164665) q[0];
rz(-2.1195124) q[1];
sx q[1];
rz(-2.4740969) q[1];
sx q[1];
rz(-0.8134841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047377) q[0];
sx q[0];
rz(-1.2053524) q[0];
sx q[0];
rz(-2.9904537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44806077) q[2];
sx q[2];
rz(-1.8984814) q[2];
sx q[2];
rz(-2.0304012) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95619394) q[1];
sx q[1];
rz(-2.6743202) q[1];
sx q[1];
rz(-2.8059792) q[1];
rz(-pi) q[2];
rz(-1.4072106) q[3];
sx q[3];
rz(-1.4814875) q[3];
sx q[3];
rz(-0.79998868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(2.2129464) q[2];
rz(-1.5422025) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(-1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20392513) q[0];
sx q[0];
rz(-1.3805905) q[0];
sx q[0];
rz(-3.0145338) q[0];
rz(-0.98310414) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(-2.3703221) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255804) q[0];
sx q[0];
rz(-0.86992747) q[0];
sx q[0];
rz(-0.46937816) q[0];
x q[1];
rz(-1.2331687) q[2];
sx q[2];
rz(-1.513483) q[2];
sx q[2];
rz(-1.415451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88395547) q[1];
sx q[1];
rz(-1.9169382) q[1];
sx q[1];
rz(-0.37948541) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0326321) q[3];
sx q[3];
rz(-1.0618883) q[3];
sx q[3];
rz(-2.5050688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9942921) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(-3.1331983) q[2];
rz(2.4781135) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(0.26507637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0405149) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(-0.44152942) q[0];
rz(2.1532374) q[1];
sx q[1];
rz(-1.1343196) q[1];
sx q[1];
rz(3.006014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0975453) q[0];
sx q[0];
rz(-1.5870461) q[0];
sx q[0];
rz(1.5784997) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1889815) q[2];
sx q[2];
rz(-2.2509607) q[2];
sx q[2];
rz(-1.176468) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8959055) q[1];
sx q[1];
rz(-2.062633) q[1];
sx q[1];
rz(0.40426429) q[1];
rz(-2.3051585) q[3];
sx q[3];
rz(-1.423701) q[3];
sx q[3];
rz(-2.7874352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2077937) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(-3.1033893) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(-2.4562522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.4811089) q[0];
sx q[0];
rz(-2.3479192) q[0];
sx q[0];
rz(-2.5307181) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.7673312) q[1];
sx q[1];
rz(-2.9023721) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0804028) q[0];
sx q[0];
rz(-1.0188658) q[0];
sx q[0];
rz(0.30776382) q[0];
rz(-pi) q[1];
rz(-1.9630505) q[2];
sx q[2];
rz(-2.0118656) q[2];
sx q[2];
rz(-2.6792996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.071453302) q[1];
sx q[1];
rz(-1.8932559) q[1];
sx q[1];
rz(-0.47688213) q[1];
x q[2];
rz(1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(0.1399006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4227582) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(-0.0017496721) q[2];
rz(0.29378978) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(-2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9002429) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(-0.84306651) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(1.6036124) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085097236) q[0];
sx q[0];
rz(-0.81696327) q[0];
sx q[0];
rz(0.2635862) q[0];
x q[1];
rz(-0.17436738) q[2];
sx q[2];
rz(-1.5182759) q[2];
sx q[2];
rz(2.8157521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5953491) q[1];
sx q[1];
rz(-1.8270632) q[1];
sx q[1];
rz(-1.3687737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6513651) q[3];
sx q[3];
rz(-1.779883) q[3];
sx q[3];
rz(0.63864708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43188492) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(-2.0951927) q[2];
rz(2.8228068) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-2.5939202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043561291) q[0];
sx q[0];
rz(-1.2579608) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(1.7421534) q[1];
sx q[1];
rz(-1.6439227) q[1];
sx q[1];
rz(1.3290149) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15547046) q[0];
sx q[0];
rz(-0.24674812) q[0];
sx q[0];
rz(-1.107649) q[0];
x q[1];
rz(1.7850661) q[2];
sx q[2];
rz(-1.8545215) q[2];
sx q[2];
rz(-1.8510712) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.15550286) q[1];
sx q[1];
rz(-0.8900607) q[1];
sx q[1];
rz(0.40386856) q[1];
x q[2];
rz(1.9372609) q[3];
sx q[3];
rz(-2.0559337) q[3];
sx q[3];
rz(-1.7462891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38137388) q[2];
sx q[2];
rz(-1.2083283) q[2];
sx q[2];
rz(1.5927429) q[2];
rz(3.0717487) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(2.2729592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0233362) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(0.94183952) q[0];
rz(2.9815004) q[1];
sx q[1];
rz(-1.6611049) q[1];
sx q[1];
rz(2.9494185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1051467) q[0];
sx q[0];
rz(-0.22686401) q[0];
sx q[0];
rz(1.731621) q[0];
x q[1];
rz(-0.040080796) q[2];
sx q[2];
rz(-2.0922086) q[2];
sx q[2];
rz(1.557795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.77171626) q[1];
sx q[1];
rz(-1.9100921) q[1];
sx q[1];
rz(-0.75668613) q[1];
rz(-pi) q[2];
rz(-2.0683769) q[3];
sx q[3];
rz(-1.5514152) q[3];
sx q[3];
rz(1.7817225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3840702) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(-2.7080217) q[2];
rz(1.5844257) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(-2.7092194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.3766377) q[0];
sx q[0];
rz(1.3442511) q[0];
rz(1.2035707) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(-1.7787836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840795) q[0];
sx q[0];
rz(-1.6005033) q[0];
sx q[0];
rz(-0.021935181) q[0];
rz(0.8261189) q[2];
sx q[2];
rz(-1.4419793) q[2];
sx q[2];
rz(1.6507738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2143933) q[1];
sx q[1];
rz(-1.3424557) q[1];
sx q[1];
rz(-1.6990927) q[1];
rz(-1.1892494) q[3];
sx q[3];
rz(-2.4182712) q[3];
sx q[3];
rz(-1.2398694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56889304) q[2];
sx q[2];
rz(-2.3685444) q[2];
sx q[2];
rz(0.9838689) q[2];
rz(0.24108663) q[3];
sx q[3];
rz(-2.2086996) q[3];
sx q[3];
rz(-2.4510395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47138658) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(-0.27467003) q[0];
rz(-1.7999016) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(-2.0955657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5902025) q[0];
sx q[0];
rz(-0.98197637) q[0];
sx q[0];
rz(-1.9639652) q[0];
x q[1];
rz(1.3730896) q[2];
sx q[2];
rz(-1.4914163) q[2];
sx q[2];
rz(-2.6153836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0053802) q[1];
sx q[1];
rz(-0.58717218) q[1];
sx q[1];
rz(0.11759742) q[1];
rz(-pi) q[2];
rz(-1.6469429) q[3];
sx q[3];
rz(-0.5042432) q[3];
sx q[3];
rz(-1.7011736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(2.7421303) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-0.96650201) q[3];
sx q[3];
rz(-2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852916) q[0];
sx q[0];
rz(-1.4194019) q[0];
sx q[0];
rz(-2.5592819) q[0];
rz(-2.9583926) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(-0.80642548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2516625) q[0];
sx q[0];
rz(-1.4668873) q[0];
sx q[0];
rz(0.088168747) q[0];
rz(2.8060032) q[2];
sx q[2];
rz(-1.3177455) q[2];
sx q[2];
rz(-2.0584681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6567071) q[1];
sx q[1];
rz(-2.2995512) q[1];
sx q[1];
rz(-1.2970112) q[1];
rz(-pi) q[2];
rz(-0.29911689) q[3];
sx q[3];
rz(-2.7795305) q[3];
sx q[3];
rz(2.5821834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9091984) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(2.3363028) q[2];
rz(-0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(-2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88937) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(-1.5994785) q[1];
sx q[1];
rz(-1.5442994) q[1];
sx q[1];
rz(-1.50179) q[1];
rz(-2.2940214) q[2];
sx q[2];
rz(-1.6568274) q[2];
sx q[2];
rz(-0.33100707) q[2];
rz(2.6977758) q[3];
sx q[3];
rz(-2.4469821) q[3];
sx q[3];
rz(-1.2367005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
