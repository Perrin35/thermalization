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
rz(2.7163765) q[0];
sx q[0];
rz(-1.334231) q[0];
sx q[0];
rz(0.38036007) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(-0.18667297) q[1];
sx q[1];
rz(0.92460728) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43792576) q[0];
sx q[0];
rz(-2.3406918) q[0];
sx q[0];
rz(0.4842224) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1670465) q[2];
sx q[2];
rz(-1.0606597) q[2];
sx q[2];
rz(1.8710006) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8421665) q[1];
sx q[1];
rz(-1.5314845) q[1];
sx q[1];
rz(-0.8419324) q[1];
rz(0.070234039) q[3];
sx q[3];
rz(-3.0133504) q[3];
sx q[3];
rz(1.2508891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0164464) q[2];
sx q[2];
rz(-2.1852198) q[2];
sx q[2];
rz(-0.98801405) q[2];
rz(0.2068578) q[3];
sx q[3];
rz(-1.4068973) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67966953) q[0];
sx q[0];
rz(-1.672687) q[0];
sx q[0];
rz(-0.2521421) q[0];
rz(-3.120046) q[1];
sx q[1];
rz(-2.4485059) q[1];
sx q[1];
rz(0.85539877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148672) q[0];
sx q[0];
rz(-1.572319) q[0];
sx q[0];
rz(-1.5437838) q[0];
x q[1];
rz(-0.62521817) q[2];
sx q[2];
rz(-2.4897235) q[2];
sx q[2];
rz(-0.21833459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.508042) q[1];
sx q[1];
rz(-1.7910069) q[1];
sx q[1];
rz(2.1160265) q[1];
rz(-1.2367494) q[3];
sx q[3];
rz(-2.2260087) q[3];
sx q[3];
rz(1.1545336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1713193) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(-0.26101905) q[2];
rz(0.039693443) q[3];
sx q[3];
rz(-1.7429765) q[3];
sx q[3];
rz(0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4259341) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(2.4470827) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(2.1515501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5607669) q[0];
sx q[0];
rz(-2.7832354) q[0];
sx q[0];
rz(2.906515) q[0];
x q[1];
rz(-0.98857356) q[2];
sx q[2];
rz(-1.7372088) q[2];
sx q[2];
rz(2.5022282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5591919) q[1];
sx q[1];
rz(-1.2123931) q[1];
sx q[1];
rz(0.47689516) q[1];
x q[2];
rz(-0.53931358) q[3];
sx q[3];
rz(-2.1997872) q[3];
sx q[3];
rz(-1.1719538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29828829) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(-2.6514371) q[2];
rz(0.35689029) q[3];
sx q[3];
rz(-0.39339742) q[3];
sx q[3];
rz(-1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056034293) q[0];
sx q[0];
rz(-1.1857251) q[0];
sx q[0];
rz(-0.84247843) q[0];
rz(2.339824) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(-0.78559771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89891035) q[0];
sx q[0];
rz(-1.8225095) q[0];
sx q[0];
rz(-0.71941914) q[0];
x q[1];
rz(-2.9688492) q[2];
sx q[2];
rz(-0.72451353) q[2];
sx q[2];
rz(2.8959664) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.468049) q[1];
sx q[1];
rz(-1.0338514) q[1];
sx q[1];
rz(1.2098486) q[1];
rz(2.9369041) q[3];
sx q[3];
rz(-0.83016268) q[3];
sx q[3];
rz(1.1799174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0164612) q[2];
sx q[2];
rz(-0.74840122) q[2];
sx q[2];
rz(-1.008519) q[2];
rz(-0.85159167) q[3];
sx q[3];
rz(-0.75751704) q[3];
sx q[3];
rz(-2.0612702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763024) q[0];
sx q[0];
rz(-2.1602614) q[0];
sx q[0];
rz(1.0705795) q[0];
rz(2.3640682) q[1];
sx q[1];
rz(-0.91064015) q[1];
sx q[1];
rz(1.0106962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3883481) q[0];
sx q[0];
rz(-1.0976315) q[0];
sx q[0];
rz(-1.1654668) q[0];
x q[1];
rz(-0.065344409) q[2];
sx q[2];
rz(-1.3962708) q[2];
sx q[2];
rz(-0.068152817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8292099) q[1];
sx q[1];
rz(-1.1151894) q[1];
sx q[1];
rz(2.4036951) q[1];
rz(-2.5794677) q[3];
sx q[3];
rz(-1.486612) q[3];
sx q[3];
rz(2.4474194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8396478) q[2];
sx q[2];
rz(-0.10493111) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(-1.7293845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3509336) q[0];
sx q[0];
rz(-2.9582294) q[0];
sx q[0];
rz(-0.80107981) q[0];
rz(-3.0084897) q[1];
sx q[1];
rz(-1.8201273) q[1];
sx q[1];
rz(2.7395693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6989845) q[0];
sx q[0];
rz(-0.8359209) q[0];
sx q[0];
rz(-2.0280272) q[0];
x q[1];
rz(-3.1394464) q[2];
sx q[2];
rz(-1.8615926) q[2];
sx q[2];
rz(0.69974297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61643663) q[1];
sx q[1];
rz(-2.0503938) q[1];
sx q[1];
rz(2.8308771) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4834255) q[3];
sx q[3];
rz(-1.6699381) q[3];
sx q[3];
rz(2.1226573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63783995) q[2];
sx q[2];
rz(-1.4053586) q[2];
sx q[2];
rz(2.3657738) q[2];
rz(1.6860298) q[3];
sx q[3];
rz(-1.173923) q[3];
sx q[3];
rz(-1.0078526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044947226) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(-2.4626379) q[0];
rz(1.8567765) q[1];
sx q[1];
rz(-2.4285451) q[1];
sx q[1];
rz(-1.7624034) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5011936) q[0];
sx q[0];
rz(-1.7441347) q[0];
sx q[0];
rz(-0.79357432) q[0];
x q[1];
rz(3.0581362) q[2];
sx q[2];
rz(-1.8493358) q[2];
sx q[2];
rz(2.4065774) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1393362) q[1];
sx q[1];
rz(-0.69076194) q[1];
sx q[1];
rz(0.73378566) q[1];
x q[2];
rz(-1.3706743) q[3];
sx q[3];
rz(-1.6255696) q[3];
sx q[3];
rz(1.9211662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8046367) q[2];
sx q[2];
rz(-0.75504428) q[2];
sx q[2];
rz(1.2389368) q[2];
rz(1.8094481) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-2.8968887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181344) q[0];
sx q[0];
rz(-1.1661538) q[0];
sx q[0];
rz(3.0031257) q[0];
rz(-0.18374099) q[1];
sx q[1];
rz(-2.467149) q[1];
sx q[1];
rz(2.8750681) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0040695) q[0];
sx q[0];
rz(-2.7752987) q[0];
sx q[0];
rz(2.0384602) q[0];
x q[1];
rz(2.1378072) q[2];
sx q[2];
rz(-1.3409753) q[2];
sx q[2];
rz(0.15574317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6147242) q[1];
sx q[1];
rz(-0.42335948) q[1];
sx q[1];
rz(-1.4335267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0098465482) q[3];
sx q[3];
rz(-0.72692633) q[3];
sx q[3];
rz(-1.3971412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0900241) q[2];
sx q[2];
rz(-1.8847621) q[2];
sx q[2];
rz(-0.23452342) q[2];
rz(1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(2.3852111) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804857) q[0];
sx q[0];
rz(-0.8661626) q[0];
sx q[0];
rz(2.3418703) q[0];
rz(-2.5526478) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(-2.934093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0711771) q[0];
sx q[0];
rz(-1.1780329) q[0];
sx q[0];
rz(1.429274) q[0];
rz(-pi) q[1];
rz(-1.2166747) q[2];
sx q[2];
rz(-2.3975044) q[2];
sx q[2];
rz(-0.085372774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60002414) q[1];
sx q[1];
rz(-0.80994058) q[1];
sx q[1];
rz(1.0066973) q[1];
rz(1.1010246) q[3];
sx q[3];
rz(-1.5964526) q[3];
sx q[3];
rz(-0.86650833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3346682) q[2];
sx q[2];
rz(-2.25756) q[2];
sx q[2];
rz(-2.7743288) q[2];
rz(1.7043097) q[3];
sx q[3];
rz(-1.1742914) q[3];
sx q[3];
rz(-1.9678763) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87840286) q[0];
sx q[0];
rz(-1.0028361) q[0];
sx q[0];
rz(-0.81018418) q[0];
rz(-2.0267678) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(2.5883163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5636175) q[0];
sx q[0];
rz(-1.5153236) q[0];
sx q[0];
rz(0.93240057) q[0];
x q[1];
rz(2.5279494) q[2];
sx q[2];
rz(-0.89362234) q[2];
sx q[2];
rz(0.26640642) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3685963) q[1];
sx q[1];
rz(-2.8008411) q[1];
sx q[1];
rz(-0.98253886) q[1];
x q[2];
rz(-0.78599591) q[3];
sx q[3];
rz(-1.3162426) q[3];
sx q[3];
rz(0.28721519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1222003) q[2];
sx q[2];
rz(-0.70851749) q[2];
sx q[2];
rz(2.9019287) q[2];
rz(-1.3278809) q[3];
sx q[3];
rz(-1.0646822) q[3];
sx q[3];
rz(-2.6740668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079530579) q[0];
sx q[0];
rz(-1.3539599) q[0];
sx q[0];
rz(0.39394105) q[0];
rz(2.486034) q[1];
sx q[1];
rz(-1.9656904) q[1];
sx q[1];
rz(0.94737731) q[1];
rz(-2.4782933) q[2];
sx q[2];
rz(-0.51021432) q[2];
sx q[2];
rz(-2.9205657) q[2];
rz(0.58047337) q[3];
sx q[3];
rz(-1.3571285) q[3];
sx q[3];
rz(0.50783689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
