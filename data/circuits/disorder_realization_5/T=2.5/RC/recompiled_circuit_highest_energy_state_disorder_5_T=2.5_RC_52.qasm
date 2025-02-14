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
rz(-0.37762168) q[0];
sx q[0];
rz(3.5701516) q[0];
sx q[0];
rz(9.188434) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(-0.16170391) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40370052) q[0];
sx q[0];
rz(-2.1503415) q[0];
sx q[0];
rz(-3.011322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19585545) q[2];
sx q[2];
rz(-1.5136949) q[2];
sx q[2];
rz(-1.7498992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.573805) q[1];
sx q[1];
rz(-1.5638788) q[1];
sx q[1];
rz(-1.597639) q[1];
rz(-pi) q[2];
rz(-0.51689619) q[3];
sx q[3];
rz(-2.7166043) q[3];
sx q[3];
rz(-1.4864511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8990495) q[2];
sx q[2];
rz(-0.87729064) q[2];
sx q[2];
rz(0.38929942) q[2];
rz(2.9111275) q[3];
sx q[3];
rz(-0.018229818) q[3];
sx q[3];
rz(-2.5267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56735754) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(1.6472598) q[0];
rz(-1.5859454) q[1];
sx q[1];
rz(-0.21265282) q[1];
sx q[1];
rz(-2.0071323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8788505) q[0];
sx q[0];
rz(-1.8244184) q[0];
sx q[0];
rz(-2.4154759) q[0];
x q[1];
rz(2.6642642) q[2];
sx q[2];
rz(-2.3466718) q[2];
sx q[2];
rz(1.8177462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33541533) q[1];
sx q[1];
rz(-1.9796438) q[1];
sx q[1];
rz(-0.93542904) q[1];
rz(-pi) q[2];
rz(-1.0026917) q[3];
sx q[3];
rz(-0.8352931) q[3];
sx q[3];
rz(-0.42772528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.050934164) q[2];
sx q[2];
rz(-2.2218573) q[2];
sx q[2];
rz(-1.301379) q[2];
rz(1.0743514) q[3];
sx q[3];
rz(-2.8315872) q[3];
sx q[3];
rz(1.5955135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12304561) q[0];
sx q[0];
rz(-2.8392241) q[0];
sx q[0];
rz(0.61755919) q[0];
rz(-2.0410208) q[1];
sx q[1];
rz(-0.01958422) q[1];
sx q[1];
rz(2.7380131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27173938) q[0];
sx q[0];
rz(-1.4560486) q[0];
sx q[0];
rz(-3.1310351) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6721272) q[2];
sx q[2];
rz(-1.0472385) q[2];
sx q[2];
rz(3.0803185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.523905) q[1];
sx q[1];
rz(-2.9343961) q[1];
sx q[1];
rz(-2.2788384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37637122) q[3];
sx q[3];
rz(-0.59089243) q[3];
sx q[3];
rz(-1.2749689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5201716) q[2];
sx q[2];
rz(-1.6941864) q[2];
sx q[2];
rz(2.3604895) q[2];
rz(2.805294) q[3];
sx q[3];
rz(-1.7016442) q[3];
sx q[3];
rz(0.70359105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4562562) q[0];
sx q[0];
rz(-0.51613134) q[0];
sx q[0];
rz(1.2180895) q[0];
rz(-1.7958027) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(-1.2340612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06935519) q[0];
sx q[0];
rz(-1.580582) q[0];
sx q[0];
rz(0.18982498) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61361112) q[2];
sx q[2];
rz(-1.3296812) q[2];
sx q[2];
rz(-0.63617731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30013903) q[1];
sx q[1];
rz(-2.6965504) q[1];
sx q[1];
rz(0.95325327) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5640573) q[3];
sx q[3];
rz(-1.5608265) q[3];
sx q[3];
rz(0.69025595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1951083) q[2];
sx q[2];
rz(-0.4230963) q[2];
sx q[2];
rz(-1.1484324) q[2];
rz(-2.3871683) q[3];
sx q[3];
rz(-1.9781338) q[3];
sx q[3];
rz(-0.26328009) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36963439) q[0];
sx q[0];
rz(-3.0535871) q[0];
sx q[0];
rz(2.5391286) q[0];
rz(2.1594436) q[1];
sx q[1];
rz(-0.0063449675) q[1];
sx q[1];
rz(-2.6314661) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0629629) q[0];
sx q[0];
rz(-1.3853243) q[0];
sx q[0];
rz(-1.3736891) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.449376) q[2];
sx q[2];
rz(-1.4842786) q[2];
sx q[2];
rz(-1.0457116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1044961) q[1];
sx q[1];
rz(-1.949661) q[1];
sx q[1];
rz(-2.563743) q[1];
x q[2];
rz(0.81437935) q[3];
sx q[3];
rz(-1.5703716) q[3];
sx q[3];
rz(-2.9546776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0977352) q[2];
sx q[2];
rz(-2.0243702) q[2];
sx q[2];
rz(2.1065693) q[2];
rz(-1.6739316) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(-2.5586832) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408836) q[0];
sx q[0];
rz(-0.98475921) q[0];
sx q[0];
rz(-2.050052) q[0];
rz(2.7409399) q[1];
sx q[1];
rz(-3.1392097) q[1];
sx q[1];
rz(-0.56652743) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7687205) q[0];
sx q[0];
rz(-2.3864288) q[0];
sx q[0];
rz(-1.7674957) q[0];
rz(-2.3344757) q[2];
sx q[2];
rz(-2.2645778) q[2];
sx q[2];
rz(-1.6859695) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57921806) q[1];
sx q[1];
rz(-1.3462344) q[1];
sx q[1];
rz(2.5689305) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6831538) q[3];
sx q[3];
rz(-1.0010202) q[3];
sx q[3];
rz(-2.7801685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94622508) q[2];
sx q[2];
rz(-0.75104284) q[2];
sx q[2];
rz(-1.6132149) q[2];
rz(1.001312) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(2.9010469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85750759) q[0];
sx q[0];
rz(-0.79428285) q[0];
sx q[0];
rz(1.9965782) q[0];
rz(-1.6192294) q[1];
sx q[1];
rz(-0.018773627) q[1];
sx q[1];
rz(1.1633263) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6451241) q[0];
sx q[0];
rz(-1.860052) q[0];
sx q[0];
rz(3.0933063) q[0];
x q[1];
rz(-1.3848035) q[2];
sx q[2];
rz(-2.459068) q[2];
sx q[2];
rz(2.1738717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8790354) q[1];
sx q[1];
rz(-1.702629) q[1];
sx q[1];
rz(-2.1798032) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9807253) q[3];
sx q[3];
rz(-0.64023521) q[3];
sx q[3];
rz(1.7750027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5160949) q[2];
sx q[2];
rz(-2.1558546) q[2];
sx q[2];
rz(-2.7682448) q[2];
rz(0.18800023) q[3];
sx q[3];
rz(-1.5865654) q[3];
sx q[3];
rz(1.9787623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2242551) q[0];
sx q[0];
rz(-0.52977109) q[0];
sx q[0];
rz(2.8944471) q[0];
rz(-0.22357926) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(-1.6252958) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14683826) q[0];
sx q[0];
rz(-0.13299599) q[0];
sx q[0];
rz(2.07703) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3749977) q[2];
sx q[2];
rz(-1.261289) q[2];
sx q[2];
rz(2.7628037) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5471955) q[1];
sx q[1];
rz(-1.5548348) q[1];
sx q[1];
rz(0.39401024) q[1];
rz(-1.5228294) q[3];
sx q[3];
rz(-0.9771416) q[3];
sx q[3];
rz(-1.1831212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94459263) q[2];
sx q[2];
rz(-0.71114117) q[2];
sx q[2];
rz(-1.5947343) q[2];
rz(-0.77557766) q[3];
sx q[3];
rz(-1.2838793) q[3];
sx q[3];
rz(-1.4625134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32770661) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(-1.1073329) q[0];
rz(-0.92512098) q[1];
sx q[1];
rz(-0.0019625891) q[1];
sx q[1];
rz(2.3855239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32996477) q[0];
sx q[0];
rz(-2.1058308) q[0];
sx q[0];
rz(-2.7828548) q[0];
rz(-pi) q[1];
rz(-1.5599361) q[2];
sx q[2];
rz(-0.53087044) q[2];
sx q[2];
rz(0.52299352) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1518095) q[1];
sx q[1];
rz(-1.0149455) q[1];
sx q[1];
rz(-1.226306) q[1];
rz(-pi) q[2];
x q[2];
rz(2.954079) q[3];
sx q[3];
rz(-2.5616259) q[3];
sx q[3];
rz(-0.43646508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6206616) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(-1.0010285) q[2];
rz(-1.777461) q[3];
sx q[3];
rz(-2.2083211) q[3];
sx q[3];
rz(-2.6076243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.119568) q[0];
sx q[0];
rz(-1.3725932) q[0];
sx q[0];
rz(0.46510988) q[0];
rz(1.8067092) q[1];
sx q[1];
rz(-2.7692134) q[1];
sx q[1];
rz(-1.575527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4441285) q[0];
sx q[0];
rz(-2.9078472) q[0];
sx q[0];
rz(1.3819737) q[0];
rz(-pi) q[1];
rz(1.0687629) q[2];
sx q[2];
rz(-0.58922592) q[2];
sx q[2];
rz(3.0844944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82619691) q[1];
sx q[1];
rz(-0.0028571833) q[1];
sx q[1];
rz(0.14551659) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0483143) q[3];
sx q[3];
rz(-0.3808379) q[3];
sx q[3];
rz(0.24254984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89432013) q[2];
sx q[2];
rz(-0.044581052) q[2];
sx q[2];
rz(-2.0559922) q[2];
rz(1.9191437) q[3];
sx q[3];
rz(-0.51210755) q[3];
sx q[3];
rz(-1.2781757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4492252) q[0];
sx q[0];
rz(-1.4341555) q[0];
sx q[0];
rz(1.7644802) q[0];
rz(-1.5832681) q[1];
sx q[1];
rz(-0.91455864) q[1];
sx q[1];
rz(0.22462489) q[1];
rz(-3.0947826) q[2];
sx q[2];
rz(-0.10086127) q[2];
sx q[2];
rz(0.33831102) q[2];
rz(-0.23595594) q[3];
sx q[3];
rz(-2.0045842) q[3];
sx q[3];
rz(1.830238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
