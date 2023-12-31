OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18626285) q[0];
sx q[0];
rz(-1.3583399) q[0];
sx q[0];
rz(2.0583378) q[0];
rz(0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(1.1320621) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.010477) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(-2.1285776) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16529103) q[3];
sx q[3];
rz(-0.2632907) q[3];
sx q[3];
rz(-2.1804682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(-2.9361172) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(-1.1516494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8173556) q[0];
sx q[0];
rz(-1.865987) q[0];
sx q[0];
rz(-0.8582219) q[0];
rz(-pi) q[1];
rz(1.0863016) q[2];
sx q[2];
rz(-2.0995579) q[2];
sx q[2];
rz(0.27258401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4558251) q[1];
sx q[1];
rz(-1.1357422) q[1];
sx q[1];
rz(2.0352092) q[1];
x q[2];
rz(-0.63668164) q[3];
sx q[3];
rz(-2.5425306) q[3];
sx q[3];
rz(0.77081313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(2.9197664) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31056988) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6037613) q[0];
sx q[0];
rz(-2.1268401) q[0];
sx q[0];
rz(-2.6105196) q[0];
rz(-0.25643202) q[2];
sx q[2];
rz(-1.4415381) q[2];
sx q[2];
rz(1.3916707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1287071) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(0.30602869) q[1];
rz(-pi) q[2];
rz(0.49719663) q[3];
sx q[3];
rz(-0.70780863) q[3];
sx q[3];
rz(-1.4149815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-0.46359584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620561) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(-0.38328538) q[0];
x q[1];
rz(-0.33275231) q[2];
sx q[2];
rz(-1.6289662) q[2];
sx q[2];
rz(-2.1259049) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5220118) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(-2.3492858) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2673244) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(2.1972426) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1221065) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(-1.7800063) q[0];
rz(-pi) q[1];
rz(1.1733426) q[2];
sx q[2];
rz(-2.3059418) q[2];
sx q[2];
rz(-0.022692516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9152865) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(-1.5096942) q[1];
rz(1.827042) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(-2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(-3.1392858) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(2.5571402) q[0];
rz(0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-0.054919682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829837) q[0];
sx q[0];
rz(-0.28124547) q[0];
sx q[0];
rz(-1.2693229) q[0];
rz(-0.018718406) q[2];
sx q[2];
rz(-1.2476377) q[2];
sx q[2];
rz(1.2013555) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1850486) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(0.59021414) q[1];
rz(0.46338007) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(0.86138553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75446689) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(-0.91032666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9757662) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(1.635701) q[0];
rz(-1.9867284) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(1.8679801) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4714204) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(2.7736204) q[1];
x q[2];
rz(2.5038239) q[3];
sx q[3];
rz(-2.31156) q[3];
sx q[3];
rz(-0.72731599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-0.30050373) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43113118) q[0];
sx q[0];
rz(-0.80695242) q[0];
sx q[0];
rz(-2.0956844) q[0];
rz(0.76086107) q[2];
sx q[2];
rz(-2.0458474) q[2];
sx q[2];
rz(-2.8617815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0771675) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(-2.6618631) q[1];
rz(-pi) q[2];
rz(0.81340202) q[3];
sx q[3];
rz(-1.8040931) q[3];
sx q[3];
rz(2.5113311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(0.78197455) q[2];
rz(-2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-0.75884563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1369143) q[0];
sx q[0];
rz(-1.5537964) q[0];
sx q[0];
rz(1.9949811) q[0];
rz(-2.8685832) q[2];
sx q[2];
rz(-0.62676478) q[2];
sx q[2];
rz(-0.11944709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89275937) q[1];
sx q[1];
rz(-1.4151238) q[1];
sx q[1];
rz(-2.4397736) q[1];
rz(-pi) q[2];
rz(1.5893448) q[3];
sx q[3];
rz(-0.66569257) q[3];
sx q[3];
rz(1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9655351) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(-2.7479991) q[0];
rz(-pi) q[1];
rz(-2.7650325) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(3.0388447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3366821) q[1];
sx q[1];
rz(-1.33178) q[1];
sx q[1];
rz(-2.3906624) q[1];
x q[2];
rz(-0.742357) q[3];
sx q[3];
rz(-1.9543813) q[3];
sx q[3];
rz(-0.63391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(0.71371901) q[2];
rz(-2.7632726) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(3.1109839) q[2];
sx q[2];
rz(-1.3749214) q[2];
sx q[2];
rz(2.2236852) q[2];
rz(1.1013423) q[3];
sx q[3];
rz(-0.51870844) q[3];
sx q[3];
rz(-1.4389256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
