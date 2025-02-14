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
rz(0.87585706) q[0];
sx q[0];
rz(2.1729204) q[0];
sx q[0];
rz(7.2090413) q[0];
rz(-2.5246188) q[1];
sx q[1];
rz(-2.4763835) q[1];
sx q[1];
rz(1.8170504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61875859) q[0];
sx q[0];
rz(-1.6346187) q[0];
sx q[0];
rz(0.60019779) q[0];
x q[1];
rz(-1.6683031) q[2];
sx q[2];
rz(-0.98819369) q[2];
sx q[2];
rz(-2.1511457) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9369389) q[1];
sx q[1];
rz(-1.1142245) q[1];
sx q[1];
rz(-0.66557933) q[1];
rz(-pi) q[2];
rz(2.2258513) q[3];
sx q[3];
rz(-2.2217882) q[3];
sx q[3];
rz(-0.30177339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2609743) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(0.79208881) q[2];
rz(-0.875862) q[3];
sx q[3];
rz(-0.10871092) q[3];
sx q[3];
rz(-2.47827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.51139128) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(2.200101) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(-0.082854465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64697826) q[0];
sx q[0];
rz(-0.74132338) q[0];
sx q[0];
rz(2.9927918) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24225927) q[2];
sx q[2];
rz(-1.2238811) q[2];
sx q[2];
rz(2.0478947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8585457) q[1];
sx q[1];
rz(-1.0632005) q[1];
sx q[1];
rz(-1.0495484) q[1];
rz(-2.4965246) q[3];
sx q[3];
rz(-2.0760025) q[3];
sx q[3];
rz(0.92407214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2443709) q[2];
sx q[2];
rz(-1.7873849) q[2];
sx q[2];
rz(-1.2574035) q[2];
rz(-0.8463549) q[3];
sx q[3];
rz(-0.1592764) q[3];
sx q[3];
rz(0.67811051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1988679) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(2.9123836) q[0];
rz(-0.21408679) q[1];
sx q[1];
rz(-0.87132088) q[1];
sx q[1];
rz(-2.7600938) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59807459) q[0];
sx q[0];
rz(-1.2208573) q[0];
sx q[0];
rz(2.475481) q[0];
x q[1];
rz(-1.3573285) q[2];
sx q[2];
rz(-0.74356438) q[2];
sx q[2];
rz(1.5876169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2974071) q[1];
sx q[1];
rz(-1.5995889) q[1];
sx q[1];
rz(-1.5178568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7945064) q[3];
sx q[3];
rz(-1.5659837) q[3];
sx q[3];
rz(-2.6322685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4182338) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(0.56824938) q[2];
rz(-1.7226284) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(-2.0645781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1133465) q[0];
sx q[0];
rz(-2.009511) q[0];
sx q[0];
rz(-0.25076732) q[0];
rz(-0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.9409404) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6687209) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(-1.0502104) q[0];
rz(-1.1293639) q[2];
sx q[2];
rz(-1.7953292) q[2];
sx q[2];
rz(-2.6972636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.975864) q[1];
sx q[1];
rz(-1.670125) q[1];
sx q[1];
rz(2.4286859) q[1];
x q[2];
rz(2.8765466) q[3];
sx q[3];
rz(-2.0816021) q[3];
sx q[3];
rz(0.068451133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67479977) q[2];
sx q[2];
rz(-2.0189221) q[2];
sx q[2];
rz(-1.8505081) q[2];
rz(0.42168266) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(-1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52294937) q[0];
sx q[0];
rz(-2.4186501) q[0];
sx q[0];
rz(-1.500754) q[0];
rz(-0.067642637) q[1];
sx q[1];
rz(-1.5539955) q[1];
sx q[1];
rz(1.1281475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2826506) q[0];
sx q[0];
rz(-2.9500486) q[0];
sx q[0];
rz(-0.19609496) q[0];
x q[1];
rz(0.35223799) q[2];
sx q[2];
rz(-1.0718498) q[2];
sx q[2];
rz(-1.7583454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.060202816) q[1];
sx q[1];
rz(-0.12559016) q[1];
sx q[1];
rz(1.1274028) q[1];
rz(-1.2913843) q[3];
sx q[3];
rz(-1.7959204) q[3];
sx q[3];
rz(-2.1373419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.03380123) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(-2.4647253) q[2];
rz(-0.90773165) q[3];
sx q[3];
rz(-1.2264484) q[3];
sx q[3];
rz(0.41456732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-0.75160471) q[0];
sx q[0];
rz(-0.76939097) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(-0.21533899) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7756336) q[0];
sx q[0];
rz(-1.7582408) q[0];
sx q[0];
rz(2.1523317) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5211283) q[2];
sx q[2];
rz(-2.81041) q[2];
sx q[2];
rz(2.0803723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29500719) q[1];
sx q[1];
rz(-0.93448105) q[1];
sx q[1];
rz(1.2656487) q[1];
x q[2];
rz(1.2329383) q[3];
sx q[3];
rz(-0.90150276) q[3];
sx q[3];
rz(-0.026642628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26663366) q[2];
sx q[2];
rz(-1.5077488) q[2];
sx q[2];
rz(-0.61666644) q[2];
rz(0.2002317) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013833372) q[0];
sx q[0];
rz(-2.5795689) q[0];
sx q[0];
rz(-0.01509893) q[0];
rz(-3.0191782) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(2.1133568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6262344) q[0];
sx q[0];
rz(-1.6263824) q[0];
sx q[0];
rz(-0.65599156) q[0];
rz(-0.21034849) q[2];
sx q[2];
rz(-0.96322434) q[2];
sx q[2];
rz(1.4841532) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15785698) q[1];
sx q[1];
rz(-0.4326371) q[1];
sx q[1];
rz(-0.067072596) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.841373) q[3];
sx q[3];
rz(-1.2660053) q[3];
sx q[3];
rz(-2.6255325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60468173) q[2];
sx q[2];
rz(-1.9313507) q[2];
sx q[2];
rz(-0.49986419) q[2];
rz(2.8258421) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(-2.1226814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71378088) q[0];
sx q[0];
rz(-1.918387) q[0];
sx q[0];
rz(2.1977303) q[0];
rz(1.4048514) q[1];
sx q[1];
rz(-1.8753139) q[1];
sx q[1];
rz(0.92165438) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.114678) q[0];
sx q[0];
rz(-1.8887547) q[0];
sx q[0];
rz(-1.5131895) q[0];
rz(-pi) q[1];
rz(-1.2723075) q[2];
sx q[2];
rz(-1.8137852) q[2];
sx q[2];
rz(1.1557494) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4186395) q[1];
sx q[1];
rz(-2.2721842) q[1];
sx q[1];
rz(-0.69461125) q[1];
rz(-1.0579487) q[3];
sx q[3];
rz(-1.5745224) q[3];
sx q[3];
rz(1.8029574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1450119) q[2];
sx q[2];
rz(-1.2715481) q[2];
sx q[2];
rz(1.5550295) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(-1.4012977) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15040511) q[0];
sx q[0];
rz(-1.0973955) q[0];
sx q[0];
rz(-0.64252585) q[0];
rz(1.7036899) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(0.14370758) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3966345) q[0];
sx q[0];
rz(-1.2681343) q[0];
sx q[0];
rz(-2.55654) q[0];
x q[1];
rz(-2.6176378) q[2];
sx q[2];
rz(-0.99159504) q[2];
sx q[2];
rz(-0.29925811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.70348008) q[1];
sx q[1];
rz(-0.92877711) q[1];
sx q[1];
rz(-1.239052) q[1];
x q[2];
rz(-1.9552574) q[3];
sx q[3];
rz(-1.1610306) q[3];
sx q[3];
rz(-2.055197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5370499) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(-2.878888) q[2];
rz(-0.25775868) q[3];
sx q[3];
rz(-0.38193211) q[3];
sx q[3];
rz(0.8530544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(-1.9211796) q[0];
rz(2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-0.89734546) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98678714) q[0];
sx q[0];
rz(-2.6098394) q[0];
sx q[0];
rz(0.24260862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4608896) q[2];
sx q[2];
rz(-1.8898003) q[2];
sx q[2];
rz(-1.6750592) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.022696115) q[1];
sx q[1];
rz(-1.5210946) q[1];
sx q[1];
rz(-1.3124052) q[1];
x q[2];
rz(3.0378689) q[3];
sx q[3];
rz(-0.96451603) q[3];
sx q[3];
rz(-1.0223573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64409488) q[2];
sx q[2];
rz(-2.2150025) q[2];
sx q[2];
rz(-3.0268055) q[2];
rz(0.74350205) q[3];
sx q[3];
rz(-1.8758352) q[3];
sx q[3];
rz(-0.44873294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153618) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(1.505898) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(2.220357) q[2];
sx q[2];
rz(-2.2050646) q[2];
sx q[2];
rz(0.83940432) q[2];
rz(-1.4033477) q[3];
sx q[3];
rz(-1.1421775) q[3];
sx q[3];
rz(0.92489064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
