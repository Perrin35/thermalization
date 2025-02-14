OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9353256) q[0];
sx q[0];
rz(4.7630035) q[0];
sx q[0];
rz(8.0061316) q[0];
rz(0.043579276) q[1];
sx q[1];
rz(-1.2830623) q[1];
sx q[1];
rz(-2.9820774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.537764) q[0];
sx q[0];
rz(-1.7366752) q[0];
sx q[0];
rz(3.1127451) q[0];
rz(2.5297727) q[2];
sx q[2];
rz(-1.3587591) q[2];
sx q[2];
rz(0.093487948) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7456513) q[1];
sx q[1];
rz(-2.6439754) q[1];
sx q[1];
rz(-1.2785589) q[1];
rz(-pi) q[2];
rz(-1.1662935) q[3];
sx q[3];
rz(-2.3784436) q[3];
sx q[3];
rz(-0.4296225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45304766) q[2];
sx q[2];
rz(-1.1587554) q[2];
sx q[2];
rz(-0.090252074) q[2];
rz(-0.074617535) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(0.41372764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2523044) q[0];
sx q[0];
rz(-1.9828718) q[0];
sx q[0];
rz(-1.4537551) q[0];
rz(0.74268913) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(0.76505605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8370253) q[0];
sx q[0];
rz(-1.6374987) q[0];
sx q[0];
rz(-2.5480258) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80185955) q[2];
sx q[2];
rz(-0.86881402) q[2];
sx q[2];
rz(-0.02846708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38854996) q[1];
sx q[1];
rz(-1.8990128) q[1];
sx q[1];
rz(-0.80372827) q[1];
rz(2.5935947) q[3];
sx q[3];
rz(-1.5685076) q[3];
sx q[3];
rz(-1.3399923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.226977) q[2];
sx q[2];
rz(-2.1299629) q[2];
sx q[2];
rz(1.4545308) q[2];
rz(-2.4948273) q[3];
sx q[3];
rz(-2.5819467) q[3];
sx q[3];
rz(1.2796992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890559) q[0];
sx q[0];
rz(-1.6707957) q[0];
sx q[0];
rz(0.045150969) q[0];
rz(1.2308944) q[1];
sx q[1];
rz(-1.3094333) q[1];
sx q[1];
rz(-2.0848354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8674744) q[0];
sx q[0];
rz(-1.8961268) q[0];
sx q[0];
rz(-2.1528457) q[0];
x q[1];
rz(2.4963794) q[2];
sx q[2];
rz(-1.1629379) q[2];
sx q[2];
rz(-2.4347664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31683982) q[1];
sx q[1];
rz(-2.415363) q[1];
sx q[1];
rz(-1.9588425) q[1];
rz(-1.2892492) q[3];
sx q[3];
rz(-0.70812884) q[3];
sx q[3];
rz(0.81048465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7354273) q[2];
sx q[2];
rz(-0.96144095) q[2];
sx q[2];
rz(2.4889634) q[2];
rz(-1.8485707) q[3];
sx q[3];
rz(-2.3725464) q[3];
sx q[3];
rz(3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3339313) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(1.7474784) q[0];
rz(1.3035424) q[1];
sx q[1];
rz(-1.1640254) q[1];
sx q[1];
rz(1.0882783) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1028839) q[0];
sx q[0];
rz(-0.81224487) q[0];
sx q[0];
rz(1.1884304) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59594229) q[2];
sx q[2];
rz(-2.1998027) q[2];
sx q[2];
rz(0.12550917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66341463) q[1];
sx q[1];
rz(-0.82163787) q[1];
sx q[1];
rz(2.4577702) q[1];
rz(-pi) q[2];
rz(-0.49934629) q[3];
sx q[3];
rz(-1.5245228) q[3];
sx q[3];
rz(-0.83016992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7465705) q[2];
sx q[2];
rz(-1.3777379) q[2];
sx q[2];
rz(0.48040381) q[2];
rz(-2.8754442) q[3];
sx q[3];
rz(-1.1689309) q[3];
sx q[3];
rz(-0.84421414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10923037) q[0];
sx q[0];
rz(-2.5205595) q[0];
sx q[0];
rz(-1.9367628) q[0];
rz(-0.042639848) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(0.94243324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401817) q[0];
sx q[0];
rz(-1.4352006) q[0];
sx q[0];
rz(1.1419161) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3187014) q[2];
sx q[2];
rz(-0.79051711) q[2];
sx q[2];
rz(2.2728777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.98676315) q[1];
sx q[1];
rz(-2.0663321) q[1];
sx q[1];
rz(-0.43191989) q[1];
rz(-pi) q[2];
rz(1.4616555) q[3];
sx q[3];
rz(-0.93459826) q[3];
sx q[3];
rz(0.075625758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11613906) q[2];
sx q[2];
rz(-2.4567273) q[2];
sx q[2];
rz(1.9174513) q[2];
rz(-0.72743607) q[3];
sx q[3];
rz(-1.6614611) q[3];
sx q[3];
rz(3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3392451) q[0];
sx q[0];
rz(-0.50528637) q[0];
sx q[0];
rz(0.25830609) q[0];
rz(0.91336617) q[1];
sx q[1];
rz(-0.86562997) q[1];
sx q[1];
rz(-2.1736653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0570404) q[0];
sx q[0];
rz(-1.0121317) q[0];
sx q[0];
rz(-1.1112122) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36138205) q[2];
sx q[2];
rz(-2.2602096) q[2];
sx q[2];
rz(2.9887082) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63686759) q[1];
sx q[1];
rz(-0.47885413) q[1];
sx q[1];
rz(0.07696159) q[1];
rz(-pi) q[2];
rz(1.0276035) q[3];
sx q[3];
rz(-2.1434727) q[3];
sx q[3];
rz(3.0610994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6496868) q[2];
sx q[2];
rz(-0.68444362) q[2];
sx q[2];
rz(-2.7223132) q[2];
rz(2.3573917) q[3];
sx q[3];
rz(-1.0201642) q[3];
sx q[3];
rz(1.3721344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.075031) q[0];
sx q[0];
rz(-0.12396585) q[0];
sx q[0];
rz(-3.1228464) q[0];
rz(-3.0491414) q[1];
sx q[1];
rz(-1.3321313) q[1];
sx q[1];
rz(-0.094955347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0844011) q[0];
sx q[0];
rz(-0.6633299) q[0];
sx q[0];
rz(-2.6003772) q[0];
x q[1];
rz(-0.59573458) q[2];
sx q[2];
rz(-1.2168665) q[2];
sx q[2];
rz(-0.46094337) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.590789) q[1];
sx q[1];
rz(-1.1491286) q[1];
sx q[1];
rz(-1.8487537) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9390887) q[3];
sx q[3];
rz(-0.42516685) q[3];
sx q[3];
rz(1.396338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2350601) q[2];
sx q[2];
rz(-1.9573213) q[2];
sx q[2];
rz(-0.85901421) q[2];
rz(2.5281995) q[3];
sx q[3];
rz(-2.9412013) q[3];
sx q[3];
rz(-1.8092417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61619806) q[0];
sx q[0];
rz(-1.8680251) q[0];
sx q[0];
rz(-1.4816351) q[0];
rz(1.067465) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(1.8225089) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1265433) q[0];
sx q[0];
rz(-1.476822) q[0];
sx q[0];
rz(0.90576147) q[0];
rz(-pi) q[1];
rz(0.51609938) q[2];
sx q[2];
rz(-1.8987499) q[2];
sx q[2];
rz(-0.60521561) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47506491) q[1];
sx q[1];
rz(-1.804978) q[1];
sx q[1];
rz(-0.98824309) q[1];
x q[2];
rz(2.8105763) q[3];
sx q[3];
rz(-0.3754456) q[3];
sx q[3];
rz(-2.7956243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8760425) q[2];
sx q[2];
rz(-0.37028131) q[2];
sx q[2];
rz(-0.041157095) q[2];
rz(-1.1228849) q[3];
sx q[3];
rz(-1.6582146) q[3];
sx q[3];
rz(-0.52367228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48162833) q[0];
sx q[0];
rz(-2.1631503) q[0];
sx q[0];
rz(1.6690669) q[0];
rz(2.4166079) q[1];
sx q[1];
rz(-1.1712733) q[1];
sx q[1];
rz(-0.57458007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32550463) q[0];
sx q[0];
rz(-2.8003152) q[0];
sx q[0];
rz(-0.14320548) q[0];
rz(0.2821183) q[2];
sx q[2];
rz(-1.9777672) q[2];
sx q[2];
rz(-0.59772432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2255937) q[1];
sx q[1];
rz(-2.288842) q[1];
sx q[1];
rz(-1.5313994) q[1];
rz(-pi) q[2];
rz(0.88317849) q[3];
sx q[3];
rz(-1.6856582) q[3];
sx q[3];
rz(0.77199793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64728111) q[2];
sx q[2];
rz(-2.8936671) q[2];
sx q[2];
rz(-0.9606804) q[2];
rz(-2.658355) q[3];
sx q[3];
rz(-1.5944696) q[3];
sx q[3];
rz(-1.5620935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0110589) q[0];
sx q[0];
rz(-2.647825) q[0];
sx q[0];
rz(0.44179398) q[0];
rz(0.68253016) q[1];
sx q[1];
rz(-1.4492757) q[1];
sx q[1];
rz(-2.0880879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5614612) q[0];
sx q[0];
rz(-0.86505689) q[0];
sx q[0];
rz(-2.8012026) q[0];
x q[1];
rz(1.9437213) q[2];
sx q[2];
rz(-1.4879543) q[2];
sx q[2];
rz(-1.6380472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23020959) q[1];
sx q[1];
rz(-2.1911096) q[1];
sx q[1];
rz(2.7432454) q[1];
x q[2];
rz(0.17591646) q[3];
sx q[3];
rz(-0.39947594) q[3];
sx q[3];
rz(-2.7586058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0423476) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(0.22259268) q[2];
rz(-1.3709204) q[3];
sx q[3];
rz(-1.7226847) q[3];
sx q[3];
rz(0.62209779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2400146) q[0];
sx q[0];
rz(-2.1299025) q[0];
sx q[0];
rz(1.0255751) q[0];
rz(-1.9758132) q[1];
sx q[1];
rz(-0.60973254) q[1];
sx q[1];
rz(2.9235074) q[1];
rz(-0.54033242) q[2];
sx q[2];
rz(-1.9987543) q[2];
sx q[2];
rz(1.0812154) q[2];
rz(1.8924186) q[3];
sx q[3];
rz(-1.5703441) q[3];
sx q[3];
rz(2.3342568) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
