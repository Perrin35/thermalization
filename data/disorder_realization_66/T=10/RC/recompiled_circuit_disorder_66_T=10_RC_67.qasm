OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(3.0193168) q[1];
sx q[1];
rz(-3.0552157) q[1];
sx q[1];
rz(3.1187305) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6883144) q[0];
sx q[0];
rz(-2.9524321) q[0];
sx q[0];
rz(2.2384089) q[0];
x q[1];
rz(-0.68732287) q[2];
sx q[2];
rz(-0.63185531) q[2];
sx q[2];
rz(2.1473715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7658246) q[1];
sx q[1];
rz(-1.4210912) q[1];
sx q[1];
rz(1.9181812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1894737) q[3];
sx q[3];
rz(-1.3631571) q[3];
sx q[3];
rz(1.729897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7227398) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(0.76618761) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-2.4092716) q[3];
sx q[3];
rz(0.58656251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3294753) q[0];
sx q[0];
rz(-2.7024039) q[0];
sx q[0];
rz(-0.086659327) q[0];
rz(0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-3.0564953) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3413977) q[0];
sx q[0];
rz(-1.5543823) q[0];
sx q[0];
rz(0.62231681) q[0];
rz(-pi) q[1];
rz(-2.000196) q[2];
sx q[2];
rz(-0.64711249) q[2];
sx q[2];
rz(-0.61691689) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9566006) q[1];
sx q[1];
rz(-1.6829832) q[1];
sx q[1];
rz(-3.0549906) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3195023) q[3];
sx q[3];
rz(-2.0204244) q[3];
sx q[3];
rz(-0.055469661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.6764199) q[2];
rz(1.1761459) q[3];
sx q[3];
rz(-0.90548235) q[3];
sx q[3];
rz(-0.24648497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48433205) q[0];
sx q[0];
rz(-3.0069139) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-0.33272818) q[1];
sx q[1];
rz(-2.2867124) q[1];
sx q[1];
rz(1.7571626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8850088) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(2.9628739) q[0];
rz(-0.39321123) q[2];
sx q[2];
rz(-2.0438497) q[2];
sx q[2];
rz(2.2676603) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3570909) q[1];
sx q[1];
rz(-2.38584) q[1];
sx q[1];
rz(0.032553629) q[1];
rz(2.1372041) q[3];
sx q[3];
rz(-2.0319788) q[3];
sx q[3];
rz(3.0571292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.3282233) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(-0.11169294) q[3];
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
rz(3.0760647) q[0];
sx q[0];
rz(-0.17834839) q[0];
sx q[0];
rz(-2.4705825) q[0];
rz(1.3440075) q[1];
sx q[1];
rz(-1.1616511) q[1];
sx q[1];
rz(2.9342594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610791) q[0];
sx q[0];
rz(-1.5797537) q[0];
sx q[0];
rz(-1.6121959) q[0];
x q[1];
rz(1.7423332) q[2];
sx q[2];
rz(-2.9651387) q[2];
sx q[2];
rz(1.8432957) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15370788) q[1];
sx q[1];
rz(-0.28370902) q[1];
sx q[1];
rz(2.1464286) q[1];
rz(-pi) q[2];
rz(-3.1158833) q[3];
sx q[3];
rz(-0.51179143) q[3];
sx q[3];
rz(-2.5185086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31071445) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-2.3332398) q[2];
rz(-2.3338142) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5508674) q[0];
sx q[0];
rz(-2.4243675) q[0];
sx q[0];
rz(-0.89548683) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-0.83509713) q[1];
sx q[1];
rz(2.6079544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1874439) q[0];
sx q[0];
rz(-1.8549518) q[0];
sx q[0];
rz(-0.28676333) q[0];
x q[1];
rz(-2.3737714) q[2];
sx q[2];
rz(-1.4066833) q[2];
sx q[2];
rz(0.62190157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4644949) q[1];
sx q[1];
rz(-0.2650731) q[1];
sx q[1];
rz(-2.5527843) q[1];
rz(-0.20547262) q[3];
sx q[3];
rz(-2.1901202) q[3];
sx q[3];
rz(1.6831786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2395893) q[2];
sx q[2];
rz(-1.5633554) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(1.2166294) q[3];
sx q[3];
rz(-2.670848) q[3];
sx q[3];
rz(2.7980428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-0.87843043) q[0];
sx q[0];
rz(1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-2.1891179) q[1];
sx q[1];
rz(-2.1369381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5414808) q[0];
sx q[0];
rz(-1.913537) q[0];
sx q[0];
rz(2.747885) q[0];
x q[1];
rz(2.9506358) q[2];
sx q[2];
rz(-2.7381884) q[2];
sx q[2];
rz(1.0952589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1130044) q[1];
sx q[1];
rz(-2.5341946) q[1];
sx q[1];
rz(2.7007156) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.011990487) q[3];
sx q[3];
rz(-1.1569835) q[3];
sx q[3];
rz(-3.0179265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(0.44815865) q[2];
rz(-2.8435977) q[3];
sx q[3];
rz(-2.4500676) q[3];
sx q[3];
rz(-2.7929849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77666831) q[0];
sx q[0];
rz(-1.3013327) q[0];
sx q[0];
rz(2.4834494) q[0];
rz(1.2843885) q[1];
sx q[1];
rz(-0.45416608) q[1];
sx q[1];
rz(2.4694494) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9854239) q[0];
sx q[0];
rz(-0.16616136) q[0];
sx q[0];
rz(-1.6155433) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1291814) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(2.0814975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0587412) q[1];
sx q[1];
rz(-0.62742678) q[1];
sx q[1];
rz(2.4861366) q[1];
rz(-pi) q[2];
rz(-0.16996202) q[3];
sx q[3];
rz(-2.4331577) q[3];
sx q[3];
rz(-1.2989312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60823524) q[2];
sx q[2];
rz(-2.50864) q[2];
sx q[2];
rz(2.7427924) q[2];
rz(0.82459015) q[3];
sx q[3];
rz(-1.7493068) q[3];
sx q[3];
rz(-2.5430172) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010715) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(0.15821247) q[0];
rz(-2.054706) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-0.80668443) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54803941) q[0];
sx q[0];
rz(-1.1386477) q[0];
sx q[0];
rz(3.0120864) q[0];
rz(2.6629692) q[2];
sx q[2];
rz(-2.9106986) q[2];
sx q[2];
rz(-0.17639562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1572691) q[1];
sx q[1];
rz(-1.2829797) q[1];
sx q[1];
rz(2.4411574) q[1];
rz(1.673647) q[3];
sx q[3];
rz(-0.30234435) q[3];
sx q[3];
rz(-0.24149382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5640101) q[2];
sx q[2];
rz(-1.6489886) q[2];
sx q[2];
rz(2.1419443) q[2];
rz(3.0380761) q[3];
sx q[3];
rz(-3.0045356) q[3];
sx q[3];
rz(-1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0209811) q[0];
sx q[0];
rz(-0.7779026) q[0];
sx q[0];
rz(-0.65752423) q[0];
rz(0.28655562) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(-2.8009169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176312) q[0];
sx q[0];
rz(-2.258856) q[0];
sx q[0];
rz(0.4171564) q[0];
x q[1];
rz(2.1358143) q[2];
sx q[2];
rz(-1.7592332) q[2];
sx q[2];
rz(-0.68067683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9790736) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(1.2390922) q[1];
rz(-pi) q[2];
rz(2.4386028) q[3];
sx q[3];
rz(-2.8841416) q[3];
sx q[3];
rz(-1.6182181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75491536) q[2];
sx q[2];
rz(-1.200054) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(-0.30161101) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(0.51914674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122445) q[0];
sx q[0];
rz(-3.0817741) q[0];
sx q[0];
rz(-2.8630032) q[0];
rz(0.57922286) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(0.07671193) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68966507) q[0];
sx q[0];
rz(-1.8849012) q[0];
sx q[0];
rz(0.12977022) q[0];
rz(-pi) q[1];
rz(1.8329343) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(-2.4208456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6381166) q[1];
sx q[1];
rz(-1.5917115) q[1];
sx q[1];
rz(2.7760844) q[1];
rz(-pi) q[2];
rz(-0.46080188) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(0.27289665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-0.45551604) q[2];
rz(-0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-0.063522696) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(-2.5554399) q[1];
sx q[1];
rz(-1.502232) q[1];
sx q[1];
rz(-1.4486817) q[1];
rz(1.099844) q[2];
sx q[2];
rz(-1.3146613) q[2];
sx q[2];
rz(1.4801499) q[2];
rz(-1.205411) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
