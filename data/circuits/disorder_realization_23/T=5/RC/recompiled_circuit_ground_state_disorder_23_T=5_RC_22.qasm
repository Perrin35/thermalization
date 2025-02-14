OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5705559) q[0];
sx q[0];
rz(-0.67222995) q[0];
sx q[0];
rz(-1.6663405) q[0];
rz(-2.1885459) q[1];
sx q[1];
rz(-0.16521984) q[1];
sx q[1];
rz(-1.2686977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3606963) q[0];
sx q[0];
rz(-1.7994045) q[0];
sx q[0];
rz(2.7318888) q[0];
x q[1];
rz(-2.2242429) q[2];
sx q[2];
rz(-1.1684061) q[2];
sx q[2];
rz(2.0584597) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7172473) q[1];
sx q[1];
rz(-1.0043036) q[1];
sx q[1];
rz(0.15498181) q[1];
x q[2];
rz(0.24390177) q[3];
sx q[3];
rz(-2.3696405) q[3];
sx q[3];
rz(3.1305135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1573726) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(0.79322195) q[2];
rz(-2.4685517) q[3];
sx q[3];
rz(-1.0854191) q[3];
sx q[3];
rz(-1.2696772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15758812) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(-0.60930914) q[0];
rz(2.9318103) q[1];
sx q[1];
rz(-2.3502626) q[1];
sx q[1];
rz(-2.1817068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5426838) q[0];
sx q[0];
rz(-3.050878) q[0];
sx q[0];
rz(0.51341052) q[0];
x q[1];
rz(2.0591956) q[2];
sx q[2];
rz(-1.5090183) q[2];
sx q[2];
rz(-0.90324963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7916225) q[1];
sx q[1];
rz(-2.4020139) q[1];
sx q[1];
rz(-1.4776358) q[1];
rz(-2.0783362) q[3];
sx q[3];
rz(-1.9572658) q[3];
sx q[3];
rz(-1.2852576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3557055) q[2];
sx q[2];
rz(-0.83928078) q[2];
sx q[2];
rz(2.6500224) q[2];
rz(-0.0056886557) q[3];
sx q[3];
rz(-1.0892884) q[3];
sx q[3];
rz(2.6277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095057644) q[0];
sx q[0];
rz(-2.6061366) q[0];
sx q[0];
rz(0.74263483) q[0];
rz(0.62478089) q[1];
sx q[1];
rz(-2.2014328) q[1];
sx q[1];
rz(-2.0754441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2258573) q[0];
sx q[0];
rz(-1.429378) q[0];
sx q[0];
rz(1.3456324) q[0];
rz(-pi) q[1];
rz(1.891648) q[2];
sx q[2];
rz(-1.3898456) q[2];
sx q[2];
rz(0.89676434) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9650363) q[1];
sx q[1];
rz(-1.4449458) q[1];
sx q[1];
rz(1.1648037) q[1];
x q[2];
rz(2.0622018) q[3];
sx q[3];
rz(-1.764445) q[3];
sx q[3];
rz(-1.0102864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8891958) q[2];
sx q[2];
rz(-0.19813457) q[2];
sx q[2];
rz(-0.55257094) q[2];
rz(0.50734723) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(2.5721917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6602537) q[0];
sx q[0];
rz(-2.5642671) q[0];
sx q[0];
rz(1.8950155) q[0];
rz(1.4056816) q[1];
sx q[1];
rz(-2.3697) q[1];
sx q[1];
rz(2.4235922) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5575422) q[0];
sx q[0];
rz(-2.1831838) q[0];
sx q[0];
rz(0.63300818) q[0];
x q[1];
rz(-1.0077072) q[2];
sx q[2];
rz(-1.3674842) q[2];
sx q[2];
rz(2.2846534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64252428) q[1];
sx q[1];
rz(-0.82083265) q[1];
sx q[1];
rz(1.7700511) q[1];
rz(-1.4425475) q[3];
sx q[3];
rz(-1.413563) q[3];
sx q[3];
rz(0.12928582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2672853) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(2.4370952) q[2];
rz(2.7590175) q[3];
sx q[3];
rz(-1.7035328) q[3];
sx q[3];
rz(0.85324919) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.660897) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(-0.30886343) q[0];
rz(-0.94912306) q[1];
sx q[1];
rz(-2.821065) q[1];
sx q[1];
rz(-0.55108756) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68252045) q[0];
sx q[0];
rz(-0.42051747) q[0];
sx q[0];
rz(-2.7608052) q[0];
rz(-pi) q[1];
rz(-0.068821235) q[2];
sx q[2];
rz(-2.2274979) q[2];
sx q[2];
rz(-2.7842885) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70610395) q[1];
sx q[1];
rz(-1.7167673) q[1];
sx q[1];
rz(-1.414243) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3651028) q[3];
sx q[3];
rz(-1.4675508) q[3];
sx q[3];
rz(-2.6342027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4264195) q[2];
sx q[2];
rz(-2.3931563) q[2];
sx q[2];
rz(0.61881649) q[2];
rz(0.62301451) q[3];
sx q[3];
rz(-0.84191936) q[3];
sx q[3];
rz(-0.4272517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039336786) q[0];
sx q[0];
rz(-2.4781041) q[0];
sx q[0];
rz(-2.6797507) q[0];
rz(0.49679187) q[1];
sx q[1];
rz(-1.8180314) q[1];
sx q[1];
rz(-1.7740446) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753938) q[0];
sx q[0];
rz(-2.4552422) q[0];
sx q[0];
rz(-1.696287) q[0];
x q[1];
rz(-1.8470079) q[2];
sx q[2];
rz(-0.51339692) q[2];
sx q[2];
rz(2.9718176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16713472) q[1];
sx q[1];
rz(-0.30213812) q[1];
sx q[1];
rz(1.4767417) q[1];
rz(-1.7943789) q[3];
sx q[3];
rz(-2.5075965) q[3];
sx q[3];
rz(0.29350212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0355012) q[2];
sx q[2];
rz(-2.0095299) q[2];
sx q[2];
rz(0.79037017) q[2];
rz(0.60449374) q[3];
sx q[3];
rz(-2.1166182) q[3];
sx q[3];
rz(0.42009556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6919493) q[0];
sx q[0];
rz(-0.88570166) q[0];
sx q[0];
rz(2.7272136) q[0];
rz(2.5212133) q[1];
sx q[1];
rz(-1.2216156) q[1];
sx q[1];
rz(1.5588123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010084933) q[0];
sx q[0];
rz(-1.4580618) q[0];
sx q[0];
rz(-1.4860064) q[0];
rz(-pi) q[1];
rz(-1.5887268) q[2];
sx q[2];
rz(-1.6758783) q[2];
sx q[2];
rz(2.004527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0225637) q[1];
sx q[1];
rz(-3.0830586) q[1];
sx q[1];
rz(-1.7332533) q[1];
x q[2];
rz(2.989678) q[3];
sx q[3];
rz(-1.9524116) q[3];
sx q[3];
rz(1.8982062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15535007) q[2];
sx q[2];
rz(-0.68779951) q[2];
sx q[2];
rz(-0.17779329) q[2];
rz(0.47884652) q[3];
sx q[3];
rz(-2.0279453) q[3];
sx q[3];
rz(3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14684045) q[0];
sx q[0];
rz(-0.9786334) q[0];
sx q[0];
rz(-3.0290208) q[0];
rz(-0.62537891) q[1];
sx q[1];
rz(-2.0140779) q[1];
sx q[1];
rz(-2.2523527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.205394) q[0];
sx q[0];
rz(-1.5279211) q[0];
sx q[0];
rz(-1.491719) q[0];
rz(-2.8577789) q[2];
sx q[2];
rz(-2.8655911) q[2];
sx q[2];
rz(-2.2021879) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1247876) q[1];
sx q[1];
rz(-1.8062002) q[1];
sx q[1];
rz(2.1076057) q[1];
x q[2];
rz(-3.0744138) q[3];
sx q[3];
rz(-1.805873) q[3];
sx q[3];
rz(2.377411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.858295) q[2];
sx q[2];
rz(-3.0136643) q[2];
sx q[2];
rz(-0.27321401) q[2];
rz(-2.9122399) q[3];
sx q[3];
rz(-1.554824) q[3];
sx q[3];
rz(1.2685512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44723085) q[0];
sx q[0];
rz(-2.2706967) q[0];
sx q[0];
rz(2.2253775) q[0];
rz(0.18237309) q[1];
sx q[1];
rz(-2.9353751) q[1];
sx q[1];
rz(1.9825541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75814509) q[0];
sx q[0];
rz(-2.8383857) q[0];
sx q[0];
rz(-0.68443735) q[0];
rz(3.0460814) q[2];
sx q[2];
rz(-1.6217124) q[2];
sx q[2];
rz(-1.7820304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9578286) q[1];
sx q[1];
rz(-1.8845585) q[1];
sx q[1];
rz(-0.16739629) q[1];
x q[2];
rz(2.4156225) q[3];
sx q[3];
rz(-1.0683224) q[3];
sx q[3];
rz(2.610725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75659043) q[2];
sx q[2];
rz(-2.341541) q[2];
sx q[2];
rz(0.315256) q[2];
rz(-0.59424019) q[3];
sx q[3];
rz(-2.2599594) q[3];
sx q[3];
rz(2.5074904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2600128) q[0];
sx q[0];
rz(-0.6686815) q[0];
sx q[0];
rz(2.9823533) q[0];
rz(-2.0533994) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(2.870627) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3581287) q[0];
sx q[0];
rz(-2.6081351) q[0];
sx q[0];
rz(1.0607884) q[0];
x q[1];
rz(0.6728716) q[2];
sx q[2];
rz(-1.4904163) q[2];
sx q[2];
rz(-2.5774235) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7587563) q[1];
sx q[1];
rz(-2.3985632) q[1];
sx q[1];
rz(0.50937517) q[1];
rz(-2.5529717) q[3];
sx q[3];
rz(-0.56097066) q[3];
sx q[3];
rz(-2.4952793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8012041) q[2];
sx q[2];
rz(-2.3694254) q[2];
sx q[2];
rz(2.8188952) q[2];
rz(-3.0835551) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(0.29840741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.4651481) q[0];
sx q[0];
rz(-1.6315176) q[0];
sx q[0];
rz(2.3254707) q[0];
rz(0.25794087) q[1];
sx q[1];
rz(-2.0057269) q[1];
sx q[1];
rz(-1.534091) q[1];
rz(2.8526974) q[2];
sx q[2];
rz(-1.8452273) q[2];
sx q[2];
rz(0.065963521) q[2];
rz(-0.54775379) q[3];
sx q[3];
rz(-2.4079101) q[3];
sx q[3];
rz(-1.1338601) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
