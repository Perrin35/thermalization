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
rz(0.63737386) q[0];
sx q[0];
rz(-2.0847991) q[0];
sx q[0];
rz(2.3699397) q[0];
rz(2.9906988) q[1];
sx q[1];
rz(-2.8291191) q[1];
sx q[1];
rz(1.6627275) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85017289) q[0];
sx q[0];
rz(-1.8628736) q[0];
sx q[0];
rz(2.7946081) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0885299) q[2];
sx q[2];
rz(-1.3446756) q[2];
sx q[2];
rz(-2.256611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2661765) q[1];
sx q[1];
rz(-2.3485942) q[1];
sx q[1];
rz(-2.7686053) q[1];
rz(-pi) q[2];
rz(-2.9247901) q[3];
sx q[3];
rz(-2.4507782) q[3];
sx q[3];
rz(1.8280713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5375157) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(0.83211952) q[2];
rz(-0.65699792) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(0.35201296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326139) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(0.64999181) q[0];
rz(0.12282148) q[1];
sx q[1];
rz(-2.717412) q[1];
sx q[1];
rz(-2.4688683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8685659) q[0];
sx q[0];
rz(-1.5590686) q[0];
sx q[0];
rz(1.5290029) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71750516) q[2];
sx q[2];
rz(-2.7360672) q[2];
sx q[2];
rz(-0.44446352) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0699464) q[1];
sx q[1];
rz(-0.59161797) q[1];
sx q[1];
rz(1.2143192) q[1];
x q[2];
rz(0.026521369) q[3];
sx q[3];
rz(-1.0511729) q[3];
sx q[3];
rz(-0.61557239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1136721) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(-0.38635722) q[2];
rz(1.1682642) q[3];
sx q[3];
rz(-1.9100274) q[3];
sx q[3];
rz(1.108235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89762178) q[0];
sx q[0];
rz(-1.6672927) q[0];
sx q[0];
rz(-3.0834055) q[0];
rz(0.30481401) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(2.286639) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2236587) q[0];
sx q[0];
rz(-3.0732647) q[0];
sx q[0];
rz(3.1177417) q[0];
x q[1];
rz(-2.0570175) q[2];
sx q[2];
rz(-2.2006319) q[2];
sx q[2];
rz(1.4818918) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.952576) q[1];
sx q[1];
rz(-1.7754648) q[1];
sx q[1];
rz(-2.3701028) q[1];
rz(-1.1596572) q[3];
sx q[3];
rz(-2.5469031) q[3];
sx q[3];
rz(0.88319544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.052224) q[2];
sx q[2];
rz(-1.8171909) q[2];
sx q[2];
rz(-1.8734056) q[2];
rz(-2.7005699) q[3];
sx q[3];
rz(-2.5206168) q[3];
sx q[3];
rz(0.84180251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12255254) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(2.4265491) q[0];
rz(2.7252281) q[1];
sx q[1];
rz(-0.48960296) q[1];
sx q[1];
rz(-0.29410902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9884777) q[0];
sx q[0];
rz(-1.5128969) q[0];
sx q[0];
rz(-0.87429509) q[0];
rz(3.0766671) q[2];
sx q[2];
rz(-1.366642) q[2];
sx q[2];
rz(0.72694187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8079559) q[1];
sx q[1];
rz(-0.2273493) q[1];
sx q[1];
rz(-0.61509404) q[1];
x q[2];
rz(2.3669621) q[3];
sx q[3];
rz(-2.3916158) q[3];
sx q[3];
rz(-0.74912723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2405582) q[2];
sx q[2];
rz(-1.7419387) q[2];
sx q[2];
rz(-0.72875363) q[2];
rz(0.48480222) q[3];
sx q[3];
rz(-2.1383643) q[3];
sx q[3];
rz(-2.9639967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6840376) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(-2.8628602) q[0];
rz(3.0740671) q[1];
sx q[1];
rz(-0.7154671) q[1];
sx q[1];
rz(-0.50594893) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8794969) q[0];
sx q[0];
rz(-2.580632) q[0];
sx q[0];
rz(0.13763388) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9437482) q[2];
sx q[2];
rz(-1.4711507) q[2];
sx q[2];
rz(1.0464335) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3619574) q[1];
sx q[1];
rz(-2.07603) q[1];
sx q[1];
rz(-2.510406) q[1];
x q[2];
rz(-1.2413195) q[3];
sx q[3];
rz(-0.70036026) q[3];
sx q[3];
rz(0.86802378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94023306) q[2];
sx q[2];
rz(-1.4505922) q[2];
sx q[2];
rz(-0.10227164) q[2];
rz(-2.8376132) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(0.48039082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893696) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(-2.4612259) q[0];
rz(0.96087372) q[1];
sx q[1];
rz(-1.5870794) q[1];
sx q[1];
rz(2.5804677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.965428) q[0];
sx q[0];
rz(-2.7211004) q[0];
sx q[0];
rz(-2.0916558) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20724736) q[2];
sx q[2];
rz(-0.85734136) q[2];
sx q[2];
rz(1.1175454) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61455124) q[1];
sx q[1];
rz(-0.30362645) q[1];
sx q[1];
rz(0.18358453) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4986491) q[3];
sx q[3];
rz(-2.3920632) q[3];
sx q[3];
rz(1.3931605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23201021) q[2];
sx q[2];
rz(-2.1640615) q[2];
sx q[2];
rz(-1.5634465) q[2];
rz(1.0163418) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(-2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4119754) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(2.6475661) q[0];
rz(1.5644851) q[1];
sx q[1];
rz(-0.99948519) q[1];
sx q[1];
rz(0.090242537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099797332) q[0];
sx q[0];
rz(-1.3372412) q[0];
sx q[0];
rz(-0.44333027) q[0];
rz(-2.162287) q[2];
sx q[2];
rz(-0.78787178) q[2];
sx q[2];
rz(2.5798714) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.75587326) q[1];
sx q[1];
rz(-1.8579646) q[1];
sx q[1];
rz(-1.1509291) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19321038) q[3];
sx q[3];
rz(-1.8487559) q[3];
sx q[3];
rz(-1.6707458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6287441) q[2];
sx q[2];
rz(-0.57400846) q[2];
sx q[2];
rz(-0.64344704) q[2];
rz(0.63055864) q[3];
sx q[3];
rz(-2.3876811) q[3];
sx q[3];
rz(3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028320463) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(-1.2233618) q[0];
rz(-2.2471097) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(-1.1462513) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664331) q[0];
sx q[0];
rz(-2.1210175) q[0];
sx q[0];
rz(2.6653637) q[0];
rz(-1.9983534) q[2];
sx q[2];
rz(-1.215286) q[2];
sx q[2];
rz(1.8380788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24487296) q[1];
sx q[1];
rz(-1.280652) q[1];
sx q[1];
rz(1.3194607) q[1];
x q[2];
rz(0.98553879) q[3];
sx q[3];
rz(-1.3511336) q[3];
sx q[3];
rz(2.1287763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2724096) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(1.6669081) q[2];
rz(2.9278751) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(-2.2739482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992582) q[0];
sx q[0];
rz(-2.5304351) q[0];
sx q[0];
rz(-2.4364731) q[0];
rz(2.5662388) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(2.5720678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25999674) q[0];
sx q[0];
rz(-1.651913) q[0];
sx q[0];
rz(-3.0641969) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51409431) q[2];
sx q[2];
rz(-2.5894458) q[2];
sx q[2];
rz(-0.29468003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25334206) q[1];
sx q[1];
rz(-2.1565795) q[1];
sx q[1];
rz(2.3164877) q[1];
rz(-pi) q[2];
rz(1.873718) q[3];
sx q[3];
rz(-2.0407131) q[3];
sx q[3];
rz(-0.63602122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9465955) q[2];
sx q[2];
rz(-2.1985998) q[2];
sx q[2];
rz(2.1627964) q[2];
rz(2.7106674) q[3];
sx q[3];
rz(-2.6543591) q[3];
sx q[3];
rz(-2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8579213) q[0];
sx q[0];
rz(-0.019275276) q[0];
sx q[0];
rz(-1.4625782) q[0];
rz(1.0390394) q[1];
sx q[1];
rz(-1.3835399) q[1];
sx q[1];
rz(1.9336112) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6504575) q[0];
sx q[0];
rz(-1.2692599) q[0];
sx q[0];
rz(1.1745321) q[0];
rz(3.1080148) q[2];
sx q[2];
rz(-1.0462282) q[2];
sx q[2];
rz(-1.7941504) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9644248) q[1];
sx q[1];
rz(-0.96098122) q[1];
sx q[1];
rz(2.9637106) q[1];
x q[2];
rz(-1.7321587) q[3];
sx q[3];
rz(-2.7415025) q[3];
sx q[3];
rz(1.0795979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.155948) q[2];
sx q[2];
rz(-1.0040823) q[2];
sx q[2];
rz(1.3204302) q[2];
rz(-0.350658) q[3];
sx q[3];
rz(-2.3242293) q[3];
sx q[3];
rz(-0.58026522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88631267) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(1.6571922) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(0.37921956) q[2];
sx q[2];
rz(-1.5046962) q[2];
sx q[2];
rz(-1.6284158) q[2];
rz(-2.8410925) q[3];
sx q[3];
rz(-0.90682744) q[3];
sx q[3];
rz(-1.711267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
