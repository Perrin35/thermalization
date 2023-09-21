OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52743071) q[0];
sx q[0];
rz(-2.3318113) q[0];
sx q[0];
rz(0.53139395) q[0];
rz(0.2375138) q[1];
sx q[1];
rz(-1.7778492) q[1];
sx q[1];
rz(-1.2030503) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7145568) q[0];
sx q[0];
rz(-1.3792975) q[0];
sx q[0];
rz(-1.7869851) q[0];
rz(-pi) q[1];
rz(2.8868944) q[2];
sx q[2];
rz(-1.7809976) q[2];
sx q[2];
rz(2.3496698) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4431475) q[1];
sx q[1];
rz(-1.5138211) q[1];
sx q[1];
rz(-2.9064102) q[1];
x q[2];
rz(3.0193127) q[3];
sx q[3];
rz(-0.29502007) q[3];
sx q[3];
rz(-2.2515841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27665859) q[2];
sx q[2];
rz(-1.1316789) q[2];
sx q[2];
rz(2.9425088) q[2];
rz(1.52786) q[3];
sx q[3];
rz(-0.49694967) q[3];
sx q[3];
rz(-2.4075107) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428225) q[0];
sx q[0];
rz(-3.0144189) q[0];
sx q[0];
rz(-2.3221827) q[0];
rz(2.8557414) q[1];
sx q[1];
rz(-0.91393036) q[1];
sx q[1];
rz(1.2664638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0985384) q[0];
sx q[0];
rz(-2.1715954) q[0];
sx q[0];
rz(-1.1061125) q[0];
rz(-pi) q[1];
rz(1.9752713) q[2];
sx q[2];
rz(-1.3504488) q[2];
sx q[2];
rz(-1.599556) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4283735) q[1];
sx q[1];
rz(-1.9698779) q[1];
sx q[1];
rz(-0.39308163) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5655079) q[3];
sx q[3];
rz(-1.1844716) q[3];
sx q[3];
rz(1.9379804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9849898) q[2];
sx q[2];
rz(-1.9282324) q[2];
sx q[2];
rz(1.161969) q[2];
rz(-3.0544288) q[3];
sx q[3];
rz(-1.5036843) q[3];
sx q[3];
rz(-3.0676837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53482985) q[0];
sx q[0];
rz(-1.1207026) q[0];
sx q[0];
rz(-0.30360046) q[0];
rz(-1.7595694) q[1];
sx q[1];
rz(-1.2761812) q[1];
sx q[1];
rz(-0.64750013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.835639) q[0];
sx q[0];
rz(-1.6157955) q[0];
sx q[0];
rz(-0.63458981) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71543773) q[2];
sx q[2];
rz(-1.9999265) q[2];
sx q[2];
rz(-1.5646343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69933575) q[1];
sx q[1];
rz(-1.2067716) q[1];
sx q[1];
rz(-1.0815716) q[1];
x q[2];
rz(2.5472574) q[3];
sx q[3];
rz(-2.1777275) q[3];
sx q[3];
rz(-2.2567574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2993762) q[2];
sx q[2];
rz(-0.78263485) q[2];
sx q[2];
rz(1.5608609) q[2];
rz(-0.98172274) q[3];
sx q[3];
rz(-1.862062) q[3];
sx q[3];
rz(-0.64341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2340387) q[0];
sx q[0];
rz(-2.3038395) q[0];
sx q[0];
rz(-2.2696944) q[0];
rz(-0.38726989) q[1];
sx q[1];
rz(-0.64278066) q[1];
sx q[1];
rz(2.8505468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5440118) q[0];
sx q[0];
rz(-0.62950069) q[0];
sx q[0];
rz(1.9660216) q[0];
rz(0.83424763) q[2];
sx q[2];
rz(-0.4220037) q[2];
sx q[2];
rz(-0.01854245) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9705022) q[1];
sx q[1];
rz(-0.88741747) q[1];
sx q[1];
rz(2.469953) q[1];
rz(-pi) q[2];
rz(-2.1471359) q[3];
sx q[3];
rz(-1.4354224) q[3];
sx q[3];
rz(-1.4801814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7724093) q[2];
sx q[2];
rz(-2.3466551) q[2];
sx q[2];
rz(-2.5975361) q[2];
rz(0.50928515) q[3];
sx q[3];
rz(-1.0638758) q[3];
sx q[3];
rz(-2.2756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13957025) q[0];
sx q[0];
rz(-1.0263356) q[0];
sx q[0];
rz(-0.64506662) q[0];
rz(0.6257261) q[1];
sx q[1];
rz(-1.2236243) q[1];
sx q[1];
rz(-0.63017875) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6117614) q[0];
sx q[0];
rz(-2.1453619) q[0];
sx q[0];
rz(1.7163506) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3216256) q[2];
sx q[2];
rz(-1.8659288) q[2];
sx q[2];
rz(0.68857869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.899257) q[1];
sx q[1];
rz(-1.3471654) q[1];
sx q[1];
rz(1.0244589) q[1];
rz(1.7230677) q[3];
sx q[3];
rz(-1.0183183) q[3];
sx q[3];
rz(0.0020364062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75382918) q[2];
sx q[2];
rz(-1.8176273) q[2];
sx q[2];
rz(1.3641664) q[2];
rz(1.8917313) q[3];
sx q[3];
rz(-1.2156237) q[3];
sx q[3];
rz(0.92145872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0548148) q[0];
sx q[0];
rz(-2.5578816) q[0];
sx q[0];
rz(-2.4247647) q[0];
rz(-0.43771276) q[1];
sx q[1];
rz(-2.4611459) q[1];
sx q[1];
rz(2.3278918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98052927) q[0];
sx q[0];
rz(-1.6335532) q[0];
sx q[0];
rz(1.3570696) q[0];
x q[1];
rz(-2.5021624) q[2];
sx q[2];
rz(-0.8317906) q[2];
sx q[2];
rz(-2.8384428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4390347) q[1];
sx q[1];
rz(-0.82779373) q[1];
sx q[1];
rz(2.031416) q[1];
rz(-pi) q[2];
rz(2.2807062) q[3];
sx q[3];
rz(-0.45300278) q[3];
sx q[3];
rz(-2.3080254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8852691) q[2];
sx q[2];
rz(-0.43748125) q[2];
sx q[2];
rz(1.1876594) q[2];
rz(-0.93196431) q[3];
sx q[3];
rz(-2.0075802) q[3];
sx q[3];
rz(2.7704172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1535783) q[0];
sx q[0];
rz(-2.1126641) q[0];
sx q[0];
rz(-2.4555092) q[0];
rz(-1.5230806) q[1];
sx q[1];
rz(-0.97421092) q[1];
sx q[1];
rz(-0.66326052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406869) q[0];
sx q[0];
rz(-1.5947043) q[0];
sx q[0];
rz(0.92047128) q[0];
x q[1];
rz(-2.8510423) q[2];
sx q[2];
rz(-0.65813375) q[2];
sx q[2];
rz(1.611329) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85147334) q[1];
sx q[1];
rz(-1.0901325) q[1];
sx q[1];
rz(2.9139247) q[1];
x q[2];
rz(2.5797964) q[3];
sx q[3];
rz(-1.1601163) q[3];
sx q[3];
rz(-0.79515776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66701165) q[2];
sx q[2];
rz(-0.8374927) q[2];
sx q[2];
rz(-0.015080301) q[2];
rz(1.0117426) q[3];
sx q[3];
rz(-2.0254617) q[3];
sx q[3];
rz(-0.57141465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886154) q[0];
sx q[0];
rz(-2.4089854) q[0];
sx q[0];
rz(-3.0342039) q[0];
rz(0.30474162) q[1];
sx q[1];
rz(-1.3184897) q[1];
sx q[1];
rz(0.4577786) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7664238) q[0];
sx q[0];
rz(-0.94978118) q[0];
sx q[0];
rz(-2.7428438) q[0];
rz(2.1019943) q[2];
sx q[2];
rz(-2.2441412) q[2];
sx q[2];
rz(-1.4199867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.01955186) q[1];
sx q[1];
rz(-1.8035839) q[1];
sx q[1];
rz(-1.2047393) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5896766) q[3];
sx q[3];
rz(-1.2593927) q[3];
sx q[3];
rz(-2.6375254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36499873) q[2];
sx q[2];
rz(-1.4564617) q[2];
sx q[2];
rz(-1.7101074) q[2];
rz(-1.7163904) q[3];
sx q[3];
rz(-0.6267572) q[3];
sx q[3];
rz(-1.2861929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872221) q[0];
sx q[0];
rz(-0.4796589) q[0];
sx q[0];
rz(-2.1881058) q[0];
rz(1.8158688) q[1];
sx q[1];
rz(-0.49294254) q[1];
sx q[1];
rz(0.58473933) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3922694) q[0];
sx q[0];
rz(-1.4587147) q[0];
sx q[0];
rz(2.6507676) q[0];
rz(-2.598001) q[2];
sx q[2];
rz(-1.2541176) q[2];
sx q[2];
rz(-1.8900332) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4084088) q[1];
sx q[1];
rz(-1.5015263) q[1];
sx q[1];
rz(0.44002156) q[1];
rz(1.0869914) q[3];
sx q[3];
rz(-0.30011794) q[3];
sx q[3];
rz(-0.3008315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3328302) q[2];
sx q[2];
rz(-2.3213883) q[2];
sx q[2];
rz(2.8653223) q[2];
rz(-2.590498) q[3];
sx q[3];
rz(-1.7421744) q[3];
sx q[3];
rz(0.38366693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.0693414) q[0];
sx q[0];
rz(-0.51420099) q[0];
sx q[0];
rz(2.4861091) q[0];
rz(1.1570702) q[1];
sx q[1];
rz(-1.3815222) q[1];
sx q[1];
rz(1.6190593) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9690543) q[0];
sx q[0];
rz(-1.2612871) q[0];
sx q[0];
rz(2.4597416) q[0];
x q[1];
rz(-1.2298146) q[2];
sx q[2];
rz(-0.91170646) q[2];
sx q[2];
rz(1.1863208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5628221) q[1];
sx q[1];
rz(-0.65014231) q[1];
sx q[1];
rz(-0.27362089) q[1];
x q[2];
rz(-0.75928648) q[3];
sx q[3];
rz(-2.5873103) q[3];
sx q[3];
rz(0.62461531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9860501) q[2];
sx q[2];
rz(-0.4077929) q[2];
sx q[2];
rz(2.299451) q[2];
rz(1.5367674) q[3];
sx q[3];
rz(-1.3804881) q[3];
sx q[3];
rz(3.1150637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355272) q[0];
sx q[0];
rz(-1.1885831) q[0];
sx q[0];
rz(2.6901235) q[0];
rz(-2.4189667) q[1];
sx q[1];
rz(-0.87651785) q[1];
sx q[1];
rz(-1.6095907) q[1];
rz(-0.026635219) q[2];
sx q[2];
rz(-1.7799108) q[2];
sx q[2];
rz(-1.5540661) q[2];
rz(1.6668672) q[3];
sx q[3];
rz(-2.4484652) q[3];
sx q[3];
rz(-0.0948003) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];