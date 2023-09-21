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
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9553298) q[0];
sx q[0];
rz(-1.7832527) q[0];
sx q[0];
rz(-2.0583378) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9287455) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(-1.1320621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1311156) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(2.1285776) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8817301) q[3];
sx q[3];
rz(-1.6136323) q[3];
sx q[3];
rz(-0.44997893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(1.1928605) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002999) q[0];
sx q[0];
rz(-2.2465758) q[0];
sx q[0];
rz(2.7594901) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58354124) q[2];
sx q[2];
rz(-1.1569287) q[2];
sx q[2];
rz(-1.5577424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0484867) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(2.6622245) q[1];
rz(-pi) q[2];
rz(1.1851951) q[3];
sx q[3];
rz(-2.0413627) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(2.9197664) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(3.085014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3765592) q[0];
sx q[0];
rz(-2.3925836) q[0];
sx q[0];
rz(0.88699938) q[0];
rz(-2.8851606) q[2];
sx q[2];
rz(-1.4415381) q[2];
sx q[2];
rz(-1.3916707) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0128855) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(0.30602869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4967381) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(-2.9063318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(2.3994989) q[0];
rz(-1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(0.46359584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0915506) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(-2.0043623) q[0];
x q[1];
rz(2.8088403) q[2];
sx q[2];
rz(-1.6289662) q[2];
sx q[2];
rz(1.0156877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78471781) q[1];
sx q[1];
rz(-2.3349635) q[1];
sx q[1];
rz(-0.23826092) q[1];
rz(-pi) q[2];
rz(-2.1257524) q[3];
sx q[3];
rz(-1.98588) q[3];
sx q[3];
rz(-2.2500028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(2.1972426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5807242) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(0.046037721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40446754) q[2];
sx q[2];
rz(-2.3239115) q[2];
sx q[2];
rz(0.58194619) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3790834) q[1];
sx q[1];
rz(-3.0804539) q[1];
sx q[1];
rz(1.5361384) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3145507) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(0.30403578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8828316) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-0.32430696) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(3.086673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024825) q[0];
sx q[0];
rz(-1.6533028) q[0];
sx q[0];
rz(1.8399747) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.9565441) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(-0.59021414) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6782126) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(0.86138553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(0.28453919) q[0];
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
rz(1.1658265) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(-1.635701) q[0];
x q[1];
rz(-1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(-1.8679801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4714204) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(2.7736204) q[1];
rz(0.72083731) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(-2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(-2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(-2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(-2.8410889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(-2.6595594) q[0];
x q[1];
rz(2.1883165) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(-0.88027871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1860736) q[1];
sx q[1];
rz(-1.7454073) q[1];
sx q[1];
rz(0.34592918) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31605966) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(-2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(2.3596181) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-0.75884563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1369143) q[0];
sx q[0];
rz(-1.5537964) q[0];
sx q[0];
rz(-1.9949811) q[0];
x q[1];
rz(0.27300948) q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(0.89275937) q[1];
sx q[1];
rz(-1.4151238) q[1];
sx q[1];
rz(-0.70181904) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5522478) q[3];
sx q[3];
rz(-0.66569257) q[3];
sx q[3];
rz(-1.1902283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(1.7193433) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(2.5316701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9655351) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(2.7479991) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62060771) q[2];
sx q[2];
rz(-2.6891516) q[2];
sx q[2];
rz(-1.1021745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48253179) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-2.798583) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(1.324211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.7667608) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(0.25272947) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
