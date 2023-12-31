OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(2.4581576) q[0];
sx q[0];
rz(12.0876) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(-1.1614769) q[1];
sx q[1];
rz(-2.5002313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37748572) q[2];
sx q[2];
rz(-1.1698206) q[2];
sx q[2];
rz(0.12667835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7176425) q[1];
sx q[1];
rz(-2.2279734) q[1];
sx q[1];
rz(-0.98492019) q[1];
x q[2];
rz(-1.1629421) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(2.4401963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59149867) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(1.452662) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801382) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(2.1226728) q[0];
rz(1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(-2.5085124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7019254) q[0];
sx q[0];
rz(-2.3283726) q[0];
sx q[0];
rz(1.6815835) q[0];
x q[1];
rz(-0.622153) q[2];
sx q[2];
rz(-0.85217798) q[2];
sx q[2];
rz(-0.46301401) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1770597) q[1];
sx q[1];
rz(-2.6033083) q[1];
sx q[1];
rz(2.2798377) q[1];
rz(1.9485336) q[3];
sx q[3];
rz(-2.4240498) q[3];
sx q[3];
rz(0.46061463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35976609) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(3.0241372) q[2];
rz(2.84058) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(-1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2837219) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(-0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(0.69082469) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3120964) q[0];
sx q[0];
rz(-1.4179686) q[0];
sx q[0];
rz(1.494708) q[0];
rz(-pi) q[1];
rz(1.3356528) q[2];
sx q[2];
rz(-2.3217839) q[2];
sx q[2];
rz(0.90607925) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9835811) q[1];
sx q[1];
rz(-0.18392715) q[1];
sx q[1];
rz(-1.4278825) q[1];
rz(-pi) q[2];
rz(1.2191804) q[3];
sx q[3];
rz(-1.2120486) q[3];
sx q[3];
rz(-0.84695942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92418015) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(3.0351191) q[2];
rz(1.4364093) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(0.52247125) q[0];
rz(0.32896313) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(2.5879588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5579917) q[0];
sx q[0];
rz(-0.86692536) q[0];
sx q[0];
rz(2.5224199) q[0];
rz(-pi) q[1];
rz(-2.372924) q[2];
sx q[2];
rz(-2.2109291) q[2];
sx q[2];
rz(1.6878355) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33448996) q[1];
sx q[1];
rz(-2.346171) q[1];
sx q[1];
rz(-2.9544178) q[1];
rz(-pi) q[2];
rz(1.5997361) q[3];
sx q[3];
rz(-0.77416285) q[3];
sx q[3];
rz(0.72284568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5840977) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-0.49475691) q[2];
rz(0.90302145) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(1.2899227) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(-2.1814573) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(0.18403149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83443123) q[0];
sx q[0];
rz(-2.8330028) q[0];
sx q[0];
rz(-0.35085268) q[0];
rz(-pi) q[1];
rz(2.5023735) q[2];
sx q[2];
rz(-0.44438617) q[2];
sx q[2];
rz(-1.5803312) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3369493) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(-0.6964535) q[1];
rz(2.9156978) q[3];
sx q[3];
rz(-2.4028006) q[3];
sx q[3];
rz(-1.6177288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(1.1408172) q[2];
rz(2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-1.1269349) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(0.88678962) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(0.14850798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69420544) q[0];
sx q[0];
rz(-2.2930657) q[0];
sx q[0];
rz(1.3556051) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4609023) q[2];
sx q[2];
rz(-1.3434778) q[2];
sx q[2];
rz(2.0109039) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.442246) q[1];
sx q[1];
rz(-1.8298501) q[1];
sx q[1];
rz(-3.0287292) q[1];
rz(-pi) q[2];
rz(-2.4466483) q[3];
sx q[3];
rz(-0.57999014) q[3];
sx q[3];
rz(-0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.285816) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-2.6532069) q[2];
rz(2.6521902) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4483036) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(-2.3235902) q[0];
rz(-1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-1.12524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97911191) q[0];
sx q[0];
rz(-1.1971803) q[0];
sx q[0];
rz(-3.0309249) q[0];
rz(-2.0411885) q[2];
sx q[2];
rz(-0.33249582) q[2];
sx q[2];
rz(-1.5645129) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4459671) q[1];
sx q[1];
rz(-1.8672767) q[1];
sx q[1];
rz(2.8545024) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4626477) q[3];
sx q[3];
rz(-0.79015398) q[3];
sx q[3];
rz(-1.8378347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0730878) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(1.4792431) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.6059426) q[0];
rz(-1.9203141) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(1.01064) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686417) q[0];
sx q[0];
rz(-0.7612345) q[0];
sx q[0];
rz(-1.3377454) q[0];
rz(0.29823093) q[2];
sx q[2];
rz(-1.094162) q[2];
sx q[2];
rz(0.44520608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38883868) q[1];
sx q[1];
rz(-2.8201582) q[1];
sx q[1];
rz(2.3366117) q[1];
x q[2];
rz(2.9270494) q[3];
sx q[3];
rz(-0.40896591) q[3];
sx q[3];
rz(1.6681125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6179787) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0555608) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(2.6468497) q[0];
rz(-2.5231979) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(-0.075597413) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48329096) q[0];
sx q[0];
rz(-2.1245983) q[0];
sx q[0];
rz(0.36853405) q[0];
rz(-1.3952903) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(2.4621071) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96502111) q[1];
sx q[1];
rz(-0.88639835) q[1];
sx q[1];
rz(-1.7040764) q[1];
x q[2];
rz(2.5209386) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(-1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33346924) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(-1.8010275) q[2];
rz(-2.8373485) q[3];
sx q[3];
rz(-2.2777568) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(2.9647968) q[0];
rz(1.2416174) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(0.44100824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4917883) q[0];
sx q[0];
rz(-1.7735964) q[0];
sx q[0];
rz(2.9359398) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6454234) q[2];
sx q[2];
rz(-1.7241038) q[2];
sx q[2];
rz(0.15664936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2446574) q[1];
sx q[1];
rz(-1.999563) q[1];
sx q[1];
rz(-2.7662617) q[1];
x q[2];
rz(1.8701843) q[3];
sx q[3];
rz(-0.96794879) q[3];
sx q[3];
rz(2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5796154) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(-0.30187541) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72398913) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(-0.026731116) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-0.84457196) q[2];
sx q[2];
rz(-1.209134) q[2];
sx q[2];
rz(-0.67187885) q[2];
rz(-1.0711014) q[3];
sx q[3];
rz(-1.0715967) q[3];
sx q[3];
rz(-1.788492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
