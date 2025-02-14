OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3616537) q[0];
sx q[0];
rz(-3.0335463) q[0];
sx q[0];
rz(7.530355) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(-1.127004) q[1];
sx q[1];
rz(-0.74581528) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0184263) q[0];
sx q[0];
rz(-1.2300876) q[0];
sx q[0];
rz(0.19827224) q[0];
rz(-3.132944) q[2];
sx q[2];
rz(-1.7315355) q[2];
sx q[2];
rz(0.5102821) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2850029) q[1];
sx q[1];
rz(-1.3599281) q[1];
sx q[1];
rz(3.0233356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6216032) q[3];
sx q[3];
rz(-1.4253748) q[3];
sx q[3];
rz(-3.036635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9790393) q[2];
sx q[2];
rz(-2.5472842) q[2];
sx q[2];
rz(-1.7621367) q[2];
rz(-2.4123794) q[3];
sx q[3];
rz(-2.3766434) q[3];
sx q[3];
rz(0.929207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5586229) q[0];
sx q[0];
rz(-1.0635149) q[0];
sx q[0];
rz(-2.9003918) q[0];
rz(-2.8855715) q[1];
sx q[1];
rz(-1.0216917) q[1];
sx q[1];
rz(2.3604438) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22750073) q[0];
sx q[0];
rz(-1.0100528) q[0];
sx q[0];
rz(-2.1888211) q[0];
rz(0.079226569) q[2];
sx q[2];
rz(-2.5753655) q[2];
sx q[2];
rz(1.5528284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9731739) q[1];
sx q[1];
rz(-0.89514625) q[1];
sx q[1];
rz(0.78898095) q[1];
rz(-pi) q[2];
rz(1.433302) q[3];
sx q[3];
rz(-1.526727) q[3];
sx q[3];
rz(-1.6018683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8476734) q[2];
sx q[2];
rz(-1.5519698) q[2];
sx q[2];
rz(-1.5576564) q[2];
rz(-0.0047575792) q[3];
sx q[3];
rz(-2.5465951) q[3];
sx q[3];
rz(2.7958272) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3874338) q[0];
sx q[0];
rz(-1.6922981) q[0];
sx q[0];
rz(0.19101983) q[0];
rz(0.35107958) q[1];
sx q[1];
rz(-2.7103238) q[1];
sx q[1];
rz(-1.692159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80181771) q[0];
sx q[0];
rz(-1.2317941) q[0];
sx q[0];
rz(-1.6247686) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75228508) q[2];
sx q[2];
rz(-2.5627081) q[2];
sx q[2];
rz(-0.22116254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23891029) q[1];
sx q[1];
rz(-2.604831) q[1];
sx q[1];
rz(-0.5384268) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6552917) q[3];
sx q[3];
rz(-2.4913906) q[3];
sx q[3];
rz(2.7335258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8423975) q[2];
sx q[2];
rz(-1.4292052) q[2];
sx q[2];
rz(-3.0136285) q[2];
rz(-1.9757102) q[3];
sx q[3];
rz(-2.2539625) q[3];
sx q[3];
rz(-0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(1.2853476) q[0];
sx q[0];
rz(-1.0839533) q[0];
sx q[0];
rz(-1.4803084) q[0];
rz(-3.0604494) q[1];
sx q[1];
rz(-2.0442043) q[1];
sx q[1];
rz(-0.88044423) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94784102) q[0];
sx q[0];
rz(-1.5955076) q[0];
sx q[0];
rz(-1.5773043) q[0];
x q[1];
rz(1.6398076) q[2];
sx q[2];
rz(-2.5216148) q[2];
sx q[2];
rz(-0.91199707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9040981) q[1];
sx q[1];
rz(-1.4597055) q[1];
sx q[1];
rz(2.6059125) q[1];
rz(-pi) q[2];
rz(-2.5838636) q[3];
sx q[3];
rz(-1.5909877) q[3];
sx q[3];
rz(-2.9021183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5389898) q[2];
sx q[2];
rz(-0.67402855) q[2];
sx q[2];
rz(1.4269786) q[2];
rz(-1.1566628) q[3];
sx q[3];
rz(-1.4256698) q[3];
sx q[3];
rz(2.9837933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2222774) q[0];
sx q[0];
rz(-0.8150402) q[0];
sx q[0];
rz(0.94006938) q[0];
rz(-2.6433511) q[1];
sx q[1];
rz(-0.91612852) q[1];
sx q[1];
rz(1.1246276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263382) q[0];
sx q[0];
rz(-0.7043411) q[0];
sx q[0];
rz(1.5096774) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0168991) q[2];
sx q[2];
rz(-2.6931525) q[2];
sx q[2];
rz(-0.61277991) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.985253) q[1];
sx q[1];
rz(-2.1731227) q[1];
sx q[1];
rz(-2.1039318) q[1];
rz(0.021291906) q[3];
sx q[3];
rz(-0.7856889) q[3];
sx q[3];
rz(-1.1710492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70790946) q[2];
sx q[2];
rz(-2.3072672) q[2];
sx q[2];
rz(1.679812) q[2];
rz(-2.8953569) q[3];
sx q[3];
rz(-2.2610531) q[3];
sx q[3];
rz(-1.0730526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4435302) q[0];
sx q[0];
rz(-1.9140697) q[0];
sx q[0];
rz(0.45502934) q[0];
rz(0.21806923) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(0.47854447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45947507) q[0];
sx q[0];
rz(-2.2322901) q[0];
sx q[0];
rz(-1.6949928) q[0];
rz(-2.5019574) q[2];
sx q[2];
rz(-1.860485) q[2];
sx q[2];
rz(1.3418496) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.45286059) q[1];
sx q[1];
rz(-1.11519) q[1];
sx q[1];
rz(-0.36700008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68799893) q[3];
sx q[3];
rz(-2.0562226) q[3];
sx q[3];
rz(1.7958876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9559418) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(-0.68598023) q[2];
rz(1.3395122) q[3];
sx q[3];
rz(-1.9728262) q[3];
sx q[3];
rz(1.3919938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5865536) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(-0.47898022) q[0];
rz(-2.2827177) q[1];
sx q[1];
rz(-0.54194599) q[1];
sx q[1];
rz(3.0368793) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0901439) q[0];
sx q[0];
rz(-1.7605172) q[0];
sx q[0];
rz(-0.39162292) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62868406) q[2];
sx q[2];
rz(-0.57839823) q[2];
sx q[2];
rz(2.4378928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4762209) q[1];
sx q[1];
rz(-1.7867861) q[1];
sx q[1];
rz(2.1588227) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0842544) q[3];
sx q[3];
rz(-2.0481589) q[3];
sx q[3];
rz(0.088378084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0571478) q[2];
sx q[2];
rz(-1.2556602) q[2];
sx q[2];
rz(-2.2596333) q[2];
rz(-3.0706578) q[3];
sx q[3];
rz(-1.8788012) q[3];
sx q[3];
rz(-0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1607745) q[0];
sx q[0];
rz(-2.6698298) q[0];
sx q[0];
rz(1.2534575) q[0];
rz(-0.71022931) q[1];
sx q[1];
rz(-1.9435147) q[1];
sx q[1];
rz(1.089878) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3364351) q[0];
sx q[0];
rz(-1.1602279) q[0];
sx q[0];
rz(-2.3768976) q[0];
rz(-pi) q[1];
rz(-2.5023191) q[2];
sx q[2];
rz(-0.54080039) q[2];
sx q[2];
rz(-1.6317473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1068258) q[1];
sx q[1];
rz(-2.2112101) q[1];
sx q[1];
rz(-2.5711715) q[1];
rz(0.57641455) q[3];
sx q[3];
rz(-1.4428328) q[3];
sx q[3];
rz(1.1157677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.39500427) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(2.11917) q[2];
rz(-2.5452781) q[3];
sx q[3];
rz(-0.78323451) q[3];
sx q[3];
rz(-2.7885126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8280867) q[0];
sx q[0];
rz(-0.49740288) q[0];
sx q[0];
rz(-1.585438) q[0];
rz(1.4900788) q[1];
sx q[1];
rz(-2.6863292) q[1];
sx q[1];
rz(-0.99686399) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8291959) q[0];
sx q[0];
rz(-1.2252843) q[0];
sx q[0];
rz(1.6457883) q[0];
rz(-pi) q[1];
rz(-3.1399675) q[2];
sx q[2];
rz(-1.7559984) q[2];
sx q[2];
rz(-2.8855326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91040033) q[1];
sx q[1];
rz(-2.2072189) q[1];
sx q[1];
rz(0.6271805) q[1];
x q[2];
rz(2.6410651) q[3];
sx q[3];
rz(-1.0382892) q[3];
sx q[3];
rz(-0.73126572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8263714) q[2];
sx q[2];
rz(-2.4085277) q[2];
sx q[2];
rz(0.23615393) q[2];
rz(-2.835623) q[3];
sx q[3];
rz(-1.9491842) q[3];
sx q[3];
rz(1.4845622) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21208256) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(2.4834852) q[0];
rz(1.2492389) q[1];
sx q[1];
rz(-2.55195) q[1];
sx q[1];
rz(0.49159893) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957053) q[0];
sx q[0];
rz(-2.108556) q[0];
sx q[0];
rz(2.0989492) q[0];
x q[1];
rz(-0.34296918) q[2];
sx q[2];
rz(-2.548647) q[2];
sx q[2];
rz(-0.85747257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9720077) q[1];
sx q[1];
rz(-1.2950724) q[1];
sx q[1];
rz(-1.5610525) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4171353) q[3];
sx q[3];
rz(-2.6479122) q[3];
sx q[3];
rz(0.069531893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3601941) q[2];
sx q[2];
rz(-0.34803826) q[2];
sx q[2];
rz(-1.6602824) q[2];
rz(-2.4252452) q[3];
sx q[3];
rz(-2.4951388) q[3];
sx q[3];
rz(-2.9334478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3502055) q[0];
sx q[0];
rz(-1.5972142) q[0];
sx q[0];
rz(0.41440339) q[0];
rz(2.1898337) q[1];
sx q[1];
rz(-1.2928243) q[1];
sx q[1];
rz(-1.181319) q[1];
rz(2.4454115) q[2];
sx q[2];
rz(-1.0715108) q[2];
sx q[2];
rz(0.61658695) q[2];
rz(2.2496102) q[3];
sx q[3];
rz(-1.4374229) q[3];
sx q[3];
rz(-2.0363804) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
