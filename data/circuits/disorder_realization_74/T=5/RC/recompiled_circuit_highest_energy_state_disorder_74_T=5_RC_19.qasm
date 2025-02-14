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
rz(2.6994045) q[0];
sx q[0];
rz(-1.9884041) q[0];
sx q[0];
rz(-3.0129504) q[0];
rz(1.6755942) q[1];
sx q[1];
rz(-2.0425551) q[1];
sx q[1];
rz(1.0817945) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32444123) q[0];
sx q[0];
rz(-1.5400346) q[0];
sx q[0];
rz(0.064662393) q[0];
rz(1.9456909) q[2];
sx q[2];
rz(-1.257295) q[2];
sx q[2];
rz(1.9913265) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0073407081) q[1];
sx q[1];
rz(-1.6414343) q[1];
sx q[1];
rz(1.7299446) q[1];
x q[2];
rz(-3.0469037) q[3];
sx q[3];
rz(-1.296954) q[3];
sx q[3];
rz(-1.7750334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0470524) q[2];
sx q[2];
rz(-2.4461942) q[2];
sx q[2];
rz(0.73235861) q[2];
rz(-3.0426466) q[3];
sx q[3];
rz(-1.0354592) q[3];
sx q[3];
rz(1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.052213) q[0];
sx q[0];
rz(-2.3591924) q[0];
sx q[0];
rz(0.041444929) q[0];
rz(-0.6913569) q[1];
sx q[1];
rz(-1.3203011) q[1];
sx q[1];
rz(0.88821214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81538686) q[0];
sx q[0];
rz(-2.0688217) q[0];
sx q[0];
rz(-0.88125687) q[0];
x q[1];
rz(-0.67866831) q[2];
sx q[2];
rz(-1.9614842) q[2];
sx q[2];
rz(-0.44579577) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8687369) q[1];
sx q[1];
rz(-1.6140198) q[1];
sx q[1];
rz(0.24300139) q[1];
x q[2];
rz(2.668374) q[3];
sx q[3];
rz(-0.56269533) q[3];
sx q[3];
rz(-2.0317047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9899675) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(-1.7612696) q[2];
rz(0.049792854) q[3];
sx q[3];
rz(-2.4536665) q[3];
sx q[3];
rz(0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77408537) q[0];
sx q[0];
rz(-2.9893576) q[0];
sx q[0];
rz(-1.9277068) q[0];
rz(-1.8269352) q[1];
sx q[1];
rz(-1.0379125) q[1];
sx q[1];
rz(-1.5705869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11078283) q[0];
sx q[0];
rz(-2.6834496) q[0];
sx q[0];
rz(-0.58266298) q[0];
x q[1];
rz(-2.7941189) q[2];
sx q[2];
rz(-0.4336401) q[2];
sx q[2];
rz(2.4739802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8705026) q[1];
sx q[1];
rz(-1.6047641) q[1];
sx q[1];
rz(-1.7351331) q[1];
rz(-pi) q[2];
rz(2.1718484) q[3];
sx q[3];
rz(-1.0625524) q[3];
sx q[3];
rz(-2.598552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7870002) q[2];
sx q[2];
rz(-2.3921693) q[2];
sx q[2];
rz(-0.3053537) q[2];
rz(0.8461771) q[3];
sx q[3];
rz(-1.927522) q[3];
sx q[3];
rz(-0.86887104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8531891) q[0];
sx q[0];
rz(-1.0193634) q[0];
sx q[0];
rz(0.97002059) q[0];
rz(0.24523973) q[1];
sx q[1];
rz(-1.0673362) q[1];
sx q[1];
rz(1.6397569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062881447) q[0];
sx q[0];
rz(-2.9406266) q[0];
sx q[0];
rz(2.4522792) q[0];
rz(-2.9059783) q[2];
sx q[2];
rz(-1.4019626) q[2];
sx q[2];
rz(1.1323358) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3093331) q[1];
sx q[1];
rz(-0.24411476) q[1];
sx q[1];
rz(0.67863087) q[1];
rz(-pi) q[2];
rz(2.4145441) q[3];
sx q[3];
rz(-1.8173976) q[3];
sx q[3];
rz(-1.5545157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5331427) q[2];
sx q[2];
rz(-1.5641944) q[2];
sx q[2];
rz(1.2987785) q[2];
rz(-2.6070144) q[3];
sx q[3];
rz(-2.5375073) q[3];
sx q[3];
rz(-2.4728313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0534441) q[0];
sx q[0];
rz(-1.3548594) q[0];
sx q[0];
rz(0.97766367) q[0];
rz(0.68556249) q[1];
sx q[1];
rz(-2.3473163) q[1];
sx q[1];
rz(-2.175144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.179305) q[0];
sx q[0];
rz(-1.1098658) q[0];
sx q[0];
rz(-2.3910752) q[0];
x q[1];
rz(0.56880237) q[2];
sx q[2];
rz(-1.2035655) q[2];
sx q[2];
rz(-2.9609704) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.492474) q[1];
sx q[1];
rz(-0.8706514) q[1];
sx q[1];
rz(-2.9597069) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55111663) q[3];
sx q[3];
rz(-1.2663906) q[3];
sx q[3];
rz(1.1486883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5886249) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(0.28263131) q[2];
rz(2.1770554) q[3];
sx q[3];
rz(-1.5805565) q[3];
sx q[3];
rz(-0.36981043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033217) q[0];
sx q[0];
rz(-2.7300457) q[0];
sx q[0];
rz(2.0630398) q[0];
rz(0.82303965) q[1];
sx q[1];
rz(-1.1230725) q[1];
sx q[1];
rz(-0.80026921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33188021) q[0];
sx q[0];
rz(-0.030411424) q[0];
sx q[0];
rz(1.2872075) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3089463) q[2];
sx q[2];
rz(-1.8712723) q[2];
sx q[2];
rz(-2.7555675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2080704) q[1];
sx q[1];
rz(-0.88729837) q[1];
sx q[1];
rz(2.3166189) q[1];
x q[2];
rz(2.0144281) q[3];
sx q[3];
rz(-2.1882957) q[3];
sx q[3];
rz(-1.8007743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1396973) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(0.099420698) q[2];
rz(0.56441489) q[3];
sx q[3];
rz(-1.5494917) q[3];
sx q[3];
rz(2.2216589) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99892202) q[0];
sx q[0];
rz(-2.8746334) q[0];
sx q[0];
rz(0.024918407) q[0];
rz(-2.0492367) q[1];
sx q[1];
rz(-1.151029) q[1];
sx q[1];
rz(2.5999462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9077907) q[0];
sx q[0];
rz(-2.4176717) q[0];
sx q[0];
rz(-1.0095327) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.815258) q[2];
sx q[2];
rz(-1.0251364) q[2];
sx q[2];
rz(-1.5424002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.4723963) q[1];
sx q[1];
rz(-1.4265359) q[1];
sx q[1];
rz(1.2394106) q[1];
rz(-pi) q[2];
rz(-2.1775167) q[3];
sx q[3];
rz(-1.7172482) q[3];
sx q[3];
rz(1.3769384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8267374) q[2];
sx q[2];
rz(-1.9172226) q[2];
sx q[2];
rz(2.0014191) q[2];
rz(-1.4402116) q[3];
sx q[3];
rz(-0.39339104) q[3];
sx q[3];
rz(2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88406554) q[0];
sx q[0];
rz(-1.1648357) q[0];
sx q[0];
rz(-0.87673941) q[0];
rz(0.17403099) q[1];
sx q[1];
rz(-1.0935676) q[1];
sx q[1];
rz(1.7736951) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6606431) q[0];
sx q[0];
rz(-0.62134734) q[0];
sx q[0];
rz(-2.0578007) q[0];
rz(-2.1166685) q[2];
sx q[2];
rz(-0.41965963) q[2];
sx q[2];
rz(0.98463204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11238449) q[1];
sx q[1];
rz(-1.9416326) q[1];
sx q[1];
rz(0.98354323) q[1];
rz(-0.075966751) q[3];
sx q[3];
rz(-1.6923762) q[3];
sx q[3];
rz(-1.0978804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4725264) q[2];
sx q[2];
rz(-1.3269227) q[2];
sx q[2];
rz(2.9533022) q[2];
rz(1.3663728) q[3];
sx q[3];
rz(-1.3683189) q[3];
sx q[3];
rz(-1.202047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5240204) q[0];
sx q[0];
rz(-3.0078648) q[0];
sx q[0];
rz(-2.4753841) q[0];
rz(2.0062402) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(1.1599783) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2238623) q[0];
sx q[0];
rz(-2.1356815) q[0];
sx q[0];
rz(0.80819545) q[0];
rz(-0.23187821) q[2];
sx q[2];
rz(-0.7674976) q[2];
sx q[2];
rz(2.9434443) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0409909) q[1];
sx q[1];
rz(-1.578465) q[1];
sx q[1];
rz(0.95485447) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0464675) q[3];
sx q[3];
rz(-0.6199277) q[3];
sx q[3];
rz(1.4038831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3975415) q[2];
sx q[2];
rz(-2.5752189) q[2];
sx q[2];
rz(2.9446824) q[2];
rz(-0.87013733) q[3];
sx q[3];
rz(-2.0700442) q[3];
sx q[3];
rz(-1.2675233) q[3];
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
rz(0.63721913) q[0];
sx q[0];
rz(-1.0727896) q[0];
sx q[0];
rz(-0.037242591) q[0];
rz(2.0540909) q[1];
sx q[1];
rz(-1.0481513) q[1];
sx q[1];
rz(2.9748999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32023889) q[0];
sx q[0];
rz(-1.8379837) q[0];
sx q[0];
rz(0.70810476) q[0];
rz(-pi) q[1];
rz(2.6564381) q[2];
sx q[2];
rz(-1.302716) q[2];
sx q[2];
rz(1.175665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6065258) q[1];
sx q[1];
rz(-1.4826315) q[1];
sx q[1];
rz(-0.29523452) q[1];
rz(0.44346614) q[3];
sx q[3];
rz(-0.16032585) q[3];
sx q[3];
rz(0.043930862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4166261) q[2];
sx q[2];
rz(-2.5278957) q[2];
sx q[2];
rz(-1.7566768) q[2];
rz(-0.40063217) q[3];
sx q[3];
rz(-1.2847565) q[3];
sx q[3];
rz(1.0111151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2597802) q[0];
sx q[0];
rz(-2.0406944) q[0];
sx q[0];
rz(2.5115321) q[0];
rz(-3.0084445) q[1];
sx q[1];
rz(-1.9825736) q[1];
sx q[1];
rz(2.0864743) q[1];
rz(1.3721977) q[2];
sx q[2];
rz(-2.9469941) q[2];
sx q[2];
rz(-0.17175248) q[2];
rz(0.12689982) q[3];
sx q[3];
rz(-0.77346934) q[3];
sx q[3];
rz(1.0897549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
