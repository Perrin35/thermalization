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
rz(1.4752969) q[0];
sx q[0];
rz(-1.2694321) q[0];
sx q[0];
rz(0.67212927) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(-2.3840005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8586977) q[0];
sx q[0];
rz(-3.0083249) q[0];
sx q[0];
rz(-0.96576502) q[0];
x q[1];
rz(-0.50067164) q[2];
sx q[2];
rz(-1.2462052) q[2];
sx q[2];
rz(-0.90510923) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.014908286) q[1];
sx q[1];
rz(-1.9796625) q[1];
sx q[1];
rz(-1.3247299) q[1];
rz(-pi) q[2];
rz(0.67528649) q[3];
sx q[3];
rz(-2.3216341) q[3];
sx q[3];
rz(-0.91778558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0932833) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(-1.3422356) q[2];
rz(-1.7736769) q[3];
sx q[3];
rz(-1.0932837) q[3];
sx q[3];
rz(-1.5526519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9369478) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(2.192705) q[0];
rz(-3.0401547) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(-2.2116275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0819433) q[0];
sx q[0];
rz(-0.23632061) q[0];
sx q[0];
rz(0.44775072) q[0];
x q[1];
rz(1.5783159) q[2];
sx q[2];
rz(-1.0656271) q[2];
sx q[2];
rz(-1.1669605) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3392433) q[1];
sx q[1];
rz(-1.4096469) q[1];
sx q[1];
rz(1.4646962) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3306777) q[3];
sx q[3];
rz(-1.0958015) q[3];
sx q[3];
rz(-2.7678633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(-0.742221) q[2];
rz(2.6774075) q[3];
sx q[3];
rz(-1.3451385) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6716229) q[0];
sx q[0];
rz(-2.9556584) q[0];
sx q[0];
rz(-1.4416913) q[0];
rz(0.16547671) q[1];
sx q[1];
rz(-0.73935699) q[1];
sx q[1];
rz(2.2742719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.353426) q[0];
sx q[0];
rz(-1.5713619) q[0];
sx q[0];
rz(-4.0325469e-05) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27133743) q[2];
sx q[2];
rz(-2.0607161) q[2];
sx q[2];
rz(-2.4957531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9359303) q[1];
sx q[1];
rz(-1.2453658) q[1];
sx q[1];
rz(1.8673351) q[1];
rz(-pi) q[2];
rz(2.1380566) q[3];
sx q[3];
rz(-2.2972882) q[3];
sx q[3];
rz(-2.9851825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9637588) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(2.3826694) q[2];
rz(2.3376076) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(-0.72833958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8300962) q[0];
sx q[0];
rz(-1.8599956) q[0];
sx q[0];
rz(2.6575644) q[0];
rz(1.6054035) q[1];
sx q[1];
rz(-1.7264629) q[1];
sx q[1];
rz(-1.7162292) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40544505) q[0];
sx q[0];
rz(-2.510294) q[0];
sx q[0];
rz(-2.5456356) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62230939) q[2];
sx q[2];
rz(-0.37823411) q[2];
sx q[2];
rz(0.44307274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0565383) q[1];
sx q[1];
rz(-0.94496545) q[1];
sx q[1];
rz(2.8686348) q[1];
rz(-pi) q[2];
rz(1.1953765) q[3];
sx q[3];
rz(-1.9945126) q[3];
sx q[3];
rz(-2.8590607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.385685) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(-0.29328406) q[2];
rz(0.048132345) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(-0.3046681) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83679477) q[0];
sx q[0];
rz(-1.7339107) q[0];
sx q[0];
rz(2.0378713) q[0];
rz(0.91148218) q[1];
sx q[1];
rz(-1.5244923) q[1];
sx q[1];
rz(-3.0326861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24507228) q[0];
sx q[0];
rz(-0.86045107) q[0];
sx q[0];
rz(2.4857387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2997024) q[2];
sx q[2];
rz(-2.5863918) q[2];
sx q[2];
rz(0.19367684) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6088691) q[1];
sx q[1];
rz(-2.2298498) q[1];
sx q[1];
rz(0.20688914) q[1];
x q[2];
rz(1.2690684) q[3];
sx q[3];
rz(-2.1273888) q[3];
sx q[3];
rz(2.3315786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4917422) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(2.8835127) q[2];
rz(-0.38170013) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(2.5842353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.7447164) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(-1.1599524) q[0];
rz(0.57811919) q[1];
sx q[1];
rz(-1.6696397) q[1];
sx q[1];
rz(0.32346183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9834546) q[0];
sx q[0];
rz(-1.1295415) q[0];
sx q[0];
rz(-2.6652314) q[0];
rz(0.44536369) q[2];
sx q[2];
rz(-2.6008743) q[2];
sx q[2];
rz(3.0034844) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9449294) q[1];
sx q[1];
rz(-0.62360901) q[1];
sx q[1];
rz(-0.39396472) q[1];
rz(-pi) q[2];
rz(-2.6810535) q[3];
sx q[3];
rz(-1.8611188) q[3];
sx q[3];
rz(-3.0825305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20299992) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(-2.2933551) q[2];
rz(1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(0.69376865) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11442014) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(0.87345901) q[0];
rz(1.2365485) q[1];
sx q[1];
rz(-0.97573391) q[1];
sx q[1];
rz(-0.083560856) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3665584) q[0];
sx q[0];
rz(-1.2511484) q[0];
sx q[0];
rz(1.699422) q[0];
x q[1];
rz(2.9026743) q[2];
sx q[2];
rz(-2.0878125) q[2];
sx q[2];
rz(0.97415249) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.04490964) q[1];
sx q[1];
rz(-2.2538141) q[1];
sx q[1];
rz(0.51814305) q[1];
rz(-pi) q[2];
rz(0.36079455) q[3];
sx q[3];
rz(-1.8550145) q[3];
sx q[3];
rz(2.7786215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43625912) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(1.7748888) q[2];
rz(2.8691835) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347539) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(2.2032264) q[0];
rz(-2.3387108) q[1];
sx q[1];
rz(-2.1736841) q[1];
sx q[1];
rz(-0.82069194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21194776) q[0];
sx q[0];
rz(-1.7045472) q[0];
sx q[0];
rz(1.9424214) q[0];
rz(1.1633881) q[2];
sx q[2];
rz(-2.738224) q[2];
sx q[2];
rz(0.58248108) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.332224) q[1];
sx q[1];
rz(-1.4648155) q[1];
sx q[1];
rz(-2.1952704) q[1];
rz(-pi) q[2];
rz(-1.197416) q[3];
sx q[3];
rz(-0.42694091) q[3];
sx q[3];
rz(0.77746848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20251033) q[2];
sx q[2];
rz(-0.15190092) q[2];
sx q[2];
rz(1.067777) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(3.06156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.985567) q[0];
sx q[0];
rz(-1.1556867) q[0];
sx q[0];
rz(0.40618968) q[0];
rz(-1.4093026) q[1];
sx q[1];
rz(-2.1235762) q[1];
sx q[1];
rz(-2.0565775) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70124088) q[0];
sx q[0];
rz(-2.1862162) q[0];
sx q[0];
rz(1.3568272) q[0];
rz(-pi) q[1];
rz(-0.0064956587) q[2];
sx q[2];
rz(-0.61031872) q[2];
sx q[2];
rz(1.4607061) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20445261) q[1];
sx q[1];
rz(-1.8984183) q[1];
sx q[1];
rz(2.3069068) q[1];
rz(-pi) q[2];
rz(-2.1187339) q[3];
sx q[3];
rz(-1.0617439) q[3];
sx q[3];
rz(-1.7804543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5558418) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(-2.2108868) q[2];
rz(1.3961004) q[3];
sx q[3];
rz(-1.7788818) q[3];
sx q[3];
rz(3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4569106) q[0];
sx q[0];
rz(-2.1645808) q[0];
sx q[0];
rz(1.3071625) q[0];
rz(-0.3859418) q[1];
sx q[1];
rz(-1.7470876) q[1];
sx q[1];
rz(1.0345667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86901122) q[0];
sx q[0];
rz(-2.9391461) q[0];
sx q[0];
rz(2.1073209) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7427086) q[2];
sx q[2];
rz(-2.6475057) q[2];
sx q[2];
rz(1.2436155) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2092691) q[1];
sx q[1];
rz(-2.0210365) q[1];
sx q[1];
rz(-0.85659124) q[1];
rz(-pi) q[2];
rz(-0.12328445) q[3];
sx q[3];
rz(-0.81466952) q[3];
sx q[3];
rz(-2.2565763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4474386) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(2.7745957) q[2];
rz(1.1191818) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601396) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(-0.057859261) q[1];
sx q[1];
rz(-1.8744938) q[1];
sx q[1];
rz(-1.4300463) q[1];
rz(2.62769) q[2];
sx q[2];
rz(-0.19795098) q[2];
sx q[2];
rz(1.2951938) q[2];
rz(1.9290942) q[3];
sx q[3];
rz(-0.90322106) q[3];
sx q[3];
rz(-1.7424104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
