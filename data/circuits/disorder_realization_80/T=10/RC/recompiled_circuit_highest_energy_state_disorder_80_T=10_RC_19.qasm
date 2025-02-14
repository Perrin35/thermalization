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
rz(-1.9266204) q[0];
sx q[0];
rz(1.5710693) q[0];
sx q[0];
rz(12.199353) q[0];
rz(3.1103599) q[1];
sx q[1];
rz(0.90489689) q[1];
sx q[1];
rz(8.7483258) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97359346) q[0];
sx q[0];
rz(-1.8449835) q[0];
sx q[0];
rz(1.4741769) q[0];
x q[1];
rz(1.1284087) q[2];
sx q[2];
rz(-2.4050875) q[2];
sx q[2];
rz(-2.1188096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2707959) q[1];
sx q[1];
rz(-2.6397986) q[1];
sx q[1];
rz(1.938799) q[1];
x q[2];
rz(-1.7368) q[3];
sx q[3];
rz(-2.2271384) q[3];
sx q[3];
rz(-2.9350077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5260432) q[2];
sx q[2];
rz(-2.575826) q[2];
sx q[2];
rz(2.2541798) q[2];
rz(-3.060107) q[3];
sx q[3];
rz(-1.2469651) q[3];
sx q[3];
rz(-0.87765774) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3990729) q[0];
sx q[0];
rz(-0.44624534) q[0];
sx q[0];
rz(1.9148069) q[0];
rz(-0.98608214) q[1];
sx q[1];
rz(-1.4619275) q[1];
sx q[1];
rz(-0.90423924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8522794) q[0];
sx q[0];
rz(-1.3903575) q[0];
sx q[0];
rz(-1.4562765) q[0];
rz(-pi) q[1];
rz(0.15707966) q[2];
sx q[2];
rz(-1.6081405) q[2];
sx q[2];
rz(3.0455644) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.095825776) q[1];
sx q[1];
rz(-2.7012022) q[1];
sx q[1];
rz(-2.7604123) q[1];
rz(-2.9504602) q[3];
sx q[3];
rz(-2.2588552) q[3];
sx q[3];
rz(-2.9241274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6282689) q[2];
sx q[2];
rz(-0.29418918) q[2];
sx q[2];
rz(0.42523709) q[2];
rz(0.50906316) q[3];
sx q[3];
rz(-2.013701) q[3];
sx q[3];
rz(1.7019255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8486479) q[0];
sx q[0];
rz(-3.0631298) q[0];
sx q[0];
rz(-1.5516094) q[0];
rz(-1.4750922) q[1];
sx q[1];
rz(-2.0547129) q[1];
sx q[1];
rz(-2.1045254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161517) q[0];
sx q[0];
rz(-2.1797928) q[0];
sx q[0];
rz(2.9854965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19107341) q[2];
sx q[2];
rz(-2.5961235) q[2];
sx q[2];
rz(-2.7857162) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.05987) q[1];
sx q[1];
rz(-1.7699426) q[1];
sx q[1];
rz(-0.14330602) q[1];
rz(-pi) q[2];
rz(-2.6821503) q[3];
sx q[3];
rz(-1.407335) q[3];
sx q[3];
rz(-1.3889988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7324149) q[2];
sx q[2];
rz(-0.38984177) q[2];
sx q[2];
rz(-1.8791492) q[2];
rz(-0.08517313) q[3];
sx q[3];
rz(-0.44819918) q[3];
sx q[3];
rz(-1.1220773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20268102) q[0];
sx q[0];
rz(-2.1124463) q[0];
sx q[0];
rz(-0.19126782) q[0];
rz(-0.8910886) q[1];
sx q[1];
rz(-2.8137408) q[1];
sx q[1];
rz(1.0506312) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3252661) q[0];
sx q[0];
rz(-2.3065595) q[0];
sx q[0];
rz(0.28006552) q[0];
rz(1.2096425) q[2];
sx q[2];
rz(-0.74336487) q[2];
sx q[2];
rz(-2.6860643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5393319) q[1];
sx q[1];
rz(-2.6289399) q[1];
sx q[1];
rz(-2.4668478) q[1];
rz(-1.0001783) q[3];
sx q[3];
rz(-1.6865391) q[3];
sx q[3];
rz(-1.9952578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8656859) q[2];
sx q[2];
rz(-2.9624717) q[2];
sx q[2];
rz(-0.56037819) q[2];
rz(-3.0748902) q[3];
sx q[3];
rz(-1.3542465) q[3];
sx q[3];
rz(-2.1393447) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0340843) q[0];
sx q[0];
rz(-1.992724) q[0];
sx q[0];
rz(-2.4565571) q[0];
rz(0.3512474) q[1];
sx q[1];
rz(-2.5076187) q[1];
sx q[1];
rz(1.3397071) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3303497) q[0];
sx q[0];
rz(-1.0817391) q[0];
sx q[0];
rz(1.5782209) q[0];
x q[1];
rz(1.9964016) q[2];
sx q[2];
rz(-2.1972924) q[2];
sx q[2];
rz(-0.76001924) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7869767) q[1];
sx q[1];
rz(-1.7324582) q[1];
sx q[1];
rz(3.1114108) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4248554) q[3];
sx q[3];
rz(-2.5261758) q[3];
sx q[3];
rz(0.36726609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54905218) q[2];
sx q[2];
rz(-1.4729187) q[2];
sx q[2];
rz(2.9259017) q[2];
rz(2.9933764) q[3];
sx q[3];
rz(-0.18311466) q[3];
sx q[3];
rz(2.3858002) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53629476) q[0];
sx q[0];
rz(-1.4301825) q[0];
sx q[0];
rz(0.17912616) q[0];
rz(2.1419549) q[1];
sx q[1];
rz(-0.78796402) q[1];
sx q[1];
rz(-2.9771908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1622551) q[0];
sx q[0];
rz(-1.3703523) q[0];
sx q[0];
rz(-0.76058273) q[0];
rz(-pi) q[1];
x q[1];
rz(1.705242) q[2];
sx q[2];
rz(-2.0847581) q[2];
sx q[2];
rz(-2.9902637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8021014) q[1];
sx q[1];
rz(-0.65971148) q[1];
sx q[1];
rz(-0.011899634) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54349676) q[3];
sx q[3];
rz(-2.7201289) q[3];
sx q[3];
rz(-1.8068562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0187692) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(1.4336047) q[2];
rz(2.150599) q[3];
sx q[3];
rz(-1.7510479) q[3];
sx q[3];
rz(-1.8979134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7551512) q[0];
sx q[0];
rz(-3.0720818) q[0];
sx q[0];
rz(-2.884927) q[0];
rz(-0.060404213) q[1];
sx q[1];
rz(-2.3237065) q[1];
sx q[1];
rz(2.7108257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11837932) q[0];
sx q[0];
rz(-2.1617365) q[0];
sx q[0];
rz(2.803938) q[0];
x q[1];
rz(2.4469023) q[2];
sx q[2];
rz(-0.81889924) q[2];
sx q[2];
rz(-1.652405) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6182729) q[1];
sx q[1];
rz(-0.91611169) q[1];
sx q[1];
rz(1.0275082) q[1];
rz(-pi) q[2];
rz(-1.3957142) q[3];
sx q[3];
rz(-0.69718593) q[3];
sx q[3];
rz(-2.8569263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7595235) q[2];
sx q[2];
rz(-0.52983317) q[2];
sx q[2];
rz(-1.0158553) q[2];
rz(1.2782512) q[3];
sx q[3];
rz(-1.2334373) q[3];
sx q[3];
rz(0.61346936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51437369) q[0];
sx q[0];
rz(-2.6173499) q[0];
sx q[0];
rz(-2.4216477) q[0];
rz(-2.4607957) q[1];
sx q[1];
rz(-2.7578208) q[1];
sx q[1];
rz(1.0926355) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1426861) q[0];
sx q[0];
rz(-0.49033316) q[0];
sx q[0];
rz(-0.057221091) q[0];
x q[1];
rz(-0.50555412) q[2];
sx q[2];
rz(-0.67086911) q[2];
sx q[2];
rz(-1.0370924) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35737513) q[1];
sx q[1];
rz(-0.89723321) q[1];
sx q[1];
rz(2.4269876) q[1];
rz(-pi) q[2];
rz(1.4482139) q[3];
sx q[3];
rz(-1.8245909) q[3];
sx q[3];
rz(-0.13939652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1277658) q[2];
sx q[2];
rz(-2.3769145) q[2];
sx q[2];
rz(-0.71316767) q[2];
rz(1.7662883) q[3];
sx q[3];
rz(-1.2262552) q[3];
sx q[3];
rz(-1.3885952) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9749198) q[0];
sx q[0];
rz(-2.9627934) q[0];
sx q[0];
rz(-1.4270225) q[0];
rz(-0.40766454) q[1];
sx q[1];
rz(-2.1033606) q[1];
sx q[1];
rz(2.8020249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996364) q[0];
sx q[0];
rz(-0.74206458) q[0];
sx q[0];
rz(-1.1789382) q[0];
rz(-3.0318674) q[2];
sx q[2];
rz(-2.7229157) q[2];
sx q[2];
rz(2.3256486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8060866) q[1];
sx q[1];
rz(-1.0644344) q[1];
sx q[1];
rz(0.36878364) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4601213) q[3];
sx q[3];
rz(-0.95451285) q[3];
sx q[3];
rz(-1.8679152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0582383) q[2];
sx q[2];
rz(-2.3149172) q[2];
sx q[2];
rz(0.44462407) q[2];
rz(0.63117635) q[3];
sx q[3];
rz(-1.3381713) q[3];
sx q[3];
rz(-1.1609424) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6889965) q[0];
sx q[0];
rz(-1.3225553) q[0];
sx q[0];
rz(3.1153862) q[0];
rz(-1.7474489) q[1];
sx q[1];
rz(-1.988764) q[1];
sx q[1];
rz(1.1411512) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8265322) q[0];
sx q[0];
rz(-0.71825832) q[0];
sx q[0];
rz(2.503977) q[0];
rz(-pi) q[1];
rz(0.43825324) q[2];
sx q[2];
rz(-1.5293596) q[2];
sx q[2];
rz(2.1626894) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2403912) q[1];
sx q[1];
rz(-1.1537885) q[1];
sx q[1];
rz(-2.1621125) q[1];
rz(-pi) q[2];
rz(0.23330063) q[3];
sx q[3];
rz(-0.58252347) q[3];
sx q[3];
rz(1.7482116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0673087) q[2];
sx q[2];
rz(-2.0412385) q[2];
sx q[2];
rz(-1.1624153) q[2];
rz(-0.31636604) q[3];
sx q[3];
rz(-1.1917944) q[3];
sx q[3];
rz(1.6587862) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3027073) q[0];
sx q[0];
rz(-0.067262983) q[0];
sx q[0];
rz(2.5433232) q[0];
rz(0.0089946714) q[1];
sx q[1];
rz(-0.79575494) q[1];
sx q[1];
rz(1.5473821) q[1];
rz(0.18890201) q[2];
sx q[2];
rz(-0.67086611) q[2];
sx q[2];
rz(-0.98626731) q[2];
rz(1.5133933) q[3];
sx q[3];
rz(-0.77779278) q[3];
sx q[3];
rz(-2.8540303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
