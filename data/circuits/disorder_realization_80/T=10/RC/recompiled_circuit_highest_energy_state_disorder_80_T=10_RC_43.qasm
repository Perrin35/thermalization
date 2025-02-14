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
rz(-0.031232746) q[1];
sx q[1];
rz(-0.90489689) q[1];
sx q[1];
rz(2.4651405) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.518153) q[0];
sx q[0];
rz(-1.4777967) q[0];
sx q[0];
rz(2.8661845) q[0];
x q[1];
rz(-2.7713369) q[2];
sx q[2];
rz(-0.91846648) q[2];
sx q[2];
rz(-0.45387646) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0259754) q[1];
sx q[1];
rz(-1.3968803) q[1];
sx q[1];
rz(-1.0976726) q[1];
rz(-pi) q[2];
rz(0.6630456) q[3];
sx q[3];
rz(-1.4395096) q[3];
sx q[3];
rz(-1.6754952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5260432) q[2];
sx q[2];
rz(-0.56576663) q[2];
sx q[2];
rz(-0.88741285) q[2];
rz(-0.081485661) q[3];
sx q[3];
rz(-1.2469651) q[3];
sx q[3];
rz(-2.2639349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(0.98608214) q[1];
sx q[1];
rz(-1.6796651) q[1];
sx q[1];
rz(2.2373534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3021224) q[0];
sx q[0];
rz(-1.6834489) q[0];
sx q[0];
rz(0.18160261) q[0];
rz(-pi) q[1];
x q[1];
rz(2.984513) q[2];
sx q[2];
rz(-1.6081405) q[2];
sx q[2];
rz(0.096028286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.820354) q[1];
sx q[1];
rz(-1.9776219) q[1];
sx q[1];
rz(-1.7443481) q[1];
rz(2.9504602) q[3];
sx q[3];
rz(-0.88273747) q[3];
sx q[3];
rz(0.21746527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6282689) q[2];
sx q[2];
rz(-2.8474035) q[2];
sx q[2];
rz(-2.7163556) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2929448) q[0];
sx q[0];
rz(-3.0631298) q[0];
sx q[0];
rz(-1.5516094) q[0];
rz(1.4750922) q[1];
sx q[1];
rz(-1.0868797) q[1];
sx q[1];
rz(1.0370673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1254409) q[0];
sx q[0];
rz(-0.96179987) q[0];
sx q[0];
rz(-2.9854965) q[0];
rz(-pi) q[1];
rz(1.6855459) q[2];
sx q[2];
rz(-2.1052573) q[2];
sx q[2];
rz(0.13335511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6810575) q[1];
sx q[1];
rz(-1.7112511) q[1];
sx q[1];
rz(1.3696425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3888168) q[3];
sx q[3];
rz(-1.1179384) q[3];
sx q[3];
rz(0.10145726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40917778) q[2];
sx q[2];
rz(-0.38984177) q[2];
sx q[2];
rz(1.8791492) q[2];
rz(-0.08517313) q[3];
sx q[3];
rz(-2.6933935) q[3];
sx q[3];
rz(1.1220773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20268102) q[0];
sx q[0];
rz(-1.0291463) q[0];
sx q[0];
rz(0.19126782) q[0];
rz(-0.8910886) q[1];
sx q[1];
rz(-2.8137408) q[1];
sx q[1];
rz(-2.0909615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7301431) q[0];
sx q[0];
rz(-2.3637584) q[0];
sx q[0];
rz(1.8671237) q[0];
x q[1];
rz(-0.3140788) q[2];
sx q[2];
rz(-2.2563612) q[2];
sx q[2];
rz(-2.2120145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6022608) q[1];
sx q[1];
rz(-0.51265279) q[1];
sx q[1];
rz(0.67474483) q[1];
rz(-pi) q[2];
rz(3.0043118) q[3];
sx q[3];
rz(-1.0044668) q[3];
sx q[3];
rz(-2.7911186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8656859) q[2];
sx q[2];
rz(-0.17912093) q[2];
sx q[2];
rz(2.5812145) q[2];
rz(3.0748902) q[3];
sx q[3];
rz(-1.3542465) q[3];
sx q[3];
rz(2.1393447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0340843) q[0];
sx q[0];
rz(-1.992724) q[0];
sx q[0];
rz(2.4565571) q[0];
rz(-2.7903453) q[1];
sx q[1];
rz(-2.5076187) q[1];
sx q[1];
rz(1.3397071) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23695859) q[0];
sx q[0];
rz(-1.5773505) q[0];
sx q[0];
rz(2.652524) q[0];
rz(1.1451911) q[2];
sx q[2];
rz(-2.1972924) q[2];
sx q[2];
rz(0.76001924) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9205528) q[1];
sx q[1];
rz(-1.5410081) q[1];
sx q[1];
rz(1.4090621) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6518526) q[3];
sx q[3];
rz(-1.9597754) q[3];
sx q[3];
rz(1.31969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54905218) q[2];
sx q[2];
rz(-1.668674) q[2];
sx q[2];
rz(-0.215691) q[2];
rz(-0.14821626) q[3];
sx q[3];
rz(-0.18311466) q[3];
sx q[3];
rz(-0.7557925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53629476) q[0];
sx q[0];
rz(-1.7114102) q[0];
sx q[0];
rz(0.17912616) q[0];
rz(0.99963775) q[1];
sx q[1];
rz(-2.3536286) q[1];
sx q[1];
rz(-2.9771908) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20227725) q[0];
sx q[0];
rz(-0.78139751) q[0];
sx q[0];
rz(2.8549749) q[0];
rz(2.6237409) q[2];
sx q[2];
rz(-1.4538063) q[2];
sx q[2];
rz(1.3530664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32443207) q[1];
sx q[1];
rz(-2.2304529) q[1];
sx q[1];
rz(-1.5800257) q[1];
rz(-pi) q[2];
rz(2.5980959) q[3];
sx q[3];
rz(-0.42146376) q[3];
sx q[3];
rz(1.3347365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1228235) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(1.4336047) q[2];
rz(2.150599) q[3];
sx q[3];
rz(-1.3905448) q[3];
sx q[3];
rz(1.8979134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3864415) q[0];
sx q[0];
rz(-0.069510892) q[0];
sx q[0];
rz(-0.25666562) q[0];
rz(-0.060404213) q[1];
sx q[1];
rz(-0.81788617) q[1];
sx q[1];
rz(0.43076691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2592436) q[0];
sx q[0];
rz(-1.8494864) q[0];
sx q[0];
rz(-0.95267991) q[0];
rz(2.45386) q[2];
sx q[2];
rz(-2.057325) q[2];
sx q[2];
rz(-0.4358049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8369981) q[1];
sx q[1];
rz(-0.82436383) q[1];
sx q[1];
rz(-0.5926822) q[1];
x q[2];
rz(-1.3957142) q[3];
sx q[3];
rz(-0.69718593) q[3];
sx q[3];
rz(0.28466636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7595235) q[2];
sx q[2];
rz(-2.6117595) q[2];
sx q[2];
rz(-2.1257373) q[2];
rz(-1.8633415) q[3];
sx q[3];
rz(-1.9081554) q[3];
sx q[3];
rz(2.5281233) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.627219) q[0];
sx q[0];
rz(-0.52424279) q[0];
sx q[0];
rz(-0.71994495) q[0];
rz(2.4607957) q[1];
sx q[1];
rz(-2.7578208) q[1];
sx q[1];
rz(2.0489571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62238111) q[0];
sx q[0];
rz(-1.5438612) q[0];
sx q[0];
rz(-0.48965298) q[0];
rz(-0.60694867) q[2];
sx q[2];
rz(-1.2649833) q[2];
sx q[2];
rz(0.94281498) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7842175) q[1];
sx q[1];
rz(-0.89723321) q[1];
sx q[1];
rz(-0.71460502) q[1];
rz(-pi) q[2];
rz(-0.25563116) q[3];
sx q[3];
rz(-1.6894332) q[3];
sx q[3];
rz(1.7411159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1277658) q[2];
sx q[2];
rz(-2.3769145) q[2];
sx q[2];
rz(-2.428425) q[2];
rz(1.3753043) q[3];
sx q[3];
rz(-1.2262552) q[3];
sx q[3];
rz(1.3885952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16667287) q[0];
sx q[0];
rz(-2.9627934) q[0];
sx q[0];
rz(1.7145702) q[0];
rz(2.7339281) q[1];
sx q[1];
rz(-1.0382321) q[1];
sx q[1];
rz(0.33956775) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9419562) q[0];
sx q[0];
rz(-2.3995281) q[0];
sx q[0];
rz(1.9626544) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5221065) q[2];
sx q[2];
rz(-1.9867989) q[2];
sx q[2];
rz(-0.69594072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.133278) q[1];
sx q[1];
rz(-0.61682683) q[1];
sx q[1];
rz(0.99442579) q[1];
x q[2];
rz(-1.6814713) q[3];
sx q[3];
rz(-0.95451285) q[3];
sx q[3];
rz(1.8679152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0582383) q[2];
sx q[2];
rz(-0.82667542) q[2];
sx q[2];
rz(0.44462407) q[2];
rz(-0.63117635) q[3];
sx q[3];
rz(-1.3381713) q[3];
sx q[3];
rz(1.1609424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6889965) q[0];
sx q[0];
rz(-1.8190374) q[0];
sx q[0];
rz(3.1153862) q[0];
rz(-1.3941437) q[1];
sx q[1];
rz(-1.1528287) q[1];
sx q[1];
rz(-2.0004415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3150605) q[0];
sx q[0];
rz(-0.71825832) q[0];
sx q[0];
rz(-0.63761561) q[0];
rz(-pi) q[1];
rz(0.097392453) q[2];
sx q[2];
rz(-0.44008128) q[2];
sx q[2];
rz(-2.4615364) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.929229) q[1];
sx q[1];
rz(-2.432669) q[1];
sx q[1];
rz(-0.89931788) q[1];
rz(0.23330063) q[3];
sx q[3];
rz(-0.58252347) q[3];
sx q[3];
rz(1.7482116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0673087) q[2];
sx q[2];
rz(-2.0412385) q[2];
sx q[2];
rz(1.9791774) q[2];
rz(-0.31636604) q[3];
sx q[3];
rz(-1.1917944) q[3];
sx q[3];
rz(1.6587862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3027073) q[0];
sx q[0];
rz(-3.0743297) q[0];
sx q[0];
rz(-0.59826941) q[0];
rz(-0.0089946714) q[1];
sx q[1];
rz(-2.3458377) q[1];
sx q[1];
rz(-1.5942106) q[1];
rz(-1.7187422) q[2];
sx q[2];
rz(-0.91397094) q[2];
sx q[2];
rz(1.915929) q[2];
rz(0.056445382) q[3];
sx q[3];
rz(-0.79462449) q[3];
sx q[3];
rz(0.36804646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
