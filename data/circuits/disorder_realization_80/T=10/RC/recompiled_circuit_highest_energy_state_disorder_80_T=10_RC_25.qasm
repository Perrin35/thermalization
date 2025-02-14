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
rz(1.2149723) q[0];
sx q[0];
rz(-1.5710693) q[0];
sx q[0];
rz(-0.36701742) q[0];
rz(3.1103599) q[1];
sx q[1];
rz(0.90489689) q[1];
sx q[1];
rz(8.7483258) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1679992) q[0];
sx q[0];
rz(-1.8449835) q[0];
sx q[0];
rz(-1.6674158) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37025578) q[2];
sx q[2];
rz(-2.2231262) q[2];
sx q[2];
rz(0.45387646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2707959) q[1];
sx q[1];
rz(-2.6397986) q[1];
sx q[1];
rz(1.938799) q[1];
rz(-pi) q[2];
rz(-1.7368) q[3];
sx q[3];
rz(-2.2271384) q[3];
sx q[3];
rz(0.20658499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5260432) q[2];
sx q[2];
rz(-2.575826) q[2];
sx q[2];
rz(-2.2541798) q[2];
rz(-3.060107) q[3];
sx q[3];
rz(-1.2469651) q[3];
sx q[3];
rz(2.2639349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74251974) q[0];
sx q[0];
rz(-0.44624534) q[0];
sx q[0];
rz(-1.2267858) q[0];
rz(-2.1555105) q[1];
sx q[1];
rz(-1.4619275) q[1];
sx q[1];
rz(0.90423924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2806633) q[0];
sx q[0];
rz(-0.21337803) q[0];
sx q[0];
rz(0.55960525) q[0];
rz(-pi) q[1];
rz(0.23444011) q[2];
sx q[2];
rz(-2.9801705) q[2];
sx q[2];
rz(1.2432673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8227742) q[1];
sx q[1];
rz(-1.4115361) q[1];
sx q[1];
rz(0.41235473) q[1];
x q[2];
rz(1.3436965) q[3];
sx q[3];
rz(-0.70992142) q[3];
sx q[3];
rz(2.6283669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51332372) q[2];
sx q[2];
rz(-2.8474035) q[2];
sx q[2];
rz(-0.42523709) q[2];
rz(2.6325295) q[3];
sx q[3];
rz(-1.1278917) q[3];
sx q[3];
rz(-1.4396671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2929448) q[0];
sx q[0];
rz(-0.078462891) q[0];
sx q[0];
rz(1.5516094) q[0];
rz(-1.6665005) q[1];
sx q[1];
rz(-1.0868797) q[1];
sx q[1];
rz(1.0370673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0161517) q[0];
sx q[0];
rz(-0.96179987) q[0];
sx q[0];
rz(-0.15609619) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6042348) q[2];
sx q[2];
rz(-1.4721057) q[2];
sx q[2];
rz(-1.3787998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4296966) q[1];
sx q[1];
rz(-0.24479228) q[1];
sx q[1];
rz(-2.1866287) q[1];
x q[2];
rz(-0.45944233) q[3];
sx q[3];
rz(-1.407335) q[3];
sx q[3];
rz(-1.7525938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7324149) q[2];
sx q[2];
rz(-2.7517509) q[2];
sx q[2];
rz(-1.2624435) q[2];
rz(-3.0564195) q[3];
sx q[3];
rz(-2.6933935) q[3];
sx q[3];
rz(-1.1220773) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9389116) q[0];
sx q[0];
rz(-1.0291463) q[0];
sx q[0];
rz(2.9503248) q[0];
rz(-0.8910886) q[1];
sx q[1];
rz(-2.8137408) q[1];
sx q[1];
rz(-2.0909615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4114496) q[0];
sx q[0];
rz(-2.3637584) q[0];
sx q[0];
rz(-1.8671237) q[0];
x q[1];
rz(2.2810535) q[2];
sx q[2];
rz(-1.8122753) q[2];
sx q[2];
rz(-0.84404404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5393319) q[1];
sx q[1];
rz(-2.6289399) q[1];
sx q[1];
rz(-0.67474483) q[1];
x q[2];
rz(-1.3587908) q[3];
sx q[3];
rz(-0.58095914) q[3];
sx q[3];
rz(0.60248366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27590671) q[2];
sx q[2];
rz(-0.17912093) q[2];
sx q[2];
rz(-2.5812145) q[2];
rz(-0.066702453) q[3];
sx q[3];
rz(-1.7873462) q[3];
sx q[3];
rz(-2.1393447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0340843) q[0];
sx q[0];
rz(-1.992724) q[0];
sx q[0];
rz(-0.68503553) q[0];
rz(2.7903453) q[1];
sx q[1];
rz(-2.5076187) q[1];
sx q[1];
rz(1.8018855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.811243) q[0];
sx q[0];
rz(-2.0598536) q[0];
sx q[0];
rz(-1.5782209) q[0];
x q[1];
rz(1.1451911) q[2];
sx q[2];
rz(-0.94430021) q[2];
sx q[2];
rz(2.3815734) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9205528) q[1];
sx q[1];
rz(-1.5410081) q[1];
sx q[1];
rz(-1.4090621) q[1];
x q[2];
rz(2.4248554) q[3];
sx q[3];
rz(-0.61541688) q[3];
sx q[3];
rz(-0.36726609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54905218) q[2];
sx q[2];
rz(-1.4729187) q[2];
sx q[2];
rz(-2.9259017) q[2];
rz(-2.9933764) q[3];
sx q[3];
rz(-0.18311466) q[3];
sx q[3];
rz(0.7557925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6052979) q[0];
sx q[0];
rz(-1.7114102) q[0];
sx q[0];
rz(-2.9624665) q[0];
rz(-0.99963775) q[1];
sx q[1];
rz(-0.78796402) q[1];
sx q[1];
rz(-2.9771908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20227725) q[0];
sx q[0];
rz(-2.3601951) q[0];
sx q[0];
rz(0.28661771) q[0];
rz(-pi) q[1];
rz(-2.6237409) q[2];
sx q[2];
rz(-1.4538063) q[2];
sx q[2];
rz(1.7885263) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32443207) q[1];
sx q[1];
rz(-2.2304529) q[1];
sx q[1];
rz(-1.5615669) q[1];
x q[2];
rz(1.3429759) q[3];
sx q[3];
rz(-1.2130672) q[3];
sx q[3];
rz(-1.9195956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0187692) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(1.7079879) q[2];
rz(-2.150599) q[3];
sx q[3];
rz(-1.7510479) q[3];
sx q[3];
rz(-1.2436793) q[3];
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
rz(-pi/2) q[3];
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
rz(1.7551512) q[0];
sx q[0];
rz(-3.0720818) q[0];
sx q[0];
rz(0.25666562) q[0];
rz(0.060404213) q[1];
sx q[1];
rz(-0.81788617) q[1];
sx q[1];
rz(2.7108257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2592436) q[0];
sx q[0];
rz(-1.8494864) q[0];
sx q[0];
rz(2.1889127) q[0];
x q[1];
rz(-0.69469037) q[2];
sx q[2];
rz(-2.3226934) q[2];
sx q[2];
rz(1.652405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3045945) q[1];
sx q[1];
rz(-2.3172288) q[1];
sx q[1];
rz(0.5926822) q[1];
rz(2.9967331) q[3];
sx q[3];
rz(-2.2552465) q[3];
sx q[3];
rz(3.0836925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3820691) q[2];
sx q[2];
rz(-0.52983317) q[2];
sx q[2];
rz(-1.0158553) q[2];
rz(1.8633415) q[3];
sx q[3];
rz(-1.2334373) q[3];
sx q[3];
rz(-0.61346936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51437369) q[0];
sx q[0];
rz(-2.6173499) q[0];
sx q[0];
rz(-0.71994495) q[0];
rz(-0.68079692) q[1];
sx q[1];
rz(-2.7578208) q[1];
sx q[1];
rz(2.0489571) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9989065) q[0];
sx q[0];
rz(-2.6512595) q[0];
sx q[0];
rz(0.057221091) q[0];
rz(0.50555412) q[2];
sx q[2];
rz(-0.67086911) q[2];
sx q[2];
rz(-2.1045002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7094575) q[1];
sx q[1];
rz(-2.1085408) q[1];
sx q[1];
rz(-0.75790578) q[1];
rz(1.4482139) q[3];
sx q[3];
rz(-1.8245909) q[3];
sx q[3];
rz(-0.13939652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0138268) q[2];
sx q[2];
rz(-0.76467815) q[2];
sx q[2];
rz(2.428425) q[2];
rz(1.7662883) q[3];
sx q[3];
rz(-1.9153374) q[3];
sx q[3];
rz(-1.7529974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9749198) q[0];
sx q[0];
rz(-2.9627934) q[0];
sx q[0];
rz(1.4270225) q[0];
rz(0.40766454) q[1];
sx q[1];
rz(-1.0382321) q[1];
sx q[1];
rz(2.8020249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996364) q[0];
sx q[0];
rz(-0.74206458) q[0];
sx q[0];
rz(1.1789382) q[0];
x q[1];
rz(-3.0318674) q[2];
sx q[2];
rz(-2.7229157) q[2];
sx q[2];
rz(2.3256486) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0083146) q[1];
sx q[1];
rz(-2.5247658) q[1];
sx q[1];
rz(-2.1471669) q[1];
rz(1.4601213) q[3];
sx q[3];
rz(-0.95451285) q[3];
sx q[3];
rz(1.8679152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0833544) q[2];
sx q[2];
rz(-2.3149172) q[2];
sx q[2];
rz(-2.6969686) q[2];
rz(-0.63117635) q[3];
sx q[3];
rz(-1.8034214) q[3];
sx q[3];
rz(-1.1609424) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6889965) q[0];
sx q[0];
rz(-1.8190374) q[0];
sx q[0];
rz(0.026206503) q[0];
rz(-1.3941437) q[1];
sx q[1];
rz(-1.1528287) q[1];
sx q[1];
rz(-2.0004415) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25309207) q[0];
sx q[0];
rz(-1.1682751) q[0];
sx q[0];
rz(-2.5293468) q[0];
rz(-pi) q[1];
rz(-3.0442002) q[2];
sx q[2];
rz(-2.7015114) q[2];
sx q[2];
rz(2.4615364) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.929229) q[1];
sx q[1];
rz(-0.7089237) q[1];
sx q[1];
rz(2.2422748) q[1];
x q[2];
rz(0.23330063) q[3];
sx q[3];
rz(-0.58252347) q[3];
sx q[3];
rz(-1.3933811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0673087) q[2];
sx q[2];
rz(-2.0412385) q[2];
sx q[2];
rz(1.9791774) q[2];
rz(-0.31636604) q[3];
sx q[3];
rz(-1.9497982) q[3];
sx q[3];
rz(-1.6587862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
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
rz(1.4228504) q[2];
sx q[2];
rz(-0.91397094) q[2];
sx q[2];
rz(1.915929) q[2];
rz(-2.347765) q[3];
sx q[3];
rz(-1.5305274) q[3];
sx q[3];
rz(-1.2423142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
