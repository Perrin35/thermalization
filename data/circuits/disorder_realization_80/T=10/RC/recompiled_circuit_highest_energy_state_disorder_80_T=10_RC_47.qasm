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
rz(-2.2366958) q[1];
sx q[1];
rz(0.67645216) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1679992) q[0];
sx q[0];
rz(-1.8449835) q[0];
sx q[0];
rz(-1.4741769) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1284087) q[2];
sx q[2];
rz(-2.4050875) q[2];
sx q[2];
rz(-1.022783) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0259754) q[1];
sx q[1];
rz(-1.7447124) q[1];
sx q[1];
rz(-1.0976726) q[1];
rz(-pi) q[2];
rz(-0.6630456) q[3];
sx q[3];
rz(-1.4395096) q[3];
sx q[3];
rz(-1.4660975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61554945) q[2];
sx q[2];
rz(-0.56576663) q[2];
sx q[2];
rz(-2.2541798) q[2];
rz(-3.060107) q[3];
sx q[3];
rz(-1.2469651) q[3];
sx q[3];
rz(-0.87765774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3990729) q[0];
sx q[0];
rz(-0.44624534) q[0];
sx q[0];
rz(-1.2267858) q[0];
rz(-0.98608214) q[1];
sx q[1];
rz(-1.4619275) q[1];
sx q[1];
rz(2.2373534) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8394703) q[0];
sx q[0];
rz(-1.6834489) q[0];
sx q[0];
rz(-2.95999) q[0];
rz(-0.15707966) q[2];
sx q[2];
rz(-1.5334522) q[2];
sx q[2];
rz(-0.096028286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0457669) q[1];
sx q[1];
rz(-2.7012022) q[1];
sx q[1];
rz(0.38118036) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2678862) q[3];
sx q[3];
rz(-1.4235157) q[3];
sx q[3];
rz(-1.2310674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6282689) q[2];
sx q[2];
rz(-2.8474035) q[2];
sx q[2];
rz(0.42523709) q[2];
rz(0.50906316) q[3];
sx q[3];
rz(-1.1278917) q[3];
sx q[3];
rz(-1.7019255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.8486479) q[0];
sx q[0];
rz(-0.078462891) q[0];
sx q[0];
rz(1.5516094) q[0];
rz(1.4750922) q[1];
sx q[1];
rz(-2.0547129) q[1];
sx q[1];
rz(2.1045254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3939115) q[0];
sx q[0];
rz(-2.5153749) q[0];
sx q[0];
rz(-1.3514723) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53735781) q[2];
sx q[2];
rz(-1.6694869) q[2];
sx q[2];
rz(-1.3787998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.711896) q[1];
sx q[1];
rz(-2.8968004) q[1];
sx q[1];
rz(0.954964) q[1];
rz(-pi) q[2];
rz(2.7855139) q[3];
sx q[3];
rz(-0.48569187) q[3];
sx q[3];
rz(2.6420223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7324149) q[2];
sx q[2];
rz(-2.7517509) q[2];
sx q[2];
rz(-1.8791492) q[2];
rz(0.08517313) q[3];
sx q[3];
rz(-2.6933935) q[3];
sx q[3];
rz(-1.1220773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20268102) q[0];
sx q[0];
rz(-2.1124463) q[0];
sx q[0];
rz(2.9503248) q[0];
rz(-2.250504) q[1];
sx q[1];
rz(-0.32785186) q[1];
sx q[1];
rz(1.0506312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8163265) q[0];
sx q[0];
rz(-2.3065595) q[0];
sx q[0];
rz(0.28006552) q[0];
x q[1];
rz(-0.3140788) q[2];
sx q[2];
rz(-0.88523141) q[2];
sx q[2];
rz(-0.92957815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85957789) q[1];
sx q[1];
rz(-1.9638464) q[1];
sx q[1];
rz(-1.9089041) q[1];
x q[2];
rz(1.0001783) q[3];
sx q[3];
rz(-1.6865391) q[3];
sx q[3];
rz(-1.1463349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27590671) q[2];
sx q[2];
rz(-0.17912093) q[2];
sx q[2];
rz(2.5812145) q[2];
rz(3.0748902) q[3];
sx q[3];
rz(-1.3542465) q[3];
sx q[3];
rz(-1.0022479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(0.68503553) q[0];
rz(-0.3512474) q[1];
sx q[1];
rz(-0.63397399) q[1];
sx q[1];
rz(-1.8018855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3303497) q[0];
sx q[0];
rz(-2.0598536) q[0];
sx q[0];
rz(-1.5782209) q[0];
rz(-0.67147246) q[2];
sx q[2];
rz(-1.9118309) q[2];
sx q[2];
rz(1.070553) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3546159) q[1];
sx q[1];
rz(-1.4091345) q[1];
sx q[1];
rz(3.1114108) q[1];
rz(-pi) q[2];
rz(2.4248554) q[3];
sx q[3];
rz(-0.61541688) q[3];
sx q[3];
rz(2.7743266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54905218) q[2];
sx q[2];
rz(-1.668674) q[2];
sx q[2];
rz(0.215691) q[2];
rz(-2.9933764) q[3];
sx q[3];
rz(-2.958478) q[3];
sx q[3];
rz(-0.7557925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6052979) q[0];
sx q[0];
rz(-1.7114102) q[0];
sx q[0];
rz(-0.17912616) q[0];
rz(2.1419549) q[1];
sx q[1];
rz(-0.78796402) q[1];
sx q[1];
rz(-2.9771908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59578204) q[0];
sx q[0];
rz(-2.3124957) q[0];
sx q[0];
rz(-1.2973644) q[0];
rz(-pi) q[1];
rz(-1.4363507) q[2];
sx q[2];
rz(-2.0847581) q[2];
sx q[2];
rz(0.15132893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8021014) q[1];
sx q[1];
rz(-2.4818812) q[1];
sx q[1];
rz(3.129693) q[1];
rz(-1.3429759) q[3];
sx q[3];
rz(-1.9285255) q[3];
sx q[3];
rz(1.2219971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0187692) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(1.4336047) q[2];
rz(-2.150599) q[3];
sx q[3];
rz(-1.3905448) q[3];
sx q[3];
rz(1.2436793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3864415) q[0];
sx q[0];
rz(-0.069510892) q[0];
sx q[0];
rz(2.884927) q[0];
rz(-3.0811884) q[1];
sx q[1];
rz(-0.81788617) q[1];
sx q[1];
rz(-0.43076691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4608972) q[0];
sx q[0];
rz(-0.67047423) q[0];
sx q[0];
rz(-1.1121502) q[0];
x q[1];
rz(2.4469023) q[2];
sx q[2];
rz(-2.3226934) q[2];
sx q[2];
rz(-1.4891877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6182729) q[1];
sx q[1];
rz(-2.225481) q[1];
sx q[1];
rz(1.0275082) q[1];
x q[2];
rz(-0.1448596) q[3];
sx q[3];
rz(-0.88634616) q[3];
sx q[3];
rz(0.057900172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.627219) q[0];
sx q[0];
rz(-0.52424279) q[0];
sx q[0];
rz(0.71994495) q[0];
rz(-0.68079692) q[1];
sx q[1];
rz(-0.38377181) q[1];
sx q[1];
rz(1.0926355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62238111) q[0];
sx q[0];
rz(-1.5977314) q[0];
sx q[0];
rz(-2.6519397) q[0];
rz(-pi) q[1];
rz(1.9377548) q[2];
sx q[2];
rz(-2.1458744) q[2];
sx q[2];
rz(-2.7197011) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7094575) q[1];
sx q[1];
rz(-1.0330519) q[1];
sx q[1];
rz(-0.75790578) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4482139) q[3];
sx q[3];
rz(-1.3170018) q[3];
sx q[3];
rz(-0.13939652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1277658) q[2];
sx q[2];
rz(-2.3769145) q[2];
sx q[2];
rz(2.428425) q[2];
rz(1.7662883) q[3];
sx q[3];
rz(-1.9153374) q[3];
sx q[3];
rz(1.3885952) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-1.0382321) q[1];
sx q[1];
rz(0.33956775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075506181) q[0];
sx q[0];
rz(-1.3097449) q[0];
sx q[0];
rz(0.86782305) q[0];
rz(1.6194862) q[2];
sx q[2];
rz(-1.9867989) q[2];
sx q[2];
rz(-2.4456519) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.133278) q[1];
sx q[1];
rz(-0.61682683) q[1];
sx q[1];
rz(0.99442579) q[1];
rz(-pi) q[2];
rz(-2.5224116) q[3];
sx q[3];
rz(-1.6610489) q[3];
sx q[3];
rz(2.7803286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0833544) q[2];
sx q[2];
rz(-2.3149172) q[2];
sx q[2];
rz(0.44462407) q[2];
rz(-2.5104163) q[3];
sx q[3];
rz(-1.8034214) q[3];
sx q[3];
rz(1.1609424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6889965) q[0];
sx q[0];
rz(-1.8190374) q[0];
sx q[0];
rz(-0.026206503) q[0];
rz(1.7474489) q[1];
sx q[1];
rz(-1.988764) q[1];
sx q[1];
rz(-1.1411512) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8265322) q[0];
sx q[0];
rz(-2.4233343) q[0];
sx q[0];
rz(2.503977) q[0];
rz(0.097392453) q[2];
sx q[2];
rz(-0.44008128) q[2];
sx q[2];
rz(-2.4615364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7375361) q[1];
sx q[1];
rz(-1.0359799) q[1];
sx q[1];
rz(-0.49015518) q[1];
rz(-pi) q[2];
rz(-0.23330063) q[3];
sx q[3];
rz(-2.5590692) q[3];
sx q[3];
rz(1.7482116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.074284) q[2];
sx q[2];
rz(-1.1003541) q[2];
sx q[2];
rz(1.1624153) q[2];
rz(-0.31636604) q[3];
sx q[3];
rz(-1.9497982) q[3];
sx q[3];
rz(1.4828064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83888533) q[0];
sx q[0];
rz(-3.0743297) q[0];
sx q[0];
rz(-0.59826941) q[0];
rz(-3.132598) q[1];
sx q[1];
rz(-0.79575494) q[1];
sx q[1];
rz(1.5473821) q[1];
rz(0.18890201) q[2];
sx q[2];
rz(-0.67086611) q[2];
sx q[2];
rz(-0.98626731) q[2];
rz(-3.0851473) q[3];
sx q[3];
rz(-0.79462449) q[3];
sx q[3];
rz(0.36804646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
