OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74237139) q[0];
sx q[0];
rz(-1.0065777) q[0];
sx q[0];
rz(0.39682927) q[0];
rz(1.3433427) q[2];
sx q[2];
rz(-1.6709575) q[2];
sx q[2];
rz(1.7286466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5035489) q[1];
sx q[1];
rz(-0.81144864) q[1];
sx q[1];
rz(-1.391295) q[1];
x q[2];
rz(1.1159775) q[3];
sx q[3];
rz(-0.67512073) q[3];
sx q[3];
rz(-1.3779373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33064476) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3246831) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(-1.0634134) q[0];
rz(0.88513199) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(0.0016454776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7749274) q[0];
sx q[0];
rz(-1.6158551) q[0];
sx q[0];
rz(-0.010618322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8029289) q[2];
sx q[2];
rz(-0.67064697) q[2];
sx q[2];
rz(-1.8530958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.590608) q[1];
sx q[1];
rz(-1.3270757) q[1];
sx q[1];
rz(-2.3477712) q[1];
x q[2];
rz(-1.8257636) q[3];
sx q[3];
rz(-1.4062738) q[3];
sx q[3];
rz(-0.69873519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(0.70811159) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-1.089383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3141146) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(-0.88721888) q[0];
rz(-pi) q[1];
rz(-1.1459648) q[2];
sx q[2];
rz(-1.0152738) q[2];
sx q[2];
rz(-2.5512763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.032420302) q[1];
sx q[1];
rz(-1.65616) q[1];
sx q[1];
rz(-1.9568155) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79077625) q[3];
sx q[3];
rz(-1.1547935) q[3];
sx q[3];
rz(-2.5115867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0638782) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(-1.1832773) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.5723715) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-0.35983905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7774178) q[0];
sx q[0];
rz(-0.67876498) q[0];
sx q[0];
rz(-1.8666408) q[0];
x q[1];
rz(-1.4570518) q[2];
sx q[2];
rz(-2.2188088) q[2];
sx q[2];
rz(-1.5485473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6872014) q[1];
sx q[1];
rz(-1.9481716) q[1];
sx q[1];
rz(-2.6881933) q[1];
rz(-1.2218277) q[3];
sx q[3];
rz(-1.8604606) q[3];
sx q[3];
rz(-1.7581913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-2.3763669) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-2.1062772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2238335) q[0];
sx q[0];
rz(-2.1609554) q[0];
sx q[0];
rz(-0.53442861) q[0];
rz(-pi) q[1];
rz(3.0683124) q[2];
sx q[2];
rz(-0.49256941) q[2];
sx q[2];
rz(-2.5626593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.5161533) q[1];
sx q[1];
rz(3.0461237) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4953793) q[3];
sx q[3];
rz(-2.8013902) q[3];
sx q[3];
rz(2.7681805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7845903) q[2];
sx q[2];
rz(-1.4214397) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.6736354) q[3];
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
rz(-2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(-0.2072269) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3375895) q[0];
sx q[0];
rz(-2.1875256) q[0];
sx q[0];
rz(1.3998652) q[0];
rz(-2.0358762) q[2];
sx q[2];
rz(-2.1207223) q[2];
sx q[2];
rz(1.20649) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44240272) q[1];
sx q[1];
rz(-1.3283722) q[1];
sx q[1];
rz(-1.3139903) q[1];
rz(-pi) q[2];
rz(-0.16222555) q[3];
sx q[3];
rz(-2.3174006) q[3];
sx q[3];
rz(-1.970286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(-1.2747814) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(-0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(0.2917372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4371944) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-2.722446) q[0];
rz(0.10642274) q[2];
sx q[2];
rz(-1.8723882) q[2];
sx q[2];
rz(2.0472722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74624324) q[1];
sx q[1];
rz(-2.8031073) q[1];
sx q[1];
rz(2.6836718) q[1];
rz(-pi) q[2];
rz(0.52280207) q[3];
sx q[3];
rz(-1.8409981) q[3];
sx q[3];
rz(-3.0772046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(-0.83731246) q[2];
rz(1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(-0.95867872) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52866919) q[0];
sx q[0];
rz(-2.0818424) q[0];
sx q[0];
rz(-2.5445166) q[0];
x q[1];
rz(1.0326951) q[2];
sx q[2];
rz(-1.8128464) q[2];
sx q[2];
rz(-0.48697105) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9740323) q[1];
sx q[1];
rz(-1.4640199) q[1];
sx q[1];
rz(1.0436996) q[1];
rz(-1.6931375) q[3];
sx q[3];
rz(-2.1515176) q[3];
sx q[3];
rz(-1.3337222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.133193) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(-2.6760496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31637329) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(0.70789106) q[0];
rz(0.12840694) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(0.94305925) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0034345) q[1];
sx q[1];
rz(-2.5596566) q[1];
sx q[1];
rz(2.483063) q[1];
rz(0.25119987) q[3];
sx q[3];
rz(-1.2217055) q[3];
sx q[3];
rz(-2.420345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7312701) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(2.5591992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1727027) q[0];
sx q[0];
rz(-0.70969289) q[0];
sx q[0];
rz(-1.2560647) q[0];
rz(-pi) q[1];
rz(1.8495314) q[2];
sx q[2];
rz(-2.6346452) q[2];
sx q[2];
rz(2.7915733) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6813587) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(1.387499) q[1];
rz(-pi) q[2];
rz(-3.0155229) q[3];
sx q[3];
rz(-2.9226916) q[3];
sx q[3];
rz(-2.253988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762155) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.8021884) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(1.8680686) q[2];
sx q[2];
rz(-2.3625629) q[2];
sx q[2];
rz(0.071803781) q[2];
rz(-2.6639832) q[3];
sx q[3];
rz(-1.1160679) q[3];
sx q[3];
rz(3.06649) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
