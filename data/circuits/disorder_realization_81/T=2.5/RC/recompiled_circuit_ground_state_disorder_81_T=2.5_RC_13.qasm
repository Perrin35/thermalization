OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17570198) q[0];
sx q[0];
rz(5.0134563) q[0];
sx q[0];
rz(11.391517) q[0];
rz(-1.084561) q[1];
sx q[1];
rz(-1.4718055) q[1];
sx q[1];
rz(0.56301277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5064237) q[0];
sx q[0];
rz(-0.20316589) q[0];
sx q[0];
rz(0.26655339) q[0];
rz(-pi) q[1];
rz(2.6649551) q[2];
sx q[2];
rz(-2.1486006) q[2];
sx q[2];
rz(-3.1057242) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35511103) q[1];
sx q[1];
rz(-1.4110089) q[1];
sx q[1];
rz(0.33336877) q[1];
rz(-pi) q[2];
rz(1.9292163) q[3];
sx q[3];
rz(-1.2120768) q[3];
sx q[3];
rz(-1.4354853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53434831) q[2];
sx q[2];
rz(-1.4899985) q[2];
sx q[2];
rz(-1.6687757) q[2];
rz(1.2502753) q[3];
sx q[3];
rz(-1.6142802) q[3];
sx q[3];
rz(0.97243029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2431353) q[0];
sx q[0];
rz(-3.0861096) q[0];
sx q[0];
rz(-2.6836416) q[0];
rz(3.0398439) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(-0.32624689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889623) q[0];
sx q[0];
rz(-1.9174308) q[0];
sx q[0];
rz(1.1244357) q[0];
x q[1];
rz(-2.5102551) q[2];
sx q[2];
rz(-0.36383087) q[2];
sx q[2];
rz(-1.199786) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65349039) q[1];
sx q[1];
rz(-2.2970169) q[1];
sx q[1];
rz(-0.86700551) q[1];
x q[2];
rz(-2.3986048) q[3];
sx q[3];
rz(-1.8077733) q[3];
sx q[3];
rz(1.5239511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.132167) q[2];
sx q[2];
rz(-0.99401179) q[2];
sx q[2];
rz(-1.8780635) q[2];
rz(-3.1215014) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(1.9115537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024071368) q[0];
sx q[0];
rz(-0.063952359) q[0];
sx q[0];
rz(-0.57620311) q[0];
rz(1.6974712) q[1];
sx q[1];
rz(-1.2497808) q[1];
sx q[1];
rz(-1.8796657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5844628) q[0];
sx q[0];
rz(-0.31079159) q[0];
sx q[0];
rz(-1.4660828) q[0];
x q[1];
rz(2.0327225) q[2];
sx q[2];
rz(-1.0903768) q[2];
sx q[2];
rz(-1.9146718) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.539725) q[1];
sx q[1];
rz(-0.66942184) q[1];
sx q[1];
rz(-0.4008534) q[1];
x q[2];
rz(-2.5188082) q[3];
sx q[3];
rz(-0.81777621) q[3];
sx q[3];
rz(-1.0914413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.55804092) q[2];
sx q[2];
rz(-0.59528196) q[2];
sx q[2];
rz(2.9884647) q[2];
rz(2.2554452) q[3];
sx q[3];
rz(-1.359442) q[3];
sx q[3];
rz(-1.2019633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.8837638) q[0];
sx q[0];
rz(-1.184329) q[0];
sx q[0];
rz(-2.9828239) q[0];
rz(-2.1067045) q[1];
sx q[1];
rz(-0.95304573) q[1];
sx q[1];
rz(-2.3854756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0466111) q[0];
sx q[0];
rz(-1.5205407) q[0];
sx q[0];
rz(0.057970164) q[0];
x q[1];
rz(0.027534187) q[2];
sx q[2];
rz(-1.3176271) q[2];
sx q[2];
rz(2.262676) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19263527) q[1];
sx q[1];
rz(-1.7757984) q[1];
sx q[1];
rz(-2.8947049) q[1];
rz(1.7773553) q[3];
sx q[3];
rz(-1.8011923) q[3];
sx q[3];
rz(1.2588866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5990344) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(0.8503882) q[2];
rz(-1.6709857) q[3];
sx q[3];
rz(-1.3266404) q[3];
sx q[3];
rz(0.82408041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4849512) q[0];
sx q[0];
rz(-1.739772) q[0];
sx q[0];
rz(-2.9630419) q[0];
rz(1.8932331) q[1];
sx q[1];
rz(-1.6105885) q[1];
sx q[1];
rz(0.53370968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2634791) q[0];
sx q[0];
rz(-1.3878667) q[0];
sx q[0];
rz(-2.7862552) q[0];
rz(-pi) q[1];
rz(-2.6421946) q[2];
sx q[2];
rz(-0.93448567) q[2];
sx q[2];
rz(0.15313521) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0438178) q[1];
sx q[1];
rz(-2.2378451) q[1];
sx q[1];
rz(-0.066867877) q[1];
rz(-pi) q[2];
rz(-1.0664135) q[3];
sx q[3];
rz(-1.4850067) q[3];
sx q[3];
rz(-2.498621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3137714) q[2];
sx q[2];
rz(-1.8905819) q[2];
sx q[2];
rz(-1.0531462) q[2];
rz(-0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(-2.7789796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5617705) q[0];
sx q[0];
rz(-2.1403911) q[0];
sx q[0];
rz(-1.9129821) q[0];
rz(-0.46663943) q[1];
sx q[1];
rz(-1.2184315) q[1];
sx q[1];
rz(-0.12900464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056696293) q[0];
sx q[0];
rz(-1.4375303) q[0];
sx q[0];
rz(0.10198089) q[0];
rz(-pi) q[1];
rz(-0.1618218) q[2];
sx q[2];
rz(-0.84748778) q[2];
sx q[2];
rz(-0.79297334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1946757) q[1];
sx q[1];
rz(-1.2763775) q[1];
sx q[1];
rz(-2.5076082) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13997649) q[3];
sx q[3];
rz(-1.8365321) q[3];
sx q[3];
rz(-1.99972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0064319) q[2];
sx q[2];
rz(-2.5401523) q[2];
sx q[2];
rz(2.637291) q[2];
rz(-2.048061) q[3];
sx q[3];
rz(-1.3239599) q[3];
sx q[3];
rz(2.1849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4920376) q[0];
sx q[0];
rz(-1.2832337) q[0];
sx q[0];
rz(-1.5851703) q[0];
rz(-0.57836142) q[1];
sx q[1];
rz(-1.2588986) q[1];
sx q[1];
rz(3.0392821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8292885) q[0];
sx q[0];
rz(-1.1745546) q[0];
sx q[0];
rz(0.90950729) q[0];
rz(-pi) q[1];
rz(1.341171) q[2];
sx q[2];
rz(-1.183325) q[2];
sx q[2];
rz(-1.6796527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7150517) q[1];
sx q[1];
rz(-0.3449966) q[1];
sx q[1];
rz(-0.47899063) q[1];
x q[2];
rz(-2.3915406) q[3];
sx q[3];
rz(-2.1382209) q[3];
sx q[3];
rz(-0.22818434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9461225) q[2];
sx q[2];
rz(-1.2637694) q[2];
sx q[2];
rz(0.61402399) q[2];
rz(-2.2139003) q[3];
sx q[3];
rz(-1.5427019) q[3];
sx q[3];
rz(-2.4641002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880661) q[0];
sx q[0];
rz(-2.8480242) q[0];
sx q[0];
rz(1.199383) q[0];
rz(-1.9623494) q[1];
sx q[1];
rz(-1.0284547) q[1];
sx q[1];
rz(-0.95338043) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3971553) q[0];
sx q[0];
rz(-2.1040475) q[0];
sx q[0];
rz(2.5313733) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.751975) q[2];
sx q[2];
rz(-2.5963514) q[2];
sx q[2];
rz(-0.94961005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2732716) q[1];
sx q[1];
rz(-1.4605323) q[1];
sx q[1];
rz(0.83302193) q[1];
rz(-pi) q[2];
rz(1.6497506) q[3];
sx q[3];
rz(-2.051226) q[3];
sx q[3];
rz(-1.3955658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.392268) q[2];
sx q[2];
rz(-0.12568036) q[2];
sx q[2];
rz(2.9607062) q[2];
rz(-1.1130029) q[3];
sx q[3];
rz(-2.1235316) q[3];
sx q[3];
rz(-0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0522633) q[0];
sx q[0];
rz(-2.2398529) q[0];
sx q[0];
rz(-0.7641291) q[0];
rz(1.0230505) q[1];
sx q[1];
rz(-1.5307531) q[1];
sx q[1];
rz(1.6463564) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7376855) q[0];
sx q[0];
rz(-1.4171184) q[0];
sx q[0];
rz(2.840848) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1713397) q[2];
sx q[2];
rz(-1.8489328) q[2];
sx q[2];
rz(-0.32626611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0937522) q[1];
sx q[1];
rz(-1.4994267) q[1];
sx q[1];
rz(-2.755295) q[1];
rz(-pi) q[2];
rz(-0.71919433) q[3];
sx q[3];
rz(-2.580504) q[3];
sx q[3];
rz(0.51855519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.19227795) q[2];
sx q[2];
rz(-2.1239026) q[2];
sx q[2];
rz(2.6940572) q[2];
rz(3.0053075) q[3];
sx q[3];
rz(-0.49301967) q[3];
sx q[3];
rz(-0.59806699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.3560781) q[0];
sx q[0];
rz(-2.1165753) q[0];
sx q[0];
rz(0.4726952) q[0];
rz(-2.7418432) q[1];
sx q[1];
rz(-2.4374297) q[1];
sx q[1];
rz(2.3230816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2855354) q[0];
sx q[0];
rz(-1.3926848) q[0];
sx q[0];
rz(2.8277603) q[0];
rz(-1.1498951) q[2];
sx q[2];
rz(-2.2031261) q[2];
sx q[2];
rz(2.9078751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1401891) q[1];
sx q[1];
rz(-1.8934939) q[1];
sx q[1];
rz(-3.138209) q[1];
rz(-pi) q[2];
rz(2.6916041) q[3];
sx q[3];
rz(-2.7594824) q[3];
sx q[3];
rz(1.1799605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5588536) q[2];
sx q[2];
rz(-1.0820729) q[2];
sx q[2];
rz(-1.2782798) q[2];
rz(-1.1996972) q[3];
sx q[3];
rz(-1.7007098) q[3];
sx q[3];
rz(-0.27967683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77547913) q[0];
sx q[0];
rz(-1.4237325) q[0];
sx q[0];
rz(1.7970418) q[0];
rz(1.4115502) q[1];
sx q[1];
rz(-2.330214) q[1];
sx q[1];
rz(-0.65705962) q[1];
rz(3.1362202) q[2];
sx q[2];
rz(-3.0006144) q[2];
sx q[2];
rz(-0.49868546) q[2];
rz(2.2349002) q[3];
sx q[3];
rz(-1.8421052) q[3];
sx q[3];
rz(2.6143034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
