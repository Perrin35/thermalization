OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8837589) q[0];
sx q[0];
rz(-0.45816445) q[0];
sx q[0];
rz(0.56295653) q[0];
rz(0.37748572) q[2];
sx q[2];
rz(-1.1698206) q[2];
sx q[2];
rz(-0.12667835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8909059) q[1];
sx q[1];
rz(-2.2911934) q[1];
sx q[1];
rz(2.5198063) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1629421) q[3];
sx q[3];
rz(-2.1601094) q[3];
sx q[3];
rz(2.4401963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.550094) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(0.47810289) q[2];
rz(-1.6889307) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(-1.0189198) q[0];
rz(-1.4787176) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(0.63308024) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86242005) q[0];
sx q[0];
rz(-2.3775568) q[0];
sx q[0];
rz(-0.11636244) q[0];
rz(-pi) q[1];
rz(-2.3929246) q[2];
sx q[2];
rz(-1.1166755) q[2];
sx q[2];
rz(1.5926966) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2411023) q[1];
sx q[1];
rz(-1.9111269) q[1];
sx q[1];
rz(1.9963005) q[1];
rz(-pi) q[2];
rz(0.31140621) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(-2.1962375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7818266) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(0.11745545) q[2];
rz(-0.30101267) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2837219) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(-0.69082469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72969499) q[0];
sx q[0];
rz(-1.6459961) q[0];
sx q[0];
rz(2.9883283) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76491852) q[2];
sx q[2];
rz(-1.3996482) q[2];
sx q[2];
rz(0.82676065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.55331992) q[1];
sx q[1];
rz(-1.5447445) q[1];
sx q[1];
rz(1.3887029) q[1];
x q[2];
rz(-1.2191804) q[3];
sx q[3];
rz(-1.2120486) q[3];
sx q[3];
rz(0.84695942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2174125) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(3.0351191) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(1.6586554) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(0.52247125) q[0];
rz(2.8126295) q[1];
sx q[1];
rz(-1.4986228) q[1];
sx q[1];
rz(-2.5879588) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5579917) q[0];
sx q[0];
rz(-2.2746673) q[0];
sx q[0];
rz(-0.61917275) q[0];
rz(-pi) q[1];
rz(2.372924) q[2];
sx q[2];
rz(-2.2109291) q[2];
sx q[2];
rz(-1.6878355) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.5987293) q[1];
sx q[1];
rz(-2.3485564) q[1];
sx q[1];
rz(-1.7584156) q[1];
x q[2];
rz(-1.5418566) q[3];
sx q[3];
rz(-2.3674298) q[3];
sx q[3];
rz(2.418747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5840977) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(-0.90302145) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(1.0850798) q[0];
rz(-0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(-0.18403149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407496) q[0];
sx q[0];
rz(-1.4662192) q[0];
sx q[0];
rz(-2.8507289) q[0];
x q[1];
rz(1.8475624) q[2];
sx q[2];
rz(-1.9230611) q[2];
sx q[2];
rz(-2.2500452) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80464333) q[1];
sx q[1];
rz(-2.3579683) q[1];
sx q[1];
rz(-0.6964535) q[1];
x q[2];
rz(-2.4155951) q[3];
sx q[3];
rz(-1.7222002) q[3];
sx q[3];
rz(0.21522537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90157834) q[2];
sx q[2];
rz(-0.34873909) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(0.42282894) q[3];
sx q[3];
rz(-1.4774277) q[3];
sx q[3];
rz(2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35247701) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-2.254803) q[0];
rz(-2.9011762) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(2.9930847) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200955) q[0];
sx q[0];
rz(-1.4098865) q[0];
sx q[0];
rz(-2.4077329) q[0];
x q[1];
rz(-0.22865061) q[2];
sx q[2];
rz(-1.6778523) q[2];
sx q[2];
rz(0.46496898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2827742) q[1];
sx q[1];
rz(-2.8595279) q[1];
sx q[1];
rz(-1.1689405) q[1];
rz(-2.4466483) q[3];
sx q[3];
rz(-0.57999014) q[3];
sx q[3];
rz(-0.65141962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85577661) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(2.6532069) q[2];
rz(-0.48940247) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69328904) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(-2.3235902) q[0];
rz(-1.2524293) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(2.0163527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1624807) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(-0.11066779) q[0];
rz(-pi) q[1];
rz(-2.9863425) q[2];
sx q[2];
rz(-1.8660188) q[2];
sx q[2];
rz(1.0709907) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.352467) q[1];
sx q[1];
rz(-1.2965634) q[1];
sx q[1];
rz(1.8791566) q[1];
rz(0.67894499) q[3];
sx q[3];
rz(-0.79015398) q[3];
sx q[3];
rz(-1.3037579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0685048) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(2.4970064) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.6059426) q[0];
rz(-1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-1.01064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7692524) q[0];
sx q[0];
rz(-1.4108037) q[0];
sx q[0];
rz(0.82323797) q[0];
rz(-pi) q[1];
rz(-1.0755195) q[2];
sx q[2];
rz(-1.3066548) q[2];
sx q[2];
rz(-1.8758945) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9606564) q[1];
sx q[1];
rz(-1.8005383) q[1];
sx q[1];
rz(-0.22682637) q[1];
rz(-0.2145433) q[3];
sx q[3];
rz(-2.7326267) q[3];
sx q[3];
rz(1.4734801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6179787) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(-0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(-0.93769658) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0555608) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6583017) q[0];
sx q[0];
rz(-2.1245983) q[0];
sx q[0];
rz(2.7730586) q[0];
rz(-pi) q[1];
rz(-1.8972626) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(1.9290123) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9675933) q[1];
sx q[1];
rz(-0.6951957) q[1];
sx q[1];
rz(2.980152) q[1];
x q[2];
rz(-0.25200744) q[3];
sx q[3];
rz(-2.5056772) q[3];
sx q[3];
rz(3.0871064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8081234) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(1.8010275) q[2];
rz(2.8373485) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2952404) q[0];
sx q[0];
rz(-1.3051935) q[0];
sx q[0];
rz(-0.17679581) q[0];
rz(-1.2416174) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(2.7005844) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947185) q[0];
sx q[0];
rz(-2.8537769) q[0];
sx q[0];
rz(-2.3527282) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15372865) q[2];
sx q[2];
rz(-1.6445465) q[2];
sx q[2];
rz(-1.7388625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1383789) q[1];
sx q[1];
rz(-2.5795476) q[1];
sx q[1];
rz(0.89488645) q[1];
x q[2];
rz(1.8701843) q[3];
sx q[3];
rz(-0.96794879) q[3];
sx q[3];
rz(-0.30764461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5796154) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(0.20475234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(-3.1148615) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(-1.0529636) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(0.7209575) q[3];
sx q[3];
rz(-2.4506035) q[3];
sx q[3];
rz(-2.639365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
