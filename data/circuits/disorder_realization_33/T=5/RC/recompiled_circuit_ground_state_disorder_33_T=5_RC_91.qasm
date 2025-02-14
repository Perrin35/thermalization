OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6991718) q[0];
sx q[0];
rz(-0.81096634) q[0];
sx q[0];
rz(0.45642689) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(-0.95093095) q[1];
sx q[1];
rz(2.9589597) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5125677) q[0];
sx q[0];
rz(-0.95172608) q[0];
sx q[0];
rz(-2.7946212) q[0];
x q[1];
rz(2.6668197) q[2];
sx q[2];
rz(-2.2150196) q[2];
sx q[2];
rz(-0.34851532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7108954) q[1];
sx q[1];
rz(-2.9646246) q[1];
sx q[1];
rz(-2.9269993) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51974082) q[3];
sx q[3];
rz(-1.2260574) q[3];
sx q[3];
rz(2.910579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7543588) q[2];
sx q[2];
rz(-1.0789472) q[2];
sx q[2];
rz(0.31164247) q[2];
rz(-1.8841057) q[3];
sx q[3];
rz(-2.8897372) q[3];
sx q[3];
rz(2.9529412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1455014) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(1.2392932) q[0];
rz(2.7481825) q[1];
sx q[1];
rz(-1.0478123) q[1];
sx q[1];
rz(2.1107103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6604583) q[0];
sx q[0];
rz(-1.0054614) q[0];
sx q[0];
rz(2.7271366) q[0];
rz(-0.71218636) q[2];
sx q[2];
rz(-2.7982494) q[2];
sx q[2];
rz(0.58420974) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46566612) q[1];
sx q[1];
rz(-0.74518004) q[1];
sx q[1];
rz(2.096677) q[1];
rz(-pi) q[2];
rz(2.5000948) q[3];
sx q[3];
rz(-2.201295) q[3];
sx q[3];
rz(1.8127831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3139412) q[2];
sx q[2];
rz(-2.2590019) q[2];
sx q[2];
rz(2.7871056) q[2];
rz(-0.83141023) q[3];
sx q[3];
rz(-0.76787132) q[3];
sx q[3];
rz(-1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77446929) q[0];
sx q[0];
rz(-0.18018436) q[0];
sx q[0];
rz(-0.48429504) q[0];
rz(0.88227415) q[1];
sx q[1];
rz(-0.87119281) q[1];
sx q[1];
rz(0.54214111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061848693) q[0];
sx q[0];
rz(-1.2164017) q[0];
sx q[0];
rz(0.3318048) q[0];
rz(-pi) q[1];
rz(-2.9464821) q[2];
sx q[2];
rz(-0.32719041) q[2];
sx q[2];
rz(1.3766152) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24607813) q[1];
sx q[1];
rz(-0.82209229) q[1];
sx q[1];
rz(0.24242927) q[1];
x q[2];
rz(0.063850689) q[3];
sx q[3];
rz(-1.8076767) q[3];
sx q[3];
rz(-0.76343918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.813628) q[2];
sx q[2];
rz(-2.4317604) q[2];
sx q[2];
rz(-1.0343118) q[2];
rz(-1.7799001) q[3];
sx q[3];
rz(-1.9618278) q[3];
sx q[3];
rz(-2.5907607) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4807602) q[0];
sx q[0];
rz(-1.7518504) q[0];
sx q[0];
rz(1.047629) q[0];
rz(-1.818559) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(0.38527647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7390201) q[0];
sx q[0];
rz(-2.077335) q[0];
sx q[0];
rz(2.6468524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4249737) q[2];
sx q[2];
rz(-1.5138549) q[2];
sx q[2];
rz(-2.5927134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2400018) q[1];
sx q[1];
rz(-1.4996254) q[1];
sx q[1];
rz(-2.932697) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0317467) q[3];
sx q[3];
rz(-2.2513736) q[3];
sx q[3];
rz(2.3504013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4177527) q[2];
sx q[2];
rz(-2.1989792) q[2];
sx q[2];
rz(2.8670132) q[2];
rz(0.56139055) q[3];
sx q[3];
rz(-2.0761108) q[3];
sx q[3];
rz(1.5076465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0059589) q[0];
sx q[0];
rz(-2.5649286) q[0];
sx q[0];
rz(2.2744001) q[0];
rz(-2.7658956) q[1];
sx q[1];
rz(-2.5521894) q[1];
sx q[1];
rz(1.9120749) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4081683) q[0];
sx q[0];
rz(-0.46981341) q[0];
sx q[0];
rz(0.3656268) q[0];
rz(-2.0887689) q[2];
sx q[2];
rz(-1.8395632) q[2];
sx q[2];
rz(-2.4160699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7244959) q[1];
sx q[1];
rz(-2.1388106) q[1];
sx q[1];
rz(2.591631) q[1];
rz(0.33809148) q[3];
sx q[3];
rz(-1.9623358) q[3];
sx q[3];
rz(-3.1126659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1162954) q[2];
sx q[2];
rz(-1.2325492) q[2];
sx q[2];
rz(-0.06289014) q[2];
rz(2.7316015) q[3];
sx q[3];
rz(-0.72628179) q[3];
sx q[3];
rz(-0.15967742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1493688) q[0];
sx q[0];
rz(-3.0719482) q[0];
sx q[0];
rz(0.83576354) q[0];
rz(1.958485) q[1];
sx q[1];
rz(-1.8712021) q[1];
sx q[1];
rz(-0.74367181) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27722699) q[0];
sx q[0];
rz(-1.0334618) q[0];
sx q[0];
rz(1.2173843) q[0];
rz(3.0871687) q[2];
sx q[2];
rz(-1.7147002) q[2];
sx q[2];
rz(-1.2765034) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3144296) q[1];
sx q[1];
rz(-1.3884228) q[1];
sx q[1];
rz(-0.86079396) q[1];
rz(0.20009508) q[3];
sx q[3];
rz(-0.26069122) q[3];
sx q[3];
rz(-0.30567769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30367294) q[2];
sx q[2];
rz(-2.643879) q[2];
sx q[2];
rz(-2.3480603) q[2];
rz(-2.7225336) q[3];
sx q[3];
rz(-1.4895118) q[3];
sx q[3];
rz(-0.52217531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1917052) q[0];
sx q[0];
rz(-0.97923034) q[0];
sx q[0];
rz(1.1631843) q[0];
rz(-0.78041068) q[1];
sx q[1];
rz(-0.33640948) q[1];
sx q[1];
rz(-1.6132678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53188092) q[0];
sx q[0];
rz(-1.6356633) q[0];
sx q[0];
rz(2.2676629) q[0];
rz(-2.0397908) q[2];
sx q[2];
rz(-0.38809478) q[2];
sx q[2];
rz(-2.3562252) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1133729) q[1];
sx q[1];
rz(-1.7689147) q[1];
sx q[1];
rz(-3.1199725) q[1];
x q[2];
rz(-2.1329857) q[3];
sx q[3];
rz(-0.48848029) q[3];
sx q[3];
rz(-0.83399663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10716001) q[2];
sx q[2];
rz(-1.9623423) q[2];
sx q[2];
rz(0.068664702) q[2];
rz(2.5717403) q[3];
sx q[3];
rz(-0.48044258) q[3];
sx q[3];
rz(0.3705875) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079085199) q[0];
sx q[0];
rz(-0.67254368) q[0];
sx q[0];
rz(-1.3154718) q[0];
rz(1.2830118) q[1];
sx q[1];
rz(-2.7048769) q[1];
sx q[1];
rz(-3.0924996) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0593215) q[0];
sx q[0];
rz(-1.4937287) q[0];
sx q[0];
rz(-1.4801816) q[0];
x q[1];
rz(0.95048381) q[2];
sx q[2];
rz(-2.0237708) q[2];
sx q[2];
rz(-0.17826232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80688991) q[1];
sx q[1];
rz(-0.8262615) q[1];
sx q[1];
rz(-2.101154) q[1];
rz(-pi) q[2];
rz(-1.4021644) q[3];
sx q[3];
rz(-2.4908427) q[3];
sx q[3];
rz(-0.73413056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38356885) q[2];
sx q[2];
rz(-2.263676) q[2];
sx q[2];
rz(0.54086584) q[2];
rz(1.0848378) q[3];
sx q[3];
rz(-2.4386051) q[3];
sx q[3];
rz(-1.7613523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.908602) q[0];
sx q[0];
rz(-0.99642307) q[0];
sx q[0];
rz(-0.10502271) q[0];
rz(-2.5999293) q[1];
sx q[1];
rz(-2.2553406) q[1];
sx q[1];
rz(-0.35596102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.134577) q[0];
sx q[0];
rz(-1.4316214) q[0];
sx q[0];
rz(0.037875847) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85552026) q[2];
sx q[2];
rz(-2.2770555) q[2];
sx q[2];
rz(0.17937961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4222504) q[1];
sx q[1];
rz(-1.0021035) q[1];
sx q[1];
rz(-0.45180068) q[1];
x q[2];
rz(-1.2508873) q[3];
sx q[3];
rz(-1.4540795) q[3];
sx q[3];
rz(3.0206783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2124704) q[2];
sx q[2];
rz(-1.4996303) q[2];
sx q[2];
rz(1.8222202) q[2];
rz(-1.8170554) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(2.1738079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.20973715) q[0];
sx q[0];
rz(-3.1071438) q[0];
sx q[0];
rz(1.4631648) q[0];
rz(-1.6819008) q[1];
sx q[1];
rz(-1.1594783) q[1];
sx q[1];
rz(2.7957338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5992085) q[0];
sx q[0];
rz(-3.0488692) q[0];
sx q[0];
rz(2.1044162) q[0];
rz(-pi) q[1];
rz(-2.9458617) q[2];
sx q[2];
rz(-2.8189341) q[2];
sx q[2];
rz(-3.0374073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1766369) q[1];
sx q[1];
rz(-2.0207016) q[1];
sx q[1];
rz(-0.69590203) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8000051) q[3];
sx q[3];
rz(-0.95529592) q[3];
sx q[3];
rz(-3.0581491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2807002) q[2];
sx q[2];
rz(-1.2714551) q[2];
sx q[2];
rz(0.26091179) q[2];
rz(-1.6711309) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(1.7117333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3093001) q[0];
sx q[0];
rz(-1.1080879) q[0];
sx q[0];
rz(2.9453887) q[0];
rz(2.590754) q[1];
sx q[1];
rz(-1.486634) q[1];
sx q[1];
rz(-2.6015729) q[1];
rz(1.3030686) q[2];
sx q[2];
rz(-2.3961551) q[2];
sx q[2];
rz(1.2533631) q[2];
rz(-2.4110386) q[3];
sx q[3];
rz(-1.3638221) q[3];
sx q[3];
rz(2.4331349) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
