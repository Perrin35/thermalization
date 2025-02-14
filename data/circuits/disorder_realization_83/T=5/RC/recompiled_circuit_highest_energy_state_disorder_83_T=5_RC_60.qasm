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
rz(-1.8981847) q[0];
sx q[0];
rz(-1.564448) q[0];
sx q[0];
rz(-0.91140437) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(3.4179847) q[1];
sx q[1];
rz(8.5951947) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4938717) q[0];
sx q[0];
rz(-0.52807489) q[0];
sx q[0];
rz(2.2654669) q[0];
x q[1];
rz(-0.24721036) q[2];
sx q[2];
rz(-1.6610378) q[2];
sx q[2];
rz(-1.5238786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0427057) q[1];
sx q[1];
rz(-1.2707842) q[1];
sx q[1];
rz(-2.3753931) q[1];
x q[2];
rz(-1.4178278) q[3];
sx q[3];
rz(-2.0024096) q[3];
sx q[3];
rz(0.6237517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36898819) q[2];
sx q[2];
rz(-0.96869865) q[2];
sx q[2];
rz(-1.9770835) q[2];
rz(-0.72689593) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(-3.0708142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0182536) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(2.8501999) q[0];
rz(0.60028752) q[1];
sx q[1];
rz(-0.95708668) q[1];
sx q[1];
rz(-0.16500638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24848973) q[0];
sx q[0];
rz(-1.4997673) q[0];
sx q[0];
rz(-1.8336589) q[0];
rz(-2.0197312) q[2];
sx q[2];
rz(-1.672272) q[2];
sx q[2];
rz(1.6368653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2875322) q[1];
sx q[1];
rz(-0.89193501) q[1];
sx q[1];
rz(1.4816585) q[1];
rz(-pi) q[2];
rz(-0.21221186) q[3];
sx q[3];
rz(-0.2538052) q[3];
sx q[3];
rz(2.5830944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0773641) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(-2.9020818) q[2];
rz(0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-0.094836205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25642446) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(1.5688131) q[0];
rz(-0.39372152) q[1];
sx q[1];
rz(-0.55501333) q[1];
sx q[1];
rz(-2.6655925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5167907) q[0];
sx q[0];
rz(-1.778852) q[0];
sx q[0];
rz(-2.1387908) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7287929) q[2];
sx q[2];
rz(-2.9523628) q[2];
sx q[2];
rz(0.10052027) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4502628) q[1];
sx q[1];
rz(-0.99585497) q[1];
sx q[1];
rz(0.49552576) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8937821) q[3];
sx q[3];
rz(-1.8988639) q[3];
sx q[3];
rz(-2.1472665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8126882) q[2];
sx q[2];
rz(-0.71949553) q[2];
sx q[2];
rz(-0.94181124) q[2];
rz(1.4224667) q[3];
sx q[3];
rz(-0.37426451) q[3];
sx q[3];
rz(-0.67371887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.5330676) q[0];
rz(-2.7948921) q[1];
sx q[1];
rz(-1.0418714) q[1];
sx q[1];
rz(0.18992058) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3118211) q[0];
sx q[0];
rz(-1.30952) q[0];
sx q[0];
rz(2.3992062) q[0];
x q[1];
rz(-1.0937017) q[2];
sx q[2];
rz(-1.3253085) q[2];
sx q[2];
rz(-2.7141822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.080481361) q[1];
sx q[1];
rz(-0.60855344) q[1];
sx q[1];
rz(1.670865) q[1];
rz(2.3307822) q[3];
sx q[3];
rz(-0.56104198) q[3];
sx q[3];
rz(1.1277792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7803663) q[2];
sx q[2];
rz(-0.85005886) q[2];
sx q[2];
rz(0.38193199) q[2];
rz(0.10968883) q[3];
sx q[3];
rz(-2.1083125) q[3];
sx q[3];
rz(0.1194574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1031951) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(2.2947445) q[0];
rz(-3.0021744) q[1];
sx q[1];
rz(-0.81652313) q[1];
sx q[1];
rz(-2.3925508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3232305) q[0];
sx q[0];
rz(-1.7182171) q[0];
sx q[0];
rz(2.4458103) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7125569) q[2];
sx q[2];
rz(-1.0899915) q[2];
sx q[2];
rz(-0.17021179) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3440123) q[1];
sx q[1];
rz(-1.7645807) q[1];
sx q[1];
rz(1.7001726) q[1];
rz(-2.0723125) q[3];
sx q[3];
rz(-1.4677047) q[3];
sx q[3];
rz(3.0018798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9355115) q[2];
sx q[2];
rz(-1.8847909) q[2];
sx q[2];
rz(-2.0743267) q[2];
rz(-0.71077985) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(-1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34273219) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(2.1730098) q[0];
rz(-2.4505278) q[1];
sx q[1];
rz(-1.3422809) q[1];
sx q[1];
rz(1.8341281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027031487) q[0];
sx q[0];
rz(-1.8252715) q[0];
sx q[0];
rz(0.64395321) q[0];
rz(0.2910462) q[2];
sx q[2];
rz(-1.7564764) q[2];
sx q[2];
rz(-3.1094375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1561574) q[1];
sx q[1];
rz(-1.5181386) q[1];
sx q[1];
rz(1.1833722) q[1];
x q[2];
rz(-0.97882459) q[3];
sx q[3];
rz(-0.32099629) q[3];
sx q[3];
rz(2.5227566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-2.9066777) q[2];
sx q[2];
rz(-1.2972181) q[2];
rz(-1.3951067) q[3];
sx q[3];
rz(-1.8078943) q[3];
sx q[3];
rz(1.5313799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7995826) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(0.94014257) q[0];
rz(2.4932585) q[1];
sx q[1];
rz(-1.2038566) q[1];
sx q[1];
rz(2.8727093) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20737442) q[0];
sx q[0];
rz(-1.8158127) q[0];
sx q[0];
rz(2.5859358) q[0];
x q[1];
rz(0.36264561) q[2];
sx q[2];
rz(-2.4055436) q[2];
sx q[2];
rz(-2.6578052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2317048) q[1];
sx q[1];
rz(-2.1429166) q[1];
sx q[1];
rz(-0.27426274) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4869789) q[3];
sx q[3];
rz(-0.5236434) q[3];
sx q[3];
rz(2.3414827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0073283422) q[2];
sx q[2];
rz(-1.2790044) q[2];
sx q[2];
rz(2.4556665) q[2];
rz(-2.4053597) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(-0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.69970423) q[0];
sx q[0];
rz(-1.2208953) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(-0.092983149) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(-2.0704796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37866805) q[0];
sx q[0];
rz(-1.0707741) q[0];
sx q[0];
rz(-2.0471694) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1125715) q[2];
sx q[2];
rz(-2.1734218) q[2];
sx q[2];
rz(-0.17652179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6735355) q[1];
sx q[1];
rz(-2.2186465) q[1];
sx q[1];
rz(-0.38783973) q[1];
rz(-pi) q[2];
rz(0.052196189) q[3];
sx q[3];
rz(-1.9813271) q[3];
sx q[3];
rz(-0.13971381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45212713) q[2];
sx q[2];
rz(-2.5048544) q[2];
sx q[2];
rz(2.9316736) q[2];
rz(3.012015) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(-1.5642222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.121948) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(-3.1280532) q[0];
rz(-0.53567046) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(-0.013669107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22937742) q[0];
sx q[0];
rz(-1.7969062) q[0];
sx q[0];
rz(-2.3588728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0837567) q[2];
sx q[2];
rz(-1.221861) q[2];
sx q[2];
rz(1.2020122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0993578) q[1];
sx q[1];
rz(-2.3546948) q[1];
sx q[1];
rz(-2.6762647) q[1];
rz(-pi) q[2];
rz(0.39710703) q[3];
sx q[3];
rz(-0.79063168) q[3];
sx q[3];
rz(1.9472029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(-1.7784485) q[2];
rz(2.3220883) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118947) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(-1.6534506) q[0];
rz(-2.7986774) q[1];
sx q[1];
rz(-1.5577134) q[1];
sx q[1];
rz(-0.75103474) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2011178) q[0];
sx q[0];
rz(-1.6342666) q[0];
sx q[0];
rz(-1.2918143) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.033634) q[2];
sx q[2];
rz(-1.7506536) q[2];
sx q[2];
rz(3.0015869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0768834) q[1];
sx q[1];
rz(-0.61184498) q[1];
sx q[1];
rz(-2.9348899) q[1];
rz(1.9225227) q[3];
sx q[3];
rz(-0.096032525) q[3];
sx q[3];
rz(2.072123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.670383) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(2.75441) q[2];
rz(0.47532982) q[3];
sx q[3];
rz(-1.1673704) q[3];
sx q[3];
rz(-0.79132426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8925856) q[0];
sx q[0];
rz(-0.96207608) q[0];
sx q[0];
rz(0.67208653) q[0];
rz(-0.36920209) q[1];
sx q[1];
rz(-2.3875356) q[1];
sx q[1];
rz(1.748132) q[1];
rz(2.4245928) q[2];
sx q[2];
rz(-1.430027) q[2];
sx q[2];
rz(1.2363557) q[2];
rz(1.9561097) q[3];
sx q[3];
rz(-0.29373071) q[3];
sx q[3];
rz(3.1391524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
