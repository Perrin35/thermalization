OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(-0.12726769) q[0];
sx q[0];
rz(0.98841086) q[0];
rz(-0.53108162) q[1];
sx q[1];
rz(3.4903033) q[1];
sx q[1];
rz(11.217584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006346) q[0];
sx q[0];
rz(-1.4057584) q[0];
sx q[0];
rz(1.6883231) q[0];
x q[1];
rz(-0.26308194) q[2];
sx q[2];
rz(-1.6533268) q[2];
sx q[2];
rz(-2.6334327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8197644) q[1];
sx q[1];
rz(-1.0294224) q[1];
sx q[1];
rz(-2.2775047) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3380321) q[3];
sx q[3];
rz(-0.27419146) q[3];
sx q[3];
rz(-2.1823332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7621883) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-2.7963426) q[2];
rz(2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(-2.1729443) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-2.7217857) q[0];
rz(0.35821113) q[1];
sx q[1];
rz(-0.63136357) q[1];
sx q[1];
rz(1.0158687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7656454) q[0];
sx q[0];
rz(-0.65212661) q[0];
sx q[0];
rz(-0.69407065) q[0];
rz(-1.3426571) q[2];
sx q[2];
rz(-1.7717138) q[2];
sx q[2];
rz(1.9659496) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8034755) q[1];
sx q[1];
rz(-0.40121597) q[1];
sx q[1];
rz(-2.7944399) q[1];
rz(0.44554168) q[3];
sx q[3];
rz(-1.1098776) q[3];
sx q[3];
rz(1.2156435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35280716) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(2.0324198) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(3.0811908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(2.2614959) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-0.59363669) q[1];
sx q[1];
rz(0.18149158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648375) q[0];
sx q[0];
rz(-2.0252882) q[0];
sx q[0];
rz(0.48671133) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3679738) q[2];
sx q[2];
rz(-0.81589375) q[2];
sx q[2];
rz(-1.1676163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5955334) q[1];
sx q[1];
rz(-1.5469578) q[1];
sx q[1];
rz(1.8438575) q[1];
rz(-0.4811901) q[3];
sx q[3];
rz(-1.7627197) q[3];
sx q[3];
rz(-0.47062518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1168388) q[2];
sx q[2];
rz(-2.0334058) q[2];
sx q[2];
rz(-1.7987569) q[2];
rz(1.3252307) q[3];
sx q[3];
rz(-0.67458761) q[3];
sx q[3];
rz(2.0480115) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6670068) q[0];
sx q[0];
rz(-2.6285567) q[0];
sx q[0];
rz(2.4530607) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-2.7929247) q[1];
sx q[1];
rz(2.9825488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57590398) q[0];
sx q[0];
rz(-0.85907798) q[0];
sx q[0];
rz(0.84337658) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37695388) q[2];
sx q[2];
rz(-1.9184343) q[2];
sx q[2];
rz(0.24573791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-0.025440865) q[1];
rz(-pi) q[2];
rz(-1.3499447) q[3];
sx q[3];
rz(-2.0851118) q[3];
sx q[3];
rz(-0.056003464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29545414) q[2];
sx q[2];
rz(-1.6590154) q[2];
sx q[2];
rz(2.7159178) q[2];
rz(0.019429026) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(0.69452906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6762125) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(0.12839578) q[0];
rz(-2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(2.593186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5655366) q[0];
sx q[0];
rz(-0.57665529) q[0];
sx q[0];
rz(1.499349) q[0];
rz(-pi) q[1];
rz(-0.19210179) q[2];
sx q[2];
rz(-0.77741277) q[2];
sx q[2];
rz(-0.32595134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7171399) q[1];
sx q[1];
rz(-1.5100749) q[1];
sx q[1];
rz(2.0408003) q[1];
rz(-0.10336419) q[3];
sx q[3];
rz(-1.169239) q[3];
sx q[3];
rz(0.38927024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8473062) q[2];
sx q[2];
rz(-0.30964482) q[2];
sx q[2];
rz(1.6748641) q[2];
rz(1.0855801) q[3];
sx q[3];
rz(-1.836136) q[3];
sx q[3];
rz(0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10228957) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(-2.8052203) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(0.99463314) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4431865) q[0];
sx q[0];
rz(-2.5289359) q[0];
sx q[0];
rz(1.1820656) q[0];
rz(3.1349796) q[2];
sx q[2];
rz(-0.50351876) q[2];
sx q[2];
rz(-2.7343482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.57833886) q[1];
sx q[1];
rz(-1.564338) q[1];
sx q[1];
rz(-1.2328813) q[1];
x q[2];
rz(-0.70049882) q[3];
sx q[3];
rz(-1.2555712) q[3];
sx q[3];
rz(-1.3198927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8424592) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(2.3343202) q[3];
sx q[3];
rz(-3.0326796) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.826236) q[0];
sx q[0];
rz(-0.5744136) q[0];
sx q[0];
rz(-2.2454967) q[0];
rz(-0.90944666) q[1];
sx q[1];
rz(-0.25032955) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5298115) q[0];
sx q[0];
rz(-1.5079594) q[0];
sx q[0];
rz(1.2114552) q[0];
rz(-pi) q[1];
rz(-0.80008614) q[2];
sx q[2];
rz(-1.1598088) q[2];
sx q[2];
rz(1.6192186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1753462) q[1];
sx q[1];
rz(-2.4840377) q[1];
sx q[1];
rz(1.46438) q[1];
x q[2];
rz(0.48052629) q[3];
sx q[3];
rz(-0.98283813) q[3];
sx q[3];
rz(-2.5939536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9748777) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(-1.2274851) q[2];
rz(3.098439) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(0.22668049) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-2.7547014) q[0];
sx q[0];
rz(2.7283227) q[0];
rz(2.5832672) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-2.4024898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055187125) q[0];
sx q[0];
rz(-2.1336745) q[0];
sx q[0];
rz(-1.6291717) q[0];
rz(3.1355255) q[2];
sx q[2];
rz(-1.1147611) q[2];
sx q[2];
rz(2.9059682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0308471) q[1];
sx q[1];
rz(-1.4111641) q[1];
sx q[1];
rz(2.2975132) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0855582) q[3];
sx q[3];
rz(-1.8348215) q[3];
sx q[3];
rz(-0.92428401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3070613) q[2];
sx q[2];
rz(-2.8246911) q[2];
sx q[2];
rz(1.1018264) q[2];
rz(2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(2.5226412) q[0];
rz(2.1221819) q[1];
sx q[1];
rz(-1.5588201) q[1];
sx q[1];
rz(0.5272665) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1934803) q[0];
sx q[0];
rz(-1.8475979) q[0];
sx q[0];
rz(0.29159082) q[0];
rz(0.97986603) q[2];
sx q[2];
rz(-0.96989378) q[2];
sx q[2];
rz(0.88496937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7098602) q[1];
sx q[1];
rz(-1.1255956) q[1];
sx q[1];
rz(0.47275895) q[1];
x q[2];
rz(-1.2654952) q[3];
sx q[3];
rz(-0.68317181) q[3];
sx q[3];
rz(2.9340569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8573389) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(-0.81106097) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5321524) q[0];
sx q[0];
rz(-1.0487707) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(2.2562064) q[1];
sx q[1];
rz(-0.65941864) q[1];
sx q[1];
rz(-0.17818174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7893716) q[0];
sx q[0];
rz(-0.86666115) q[0];
sx q[0];
rz(0.88130086) q[0];
rz(0.41583305) q[2];
sx q[2];
rz(-2.4095979) q[2];
sx q[2];
rz(0.40643613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0308932) q[1];
sx q[1];
rz(-0.57037121) q[1];
sx q[1];
rz(-1.13899) q[1];
x q[2];
rz(1.448477) q[3];
sx q[3];
rz(-1.9983851) q[3];
sx q[3];
rz(1.7784255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(-0.82328063) q[2];
rz(-0.12100425) q[3];
sx q[3];
rz(-2.3779317) q[3];
sx q[3];
rz(-2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.9085893) q[0];
sx q[0];
rz(-1.0738666) q[0];
sx q[0];
rz(-0.67169541) q[0];
rz(1.2188777) q[1];
sx q[1];
rz(-1.3296483) q[1];
sx q[1];
rz(1.7760361) q[1];
rz(-1.6505966) q[2];
sx q[2];
rz(-2.7814293) q[2];
sx q[2];
rz(-0.89520892) q[2];
rz(2.8015295) q[3];
sx q[3];
rz(-2.1652514) q[3];
sx q[3];
rz(1.3596331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
