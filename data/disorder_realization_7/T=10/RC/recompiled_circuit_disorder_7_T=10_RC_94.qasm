OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1768567) q[0];
sx q[0];
rz(-0.86369792) q[0];
sx q[0];
rz(0.20110826) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62030828) q[0];
sx q[0];
rz(-1.4518019) q[0];
sx q[0];
rz(1.3814397) q[0];
rz(0.5958545) q[2];
sx q[2];
rz(-2.7239954) q[2];
sx q[2];
rz(0.28996224) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98102346) q[1];
sx q[1];
rz(-1.8796646) q[1];
sx q[1];
rz(-0.94998756) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0073644) q[3];
sx q[3];
rz(-0.3218739) q[3];
sx q[3];
rz(0.25062996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25519249) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24793808) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(0.98518103) q[0];
rz(-2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(3.0342297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8132919) q[0];
sx q[0];
rz(-1.0257226) q[0];
sx q[0];
rz(-2.6585447) q[0];
rz(-2.7254843) q[2];
sx q[2];
rz(-2.5118718) q[2];
sx q[2];
rz(-2.1925418) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7992236) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(2.6745822) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.739005) q[3];
sx q[3];
rz(-1.9544) q[3];
sx q[3];
rz(-0.83354359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(-2.4798933) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(-2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(-0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-1.0823762) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530281) q[0];
sx q[0];
rz(-1.0596501) q[0];
sx q[0];
rz(1.9938064) q[0];
rz(0.15180397) q[2];
sx q[2];
rz(-1.2200583) q[2];
sx q[2];
rz(1.1317859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(-0.2281245) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8941746) q[3];
sx q[3];
rz(-1.7824899) q[3];
sx q[3];
rz(0.64228234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(-0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-2.6534973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49022084) q[0];
sx q[0];
rz(-1.106456) q[0];
sx q[0];
rz(-0.086829348) q[0];
x q[1];
rz(0.79779063) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(0.70300245) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0636053) q[1];
sx q[1];
rz(-1.6997153) q[1];
sx q[1];
rz(1.414826) q[1];
rz(-pi) q[2];
rz(-2.5591764) q[3];
sx q[3];
rz(-1.1629259) q[3];
sx q[3];
rz(-1.8862612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35486832) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.2549843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5421211) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(2.945074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.927658) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(1.2769804) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1097124) q[2];
sx q[2];
rz(-1.5421252) q[2];
sx q[2];
rz(1.920514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6477485) q[1];
sx q[1];
rz(-1.3953679) q[1];
sx q[1];
rz(2.8791048) q[1];
x q[2];
rz(-2.8653141) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1281517) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(0.33368567) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(-1.3938168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48702792) q[0];
sx q[0];
rz(-1.2396887) q[0];
sx q[0];
rz(1.6513255) q[0];
rz(-pi) q[1];
rz(-0.11153531) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(0.44484777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85642568) q[1];
sx q[1];
rz(-1.835583) q[1];
sx q[1];
rz(1.1681639) q[1];
rz(0.49683797) q[3];
sx q[3];
rz(-0.79998575) q[3];
sx q[3];
rz(2.5132688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(1.6585763) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(0.21025118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1649095) q[0];
sx q[0];
rz(-1.7047791) q[0];
sx q[0];
rz(-0.066819013) q[0];
rz(-pi) q[1];
rz(-2.9779151) q[2];
sx q[2];
rz(-1.6357161) q[2];
sx q[2];
rz(0.76771969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5308295) q[1];
sx q[1];
rz(-1.3794583) q[1];
sx q[1];
rz(2.8192239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93123575) q[3];
sx q[3];
rz(-0.82914549) q[3];
sx q[3];
rz(-2.2614546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.64637) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(-2.1648724) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8551089) q[0];
sx q[0];
rz(-0.47874641) q[0];
sx q[0];
rz(0.32740645) q[0];
rz(-pi) q[1];
rz(1.1114242) q[2];
sx q[2];
rz(-0.18006912) q[2];
sx q[2];
rz(-2.0153231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2663914) q[1];
sx q[1];
rz(-0.57386639) q[1];
sx q[1];
rz(-2.2465474) q[1];
rz(2.1920491) q[3];
sx q[3];
rz(-0.81406677) q[3];
sx q[3];
rz(2.4809588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(0.89909536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16222787) q[0];
sx q[0];
rz(-0.75408903) q[0];
sx q[0];
rz(2.2158932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0752418) q[2];
sx q[2];
rz(-2.5214508) q[2];
sx q[2];
rz(-1.3401741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60914674) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(1.2389517) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5702463) q[3];
sx q[3];
rz(-2.7714202) q[3];
sx q[3];
rz(-2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(-2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(2.0237645) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-0.39696473) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871646) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(-2.5973395) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30585395) q[2];
sx q[2];
rz(-1.9216188) q[2];
sx q[2];
rz(0.84301126) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6798903) q[1];
sx q[1];
rz(-1.2352518) q[1];
sx q[1];
rz(3.0479792) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40542094) q[3];
sx q[3];
rz(-2.0339031) q[3];
sx q[3];
rz(2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.94840702) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(3.0449384) q[2];
rz(1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5861355) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(0.10791049) q[2];
sx q[2];
rz(-1.4866598) q[2];
sx q[2];
rz(1.286081) q[2];
rz(1.312064) q[3];
sx q[3];
rz(-1.4553634) q[3];
sx q[3];
rz(2.8298557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];