OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(2.9404844) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1683567) q[0];
sx q[0];
rz(-1.7587979) q[0];
sx q[0];
rz(0.12113916) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7896499) q[2];
sx q[2];
rz(-1.3411739) q[2];
sx q[2];
rz(-2.4156092) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98102346) q[1];
sx q[1];
rz(-1.261928) q[1];
sx q[1];
rz(0.94998756) q[1];
x q[2];
rz(-2.1342282) q[3];
sx q[3];
rz(-2.8197188) q[3];
sx q[3];
rz(-0.25062996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(1.1039929) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(2.7424116) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023022368) q[0];
sx q[0];
rz(-1.9792299) q[0];
sx q[0];
rz(-2.1711711) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41610833) q[2];
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
rz(0.7992236) q[1];
sx q[1];
rz(-0.93042004) q[1];
sx q[1];
rz(0.46701042) q[1];
x q[2];
rz(-2.7530305) q[3];
sx q[3];
rz(-1.7266759) q[3];
sx q[3];
rz(0.673783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5169516) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(2.4798933) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(-3.1012428) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28856453) q[0];
sx q[0];
rz(-2.0819426) q[0];
sx q[0];
rz(-1.1477863) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15180397) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(-1.1317859) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97835097) q[1];
sx q[1];
rz(-2.8132073) q[1];
sx q[1];
rz(-2.3204625) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8941746) q[3];
sx q[3];
rz(-1.3591027) q[3];
sx q[3];
rz(-0.64228234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-2.5877) q[2];
rz(1.6484377) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4937113) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(0.035382263) q[0];
rz(-2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(0.48809537) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49022084) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(3.0547633) q[0];
rz(-pi) q[1];
rz(-2.5627665) q[2];
sx q[2];
rz(-1.0597214) q[2];
sx q[2];
rz(1.6312815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6285703) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(3.0111074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.475708) q[3];
sx q[3];
rz(-0.69722359) q[3];
sx q[3];
rz(-0.22660412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
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
rz(-1.5636684) q[1];
sx q[1];
rz(-2.945074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.927658) q[0];
sx q[0];
rz(-1.8727344) q[0];
sx q[0];
rz(1.2769804) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6351661) q[2];
sx q[2];
rz(-0.46191051) q[2];
sx q[2];
rz(-0.40735746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6477485) q[1];
sx q[1];
rz(-1.3953679) q[1];
sx q[1];
rz(-0.26248787) q[1];
rz(-pi) q[2];
rz(2.5971562) q[3];
sx q[3];
rz(-0.32005537) q[3];
sx q[3];
rz(-2.4879932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(-2.807907) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.82814) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(2.939558) q[0];
rz(-1.0728041) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(-1.7477759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48702792) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(1.4902671) q[0];
rz(-pi) q[1];
rz(-3.0300573) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(-0.44484777) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2653633) q[1];
sx q[1];
rz(-2.6637042) q[1];
sx q[1];
rz(-2.1761314) q[1];
rz(-pi) q[2];
rz(2.4059541) q[3];
sx q[3];
rz(-1.2218352) q[3];
sx q[3];
rz(2.5603106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(0.8423155) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(2.1202309) q[0];
rz(-2.0102603) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(0.21025118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410941) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(1.11087) q[0];
x q[1];
rz(-2.7619751) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10340313) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(1.3693387) q[1];
rz(2.5641714) q[3];
sx q[3];
rz(-0.93772674) q[3];
sx q[3];
rz(-3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8551089) q[0];
sx q[0];
rz(-0.47874641) q[0];
sx q[0];
rz(-0.32740645) q[0];
rz(-pi) q[1];
rz(2.0301684) q[2];
sx q[2];
rz(-0.18006912) q[2];
sx q[2];
rz(-1.1262696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8969438) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(-2.0379373) q[1];
rz(2.2817357) q[3];
sx q[3];
rz(-1.1338187) q[3];
sx q[3];
rz(-2.6881998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(0.55244279) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66822806) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(-0.92064944) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(2.2424973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9793648) q[0];
sx q[0];
rz(-0.75408903) q[0];
sx q[0];
rz(0.92569949) q[0];
rz(2.0752418) q[2];
sx q[2];
rz(-2.5214508) q[2];
sx q[2];
rz(-1.3401741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5324459) q[1];
sx q[1];
rz(-1.6858613) q[1];
sx q[1];
rz(-1.2389517) q[1];
rz(-0.00021342834) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(-0.40094024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(2.7446279) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3544281) q[0];
sx q[0];
rz(-1.0930601) q[0];
sx q[0];
rz(0.54425311) q[0];
rz(-pi) q[1];
rz(-0.88231477) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(0.099909401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18394477) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(1.8326879) q[1];
rz(-pi) q[2];
rz(0.90225737) q[3];
sx q[3];
rz(-0.60562953) q[3];
sx q[3];
rz(2.827364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94840702) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(0.096654264) q[2];
rz(1.7276673) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(-1.2697521) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(-2.4772714) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(3.0222223) q[3];
sx q[3];
rz(-1.8277677) q[3];
sx q[3];
rz(1.2895332) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];