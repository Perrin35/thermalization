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
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(2.5147658) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39591046) q[0];
sx q[0];
rz(-0.22326176) q[0];
sx q[0];
rz(2.1366871) q[0];
x q[1];
rz(2.7896499) q[2];
sx q[2];
rz(-1.3411739) q[2];
sx q[2];
rz(-0.72598347) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1498897) q[1];
sx q[1];
rz(-0.68421364) q[1];
sx q[1];
rz(-2.0725155) q[1];
rz(-2.1342282) q[3];
sx q[3];
rz(-0.3218739) q[3];
sx q[3];
rz(0.25062996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(-1.8480776) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-0.75479341) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(3.0342297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023022368) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(-0.9704216) q[0];
rz(-0.58789247) q[2];
sx q[2];
rz(-1.8111472) q[2];
sx q[2];
rz(-2.8628778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0641891) q[1];
sx q[1];
rz(-1.940155) q[1];
sx q[1];
rz(0.87537745) q[1];
x q[2];
rz(1.739005) q[3];
sx q[3];
rz(-1.9544) q[3];
sx q[3];
rz(0.83354359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5169516) q[2];
sx q[2];
rz(-1.3890283) q[2];
sx q[2];
rz(-0.66169935) q[2];
rz(3.0858357) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(2.8957446) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-0.89825392) q[1];
sx q[1];
rz(-1.0823762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6862415) q[0];
sx q[0];
rz(-0.65128122) q[0];
sx q[0];
rz(0.63182802) q[0];
rz(1.2163095) q[2];
sx q[2];
rz(-1.4282994) q[2];
sx q[2];
rz(2.6500677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.97835097) q[1];
sx q[1];
rz(-0.32838531) q[1];
sx q[1];
rz(-2.3204625) q[1];
rz(-2.1654487) q[3];
sx q[3];
rz(-0.38446063) q[3];
sx q[3];
rz(0.36851766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(-1.4931549) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-3.1062104) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(2.6534973) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513718) q[0];
sx q[0];
rz(-1.106456) q[0];
sx q[0];
rz(-0.086829348) q[0];
rz(-0.79779063) q[2];
sx q[2];
rz(-2.3893223) q[2];
sx q[2];
rz(-2.4385902) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51302233) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(-3.0111074) q[1];
rz(-pi) q[2];
rz(0.66588464) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(2.9149885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-0.50746894) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(2.7713293) q[0];
rz(-2.8778991) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(2.945074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58013433) q[0];
sx q[0];
rz(-0.41813865) q[0];
sx q[0];
rz(-2.3925376) q[0];
rz(1.5064266) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49911753) q[1];
sx q[1];
rz(-2.827008) q[1];
sx q[1];
rz(0.59928526) q[1];
rz(2.5971562) q[3];
sx q[3];
rz(-2.8215373) q[3];
sx q[3];
rz(2.4879932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(-0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.82814) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(-1.0728041) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(1.3938168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545647) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(-1.4902671) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11153531) q[2];
sx q[2];
rz(-1.9431912) q[2];
sx q[2];
rz(-2.6967449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.285167) q[1];
sx q[1];
rz(-1.835583) q[1];
sx q[1];
rz(1.1681639) q[1];
x q[2];
rz(-0.49683797) q[3];
sx q[3];
rz(-2.3416069) q[3];
sx q[3];
rz(-0.62832384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(-1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(-2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-0.21025118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1649095) q[0];
sx q[0];
rz(-1.7047791) q[0];
sx q[0];
rz(3.0747736) q[0];
x q[1];
rz(-2.7619751) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(1.9642252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0381895) q[1];
sx q[1];
rz(-1.8870755) q[1];
sx q[1];
rz(-1.7722539) q[1];
rz(-pi) q[2];
rz(2.290091) q[3];
sx q[3];
rz(-1.1151033) q[3];
sx q[3];
rz(-1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.3898536) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.64637) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(-2.5543509) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(0.97672021) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644972) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(0.45678267) q[0];
rz(-1.1114242) q[2];
sx q[2];
rz(-2.9615235) q[2];
sx q[2];
rz(1.1262696) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2446489) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(1.1036554) q[1];
rz(2.5891853) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.4668902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-1.027511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(-2.2209432) q[0];
rz(3.0905511) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(-0.89909536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16222787) q[0];
sx q[0];
rz(-2.3875036) q[0];
sx q[0];
rz(0.92569949) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8092438) q[2];
sx q[2];
rz(-2.1045448) q[2];
sx q[2];
rz(-0.74408434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60914674) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(1.902641) q[1];
rz(-1.2006239) q[3];
sx q[3];
rz(-1.5709953) q[3];
sx q[3];
rz(1.9716594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(-0.71240187) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(-2.7446279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3544281) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(0.54425311) q[0];
rz(-pi) q[1];
rz(-2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(0.099909401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0015928) q[1];
sx q[1];
rz(-1.4824176) q[1];
sx q[1];
rz(1.9077076) q[1];
rz(-pi) q[2];
rz(2.2393353) q[3];
sx q[3];
rz(-0.60562953) q[3];
sx q[3];
rz(-2.827364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-3.0449384) q[2];
rz(-1.7276673) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(1.2697521) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-2.4772714) q[2];
sx q[2];
rz(-0.13673377) q[2];
sx q[2];
rz(-2.7665334) q[2];
rz(-1.9962911) q[3];
sx q[3];
rz(-0.28278657) q[3];
sx q[3];
rz(0.84859802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
