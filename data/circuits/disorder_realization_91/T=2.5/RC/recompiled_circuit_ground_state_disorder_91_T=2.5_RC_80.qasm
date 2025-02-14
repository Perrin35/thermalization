OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.33136) q[0];
sx q[0];
rz(-1.2547837) q[0];
sx q[0];
rz(-2.4956508) q[0];
rz(1.6826001) q[1];
sx q[1];
rz(-3.0134633) q[1];
sx q[1];
rz(0.66817862) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2006045) q[0];
sx q[0];
rz(-1.4580237) q[0];
sx q[0];
rz(-0.24399816) q[0];
rz(-pi) q[1];
rz(1.1485841) q[2];
sx q[2];
rz(-2.4496884) q[2];
sx q[2];
rz(-2.1466604) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9757802) q[1];
sx q[1];
rz(-2.4620612) q[1];
sx q[1];
rz(0.91109101) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2916037) q[3];
sx q[3];
rz(-1.9487983) q[3];
sx q[3];
rz(-1.8545618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51840034) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(1.0830967) q[2];
rz(2.8862503) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1622247) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-3.1080143) q[0];
rz(-2.0032739) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(-0.23695645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36719272) q[0];
sx q[0];
rz(-2.8731308) q[0];
sx q[0];
rz(2.5949097) q[0];
rz(2.9541624) q[2];
sx q[2];
rz(-1.1016365) q[2];
sx q[2];
rz(-2.3694866) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2365429) q[1];
sx q[1];
rz(-1.7631751) q[1];
sx q[1];
rz(-3.1156999) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99144793) q[3];
sx q[3];
rz(-1.1889403) q[3];
sx q[3];
rz(-0.53296158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61243764) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(0.1203514) q[2];
rz(-2.555661) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(1.2683292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878207) q[0];
sx q[0];
rz(-2.1757941) q[0];
sx q[0];
rz(-0.98814386) q[0];
rz(-1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(1.5536701) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33086209) q[0];
sx q[0];
rz(-1.5867044) q[0];
sx q[0];
rz(1.5718979) q[0];
rz(-pi) q[1];
rz(-2.6695197) q[2];
sx q[2];
rz(-0.1506714) q[2];
sx q[2];
rz(0.27418384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7657679) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(3.1118666) q[1];
rz(-1.1497028) q[3];
sx q[3];
rz(-1.0437056) q[3];
sx q[3];
rz(-2.2452109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3042018) q[2];
sx q[2];
rz(-1.7981671) q[2];
sx q[2];
rz(-0.070579441) q[2];
rz(-2.6523759) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(0.14744559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1864784) q[0];
sx q[0];
rz(-0.69545737) q[0];
sx q[0];
rz(1.4008993) q[0];
rz(2.9715624) q[1];
sx q[1];
rz(-1.731512) q[1];
sx q[1];
rz(-1.9084575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2182541) q[0];
sx q[0];
rz(-1.6955351) q[0];
sx q[0];
rz(-1.9943004) q[0];
rz(-0.71738404) q[2];
sx q[2];
rz(-2.0307433) q[2];
sx q[2];
rz(1.9460974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68076949) q[1];
sx q[1];
rz(-1.5857547) q[1];
sx q[1];
rz(1.8504924) q[1];
rz(-pi) q[2];
rz(-1.5587658) q[3];
sx q[3];
rz(-2.0352484) q[3];
sx q[3];
rz(2.8248276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79973811) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(-0.62961659) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(1.4153882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.8585696) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(-0.10230219) q[0];
rz(1.6434068) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(2.7844875) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5095561) q[0];
sx q[0];
rz(-2.8297544) q[0];
sx q[0];
rz(0.76550342) q[0];
rz(-pi) q[1];
rz(-0.12320766) q[2];
sx q[2];
rz(-1.5977309) q[2];
sx q[2];
rz(-1.7719442) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2374463) q[1];
sx q[1];
rz(-2.0761298) q[1];
sx q[1];
rz(-1.2282441) q[1];
rz(-pi) q[2];
rz(-0.83729845) q[3];
sx q[3];
rz(-2.5243763) q[3];
sx q[3];
rz(-1.069998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6685278) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(1.3225887) q[2];
rz(1.708301) q[3];
sx q[3];
rz(-1.3168443) q[3];
sx q[3];
rz(-0.1639666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.598031) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(1.7033956) q[0];
rz(-1.2403129) q[1];
sx q[1];
rz(-2.4867609) q[1];
sx q[1];
rz(-1.3528489) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.529218) q[0];
sx q[0];
rz(-1.1657526) q[0];
sx q[0];
rz(3.0390826) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4022302) q[2];
sx q[2];
rz(-1.050808) q[2];
sx q[2];
rz(1.2181768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53735087) q[1];
sx q[1];
rz(-0.64879319) q[1];
sx q[1];
rz(0.93334812) q[1];
x q[2];
rz(0.26627906) q[3];
sx q[3];
rz(-2.7334573) q[3];
sx q[3];
rz(0.88312393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9271586) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(0.39043179) q[2];
rz(-1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(-0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.696233) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(-1.4129289) q[0];
rz(1.4631118) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(-0.63953343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4832925) q[0];
sx q[0];
rz(-0.72972238) q[0];
sx q[0];
rz(0.99131063) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15070559) q[2];
sx q[2];
rz(-1.6455212) q[2];
sx q[2];
rz(-0.13241235) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29296103) q[1];
sx q[1];
rz(-2.1994414) q[1];
sx q[1];
rz(-1.2413003) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1763379) q[3];
sx q[3];
rz(-1.6812075) q[3];
sx q[3];
rz(2.4019456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9245236) q[2];
sx q[2];
rz(-1.3288493) q[2];
sx q[2];
rz(-2.6893943) q[2];
rz(-2.9186115) q[3];
sx q[3];
rz(-2.860234) q[3];
sx q[3];
rz(-2.9862459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13101354) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(2.9796694) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(0.23439342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18412515) q[0];
sx q[0];
rz(-2.1229366) q[0];
sx q[0];
rz(-0.87463057) q[0];
x q[1];
rz(-0.35104819) q[2];
sx q[2];
rz(-0.75134885) q[2];
sx q[2];
rz(-2.4780432) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3752611) q[1];
sx q[1];
rz(-2.2041781) q[1];
sx q[1];
rz(-2.9087025) q[1];
rz(-pi) q[2];
rz(-2.864847) q[3];
sx q[3];
rz(-1.6349155) q[3];
sx q[3];
rz(-1.8217877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8552385) q[2];
sx q[2];
rz(-1.0711292) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(-0.74328077) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(-2.4519517) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8648935) q[0];
sx q[0];
rz(-1.6962637) q[0];
sx q[0];
rz(-2.0546761) q[0];
rz(-1.9089606) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(-1.8392275) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24827458) q[0];
sx q[0];
rz(-2.9030307) q[0];
sx q[0];
rz(2.6435411) q[0];
x q[1];
rz(-2.7872269) q[2];
sx q[2];
rz(-1.003689) q[2];
sx q[2];
rz(1.815955) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2375396) q[1];
sx q[1];
rz(-1.9642648) q[1];
sx q[1];
rz(-0.57668011) q[1];
rz(-pi) q[2];
rz(-1.5518673) q[3];
sx q[3];
rz(-0.94784465) q[3];
sx q[3];
rz(1.2894693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46633259) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(2.7272398) q[2];
rz(0.77053344) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(2.2079302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122124) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(-0.46863753) q[0];
rz(-1.8679484) q[1];
sx q[1];
rz(-1.8654774) q[1];
sx q[1];
rz(2.830107) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5632919) q[0];
sx q[0];
rz(-1.5130763) q[0];
sx q[0];
rz(-2.2895534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6577254) q[2];
sx q[2];
rz(-1.4540028) q[2];
sx q[2];
rz(0.96125666) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7817318) q[1];
sx q[1];
rz(-1.5241429) q[1];
sx q[1];
rz(2.0585039) q[1];
rz(2.0592779) q[3];
sx q[3];
rz(-1.0955819) q[3];
sx q[3];
rz(2.7431874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0353388) q[2];
sx q[2];
rz(-1.9474578) q[2];
sx q[2];
rz(1.3930456) q[2];
rz(2.2140908) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(2.2813796) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504234) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(2.9604079) q[1];
sx q[1];
rz(-2.1962427) q[1];
sx q[1];
rz(-0.93999351) q[1];
rz(-2.9857582) q[2];
sx q[2];
rz(-1.0424141) q[2];
sx q[2];
rz(-0.14137683) q[2];
rz(-0.60605449) q[3];
sx q[3];
rz(-0.99679508) q[3];
sx q[3];
rz(1.5011277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
