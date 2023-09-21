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
rz(-0.20110826) q[0];
rz(1.7445298) q[1];
sx q[1];
rz(-1.3367329) q[1];
sx q[1];
rz(0.62682682) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(1.0049055) q[0];
rz(-pi) q[1];
rz(-2.5457382) q[2];
sx q[2];
rz(-0.41759727) q[2];
sx q[2];
rz(2.8516304) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1605692) q[1];
sx q[1];
rz(-1.8796646) q[1];
sx q[1];
rz(0.94998756) q[1];
rz(1.8455891) q[3];
sx q[3];
rz(-1.7405675) q[3];
sx q[3];
rz(2.3613289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(-1.2935151) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24793808) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(2.1564116) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-3.0342297) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1194612) q[0];
sx q[0];
rz(-2.4298926) q[0];
sx q[0];
rz(2.2244781) q[0];
rz(-pi) q[1];
rz(-2.7254843) q[2];
sx q[2];
rz(-0.62972087) q[2];
sx q[2];
rz(2.1925418) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.7992236) q[1];
sx q[1];
rz(-0.93042004) q[1];
sx q[1];
rz(-2.6745822) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38856216) q[3];
sx q[3];
rz(-1.4149168) q[3];
sx q[3];
rz(0.673783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(2.4798933) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67702883) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(-3.1012428) q[0];
rz(0.57998747) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-2.0592164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28856453) q[0];
sx q[0];
rz(-2.0819426) q[0];
sx q[0];
rz(1.9938064) q[0];
x q[1];
rz(-0.15180397) q[2];
sx q[2];
rz(-1.9215343) q[2];
sx q[2];
rz(1.1317859) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.3324932) q[1];
sx q[1];
rz(-2.9134681) q[1];
rz(-pi) q[2];
rz(2.9186967) q[3];
sx q[3];
rz(-1.2548903) q[3];
sx q[3];
rz(2.1427758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(-0.55389261) q[2];
rz(1.4931549) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-2.557502) q[0];
sx q[0];
rz(-3.1062104) q[0];
rz(0.9961876) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-0.48809537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0220538) q[0];
sx q[0];
rz(-1.6484123) q[0];
sx q[0];
rz(1.104943) q[0];
x q[1];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0779873) q[1];
sx q[1];
rz(-1.6997153) q[1];
sx q[1];
rz(-1.414826) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.048269) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(-2.570591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35486832) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(1.2549843) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-0.71389714) q[0];
sx q[0];
rz(0.3702634) q[0];
rz(2.8778991) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-2.945074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21393468) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(1.2769804) q[0];
rz(-2.0318803) q[2];
sx q[2];
rz(-1.5994674) q[2];
sx q[2];
rz(1.920514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6477485) q[1];
sx q[1];
rz(-1.7462247) q[1];
sx q[1];
rz(0.26248787) q[1];
rz(1.4007832) q[3];
sx q[3];
rz(-1.8432872) q[3];
sx q[3];
rz(1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(0.91709843) q[2];
rz(0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-0.33368567) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31345263) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(-2.939558) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(1.7477759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73035875) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(2.911724) q[0];
x q[1];
rz(1.1962842) q[2];
sx q[2];
rz(-1.4669344) q[2];
sx q[2];
rz(-1.166677) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2653633) q[1];
sx q[1];
rz(-0.47788844) q[1];
sx q[1];
rz(-2.1761314) q[1];
x q[2];
rz(-1.1145669) q[3];
sx q[3];
rz(-2.2531415) q[3];
sx q[3];
rz(1.8519459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(-2.2992772) q[2];
rz(-2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.9313415) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410941) q[0];
sx q[0];
rz(-2.9919618) q[0];
sx q[0];
rz(2.0307226) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5049997) q[2];
sx q[2];
rz(-1.4074667) q[2];
sx q[2];
rz(2.3492299) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5308295) q[1];
sx q[1];
rz(-1.3794583) q[1];
sx q[1];
rz(-2.8192239) q[1];
rz(-2.5641714) q[3];
sx q[3];
rz(-0.93772674) q[3];
sx q[3];
rz(3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-0.41346082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64637) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(2.5543509) q[0];
rz(-2.6953221) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(-0.97672021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57709549) q[0];
sx q[0];
rz(-1.7194887) q[0];
sx q[0];
rz(0.45678267) q[0];
rz(1.1114242) q[2];
sx q[2];
rz(-0.18006912) q[2];
sx q[2];
rz(-2.0153231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2446489) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(2.0379373) q[1];
rz(-pi) q[2];
rz(-0.94954357) q[3];
sx q[3];
rz(-2.3275259) q[3];
sx q[3];
rz(-2.4809588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(-0.55244279) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(-0.89909536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102525) q[0];
sx q[0];
rz(-1.9950584) q[0];
sx q[0];
rz(2.2146241) q[0];
rz(1.012072) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(0.65288359) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8583262) q[1];
sx q[1];
rz(-2.7910633) q[1];
sx q[1];
rz(1.9117029) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1413792) q[3];
sx q[3];
rz(-1.2006239) q[3];
sx q[3];
rz(-2.7406524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(-2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(2.2588363) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1966232) q[0];
sx q[0];
rz(-1.0930177) q[0];
sx q[0];
rz(-2.1150132) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2043578) q[2];
sx q[2];
rz(-1.2841184) q[2];
sx q[2];
rz(2.3057126) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18394477) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(-1.8326879) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7361717) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(3.0449384) q[2];
rz(-1.7276673) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(1.1269425) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5861355) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(-1.8718406) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-0.66432129) q[2];
sx q[2];
rz(-3.0048589) q[2];
sx q[2];
rz(0.37505925) q[2];
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
