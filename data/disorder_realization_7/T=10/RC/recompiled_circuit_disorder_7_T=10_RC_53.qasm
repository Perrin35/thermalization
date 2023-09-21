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
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(-0.62682682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97323595) q[0];
sx q[0];
rz(-1.7587979) q[0];
sx q[0];
rz(3.0204535) q[0];
x q[1];
rz(-1.3267514) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(-0.92820456) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1498897) q[1];
sx q[1];
rz(-2.457379) q[1];
sx q[1];
rz(-2.0725155) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2960035) q[3];
sx q[3];
rz(-1.7405675) q[3];
sx q[3];
rz(-0.78026375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24793808) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(-0.10736297) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1185703) q[0];
sx q[0];
rz(-1.9792299) q[0];
sx q[0];
rz(-0.9704216) q[0];
rz(-0.41610833) q[2];
sx q[2];
rz(-0.62972087) q[2];
sx q[2];
rz(-2.1925418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0641891) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(0.87537745) q[1];
rz(-pi) q[2];
rz(-1.739005) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(0.83354359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-0.055756904) q[3];
sx q[3];
rz(-2.7825833) q[3];
sx q[3];
rz(-0.24584809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702883) q[0];
sx q[0];
rz(-0.55389535) q[0];
sx q[0];
rz(-0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(1.0823762) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530281) q[0];
sx q[0];
rz(-2.0819426) q[0];
sx q[0];
rz(1.1477863) q[0];
rz(1.9627377) q[2];
sx q[2];
rz(-2.760663) q[2];
sx q[2];
rz(-2.4286963) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3861474) q[1];
sx q[1];
rz(-1.3324932) q[1];
sx q[1];
rz(-0.2281245) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2228959) q[3];
sx q[3];
rz(-1.8867023) q[3];
sx q[3];
rz(-0.99881682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13517705) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(-1.6484377) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(-3.1062104) q[0];
rz(-0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(-0.48809537) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6513718) q[0];
sx q[0];
rz(-1.106456) q[0];
sx q[0];
rz(3.0547633) q[0];
rz(0.79779063) q[2];
sx q[2];
rz(-0.75227037) q[2];
sx q[2];
rz(-2.4385902) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51302233) q[1];
sx q[1];
rz(-1.4161308) q[1];
sx q[1];
rz(-0.13048529) q[1];
x q[2];
rz(2.048269) q[3];
sx q[3];
rz(-1.0415047) q[3];
sx q[3];
rz(0.57100163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35486832) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(2.6341237) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.8488665) q[3];
sx q[3];
rz(-1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(0.19651861) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21393468) q[0];
sx q[0];
rz(-1.2688583) q[0];
sx q[0];
rz(1.8646122) q[0];
x q[1];
rz(-3.1095805) q[2];
sx q[2];
rz(-1.1099166) q[2];
sx q[2];
rz(-0.33547685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49911753) q[1];
sx q[1];
rz(-0.31458464) q[1];
sx q[1];
rz(2.5423074) q[1];
x q[2];
rz(-1.7408095) q[3];
sx q[3];
rz(-1.2983054) q[3];
sx q[3];
rz(-1.221399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(2.2244942) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(-2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.82814) q[0];
sx q[0];
rz(-1.7001292) q[0];
sx q[0];
rz(2.939558) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(-1.3938168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4112339) q[0];
sx q[0];
rz(-2.8011836) q[0];
sx q[0];
rz(2.911724) q[0];
rz(-pi) q[1];
rz(1.1962842) q[2];
sx q[2];
rz(-1.6746582) q[2];
sx q[2];
rz(1.166677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.285167) q[1];
sx q[1];
rz(-1.835583) q[1];
sx q[1];
rz(-1.1681639) q[1];
rz(-1.1145669) q[3];
sx q[3];
rz(-0.8884512) q[3];
sx q[3];
rz(1.2896468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.6615901) q[2];
sx q[2];
rz(2.2992772) q[2];
rz(-1.0874282) q[3];
sx q[3];
rz(-0.73173404) q[3];
sx q[3];
rz(1.4830164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-2.1202309) q[0];
rz(-2.0102603) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(2.9313415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41482571) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(-1.7050752) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9779151) q[2];
sx q[2];
rz(-1.6357161) q[2];
sx q[2];
rz(0.76771969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10340313) q[1];
sx q[1];
rz(-1.8870755) q[1];
sx q[1];
rz(1.3693387) q[1];
rz(-2.5641714) q[3];
sx q[3];
rz(-2.2038659) q[3];
sx q[3];
rz(-3.0949743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(1.7604527) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-2.7281318) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(0.97672021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2205692) q[0];
sx q[0];
rz(-2.0221634) q[0];
sx q[0];
rz(1.736182) q[0];
x q[1];
rz(1.1114242) q[2];
sx q[2];
rz(-0.18006912) q[2];
sx q[2];
rz(1.1262696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8969438) q[1];
sx q[1];
rz(-1.9172501) q[1];
sx q[1];
rz(-2.0379373) q[1];
x q[2];
rz(2.5891853) q[3];
sx q[3];
rz(-0.93821412) q[3];
sx q[3];
rz(1.6747024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-0.70827168) q[2];
sx q[2];
rz(2.5891499) q[2];
rz(0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(-2.1140816) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4733646) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(-2.2209432) q[0];
rz(0.051041516) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(0.89909536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5024922) q[0];
sx q[0];
rz(-2.149625) q[0];
sx q[0];
rz(-0.51410189) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1295206) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(-0.65288359) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60914674) q[1];
sx q[1];
rz(-1.4557314) q[1];
sx q[1];
rz(-1.2389517) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2006239) q[3];
sx q[3];
rz(-1.5705974) q[3];
sx q[3];
rz(-1.1699333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-0.070572704) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(-2.3032522) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-2.7446279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94496942) q[0];
sx q[0];
rz(-2.0485749) q[0];
sx q[0];
rz(-1.0265795) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2592779) q[2];
sx q[2];
rz(-2.6803662) q[2];
sx q[2];
rz(-3.0416833) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.18394477) q[1];
sx q[1];
rz(-0.34788222) q[1];
sx q[1];
rz(-1.8326879) q[1];
x q[2];
rz(-2.2393353) q[3];
sx q[3];
rz(-0.60562953) q[3];
sx q[3];
rz(-0.31422869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-1.1380514) q[3];
sx q[3];
rz(-1.1269425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(0.10791049) q[2];
sx q[2];
rz(-1.4866598) q[2];
sx q[2];
rz(1.286081) q[2];
rz(0.11937033) q[3];
sx q[3];
rz(-1.3138249) q[3];
sx q[3];
rz(-1.8520595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
