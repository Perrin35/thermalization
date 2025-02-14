OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(-1.1046326) q[0];
rz(-1.8072577) q[1];
sx q[1];
rz(5.5601064) q[1];
sx q[1];
rz(13.830166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1570779) q[0];
sx q[0];
rz(-2.4369168) q[0];
sx q[0];
rz(-1.9542171) q[0];
rz(-0.96226389) q[2];
sx q[2];
rz(-0.66676408) q[2];
sx q[2];
rz(-0.61362574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0824737) q[1];
sx q[1];
rz(-2.8063941) q[1];
sx q[1];
rz(-1.3223865) q[1];
rz(-2.8975211) q[3];
sx q[3];
rz(-1.1435978) q[3];
sx q[3];
rz(-1.106316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15452142) q[2];
sx q[2];
rz(-2.0646586) q[2];
sx q[2];
rz(1.7729574) q[2];
rz(-0.23855071) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(-0.11428782) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7411165) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(1.9912632) q[0];
rz(0.087609619) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(1.5709343) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8920736) q[0];
sx q[0];
rz(-2.6457835) q[0];
sx q[0];
rz(0.22479381) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17878647) q[2];
sx q[2];
rz(-1.6331722) q[2];
sx q[2];
rz(1.5966001) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9894692) q[1];
sx q[1];
rz(-1.996576) q[1];
sx q[1];
rz(-2.2683737) q[1];
rz(-pi) q[2];
rz(-0.47111311) q[3];
sx q[3];
rz(-0.84976746) q[3];
sx q[3];
rz(-1.6278933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7972083) q[2];
sx q[2];
rz(-2.4281561) q[2];
sx q[2];
rz(2.9551282) q[2];
rz(2.3421085) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(0.98108393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28476533) q[0];
sx q[0];
rz(-1.7600049) q[0];
sx q[0];
rz(0.10936603) q[0];
rz(-2.5378387) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(-2.107479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4277735) q[0];
sx q[0];
rz(-1.7805011) q[0];
sx q[0];
rz(-3.0056277) q[0];
rz(-pi) q[1];
rz(-0.66884235) q[2];
sx q[2];
rz(-0.1136264) q[2];
sx q[2];
rz(-2.4955028) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2455403) q[1];
sx q[1];
rz(-0.29948343) q[1];
sx q[1];
rz(-2.7761032) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3699371) q[3];
sx q[3];
rz(-1.6724186) q[3];
sx q[3];
rz(-1.9664498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6094531) q[2];
sx q[2];
rz(-2.1707462) q[2];
sx q[2];
rz(2.5965221) q[2];
rz(-0.80101454) q[3];
sx q[3];
rz(-2.2478734) q[3];
sx q[3];
rz(-3.0494704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545749) q[0];
sx q[0];
rz(-1.6681404) q[0];
sx q[0];
rz(-0.45904485) q[0];
rz(1.1162988) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(-2.3150516) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61454489) q[0];
sx q[0];
rz(-1.9226719) q[0];
sx q[0];
rz(1.1220758) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0439348) q[2];
sx q[2];
rz(-1.6689166) q[2];
sx q[2];
rz(-2.1922534) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0361358) q[1];
sx q[1];
rz(-1.030827) q[1];
sx q[1];
rz(-0.2293052) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6891805) q[3];
sx q[3];
rz(-1.2384129) q[3];
sx q[3];
rz(-0.82086241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35215846) q[2];
sx q[2];
rz(-0.23098478) q[2];
sx q[2];
rz(2.3918772) q[2];
rz(2.0189144) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(2.1813755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66626755) q[0];
sx q[0];
rz(-2.1551977) q[0];
sx q[0];
rz(0.112003) q[0];
rz(0.22459596) q[1];
sx q[1];
rz(-1.9738395) q[1];
sx q[1];
rz(-2.3602643) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8664426) q[0];
sx q[0];
rz(-2.7366182) q[0];
sx q[0];
rz(1.6385796) q[0];
x q[1];
rz(-1.6649328) q[2];
sx q[2];
rz(-1.899666) q[2];
sx q[2];
rz(0.48966792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0722149) q[1];
sx q[1];
rz(-2.2615914) q[1];
sx q[1];
rz(2.1363972) q[1];
rz(-pi) q[2];
rz(-2.1634401) q[3];
sx q[3];
rz(-1.4362122) q[3];
sx q[3];
rz(-2.5434567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3096699) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(0.93264467) q[2];
rz(-2.0617088) q[3];
sx q[3];
rz(-0.44411689) q[3];
sx q[3];
rz(0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90750736) q[0];
sx q[0];
rz(-2.0103173) q[0];
sx q[0];
rz(-1.750741) q[0];
rz(-2.8098409) q[1];
sx q[1];
rz(-1.876588) q[1];
sx q[1];
rz(1.5914241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0858147) q[0];
sx q[0];
rz(-1.6059438) q[0];
sx q[0];
rz(2.2136392) q[0];
rz(-2.6820002) q[2];
sx q[2];
rz(-1.6038648) q[2];
sx q[2];
rz(3.1350373) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6350123) q[1];
sx q[1];
rz(-1.6939347) q[1];
sx q[1];
rz(1.7282196) q[1];
x q[2];
rz(0.8121536) q[3];
sx q[3];
rz(-2.6366173) q[3];
sx q[3];
rz(-2.2157089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3127689) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(-1.5213607) q[2];
rz(0.96794266) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(-0.91606417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5172326) q[0];
sx q[0];
rz(-0.42651287) q[0];
sx q[0];
rz(-2.970001) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.2587222) q[1];
sx q[1];
rz(1.4412057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2644669) q[0];
sx q[0];
rz(-1.4466337) q[0];
sx q[0];
rz(-0.32372253) q[0];
rz(-pi) q[1];
rz(1.4588548) q[2];
sx q[2];
rz(-2.0840692) q[2];
sx q[2];
rz(1.3789195) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4644147) q[1];
sx q[1];
rz(-0.59301335) q[1];
sx q[1];
rz(-0.89927425) q[1];
x q[2];
rz(3.0002241) q[3];
sx q[3];
rz(-2.1297567) q[3];
sx q[3];
rz(-1.5822922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8695996) q[2];
sx q[2];
rz(-1.4078434) q[2];
sx q[2];
rz(-0.54455152) q[2];
rz(0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(-0.47206363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8498103) q[0];
sx q[0];
rz(-2.6047459) q[0];
sx q[0];
rz(0.52919069) q[0];
rz(-0.57580194) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(1.1189438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42903439) q[0];
sx q[0];
rz(-0.055744113) q[0];
sx q[0];
rz(-0.38614614) q[0];
rz(-pi) q[1];
x q[1];
rz(2.697102) q[2];
sx q[2];
rz(-2.7348282) q[2];
sx q[2];
rz(0.14061804) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8994979) q[1];
sx q[1];
rz(-0.78259727) q[1];
sx q[1];
rz(2.461754) q[1];
rz(-pi) q[2];
rz(-0.20967926) q[3];
sx q[3];
rz(-2.0884313) q[3];
sx q[3];
rz(-2.8594494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90073663) q[2];
sx q[2];
rz(-0.46515981) q[2];
sx q[2];
rz(-2.4294803) q[2];
rz(-1.5322878) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(1.3473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3466472) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(1.0711063) q[0];
rz(-1.935293) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(-1.4869022) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95222118) q[0];
sx q[0];
rz(-0.69714386) q[0];
sx q[0];
rz(2.8482262) q[0];
x q[1];
rz(2.1558482) q[2];
sx q[2];
rz(-2.2170456) q[2];
sx q[2];
rz(1.0111292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35569977) q[1];
sx q[1];
rz(-1.744387) q[1];
sx q[1];
rz(-2.8462571) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28122357) q[3];
sx q[3];
rz(-1.1373925) q[3];
sx q[3];
rz(3.0126257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.508076) q[2];
sx q[2];
rz(-2.241892) q[2];
sx q[2];
rz(3.0637975) q[2];
rz(-0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(2.0859065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28854293) q[0];
sx q[0];
rz(-2.396614) q[0];
sx q[0];
rz(-2.7689834) q[0];
rz(2.695072) q[1];
sx q[1];
rz(-1.6852854) q[1];
sx q[1];
rz(2.2850697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3055182) q[0];
sx q[0];
rz(-1.5139607) q[0];
sx q[0];
rz(-1.5268121) q[0];
rz(-0.90214731) q[2];
sx q[2];
rz(-2.9716316) q[2];
sx q[2];
rz(-0.68357498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0622934) q[1];
sx q[1];
rz(-0.41626272) q[1];
sx q[1];
rz(1.8658616) q[1];
x q[2];
rz(2.4803745) q[3];
sx q[3];
rz(-2.2228129) q[3];
sx q[3];
rz(-2.2673502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38241688) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(2.5860533) q[2];
rz(0.38604745) q[3];
sx q[3];
rz(-1.4945533) q[3];
sx q[3];
rz(2.4773795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1356708) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(0.80264965) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(-2.7886709) q[2];
sx q[2];
rz(-1.4581994) q[2];
sx q[2];
rz(-1.1117473) q[2];
rz(-1.4925719) q[3];
sx q[3];
rz(-0.73220069) q[3];
sx q[3];
rz(-2.2673741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
