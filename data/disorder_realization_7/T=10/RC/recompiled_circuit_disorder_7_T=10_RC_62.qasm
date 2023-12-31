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
rz(-0.62682682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7456822) q[0];
sx q[0];
rz(-2.9183309) q[0];
sx q[0];
rz(-2.1366871) q[0];
rz(-pi) q[1];
rz(1.3267514) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(0.92820456) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37571733) q[1];
sx q[1];
rz(-2.1581576) q[1];
sx q[1];
rz(-2.767763) q[1];
x q[2];
rz(-0.17625531) q[3];
sx q[3];
rz(-1.8415383) q[3];
sx q[3];
rz(-2.3034629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(-1.2935151) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(-0.75479341) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-2.0825443) q[0];
sx q[0];
rz(-2.1564116) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023022368) q[0];
sx q[0];
rz(-1.1623628) q[0];
sx q[0];
rz(-0.9704216) q[0];
rz(-pi) q[1];
rz(1.2843578) q[2];
sx q[2];
rz(-2.1396494) q[2];
sx q[2];
rz(1.4494277) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0774035) q[1];
sx q[1];
rz(-1.940155) q[1];
sx q[1];
rz(-2.2662152) q[1];
x q[2];
rz(-1.4025877) q[3];
sx q[3];
rz(-1.1871927) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5169516) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(0.66169935) q[2];
rz(-3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530281) q[0];
sx q[0];
rz(-1.0596501) q[0];
sx q[0];
rz(-1.9938064) q[0];
rz(-pi) q[1];
rz(-0.15180397) q[2];
sx q[2];
rz(-1.2200583) q[2];
sx q[2];
rz(2.0098067) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3861474) q[1];
sx q[1];
rz(-1.8090994) q[1];
sx q[1];
rz(2.9134681) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1654487) q[3];
sx q[3];
rz(-0.38446063) q[3];
sx q[3];
rz(2.773075) q[3];
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
rz(-1.6484377) q[3];
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
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4937113) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(0.9961876) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(0.48809537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195388) q[0];
sx q[0];
rz(-1.6484123) q[0];
sx q[0];
rz(1.104943) q[0];
x q[1];
rz(-2.5627665) q[2];
sx q[2];
rz(-1.0597214) q[2];
sx q[2];
rz(-1.5103112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51302233) q[1];
sx q[1];
rz(-1.7254618) q[1];
sx q[1];
rz(0.13048529) q[1];
x q[2];
rz(-2.475708) q[3];
sx q[3];
rz(-2.4443691) q[3];
sx q[3];
rz(-0.22660412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-2.6341237) q[2];
rz(1.1821702) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5421211) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(2.8778991) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(0.19651861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21393468) q[0];
sx q[0];
rz(-1.8727344) q[0];
sx q[0];
rz(1.2769804) q[0];
rz(-pi) q[1];
rz(3.1095805) q[2];
sx q[2];
rz(-2.0316761) q[2];
sx q[2];
rz(2.8061158) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49911753) q[1];
sx q[1];
rz(-0.31458464) q[1];
sx q[1];
rz(-2.5423074) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5444364) q[3];
sx q[3];
rz(-0.32005537) q[3];
sx q[3];
rz(2.4879932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-0.91709843) q[2];
rz(2.3069978) q[3];
sx q[3];
rz(-1.83225) q[3];
sx q[3];
rz(2.807907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(-2.939558) q[0];
rz(2.0687885) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(1.3938168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4112339) q[0];
sx q[0];
rz(-2.8011836) q[0];
sx q[0];
rz(-0.22986869) q[0];
rz(-1.9453085) q[2];
sx q[2];
rz(-1.4669344) q[2];
sx q[2];
rz(1.9749157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2653633) q[1];
sx q[1];
rz(-2.6637042) q[1];
sx q[1];
rz(2.1761314) q[1];
rz(-2.6447547) q[3];
sx q[3];
rz(-2.3416069) q[3];
sx q[3];
rz(0.62832384) q[3];
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
rz(2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52952805) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(1.0213617) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-2.1919057) q[1];
sx q[1];
rz(-2.9313415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049858) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(-1.11087) q[0];
rz(-pi) q[1];
rz(0.37961752) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(-1.1773674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5308295) q[1];
sx q[1];
rz(-1.3794583) q[1];
sx q[1];
rz(-0.32236871) q[1];
x q[2];
rz(2.2103569) q[3];
sx q[3];
rz(-0.82914549) q[3];
sx q[3];
rz(0.88013807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82688275) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(-1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(-2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49522266) q[0];
sx q[0];
rz(-0.78283993) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(2.6953221) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(0.97672021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2864838) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(-2.8141862) q[0];
rz(0.080539695) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(-1.5921519) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8969438) q[1];
sx q[1];
rz(-1.2243425) q[1];
sx q[1];
rz(2.0379373) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55240734) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.4668902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4010767) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(-2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(-2.2209432) q[0];
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
rz(-0.63910045) q[0];
sx q[0];
rz(-0.99196767) q[0];
sx q[0];
rz(-0.51410189) q[0];
rz(2.1295206) q[2];
sx q[2];
rz(-1.2860824) q[2];
sx q[2];
rz(2.4887091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2832665) q[1];
sx q[1];
rz(-0.35052931) q[1];
sx q[1];
rz(1.9117029) q[1];
rz(1.5713463) q[3];
sx q[3];
rz(-2.7714202) q[3];
sx q[3];
rz(-2.7412424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(0.88275638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.33912441) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(-1.1178281) q[0];
rz(-2.4291908) q[1];
sx q[1];
rz(-2.5103705) q[1];
sx q[1];
rz(2.7446279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3544281) q[0];
sx q[0];
rz(-2.0485326) q[0];
sx q[0];
rz(-2.5973395) q[0];
x q[1];
rz(2.8357387) q[2];
sx q[2];
rz(-1.2199739) q[2];
sx q[2];
rz(-0.84301126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18394477) q[1];
sx q[1];
rz(-2.7937104) q[1];
sx q[1];
rz(-1.3089048) q[1];
rz(-pi) q[2];
rz(-2.0685365) q[3];
sx q[3];
rz(-1.2101678) q[3];
sx q[3];
rz(-1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(0.096654264) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.5554572) q[0];
sx q[0];
rz(-1.1145034) q[0];
sx q[0];
rz(0.282423) q[0];
rz(1.2697521) q[1];
sx q[1];
rz(-1.7813663) q[1];
sx q[1];
rz(-2.355994) q[1];
rz(-1.6554228) q[2];
sx q[2];
rz(-1.463269) q[2];
sx q[2];
rz(-0.29381897) q[2];
rz(-1.1453015) q[3];
sx q[3];
rz(-2.8588061) q[3];
sx q[3];
rz(-2.2929946) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
