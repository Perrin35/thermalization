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
rz(4.9464524) q[1];
sx q[1];
rz(10.051605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1683567) q[0];
sx q[0];
rz(-1.7587979) q[0];
sx q[0];
rz(0.12113916) q[0];
rz(-1.8148412) q[2];
sx q[2];
rz(-1.913117) q[2];
sx q[2];
rz(0.92820456) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7658753) q[1];
sx q[1];
rz(-0.98343508) q[1];
sx q[1];
rz(0.37382965) q[1];
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
rz(-2.8864002) q[2];
sx q[2];
rz(-0.93078405) q[2];
sx q[2];
rz(1.2935151) q[2];
rz(-1.1039929) q[3];
sx q[3];
rz(-2.2938426) q[3];
sx q[3];
rz(-2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(-0.39918104) q[1];
sx q[1];
rz(-1.2626516) q[1];
sx q[1];
rz(0.10736297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3283008) q[0];
sx q[0];
rz(-2.11587) q[0];
sx q[0];
rz(-0.48304793) q[0];
x q[1];
rz(2.7254843) q[2];
sx q[2];
rz(-2.5118718) q[2];
sx q[2];
rz(2.1925418) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3423691) q[1];
sx q[1];
rz(-2.2111726) q[1];
sx q[1];
rz(0.46701042) q[1];
rz(-pi) q[2];
rz(1.4025877) q[3];
sx q[3];
rz(-1.9544) q[3];
sx q[3];
rz(2.3080491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5169516) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645638) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(3.1012428) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-2.0592164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6862415) q[0];
sx q[0];
rz(-0.65128122) q[0];
sx q[0];
rz(2.5097646) q[0];
rz(-pi) q[1];
rz(1.9252831) q[2];
sx q[2];
rz(-1.7132932) q[2];
sx q[2];
rz(2.6500677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97835097) q[1];
sx q[1];
rz(-0.32838531) q[1];
sx q[1];
rz(2.3204625) q[1];
rz(-0.97614395) q[3];
sx q[3];
rz(-2.757132) q[3];
sx q[3];
rz(-2.773075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0064156) q[2];
sx q[2];
rz(-1.8879031) q[2];
sx q[2];
rz(2.5877) q[2];
rz(1.6484377) q[3];
sx q[3];
rz(-2.0886383) q[3];
sx q[3];
rz(-0.55251399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(3.1062104) q[0];
rz(0.9961876) q[1];
sx q[1];
rz(-1.6085947) q[1];
sx q[1];
rz(-0.48809537) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49022084) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(-0.086829348) q[0];
rz(2.1610291) q[2];
sx q[2];
rz(-2.0681941) q[2];
sx q[2];
rz(2.892708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0636053) q[1];
sx q[1];
rz(-1.6997153) q[1];
sx q[1];
rz(-1.7267666) q[1];
rz(-pi) q[2];
rz(2.5591764) q[3];
sx q[3];
rz(-1.1629259) q[3];
sx q[3];
rz(-1.2553314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7867243) q[2];
sx q[2];
rz(-2.1264666) q[2];
sx q[2];
rz(-0.50746894) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(1.8866084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(-0.3702634) q[0];
rz(0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(2.945074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4465966) q[0];
sx q[0];
rz(-1.2906404) q[0];
sx q[0];
rz(-2.8269935) q[0];
x q[1];
rz(1.1097124) q[2];
sx q[2];
rz(-1.5994674) q[2];
sx q[2];
rz(-1.2210786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6477485) q[1];
sx q[1];
rz(-1.3953679) q[1];
sx q[1];
rz(2.8791048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8653141) q[3];
sx q[3];
rz(-1.4071138) q[3];
sx q[3];
rz(-0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.013441) q[2];
sx q[2];
rz(-0.95169607) q[2];
sx q[2];
rz(-0.91709843) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31345263) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(2.939558) q[0];
rz(-2.0687885) q[1];
sx q[1];
rz(-1.6348811) q[1];
sx q[1];
rz(-1.3938168) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545647) q[0];
sx q[0];
rz(-1.9019039) q[0];
sx q[0];
rz(1.6513255) q[0];
rz(-pi) q[1];
rz(3.0300573) q[2];
sx q[2];
rz(-1.1984014) q[2];
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
rz(1.8762293) q[1];
sx q[1];
rz(-2.6637042) q[1];
sx q[1];
rz(2.1761314) q[1];
rz(0.49683797) q[3];
sx q[3];
rz(-0.79998575) q[3];
sx q[3];
rz(2.5132688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94971925) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(-0.8423155) q[2];
rz(1.0874282) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(-1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6120646) q[0];
sx q[0];
rz(-1.2696711) q[0];
sx q[0];
rz(-1.0213617) q[0];
rz(-1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(2.9313415) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41482571) q[0];
sx q[0];
rz(-1.504577) q[0];
sx q[0];
rz(1.7050752) q[0];
x q[1];
rz(0.16367754) q[2];
sx q[2];
rz(-1.5058766) q[2];
sx q[2];
rz(2.373873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0381895) q[1];
sx q[1];
rz(-1.2545171) q[1];
sx q[1];
rz(-1.7722539) q[1];
x q[2];
rz(-2.5641714) q[3];
sx q[3];
rz(-2.2038659) q[3];
sx q[3];
rz(0.046618332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82688275) q[2];
sx q[2];
rz(-1.3898536) q[2];
sx q[2];
rz(1.38114) q[2];
rz(-0.76210493) q[3];
sx q[3];
rz(-2.6657181) q[3];
sx q[3];
rz(-0.41346082) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.64637) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(-2.6953221) q[1];
sx q[1];
rz(-1.8621567) q[1];
sx q[1];
rz(-2.1648724) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8551089) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(-2.8141862) q[0];
x q[1];
rz(-3.061053) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(-1.5921519) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2446489) q[1];
sx q[1];
rz(-1.9172501) q[1];
sx q[1];
rz(2.0379373) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2817357) q[3];
sx q[3];
rz(-1.1338187) q[3];
sx q[3];
rz(-2.6881998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4010767) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(2.5891499) q[2];
rz(-2.6756514) q[3];
sx q[3];
rz(-1.7539141) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66822806) q[0];
sx q[0];
rz(-0.79376525) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(3.0905511) q[1];
sx q[1];
rz(-1.5850681) q[1];
sx q[1];
rz(2.2424973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102525) q[0];
sx q[0];
rz(-1.9950584) q[0];
sx q[0];
rz(-0.92696855) q[0];
rz(1.0663509) q[2];
sx q[2];
rz(-2.5214508) q[2];
sx q[2];
rz(1.3401741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92210772) q[1];
sx q[1];
rz(-1.2412294) q[1];
sx q[1];
rz(-0.12164128) q[1];
rz(-1.9409688) q[3];
sx q[3];
rz(-1.5709953) q[3];
sx q[3];
rz(-1.9716594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94545323) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(-3.0370039) q[2];
rz(2.3032522) q[3];
sx q[3];
rz(-2.392277) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33912441) q[0];
sx q[0];
rz(-0.81766468) q[0];
sx q[0];
rz(-2.0237645) q[0];
rz(0.71240187) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-2.7446279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3544281) q[0];
sx q[0];
rz(-1.0930601) q[0];
sx q[0];
rz(2.5973395) q[0];
rz(-0.30585395) q[2];
sx q[2];
rz(-1.9216188) q[2];
sx q[2];
rz(0.84301126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9576479) q[1];
sx q[1];
rz(-0.34788222) q[1];
sx q[1];
rz(-1.3089048) q[1];
rz(-1.0730562) q[3];
sx q[3];
rz(-1.2101678) q[3];
sx q[3];
rz(1.3090759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94840702) q[2];
sx q[2];
rz(-0.90017515) q[2];
sx q[2];
rz(3.0449384) q[2];
rz(1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-0.10791049) q[2];
sx q[2];
rz(-1.6549329) q[2];
sx q[2];
rz(-1.8555117) q[2];
rz(-1.8295287) q[3];
sx q[3];
rz(-1.4553634) q[3];
sx q[3];
rz(2.8298557) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
