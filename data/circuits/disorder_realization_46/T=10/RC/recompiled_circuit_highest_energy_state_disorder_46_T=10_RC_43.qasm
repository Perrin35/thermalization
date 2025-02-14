OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.96002785) q[0];
sx q[0];
rz(-2.5224944) q[0];
sx q[0];
rz(1.7382789) q[0];
rz(2.4461441) q[1];
sx q[1];
rz(-0.56013501) q[1];
sx q[1];
rz(-2.4897895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90690985) q[0];
sx q[0];
rz(-1.9056589) q[0];
sx q[0];
rz(-0.85483179) q[0];
rz(-0.02044631) q[2];
sx q[2];
rz(-2.8600577) q[2];
sx q[2];
rz(-2.881518) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6095602) q[1];
sx q[1];
rz(-1.3390307) q[1];
sx q[1];
rz(1.5941117) q[1];
rz(-pi) q[2];
rz(-1.2816209) q[3];
sx q[3];
rz(-0.77176982) q[3];
sx q[3];
rz(-1.650633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2191676) q[2];
sx q[2];
rz(-1.3895915) q[2];
sx q[2];
rz(0.31828848) q[2];
rz(-0.90218941) q[3];
sx q[3];
rz(-2.2212494) q[3];
sx q[3];
rz(-0.86377803) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0133463) q[0];
sx q[0];
rz(-2.8885169) q[0];
sx q[0];
rz(2.6213562) q[0];
rz(-3.0546313) q[1];
sx q[1];
rz(-2.088701) q[1];
sx q[1];
rz(-0.797995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8142321) q[0];
sx q[0];
rz(-1.8810187) q[0];
sx q[0];
rz(-1.7807175) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9676759) q[2];
sx q[2];
rz(-2.3905281) q[2];
sx q[2];
rz(-2.8349595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50988156) q[1];
sx q[1];
rz(-2.7260511) q[1];
sx q[1];
rz(2.7870537) q[1];
x q[2];
rz(-1.9738484) q[3];
sx q[3];
rz(-1.9857996) q[3];
sx q[3];
rz(0.61674207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.308455) q[2];
sx q[2];
rz(-1.986958) q[2];
sx q[2];
rz(-0.40033611) q[2];
rz(-0.18276246) q[3];
sx q[3];
rz(-2.5623613) q[3];
sx q[3];
rz(-1.8496752) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12277814) q[0];
sx q[0];
rz(-3.0570539) q[0];
sx q[0];
rz(-1.9345181) q[0];
rz(1.5030376) q[1];
sx q[1];
rz(-1.8820347) q[1];
sx q[1];
rz(-0.17301339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2234576) q[0];
sx q[0];
rz(-1.6561145) q[0];
sx q[0];
rz(-1.9568018) q[0];
rz(3.0315057) q[2];
sx q[2];
rz(-2.1497257) q[2];
sx q[2];
rz(0.26012173) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.26721482) q[1];
sx q[1];
rz(-2.3372531) q[1];
sx q[1];
rz(0.72093236) q[1];
x q[2];
rz(0.2263308) q[3];
sx q[3];
rz(-2.3809098) q[3];
sx q[3];
rz(0.62225354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59306899) q[2];
sx q[2];
rz(-2.5544281) q[2];
sx q[2];
rz(-2.1068088) q[2];
rz(-2.2506524) q[3];
sx q[3];
rz(-0.8689298) q[3];
sx q[3];
rz(-1.2263891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0767176) q[0];
sx q[0];
rz(-2.1893976) q[0];
sx q[0];
rz(-0.12884831) q[0];
rz(2.7128362) q[1];
sx q[1];
rz(-1.4219475) q[1];
sx q[1];
rz(-1.6645924) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7744666) q[0];
sx q[0];
rz(-1.9156349) q[0];
sx q[0];
rz(-1.4982629) q[0];
rz(-pi) q[1];
rz(2.693839) q[2];
sx q[2];
rz(-1.9734255) q[2];
sx q[2];
rz(-1.3626217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4739371) q[1];
sx q[1];
rz(-1.9868816) q[1];
sx q[1];
rz(0.23148424) q[1];
rz(3.1404895) q[3];
sx q[3];
rz(-2.1620521) q[3];
sx q[3];
rz(1.3659034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9513272) q[2];
sx q[2];
rz(-1.724388) q[2];
sx q[2];
rz(-1.3235271) q[2];
rz(1.4642814) q[3];
sx q[3];
rz(-2.5530294) q[3];
sx q[3];
rz(2.467449) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4071963) q[0];
sx q[0];
rz(-2.7792327) q[0];
sx q[0];
rz(-1.4417484) q[0];
rz(-0.4231407) q[1];
sx q[1];
rz(-0.95760456) q[1];
sx q[1];
rz(0.29022455) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26849982) q[0];
sx q[0];
rz(-2.1999176) q[0];
sx q[0];
rz(1.9353348) q[0];
x q[1];
rz(-2.5686479) q[2];
sx q[2];
rz(-2.101295) q[2];
sx q[2];
rz(0.45924074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6001989) q[1];
sx q[1];
rz(-3.0046084) q[1];
sx q[1];
rz(2.7776657) q[1];
x q[2];
rz(-2.8338088) q[3];
sx q[3];
rz(-1.8202542) q[3];
sx q[3];
rz(1.4741999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4238802) q[2];
sx q[2];
rz(-0.78705698) q[2];
sx q[2];
rz(-2.1948953) q[2];
rz(0.26149073) q[3];
sx q[3];
rz(-1.3638834) q[3];
sx q[3];
rz(-0.3869032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70447266) q[0];
sx q[0];
rz(-1.8081212) q[0];
sx q[0];
rz(-2.3719846) q[0];
rz(2.9656124) q[1];
sx q[1];
rz(-1.8006005) q[1];
sx q[1];
rz(-1.5699068) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.044445) q[0];
sx q[0];
rz(-2.189069) q[0];
sx q[0];
rz(-0.78938214) q[0];
x q[1];
rz(-1.6506344) q[2];
sx q[2];
rz(-0.72070272) q[2];
sx q[2];
rz(2.0699163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1842725) q[1];
sx q[1];
rz(-1.2046923) q[1];
sx q[1];
rz(-1.7680607) q[1];
x q[2];
rz(-2.3442101) q[3];
sx q[3];
rz(-2.4719072) q[3];
sx q[3];
rz(2.0978239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1311329) q[2];
sx q[2];
rz(-2.185148) q[2];
sx q[2];
rz(1.3156816) q[2];
rz(-1.2223505) q[3];
sx q[3];
rz(-2.0171916) q[3];
sx q[3];
rz(2.840672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-2.9199801) q[0];
sx q[0];
rz(-1.6918809) q[0];
sx q[0];
rz(-1.8020887) q[0];
rz(-1.1969396) q[1];
sx q[1];
rz(-2.0868128) q[1];
sx q[1];
rz(-2.1766677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366936) q[0];
sx q[0];
rz(-1.6654547) q[0];
sx q[0];
rz(2.6306334) q[0];
x q[1];
rz(0.2946202) q[2];
sx q[2];
rz(-1.9082532) q[2];
sx q[2];
rz(1.9532415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5894789) q[1];
sx q[1];
rz(-0.2355816) q[1];
sx q[1];
rz(2.7757591) q[1];
rz(-2.6581702) q[3];
sx q[3];
rz(-1.7258378) q[3];
sx q[3];
rz(0.29696908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40184608) q[2];
sx q[2];
rz(-1.9332989) q[2];
sx q[2];
rz(1.9604663) q[2];
rz(-3.127023) q[3];
sx q[3];
rz(-2.7404692) q[3];
sx q[3];
rz(-1.9692263) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75083098) q[0];
sx q[0];
rz(-1.2068692) q[0];
sx q[0];
rz(-2.3857351) q[0];
rz(0.66468704) q[1];
sx q[1];
rz(-1.942778) q[1];
sx q[1];
rz(-1.8519148) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15890789) q[0];
sx q[0];
rz(-2.3450627) q[0];
sx q[0];
rz(0.089287306) q[0];
rz(-pi) q[1];
rz(-0.47814887) q[2];
sx q[2];
rz(-0.59403803) q[2];
sx q[2];
rz(-0.65982925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7482704) q[1];
sx q[1];
rz(-0.38667187) q[1];
sx q[1];
rz(2.3063117) q[1];
rz(-2.8347551) q[3];
sx q[3];
rz(-2.8182013) q[3];
sx q[3];
rz(-2.8093114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0083302) q[2];
sx q[2];
rz(-1.3349168) q[2];
sx q[2];
rz(1.5745715) q[2];
rz(1.3957006) q[3];
sx q[3];
rz(-1.7392266) q[3];
sx q[3];
rz(-0.76656669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.72346) q[0];
sx q[0];
rz(-2.8657931) q[0];
sx q[0];
rz(2.7511399) q[0];
rz(-2.0402724) q[1];
sx q[1];
rz(-1.3464255) q[1];
sx q[1];
rz(-0.93213814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2208122) q[0];
sx q[0];
rz(-1.5340065) q[0];
sx q[0];
rz(3.1139741) q[0];
rz(-pi) q[1];
rz(-1.522077) q[2];
sx q[2];
rz(-2.1095968) q[2];
sx q[2];
rz(-0.032291238) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.336672) q[1];
sx q[1];
rz(-1.3778566) q[1];
sx q[1];
rz(2.8590917) q[1];
rz(-0.45472179) q[3];
sx q[3];
rz(-0.51815663) q[3];
sx q[3];
rz(-1.3738812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0497047) q[2];
sx q[2];
rz(-2.4447417) q[2];
sx q[2];
rz(1.6799124) q[2];
rz(-1.0852496) q[3];
sx q[3];
rz(-2.0094252) q[3];
sx q[3];
rz(0.6404883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.853448) q[0];
sx q[0];
rz(-0.65971056) q[0];
sx q[0];
rz(-1.1261384) q[0];
rz(2.0938734) q[1];
sx q[1];
rz(-0.62647096) q[1];
sx q[1];
rz(-1.610021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60309339) q[0];
sx q[0];
rz(-1.5694008) q[0];
sx q[0];
rz(3.1183232) q[0];
rz(-1.3381027) q[2];
sx q[2];
rz(-2.0183211) q[2];
sx q[2];
rz(1.9767424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23352392) q[1];
sx q[1];
rz(-1.6485813) q[1];
sx q[1];
rz(0.064957182) q[1];
x q[2];
rz(-1.6592386) q[3];
sx q[3];
rz(-1.0103265) q[3];
sx q[3];
rz(1.1286398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5213726) q[2];
sx q[2];
rz(-0.81338254) q[2];
sx q[2];
rz(0.20624557) q[2];
rz(-0.5591875) q[3];
sx q[3];
rz(-1.5234448) q[3];
sx q[3];
rz(-2.0172113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6013721) q[0];
sx q[0];
rz(-1.843353) q[0];
sx q[0];
rz(-2.699615) q[0];
rz(-1.9961009) q[1];
sx q[1];
rz(-1.3064697) q[1];
sx q[1];
rz(-1.3189955) q[1];
rz(2.8679642) q[2];
sx q[2];
rz(-1.3009334) q[2];
sx q[2];
rz(1.0197659) q[2];
rz(2.4174684) q[3];
sx q[3];
rz(-1.4403314) q[3];
sx q[3];
rz(-2.452313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
