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
rz(1.3069557) q[0];
sx q[0];
rz(-2.262158) q[0];
sx q[0];
rz(0.49361324) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8956229) q[0];
sx q[0];
rz(-2.1969271) q[0];
sx q[0];
rz(2.5139721) q[0];
rz(-pi) q[1];
rz(0.25716146) q[2];
sx q[2];
rz(-2.2182121) q[2];
sx q[2];
rz(-2.0055298) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8779187) q[1];
sx q[1];
rz(-1.3527737) q[1];
sx q[1];
rz(2.430357) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93500626) q[3];
sx q[3];
rz(-1.8719487) q[3];
sx q[3];
rz(3.0060535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1349858) q[2];
sx q[2];
rz(-1.1400433) q[2];
sx q[2];
rz(-2.0331649) q[2];
rz(1.8290352) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(0.31497064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38158622) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(-0.015856892) q[0];
rz(-0.99041692) q[1];
sx q[1];
rz(-0.76737338) q[1];
sx q[1];
rz(-2.2784065) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8548633) q[0];
sx q[0];
rz(-2.3026288) q[0];
sx q[0];
rz(1.7006204) q[0];
rz(-pi) q[1];
rz(1.8998186) q[2];
sx q[2];
rz(-1.8890427) q[2];
sx q[2];
rz(-1.9538188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3962209) q[1];
sx q[1];
rz(-2.6156151) q[1];
sx q[1];
rz(0.35244103) q[1];
rz(2.2926278) q[3];
sx q[3];
rz(-1.6833653) q[3];
sx q[3];
rz(1.4813042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70011675) q[2];
sx q[2];
rz(-2.8642544) q[2];
sx q[2];
rz(0.46879834) q[2];
rz(-0.68764728) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.1402635) q[0];
sx q[0];
rz(-0.92324531) q[0];
sx q[0];
rz(-2.4053251) q[0];
rz(0.34321347) q[1];
sx q[1];
rz(-2.2540269) q[1];
sx q[1];
rz(0.18377486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4030754) q[0];
sx q[0];
rz(-2.3949497) q[0];
sx q[0];
rz(0.097642032) q[0];
rz(-pi) q[1];
rz(0.69664069) q[2];
sx q[2];
rz(-1.9844779) q[2];
sx q[2];
rz(2.6282642) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50810516) q[1];
sx q[1];
rz(-2.268666) q[1];
sx q[1];
rz(1.9644323) q[1];
x q[2];
rz(1.0314437) q[3];
sx q[3];
rz(-1.8832409) q[3];
sx q[3];
rz(2.7613918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3595714) q[2];
sx q[2];
rz(-2.6247793) q[2];
sx q[2];
rz(-0.51592958) q[2];
rz(2.0333717) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(1.9903323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9786523) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(0.51280713) q[0];
rz(-0.45979744) q[1];
sx q[1];
rz(-1.5492946) q[1];
sx q[1];
rz(-0.76382452) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490246) q[0];
sx q[0];
rz(-1.4563174) q[0];
sx q[0];
rz(-3.0106198) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2470719) q[2];
sx q[2];
rz(-1.9678332) q[2];
sx q[2];
rz(-2.4415602) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0415548) q[1];
sx q[1];
rz(-0.44719346) q[1];
sx q[1];
rz(0.87765043) q[1];
x q[2];
rz(-2.7598792) q[3];
sx q[3];
rz(-0.36892051) q[3];
sx q[3];
rz(0.58894181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0189455) q[2];
sx q[2];
rz(-0.15572369) q[2];
sx q[2];
rz(0.32749185) q[2];
rz(0.28199768) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(-1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649488) q[0];
sx q[0];
rz(-0.38589859) q[0];
sx q[0];
rz(0.97736812) q[0];
rz(0.36879677) q[1];
sx q[1];
rz(-0.86124033) q[1];
sx q[1];
rz(1.7288953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1288797) q[0];
sx q[0];
rz(-2.0600208) q[0];
sx q[0];
rz(2.9262744) q[0];
x q[1];
rz(1.5951889) q[2];
sx q[2];
rz(-2.6246723) q[2];
sx q[2];
rz(0.90205294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2569132) q[1];
sx q[1];
rz(-1.8533145) q[1];
sx q[1];
rz(0.032508548) q[1];
x q[2];
rz(-2.5630412) q[3];
sx q[3];
rz(-1.2935841) q[3];
sx q[3];
rz(-1.8892242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8480924) q[2];
sx q[2];
rz(-0.81587452) q[2];
sx q[2];
rz(2.9917742) q[2];
rz(2.7430429) q[3];
sx q[3];
rz(-1.6179251) q[3];
sx q[3];
rz(-2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3715816) q[0];
sx q[0];
rz(-0.12938975) q[0];
sx q[0];
rz(-2.5883801) q[0];
rz(2.5999056) q[1];
sx q[1];
rz(-2.0952416) q[1];
sx q[1];
rz(-0.4479301) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.124156) q[0];
sx q[0];
rz(-1.4261591) q[0];
sx q[0];
rz(2.7516624) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42945736) q[2];
sx q[2];
rz(-1.9782559) q[2];
sx q[2];
rz(0.40312672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6363931) q[1];
sx q[1];
rz(-1.4318887) q[1];
sx q[1];
rz(1.9646775) q[1];
x q[2];
rz(-0.33282307) q[3];
sx q[3];
rz(-1.5170005) q[3];
sx q[3];
rz(-1.7032361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5003659) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(2.3424303) q[2];
rz(2.7568119) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(2.1414115) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7810829) q[0];
sx q[0];
rz(-2.0979083) q[0];
sx q[0];
rz(0.59180301) q[0];
rz(2.7599755) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(-0.15242481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7397933) q[0];
sx q[0];
rz(-1.8795947) q[0];
sx q[0];
rz(1.0303506) q[0];
x q[1];
rz(2.6640058) q[2];
sx q[2];
rz(-1.7268333) q[2];
sx q[2];
rz(-2.0618771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.287474) q[1];
sx q[1];
rz(-1.2240613) q[1];
sx q[1];
rz(1.5902014) q[1];
rz(-2.3197738) q[3];
sx q[3];
rz(-0.83335175) q[3];
sx q[3];
rz(-0.59274537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3064208) q[2];
sx q[2];
rz(-2.1495337) q[2];
sx q[2];
rz(1.6888118) q[2];
rz(-2.6460904) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(-2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.751916) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(0.16146846) q[0];
rz(0.031919315) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(-1.2772824) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.204435) q[0];
sx q[0];
rz(-1.3124518) q[0];
sx q[0];
rz(-2.4169371) q[0];
x q[1];
rz(2.0525682) q[2];
sx q[2];
rz(-2.6725997) q[2];
sx q[2];
rz(-2.2828141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1254221) q[1];
sx q[1];
rz(-1.6787663) q[1];
sx q[1];
rz(0.30534621) q[1];
x q[2];
rz(1.6971436) q[3];
sx q[3];
rz(-0.29623162) q[3];
sx q[3];
rz(1.1100696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6892467) q[2];
sx q[2];
rz(-0.9114868) q[2];
sx q[2];
rz(2.7455043) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(2.2649435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24092995) q[0];
sx q[0];
rz(-3.1234968) q[0];
sx q[0];
rz(-0.15116365) q[0];
rz(2.1647272) q[1];
sx q[1];
rz(-0.38115373) q[1];
sx q[1];
rz(-2.3968598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3522986) q[0];
sx q[0];
rz(-1.2937102) q[0];
sx q[0];
rz(0.21488551) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2538386) q[2];
sx q[2];
rz(-2.709667) q[2];
sx q[2];
rz(-0.077684075) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7691466) q[1];
sx q[1];
rz(-0.92248864) q[1];
sx q[1];
rz(1.5725971) q[1];
x q[2];
rz(-1.3296534) q[3];
sx q[3];
rz(-1.9134812) q[3];
sx q[3];
rz(1.2872651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(2.1935479) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(2.0675366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093216151) q[0];
sx q[0];
rz(-2.8792448) q[0];
sx q[0];
rz(-2.9076305) q[0];
rz(-1.1853064) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(-1.8918461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3435182) q[0];
sx q[0];
rz(-1.8600704) q[0];
sx q[0];
rz(-0.57855655) q[0];
x q[1];
rz(0.62057497) q[2];
sx q[2];
rz(-1.7399551) q[2];
sx q[2];
rz(-2.754654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1703012) q[1];
sx q[1];
rz(-0.93496041) q[1];
sx q[1];
rz(0.88199349) q[1];
rz(-pi) q[2];
rz(1.7464645) q[3];
sx q[3];
rz(-0.89699706) q[3];
sx q[3];
rz(1.2487981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11448161) q[2];
sx q[2];
rz(-1.9805084) q[2];
sx q[2];
rz(-1.4410045) q[2];
rz(-0.13127413) q[3];
sx q[3];
rz(-2.6797397) q[3];
sx q[3];
rz(2.89768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048653614) q[0];
sx q[0];
rz(-2.1760512) q[0];
sx q[0];
rz(-2.0078134) q[0];
rz(1.0406021) q[1];
sx q[1];
rz(-0.84747172) q[1];
sx q[1];
rz(-0.52451959) q[1];
rz(1.9807627) q[2];
sx q[2];
rz(-0.66902918) q[2];
sx q[2];
rz(1.4794028) q[2];
rz(1.7215988) q[3];
sx q[3];
rz(-0.94298132) q[3];
sx q[3];
rz(1.0366221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
