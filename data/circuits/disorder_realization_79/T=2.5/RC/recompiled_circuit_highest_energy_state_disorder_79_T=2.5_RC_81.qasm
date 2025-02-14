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
rz(-2.4920298) q[0];
sx q[0];
rz(-2.6464638) q[0];
sx q[0];
rz(2.8861217) q[0];
rz(-0.35222346) q[1];
sx q[1];
rz(4.5402266) q[1];
sx q[1];
rz(11.160688) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2737925) q[0];
sx q[0];
rz(-1.5711014) q[0];
sx q[0];
rz(3.1381232) q[0];
x q[1];
rz(-0.5159401) q[2];
sx q[2];
rz(-1.6637515) q[2];
sx q[2];
rz(-3.1414714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1683301) q[1];
sx q[1];
rz(-1.5673866) q[1];
sx q[1];
rz(1.9539321) q[1];
rz(-pi) q[2];
x q[2];
rz(0.051298012) q[3];
sx q[3];
rz(-2.5729532) q[3];
sx q[3];
rz(0.37686791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5070255) q[2];
sx q[2];
rz(-3.1105803) q[2];
sx q[2];
rz(-0.8325038) q[2];
rz(-0.80820525) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(-2.8830849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38043624) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(-0.1880745) q[0];
rz(-0.07218083) q[1];
sx q[1];
rz(-1.0344104) q[1];
sx q[1];
rz(1.5123051) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10764194) q[0];
sx q[0];
rz(-2.788262) q[0];
sx q[0];
rz(-2.8046542) q[0];
rz(1.5248271) q[2];
sx q[2];
rz(-1.6013923) q[2];
sx q[2];
rz(-2.033288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9058414) q[1];
sx q[1];
rz(-3.0260147) q[1];
sx q[1];
rz(0.53821941) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80792369) q[3];
sx q[3];
rz(-1.9821385) q[3];
sx q[3];
rz(1.5613256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9709116) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(-1.8485273) q[2];
rz(2.7944148) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(-3.0612883) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8186571) q[0];
sx q[0];
rz(-1.289239) q[0];
sx q[0];
rz(0.47054189) q[0];
rz(-1.200354) q[1];
sx q[1];
rz(-0.73089868) q[1];
sx q[1];
rz(-1.2568731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6171489) q[0];
sx q[0];
rz(-1.9223419) q[0];
sx q[0];
rz(2.1407513) q[0];
x q[1];
rz(0.09229163) q[2];
sx q[2];
rz(-1.7053635) q[2];
sx q[2];
rz(-1.2876266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4670406) q[1];
sx q[1];
rz(-1.5853154) q[1];
sx q[1];
rz(3.1402223) q[1];
rz(-pi) q[2];
rz(1.1670349) q[3];
sx q[3];
rz(-2.4903273) q[3];
sx q[3];
rz(2.2869956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54258004) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(-0.92213321) q[2];
rz(-1.3433836) q[3];
sx q[3];
rz(-2.193439) q[3];
sx q[3];
rz(-1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682997) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(0.44019765) q[0];
rz(-1.531155) q[1];
sx q[1];
rz(-1.4819769) q[1];
sx q[1];
rz(0.24756113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.089624) q[0];
sx q[0];
rz(-2.0735093) q[0];
sx q[0];
rz(-0.49552709) q[0];
rz(1.6610959) q[2];
sx q[2];
rz(-1.4019499) q[2];
sx q[2];
rz(-0.24659469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2046856) q[1];
sx q[1];
rz(-1.2720894) q[1];
sx q[1];
rz(-0.077110962) q[1];
x q[2];
rz(-0.04051716) q[3];
sx q[3];
rz(-0.64896331) q[3];
sx q[3];
rz(1.9891091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89746785) q[2];
sx q[2];
rz(-1.030587) q[2];
sx q[2];
rz(1.6991276) q[2];
rz(0.20081946) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(-1.1403181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15631256) q[0];
sx q[0];
rz(-2.9887178) q[0];
sx q[0];
rz(-0.54541624) q[0];
rz(0.044005752) q[1];
sx q[1];
rz(-0.017887201) q[1];
sx q[1];
rz(0.47637475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2891083) q[0];
sx q[0];
rz(-1.327564) q[0];
sx q[0];
rz(1.4840675) q[0];
rz(-pi) q[1];
rz(-0.77875359) q[2];
sx q[2];
rz(-0.54355946) q[2];
sx q[2];
rz(-1.5875848) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.093845) q[1];
sx q[1];
rz(-0.72870164) q[1];
sx q[1];
rz(-1.8456259) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3333605) q[3];
sx q[3];
rz(-2.6202706) q[3];
sx q[3];
rz(0.1791187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62006092) q[2];
sx q[2];
rz(-2.4510577) q[2];
sx q[2];
rz(0.36038348) q[2];
rz(-0.22200577) q[3];
sx q[3];
rz(-1.5304151) q[3];
sx q[3];
rz(-1.7300026) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5510657) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(1.5850413) q[0];
rz(-0.12697728) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(0.089687673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28676867) q[0];
sx q[0];
rz(-0.78203219) q[0];
sx q[0];
rz(0.8013335) q[0];
rz(-pi) q[1];
rz(0.96025567) q[2];
sx q[2];
rz(-0.88443236) q[2];
sx q[2];
rz(-0.45623623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8627734) q[1];
sx q[1];
rz(-1.8012281) q[1];
sx q[1];
rz(2.6342175) q[1];
rz(-pi) q[2];
rz(0.50181194) q[3];
sx q[3];
rz(-0.66365051) q[3];
sx q[3];
rz(0.36728951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.108532) q[2];
sx q[2];
rz(-0.33691275) q[2];
sx q[2];
rz(2.8899371) q[2];
rz(0.039994914) q[3];
sx q[3];
rz(-2.7997041) q[3];
sx q[3];
rz(-2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104178) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(-0.41496667) q[0];
rz(1.4893432) q[1];
sx q[1];
rz(-3.0840315) q[1];
sx q[1];
rz(-2.8156978) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0123169) q[0];
sx q[0];
rz(-1.1291478) q[0];
sx q[0];
rz(2.737816) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.089881534) q[2];
sx q[2];
rz(-1.0542068) q[2];
sx q[2];
rz(-2.3013249) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.53908097) q[1];
sx q[1];
rz(-1.8274787) q[1];
sx q[1];
rz(-2.6246482) q[1];
x q[2];
rz(1.8091466) q[3];
sx q[3];
rz(-0.311626) q[3];
sx q[3];
rz(0.99743227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2275647) q[2];
sx q[2];
rz(-1.808337) q[2];
sx q[2];
rz(-2.0951159) q[2];
rz(2.922831) q[3];
sx q[3];
rz(-2.1121787) q[3];
sx q[3];
rz(1.7441162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47531146) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(2.6783491) q[0];
rz(2.8766368) q[1];
sx q[1];
rz(-3.1400561) q[1];
sx q[1];
rz(-1.6418246) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87141052) q[0];
sx q[0];
rz(-2.0034916) q[0];
sx q[0];
rz(-0.020633296) q[0];
x q[1];
rz(1.5157264) q[2];
sx q[2];
rz(-1.1257505) q[2];
sx q[2];
rz(2.9663756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25019885) q[1];
sx q[1];
rz(-1.3490744) q[1];
sx q[1];
rz(-3.1117597) q[1];
rz(-pi) q[2];
rz(1.8345233) q[3];
sx q[3];
rz(-2.8805507) q[3];
sx q[3];
rz(1.8393976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2396592) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(-1.8233914) q[2];
rz(-0.0031331172) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(-1.1656632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5636633) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(-2.4270015) q[0];
rz(2.928012) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.2108796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058556341) q[0];
sx q[0];
rz(-2.8175111) q[0];
sx q[0];
rz(1.9084318) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7316852) q[2];
sx q[2];
rz(-2.4033467) q[2];
sx q[2];
rz(1.8869836) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80187243) q[1];
sx q[1];
rz(-1.045671) q[1];
sx q[1];
rz(2.8169291) q[1];
x q[2];
rz(-2.9835002) q[3];
sx q[3];
rz(-1.3835581) q[3];
sx q[3];
rz(-0.15915376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96656281) q[2];
sx q[2];
rz(-0.021447072) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(0.35061947) q[3];
sx q[3];
rz(-1.5712761) q[3];
sx q[3];
rz(1.143379) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458934) q[0];
sx q[0];
rz(-0.92246711) q[0];
sx q[0];
rz(-2.5272227) q[0];
rz(-3.0820097) q[1];
sx q[1];
rz(-0.091484286) q[1];
sx q[1];
rz(1.7170067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88286582) q[0];
sx q[0];
rz(-0.75260163) q[0];
sx q[0];
rz(-1.9310139) q[0];
rz(-pi) q[1];
rz(-0.015083357) q[2];
sx q[2];
rz(-3.0102979) q[2];
sx q[2];
rz(-0.35485902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.0265621) q[1];
sx q[1];
rz(-1.4938338) q[1];
sx q[1];
rz(1.413024) q[1];
rz(-pi) q[2];
rz(2.0395117) q[3];
sx q[3];
rz(-0.45022435) q[3];
sx q[3];
rz(-2.1463822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0397348) q[2];
sx q[2];
rz(-2.8701344) q[2];
sx q[2];
rz(-2.6177935) q[2];
rz(-1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(0.4642134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6687748) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(1.7300023) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(0.642943) q[2];
sx q[2];
rz(-1.8131154) q[2];
sx q[2];
rz(-2.9396069) q[2];
rz(2.7649391) q[3];
sx q[3];
rz(-0.49643825) q[3];
sx q[3];
rz(1.7986922) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
