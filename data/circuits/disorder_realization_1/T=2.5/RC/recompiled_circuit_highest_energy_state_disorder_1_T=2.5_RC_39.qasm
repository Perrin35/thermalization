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
rz(0.83956194) q[0];
sx q[0];
rz(3.8008939) q[0];
sx q[0];
rz(10.492926) q[0];
rz(-2.2358492) q[1];
sx q[1];
rz(-1.716776) q[1];
sx q[1];
rz(-2.4832771) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29355958) q[0];
sx q[0];
rz(-1.2633889) q[0];
sx q[0];
rz(-2.2821167) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80950244) q[2];
sx q[2];
rz(-2.2368119) q[2];
sx q[2];
rz(-1.314338) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37216045) q[1];
sx q[1];
rz(-2.0733412) q[1];
sx q[1];
rz(-2.2519396) q[1];
x q[2];
rz(2.108889) q[3];
sx q[3];
rz(-1.3815839) q[3];
sx q[3];
rz(-3.12836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98755223) q[2];
sx q[2];
rz(-2.0200358) q[2];
sx q[2];
rz(-1.424632) q[2];
rz(-2.3512261) q[3];
sx q[3];
rz(-0.75420403) q[3];
sx q[3];
rz(-2.8823901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454999) q[0];
sx q[0];
rz(-2.6232965) q[0];
sx q[0];
rz(-1.8081283) q[0];
rz(1.2893527) q[1];
sx q[1];
rz(-2.5080296) q[1];
sx q[1];
rz(-3.0234171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8991535) q[0];
sx q[0];
rz(-0.94510539) q[0];
sx q[0];
rz(-3.1345128) q[0];
x q[1];
rz(-1.0071272) q[2];
sx q[2];
rz(-1.360513) q[2];
sx q[2];
rz(2.0521833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0314583) q[1];
sx q[1];
rz(-2.2074568) q[1];
sx q[1];
rz(-2.6183271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2977734) q[3];
sx q[3];
rz(-1.9007933) q[3];
sx q[3];
rz(1.291549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.75225836) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(2.3040859) q[2];
rz(2.5181455) q[3];
sx q[3];
rz(-2.0313171) q[3];
sx q[3];
rz(0.71201223) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74152827) q[0];
sx q[0];
rz(-2.9900804) q[0];
sx q[0];
rz(-2.9140299) q[0];
rz(-0.96861068) q[1];
sx q[1];
rz(-0.89185762) q[1];
sx q[1];
rz(1.6516986) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9609531) q[0];
sx q[0];
rz(-1.102987) q[0];
sx q[0];
rz(2.803928) q[0];
rz(1.3883115) q[2];
sx q[2];
rz(-2.9767155) q[2];
sx q[2];
rz(-0.08872513) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3453133) q[1];
sx q[1];
rz(-1.5070762) q[1];
sx q[1];
rz(-2.696852) q[1];
rz(-pi) q[2];
rz(-1.8996542) q[3];
sx q[3];
rz(-0.76995459) q[3];
sx q[3];
rz(2.9640523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0554793) q[2];
sx q[2];
rz(-1.5986634) q[2];
sx q[2];
rz(2.1231269) q[2];
rz(0.57824072) q[3];
sx q[3];
rz(-0.48931229) q[3];
sx q[3];
rz(-1.9448634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8483491) q[0];
sx q[0];
rz(-1.542955) q[0];
sx q[0];
rz(-1.4620713) q[0];
rz(-0.86817414) q[1];
sx q[1];
rz(-1.9055007) q[1];
sx q[1];
rz(-2.6624534) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7650267) q[0];
sx q[0];
rz(-2.2181795) q[0];
sx q[0];
rz(1.8996973) q[0];
rz(1.2422945) q[2];
sx q[2];
rz(-1.1525894) q[2];
sx q[2];
rz(-0.70786661) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6304178) q[1];
sx q[1];
rz(-1.5197521) q[1];
sx q[1];
rz(-0.45089727) q[1];
x q[2];
rz(-1.8310184) q[3];
sx q[3];
rz(-1.3129915) q[3];
sx q[3];
rz(-1.5618522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9118871) q[2];
sx q[2];
rz(-0.90596002) q[2];
sx q[2];
rz(-2.4521496) q[2];
rz(0.42129579) q[3];
sx q[3];
rz(-2.3861986) q[3];
sx q[3];
rz(-1.1675444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1223849) q[0];
sx q[0];
rz(-0.27007073) q[0];
sx q[0];
rz(-0.58575678) q[0];
rz(-2.6981804) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(-0.74055368) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6950615) q[0];
sx q[0];
rz(-2.2389445) q[0];
sx q[0];
rz(0.88094385) q[0];
x q[1];
rz(-1.9890063) q[2];
sx q[2];
rz(-2.4997992) q[2];
sx q[2];
rz(0.60371232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49879229) q[1];
sx q[1];
rz(-1.8329002) q[1];
sx q[1];
rz(-2.8828055) q[1];
x q[2];
rz(-0.70526975) q[3];
sx q[3];
rz(-0.37644769) q[3];
sx q[3];
rz(1.173526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.022078557) q[2];
sx q[2];
rz(-1.4351279) q[2];
sx q[2];
rz(-0.73149991) q[2];
rz(0.59705192) q[3];
sx q[3];
rz(-2.4686333) q[3];
sx q[3];
rz(1.8895684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.6325833) q[0];
sx q[0];
rz(-0.80061299) q[0];
sx q[0];
rz(1.5923694) q[0];
rz(1.0410615) q[1];
sx q[1];
rz(-0.58716232) q[1];
sx q[1];
rz(0.36137533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3665112) q[0];
sx q[0];
rz(-1.4517227) q[0];
sx q[0];
rz(-2.4111758) q[0];
rz(-3.1388624) q[2];
sx q[2];
rz(-1.3133674) q[2];
sx q[2];
rz(-0.67104679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92897086) q[1];
sx q[1];
rz(-2.5849301) q[1];
sx q[1];
rz(0.61584872) q[1];
x q[2];
rz(-1.2986211) q[3];
sx q[3];
rz(-2.5243773) q[3];
sx q[3];
rz(-0.77780741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3352215) q[2];
sx q[2];
rz(-2.1101895) q[2];
sx q[2];
rz(-1.5058937) q[2];
rz(2.7527572) q[3];
sx q[3];
rz(-2.5917228) q[3];
sx q[3];
rz(2.4839362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1049221) q[0];
sx q[0];
rz(-2.4668283) q[0];
sx q[0];
rz(-2.2210333) q[0];
rz(-2.8432644) q[1];
sx q[1];
rz(-2.0896656) q[1];
sx q[1];
rz(1.1632318) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0981623) q[0];
sx q[0];
rz(-2.8170108) q[0];
sx q[0];
rz(0.030822072) q[0];
rz(-pi) q[1];
rz(0.40980659) q[2];
sx q[2];
rz(-2.3075477) q[2];
sx q[2];
rz(0.17716874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5015848) q[1];
sx q[1];
rz(-2.5951276) q[1];
sx q[1];
rz(0.0025005977) q[1];
x q[2];
rz(-3.0024192) q[3];
sx q[3];
rz(-1.2746029) q[3];
sx q[3];
rz(-2.938066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3433015) q[2];
sx q[2];
rz(-1.2275262) q[2];
sx q[2];
rz(0.11246559) q[2];
rz(0.28137842) q[3];
sx q[3];
rz(-0.76203841) q[3];
sx q[3];
rz(-1.5587829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6211468) q[0];
sx q[0];
rz(-2.1805094) q[0];
sx q[0];
rz(-3.0302826) q[0];
rz(-0.81980199) q[1];
sx q[1];
rz(-2.3078121) q[1];
sx q[1];
rz(0.05096635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539885) q[0];
sx q[0];
rz(-1.4021356) q[0];
sx q[0];
rz(-1.8722562) q[0];
rz(-pi) q[1];
rz(2.5185263) q[2];
sx q[2];
rz(-1.3871117) q[2];
sx q[2];
rz(-0.82090917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2735034) q[1];
sx q[1];
rz(-3.00666) q[1];
sx q[1];
rz(-1.8084303) q[1];
x q[2];
rz(-2.2045129) q[3];
sx q[3];
rz(-1.9670319) q[3];
sx q[3];
rz(1.8609488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.082197949) q[2];
sx q[2];
rz(-1.1092564) q[2];
sx q[2];
rz(1.8471897) q[2];
rz(0.95862499) q[3];
sx q[3];
rz(-0.38161033) q[3];
sx q[3];
rz(-0.81058782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628919) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(-2.1953951) q[0];
rz(-1.2092489) q[1];
sx q[1];
rz(-0.47018662) q[1];
sx q[1];
rz(-1.5399923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9935051) q[0];
sx q[0];
rz(-0.74310857) q[0];
sx q[0];
rz(-1.1412727) q[0];
rz(2.7230324) q[2];
sx q[2];
rz(-2.6113434) q[2];
sx q[2];
rz(-2.6061932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29344236) q[1];
sx q[1];
rz(-0.60876011) q[1];
sx q[1];
rz(0.3984299) q[1];
x q[2];
rz(2.5862365) q[3];
sx q[3];
rz(-1.7512519) q[3];
sx q[3];
rz(0.79977913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20751791) q[2];
sx q[2];
rz(-1.3413651) q[2];
sx q[2];
rz(-1.6259469) q[2];
rz(2.3385284) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(1.1481185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.76252) q[0];
sx q[0];
rz(-1.2305434) q[0];
sx q[0];
rz(-1.8050964) q[0];
rz(1.9876009) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(-1.1755099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0031457) q[0];
sx q[0];
rz(-1.8120736) q[0];
sx q[0];
rz(2.3958415) q[0];
x q[1];
rz(1.8075502) q[2];
sx q[2];
rz(-2.3375698) q[2];
sx q[2];
rz(1.5244816) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.43688135) q[1];
sx q[1];
rz(-2.993104) q[1];
sx q[1];
rz(-0.36098934) q[1];
rz(-pi) q[2];
rz(1.2639579) q[3];
sx q[3];
rz(-0.995702) q[3];
sx q[3];
rz(1.3098615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0846348) q[2];
sx q[2];
rz(-1.4738169) q[2];
sx q[2];
rz(-1.8353621) q[2];
rz(-0.42825395) q[3];
sx q[3];
rz(-1.3166602) q[3];
sx q[3];
rz(0.38124198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8214977) q[0];
sx q[0];
rz(-1.7392673) q[0];
sx q[0];
rz(2.0489954) q[0];
rz(2.8749657) q[1];
sx q[1];
rz(-1.0933924) q[1];
sx q[1];
rz(-0.46270121) q[1];
rz(-2.5960283) q[2];
sx q[2];
rz(-2.5013899) q[2];
sx q[2];
rz(-1.987006) q[2];
rz(-1.7206031) q[3];
sx q[3];
rz(-1.5021742) q[3];
sx q[3];
rz(-1.4175137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
