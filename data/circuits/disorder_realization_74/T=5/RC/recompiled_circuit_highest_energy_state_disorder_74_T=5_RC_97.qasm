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
rz(-0.44218818) q[0];
sx q[0];
rz(-1.1531885) q[0];
sx q[0];
rz(3.0129504) q[0];
rz(1.6755942) q[1];
sx q[1];
rz(-2.0425551) q[1];
sx q[1];
rz(1.0817945) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32444123) q[0];
sx q[0];
rz(-1.6015581) q[0];
sx q[0];
rz(0.064662393) q[0];
x q[1];
rz(2.2954582) q[2];
sx q[2];
rz(-2.6577009) q[2];
sx q[2];
rz(2.0563375) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5521277) q[1];
sx q[1];
rz(-1.4120483) q[1];
sx q[1];
rz(0.071539025) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2957828) q[3];
sx q[3];
rz(-1.661947) q[3];
sx q[3];
rz(0.22991523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0470524) q[2];
sx q[2];
rz(-0.69539842) q[2];
sx q[2];
rz(2.409234) q[2];
rz(-0.098946027) q[3];
sx q[3];
rz(-1.0354592) q[3];
sx q[3];
rz(-1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089379646) q[0];
sx q[0];
rz(-0.78240028) q[0];
sx q[0];
rz(0.041444929) q[0];
rz(-0.6913569) q[1];
sx q[1];
rz(-1.8212916) q[1];
sx q[1];
rz(2.2533805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81538686) q[0];
sx q[0];
rz(-2.0688217) q[0];
sx q[0];
rz(-2.2603358) q[0];
rz(-2.5609547) q[2];
sx q[2];
rz(-2.3742445) q[2];
sx q[2];
rz(-1.5661756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2728558) q[1];
sx q[1];
rz(-1.6140198) q[1];
sx q[1];
rz(-0.24300139) q[1];
rz(-pi) q[2];
rz(-0.51155297) q[3];
sx q[3];
rz(-1.8163876) q[3];
sx q[3];
rz(0.052215271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1516252) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(-1.3803231) q[2];
rz(-0.049792854) q[3];
sx q[3];
rz(-0.68792611) q[3];
sx q[3];
rz(0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3675073) q[0];
sx q[0];
rz(-2.9893576) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(-1.3146575) q[1];
sx q[1];
rz(-1.0379125) q[1];
sx q[1];
rz(1.5705869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0308098) q[0];
sx q[0];
rz(-0.45814308) q[0];
sx q[0];
rz(-0.58266298) q[0];
rz(1.7271829) q[2];
sx q[2];
rz(-1.9769353) q[2];
sx q[2];
rz(2.0942424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.27109001) q[1];
sx q[1];
rz(-1.5368286) q[1];
sx q[1];
rz(1.7351331) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3486628) q[3];
sx q[3];
rz(-0.76628162) q[3];
sx q[3];
rz(-2.7310581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7870002) q[2];
sx q[2];
rz(-2.3921693) q[2];
sx q[2];
rz(-0.3053537) q[2];
rz(-0.8461771) q[3];
sx q[3];
rz(-1.2140707) q[3];
sx q[3];
rz(2.2727216) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531891) q[0];
sx q[0];
rz(-2.1222293) q[0];
sx q[0];
rz(0.97002059) q[0];
rz(-2.8963529) q[1];
sx q[1];
rz(-1.0673362) q[1];
sx q[1];
rz(-1.5018357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0787112) q[0];
sx q[0];
rz(-2.9406266) q[0];
sx q[0];
rz(-0.68931343) q[0];
rz(-1.7443329) q[2];
sx q[2];
rz(-1.8029986) q[2];
sx q[2];
rz(-2.6628138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1388772) q[1];
sx q[1];
rz(-1.7600696) q[1];
sx q[1];
rz(1.4156962) q[1];
x q[2];
rz(-0.72704859) q[3];
sx q[3];
rz(-1.3241951) q[3];
sx q[3];
rz(1.5545157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60845) q[2];
sx q[2];
rz(-1.5641944) q[2];
sx q[2];
rz(-1.2987785) q[2];
rz(-0.53457824) q[3];
sx q[3];
rz(-2.5375073) q[3];
sx q[3];
rz(2.4728313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088148549) q[0];
sx q[0];
rz(-1.3548594) q[0];
sx q[0];
rz(2.163929) q[0];
rz(2.4560302) q[1];
sx q[1];
rz(-2.3473163) q[1];
sx q[1];
rz(2.175144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880497) q[0];
sx q[0];
rz(-0.85659638) q[0];
sx q[0];
rz(0.62936737) q[0];
x q[1];
rz(1.9990997) q[2];
sx q[2];
rz(-2.097528) q[2];
sx q[2];
rz(1.525804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2144262) q[1];
sx q[1];
rz(-2.4220903) q[1];
sx q[1];
rz(1.359316) q[1];
rz(1.9240941) q[3];
sx q[3];
rz(-1.0477001) q[3];
sx q[3];
rz(2.9016665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5886249) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(2.8589613) q[2];
rz(2.1770554) q[3];
sx q[3];
rz(-1.5610361) q[3];
sx q[3];
rz(0.36981043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.033217) q[0];
sx q[0];
rz(-0.41154698) q[0];
sx q[0];
rz(-2.0630398) q[0];
rz(-0.82303965) q[1];
sx q[1];
rz(-2.0185202) q[1];
sx q[1];
rz(-0.80026921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5259994) q[0];
sx q[0];
rz(-1.5999927) q[0];
sx q[0];
rz(0.008511624) q[0];
rz(1.2563317) q[2];
sx q[2];
rz(-1.8654746) q[2];
sx q[2];
rz(-1.8626358) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89085244) q[1];
sx q[1];
rz(-1.0167334) q[1];
sx q[1];
rz(2.3046369) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5979302) q[3];
sx q[3];
rz(-2.398536) q[3];
sx q[3];
rz(-2.0280968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1396973) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(0.099420698) q[2];
rz(2.5771778) q[3];
sx q[3];
rz(-1.5921009) q[3];
sx q[3];
rz(2.2216589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99892202) q[0];
sx q[0];
rz(-0.26695928) q[0];
sx q[0];
rz(3.1166742) q[0];
rz(2.0492367) q[1];
sx q[1];
rz(-1.151029) q[1];
sx q[1];
rz(0.54164642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3643862) q[0];
sx q[0];
rz(-1.2105252) q[0];
sx q[0];
rz(-2.2132753) q[0];
rz(-0.32633467) q[2];
sx q[2];
rz(-1.0251364) q[2];
sx q[2];
rz(1.5424002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6691963) q[1];
sx q[1];
rz(-1.4265359) q[1];
sx q[1];
rz(1.902182) q[1];
rz(-2.9639313) q[3];
sx q[3];
rz(-2.170106) q[3];
sx q[3];
rz(2.8467972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3148552) q[2];
sx q[2];
rz(-1.9172226) q[2];
sx q[2];
rz(-2.0014191) q[2];
rz(1.7013811) q[3];
sx q[3];
rz(-0.39339104) q[3];
sx q[3];
rz(2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88406554) q[0];
sx q[0];
rz(-1.1648357) q[0];
sx q[0];
rz(2.2648532) q[0];
rz(0.17403099) q[1];
sx q[1];
rz(-2.048025) q[1];
sx q[1];
rz(-1.7736951) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48094952) q[0];
sx q[0];
rz(-2.5202453) q[0];
sx q[0];
rz(-2.0578007) q[0];
rz(-1.935101) q[2];
sx q[2];
rz(-1.7839396) q[2];
sx q[2];
rz(-3.0618966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9570043) q[1];
sx q[1];
rz(-0.68268973) q[1];
sx q[1];
rz(-0.95889389) q[1];
rz(-pi) q[2];
rz(-2.1266227) q[3];
sx q[3];
rz(-0.14326247) q[3];
sx q[3];
rz(1.48326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4725264) q[2];
sx q[2];
rz(-1.8146699) q[2];
sx q[2];
rz(2.9533022) q[2];
rz(1.3663728) q[3];
sx q[3];
rz(-1.7732737) q[3];
sx q[3];
rz(1.202047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61757225) q[0];
sx q[0];
rz(-3.0078648) q[0];
sx q[0];
rz(2.4753841) q[0];
rz(-1.1353525) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(1.1599783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839417) q[0];
sx q[0];
rz(-0.91383024) q[0];
sx q[0];
rz(-2.3131671) q[0];
rz(-1.3526037) q[2];
sx q[2];
rz(-0.82880965) q[2];
sx q[2];
rz(0.11889501) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6659703) q[1];
sx q[1];
rz(-0.95487528) q[1];
sx q[1];
rz(-0.0093950942) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0464675) q[3];
sx q[3];
rz(-2.5216649) q[3];
sx q[3];
rz(1.4038831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3975415) q[2];
sx q[2];
rz(-0.56637374) q[2];
sx q[2];
rz(0.19691021) q[2];
rz(-0.87013733) q[3];
sx q[3];
rz(-2.0700442) q[3];
sx q[3];
rz(-1.2675233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5043735) q[0];
sx q[0];
rz(-2.0688031) q[0];
sx q[0];
rz(-3.1043501) q[0];
rz(1.0875018) q[1];
sx q[1];
rz(-1.0481513) q[1];
sx q[1];
rz(-2.9748999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.113362) q[0];
sx q[0];
rz(-0.89269243) q[0];
sx q[0];
rz(1.9166758) q[0];
rz(2.6092763) q[2];
sx q[2];
rz(-2.5925143) q[2];
sx q[2];
rz(-2.2808569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5350668) q[1];
sx q[1];
rz(-1.4826315) q[1];
sx q[1];
rz(2.8463581) q[1];
x q[2];
rz(1.6400723) q[3];
sx q[3];
rz(-1.426094) q[3];
sx q[3];
rz(-2.7370344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4166261) q[2];
sx q[2];
rz(-0.61369696) q[2];
sx q[2];
rz(-1.7566768) q[2];
rz(-2.7409605) q[3];
sx q[3];
rz(-1.8568361) q[3];
sx q[3];
rz(1.0111151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8818125) q[0];
sx q[0];
rz(-1.1008982) q[0];
sx q[0];
rz(-0.63006054) q[0];
rz(0.13314816) q[1];
sx q[1];
rz(-1.9825736) q[1];
sx q[1];
rz(2.0864743) q[1];
rz(-3.1027267) q[2];
sx q[2];
rz(-1.7615223) q[2];
sx q[2];
rz(2.7675235) q[2];
rz(0.12689982) q[3];
sx q[3];
rz(-0.77346934) q[3];
sx q[3];
rz(1.0897549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
