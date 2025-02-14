OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8920933) q[0];
sx q[0];
rz(-1.1846932) q[0];
sx q[0];
rz(2.0062334) q[0];
rz(-0.42179498) q[1];
sx q[1];
rz(-2.3611958) q[1];
sx q[1];
rz(2.5193522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6506158) q[0];
sx q[0];
rz(-2.9497006) q[0];
sx q[0];
rz(-2.2723115) q[0];
rz(-pi) q[1];
rz(0.038921629) q[2];
sx q[2];
rz(-2.0177671) q[2];
sx q[2];
rz(-2.1291422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65917817) q[1];
sx q[1];
rz(-0.59450048) q[1];
sx q[1];
rz(-3.1207419) q[1];
rz(-0.99986003) q[3];
sx q[3];
rz(-1.3821162) q[3];
sx q[3];
rz(-0.41485559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3368095) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(-0.08610227) q[2];
rz(0.51186776) q[3];
sx q[3];
rz(-1.0443001) q[3];
sx q[3];
rz(-0.69387236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43857273) q[0];
sx q[0];
rz(-2.1874671) q[0];
sx q[0];
rz(2.2254206) q[0];
rz(0.29197261) q[1];
sx q[1];
rz(-0.59140721) q[1];
sx q[1];
rz(2.3672262) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2949897) q[0];
sx q[0];
rz(-1.2520391) q[0];
sx q[0];
rz(1.9413207) q[0];
x q[1];
rz(-1.5541114) q[2];
sx q[2];
rz(-0.34539255) q[2];
sx q[2];
rz(-2.2881803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91097478) q[1];
sx q[1];
rz(-2.9819192) q[1];
sx q[1];
rz(0.69024872) q[1];
rz(-pi) q[2];
rz(-0.99454576) q[3];
sx q[3];
rz(-0.835383) q[3];
sx q[3];
rz(-1.6540838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3463717) q[2];
sx q[2];
rz(-0.56196153) q[2];
sx q[2];
rz(-0.58491659) q[2];
rz(-1.6631205) q[3];
sx q[3];
rz(-1.9641179) q[3];
sx q[3];
rz(1.6519206) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65603489) q[0];
sx q[0];
rz(-2.1182883) q[0];
sx q[0];
rz(-2.5265332) q[0];
rz(-0.3282322) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(-2.097791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8449677) q[0];
sx q[0];
rz(-0.052694885) q[0];
sx q[0];
rz(3.0449163) q[0];
rz(-pi) q[1];
rz(1.5659901) q[2];
sx q[2];
rz(-1.9991181) q[2];
sx q[2];
rz(2.7132963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17576007) q[1];
sx q[1];
rz(-1.1783667) q[1];
sx q[1];
rz(-0.85939851) q[1];
rz(-0.7316771) q[3];
sx q[3];
rz(-0.51245514) q[3];
sx q[3];
rz(0.75430124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.388776) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(0.30889312) q[2];
rz(-2.654352) q[3];
sx q[3];
rz(-1.7083218) q[3];
sx q[3];
rz(0.88700956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76172817) q[0];
sx q[0];
rz(-2.3966615) q[0];
sx q[0];
rz(0.47759011) q[0];
rz(1.2848805) q[1];
sx q[1];
rz(-1.711859) q[1];
sx q[1];
rz(0.46599785) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6687688) q[0];
sx q[0];
rz(-1.6472938) q[0];
sx q[0];
rz(2.5868505) q[0];
rz(-2.3498769) q[2];
sx q[2];
rz(-0.60388619) q[2];
sx q[2];
rz(1.1271267) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5090967) q[1];
sx q[1];
rz(-1.6057311) q[1];
sx q[1];
rz(-0.81991244) q[1];
x q[2];
rz(-2.1804564) q[3];
sx q[3];
rz(-2.1603322) q[3];
sx q[3];
rz(1.846755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8183225) q[2];
sx q[2];
rz(-2.6939836) q[2];
sx q[2];
rz(-1.4480048) q[2];
rz(1.4130392) q[3];
sx q[3];
rz(-1.6084684) q[3];
sx q[3];
rz(1.1477227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2959761) q[0];
sx q[0];
rz(-0.71346658) q[0];
sx q[0];
rz(-0.071685858) q[0];
rz(1.4777615) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(0.62612265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6090995) q[0];
sx q[0];
rz(-0.83864826) q[0];
sx q[0];
rz(0.12153836) q[0];
x q[1];
rz(3.0813498) q[2];
sx q[2];
rz(-0.92298302) q[2];
sx q[2];
rz(2.197165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4461527) q[1];
sx q[1];
rz(-1.3976516) q[1];
sx q[1];
rz(2.3320564) q[1];
rz(2.7887129) q[3];
sx q[3];
rz(-1.7418752) q[3];
sx q[3];
rz(-1.4346348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6757297) q[2];
sx q[2];
rz(-2.6079123) q[2];
sx q[2];
rz(1.6443058) q[2];
rz(0.40694445) q[3];
sx q[3];
rz(-0.99138433) q[3];
sx q[3];
rz(2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1740455) q[0];
sx q[0];
rz(-1.8438735) q[0];
sx q[0];
rz(-2.8795854) q[0];
rz(1.2213446) q[1];
sx q[1];
rz(-1.6294934) q[1];
sx q[1];
rz(1.2006203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1977494) q[0];
sx q[0];
rz(-0.89559865) q[0];
sx q[0];
rz(0.64267107) q[0];
rz(-pi) q[1];
rz(-2.1706743) q[2];
sx q[2];
rz(-2.3339503) q[2];
sx q[2];
rz(1.3563479) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3221601) q[1];
sx q[1];
rz(-1.3272078) q[1];
sx q[1];
rz(1.0526245) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1061078) q[3];
sx q[3];
rz(-1.662613) q[3];
sx q[3];
rz(-2.6694856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0210586) q[2];
sx q[2];
rz(-2.6594682) q[2];
sx q[2];
rz(-2.6684707) q[2];
rz(2.1038697) q[3];
sx q[3];
rz(-0.40268746) q[3];
sx q[3];
rz(1.13824) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6556019) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(-2.3714491) q[0];
rz(-2.7262402) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(-1.279668) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3949385) q[0];
sx q[0];
rz(-1.5459367) q[0];
sx q[0];
rz(2.3034873) q[0];
rz(-pi) q[1];
rz(-1.2998313) q[2];
sx q[2];
rz(-0.44579664) q[2];
sx q[2];
rz(-0.4933756) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34977594) q[1];
sx q[1];
rz(-0.27128759) q[1];
sx q[1];
rz(-0.94601814) q[1];
rz(-pi) q[2];
rz(-1.5611675) q[3];
sx q[3];
rz(-1.0917679) q[3];
sx q[3];
rz(-1.1232532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27560774) q[2];
sx q[2];
rz(-1.7476247) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(-1.154493) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.481021) q[0];
sx q[0];
rz(-0.72089973) q[0];
sx q[0];
rz(0.50203669) q[0];
rz(-1.0188811) q[1];
sx q[1];
rz(-1.4733682) q[1];
sx q[1];
rz(-0.92207164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88387892) q[0];
sx q[0];
rz(-1.3259616) q[0];
sx q[0];
rz(-1.4725982) q[0];
x q[1];
rz(1.3241295) q[2];
sx q[2];
rz(-2.1804159) q[2];
sx q[2];
rz(2.7423046) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.079262924) q[1];
sx q[1];
rz(-3.0714066) q[1];
sx q[1];
rz(0.90026469) q[1];
rz(1.5899379) q[3];
sx q[3];
rz(-0.12775207) q[3];
sx q[3];
rz(-2.9763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8036538) q[2];
sx q[2];
rz(-1.2793845) q[2];
sx q[2];
rz(0.12019084) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-2.7836697) q[3];
sx q[3];
rz(3.1357989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040076479) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(1.6325604) q[0];
rz(0.82935968) q[1];
sx q[1];
rz(-2.6396535) q[1];
sx q[1];
rz(2.0109743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42989995) q[0];
sx q[0];
rz(-1.2303174) q[0];
sx q[0];
rz(-2.8116624) q[0];
rz(2.0810602) q[2];
sx q[2];
rz(-2.2203373) q[2];
sx q[2];
rz(2.5819957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.045550195) q[1];
sx q[1];
rz(-0.99771959) q[1];
sx q[1];
rz(2.6144774) q[1];
x q[2];
rz(-2.8895414) q[3];
sx q[3];
rz(-1.0672616) q[3];
sx q[3];
rz(-1.9084712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.074038) q[2];
sx q[2];
rz(-2.270547) q[2];
sx q[2];
rz(-1.7769495) q[2];
rz(3.126295) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(-1.8196222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.44052723) q[0];
sx q[0];
rz(-0.65182132) q[0];
sx q[0];
rz(0.1524674) q[0];
rz(1.3854148) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(-0.32858953) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5009896) q[0];
sx q[0];
rz(-3.1251934) q[0];
sx q[0];
rz(-2.5522405) q[0];
x q[1];
rz(-1.3297021) q[2];
sx q[2];
rz(-0.82359353) q[2];
sx q[2];
rz(1.6833978) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8161502) q[1];
sx q[1];
rz(-2.901536) q[1];
sx q[1];
rz(1.6961967) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64528193) q[3];
sx q[3];
rz(-1.2835763) q[3];
sx q[3];
rz(-0.20791277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3370328) q[2];
sx q[2];
rz(-2.0064662) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(-2.1870901) q[3];
sx q[3];
rz(-1.4408305) q[3];
sx q[3];
rz(-2.8452828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8148707) q[0];
sx q[0];
rz(-2.1537415) q[0];
sx q[0];
rz(1.4065773) q[0];
rz(-1.8234491) q[1];
sx q[1];
rz(-1.5144009) q[1];
sx q[1];
rz(-2.3419793) q[1];
rz(0.90262765) q[2];
sx q[2];
rz(-3.0716574) q[2];
sx q[2];
rz(2.8989094) q[2];
rz(1.3621422) q[3];
sx q[3];
rz(-0.91051523) q[3];
sx q[3];
rz(-0.12403535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
