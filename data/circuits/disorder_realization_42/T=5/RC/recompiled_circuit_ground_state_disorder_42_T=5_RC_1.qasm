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
rz(-1.1353593) q[0];
rz(2.7197977) q[1];
sx q[1];
rz(-0.78039688) q[1];
sx q[1];
rz(0.62224046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6506158) q[0];
sx q[0];
rz(-2.9497006) q[0];
sx q[0];
rz(2.2723115) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.102671) q[2];
sx q[2];
rz(-1.1238255) q[2];
sx q[2];
rz(-1.0124504) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4824145) q[1];
sx q[1];
rz(-0.59450048) q[1];
sx q[1];
rz(3.1207419) q[1];
rz(-pi) q[2];
rz(1.9104426) q[3];
sx q[3];
rz(-2.5435919) q[3];
sx q[3];
rz(-1.4400618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8047831) q[2];
sx q[2];
rz(-1.6853036) q[2];
sx q[2];
rz(0.08610227) q[2];
rz(2.6297249) q[3];
sx q[3];
rz(-1.0443001) q[3];
sx q[3];
rz(-2.4477203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43857273) q[0];
sx q[0];
rz(-2.1874671) q[0];
sx q[0];
rz(-2.2254206) q[0];
rz(0.29197261) q[1];
sx q[1];
rz(-2.5501854) q[1];
sx q[1];
rz(-2.3672262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5385434) q[0];
sx q[0];
rz(-1.2197681) q[0];
sx q[0];
rz(2.8013264) q[0];
x q[1];
rz(1.2254481) q[2];
sx q[2];
rz(-1.5651476) q[2];
sx q[2];
rz(-2.408509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2306179) q[1];
sx q[1];
rz(-2.9819192) q[1];
sx q[1];
rz(-2.4513439) q[1];
rz(0.82335681) q[3];
sx q[3];
rz(-1.1548448) q[3];
sx q[3];
rz(0.49440629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3463717) q[2];
sx q[2];
rz(-2.5796311) q[2];
sx q[2];
rz(0.58491659) q[2];
rz(1.6631205) q[3];
sx q[3];
rz(-1.9641179) q[3];
sx q[3];
rz(1.4896721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(0.65603489) q[0];
sx q[0];
rz(-2.1182883) q[0];
sx q[0];
rz(0.61505944) q[0];
rz(0.3282322) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(-1.0438017) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39343483) q[0];
sx q[0];
rz(-1.5183477) q[0];
sx q[0];
rz(1.5657052) q[0];
rz(1.5756025) q[2];
sx q[2];
rz(-1.1424746) q[2];
sx q[2];
rz(-0.42829633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4281299) q[1];
sx q[1];
rz(-0.9232115) q[1];
sx q[1];
rz(0.50010276) q[1];
x q[2];
rz(2.4099156) q[3];
sx q[3];
rz(-2.6291375) q[3];
sx q[3];
rz(-0.75430124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.388776) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(-2.8326995) q[2];
rz(-0.4872407) q[3];
sx q[3];
rz(-1.7083218) q[3];
sx q[3];
rz(-0.88700956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76172817) q[0];
sx q[0];
rz(-0.74493113) q[0];
sx q[0];
rz(2.6640025) q[0];
rz(1.2848805) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(-0.46599785) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024740273) q[0];
sx q[0];
rz(-0.55944397) q[0];
sx q[0];
rz(-2.9970905) q[0];
x q[1];
rz(-2.3498769) q[2];
sx q[2];
rz(-0.60388619) q[2];
sx q[2];
rz(-2.0144659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5090967) q[1];
sx q[1];
rz(-1.5358616) q[1];
sx q[1];
rz(0.81991244) q[1];
rz(0.9611363) q[3];
sx q[3];
rz(-0.98126047) q[3];
sx q[3];
rz(-1.846755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3232702) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(-1.6935879) q[2];
rz(-1.4130392) q[3];
sx q[3];
rz(-1.5331242) q[3];
sx q[3];
rz(1.1477227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2959761) q[0];
sx q[0];
rz(-0.71346658) q[0];
sx q[0];
rz(-3.0699068) q[0];
rz(1.6638311) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(-0.62612265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4283765) q[0];
sx q[0];
rz(-2.4012743) q[0];
sx q[0];
rz(-1.4367144) q[0];
rz(-pi) q[1];
rz(3.0813498) q[2];
sx q[2];
rz(-2.2186096) q[2];
sx q[2];
rz(0.94442764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0542293) q[1];
sx q[1];
rz(-0.77682553) q[1];
sx q[1];
rz(1.3224949) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3887229) q[3];
sx q[3];
rz(-1.2232878) q[3];
sx q[3];
rz(0.073542882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46586299) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(-1.6443058) q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-2.1502083) q[3];
sx q[3];
rz(-2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96754718) q[0];
sx q[0];
rz(-1.2977192) q[0];
sx q[0];
rz(-0.26200727) q[0];
rz(1.920248) q[1];
sx q[1];
rz(-1.6294934) q[1];
sx q[1];
rz(-1.2006203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81075087) q[0];
sx q[0];
rz(-2.057632) q[0];
sx q[0];
rz(2.3563516) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97091834) q[2];
sx q[2];
rz(-2.3339503) q[2];
sx q[2];
rz(1.7852448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38528864) q[1];
sx q[1];
rz(-2.0722162) q[1];
sx q[1];
rz(0.27863592) q[1];
rz(-pi) q[2];
rz(-1.2030199) q[3];
sx q[3];
rz(-3.0431755) q[3];
sx q[3];
rz(-0.84151387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0210586) q[2];
sx q[2];
rz(-0.48212442) q[2];
sx q[2];
rz(-2.6684707) q[2];
rz(1.037723) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(-2.0033526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6556019) q[0];
sx q[0];
rz(-2.7779873) q[0];
sx q[0];
rz(2.3714491) q[0];
rz(-2.7262402) q[1];
sx q[1];
rz(-1.5130006) q[1];
sx q[1];
rz(1.279668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19821985) q[0];
sx q[0];
rz(-0.83838338) q[0];
sx q[0];
rz(3.1081568) q[0];
rz(3.0143731) q[2];
sx q[2];
rz(-1.9992277) q[2];
sx q[2];
rz(0.7920533) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7918167) q[1];
sx q[1];
rz(-0.27128759) q[1];
sx q[1];
rz(-0.94601814) q[1];
x q[2];
rz(3.123056) q[3];
sx q[3];
rz(-0.47911767) q[3];
sx q[3];
rz(-1.1441413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27560774) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(-2.7371791) q[2];
rz(-1.154493) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(-1.1852468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6605717) q[0];
sx q[0];
rz(-0.72089973) q[0];
sx q[0];
rz(2.639556) q[0];
rz(1.0188811) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(2.219521) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88387892) q[0];
sx q[0];
rz(-1.8156311) q[0];
sx q[0];
rz(-1.4725982) q[0];
rz(-pi) q[1];
rz(0.62412213) q[2];
sx q[2];
rz(-1.7723473) q[2];
sx q[2];
rz(1.0283284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0623297) q[1];
sx q[1];
rz(-0.070186071) q[1];
sx q[1];
rz(-2.241328) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6985252) q[3];
sx q[3];
rz(-1.5683577) q[3];
sx q[3];
rz(1.7170441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33793882) q[2];
sx q[2];
rz(-1.2793845) q[2];
sx q[2];
rz(3.0214018) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-2.7836697) q[3];
sx q[3];
rz(3.1357989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1015162) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(1.5090322) q[0];
rz(0.82935968) q[1];
sx q[1];
rz(-2.6396535) q[1];
sx q[1];
rz(-1.1306184) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7735159) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(0.83004029) q[0];
rz(1.0605325) q[2];
sx q[2];
rz(-0.9212554) q[2];
sx q[2];
rz(-0.559597) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77525951) q[1];
sx q[1];
rz(-0.75804064) q[1];
sx q[1];
rz(-0.90866462) q[1];
rz(1.1457024) q[3];
sx q[3];
rz(-2.5833874) q[3];
sx q[3];
rz(2.3987215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.074038) q[2];
sx q[2];
rz(-0.87104565) q[2];
sx q[2];
rz(-1.3646431) q[2];
rz(-3.126295) q[3];
sx q[3];
rz(-0.8316032) q[3];
sx q[3];
rz(-1.8196222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052723) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(0.1524674) q[0];
rz(1.3854148) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(-0.32858953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91157527) q[0];
sx q[0];
rz(-1.5844288) q[0];
sx q[0];
rz(1.5616807) q[0];
rz(1.3297021) q[2];
sx q[2];
rz(-2.3179991) q[2];
sx q[2];
rz(-1.4581949) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45450452) q[1];
sx q[1];
rz(-1.8089314) q[1];
sx q[1];
rz(0.030605153) q[1];
rz(-0.45654065) q[3];
sx q[3];
rz(-2.4437063) q[3];
sx q[3];
rz(2.1386351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3370328) q[2];
sx q[2];
rz(-2.0064662) q[2];
sx q[2];
rz(-2.0648773) q[2];
rz(-0.95450258) q[3];
sx q[3];
rz(-1.4408305) q[3];
sx q[3];
rz(-0.29630989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3267219) q[0];
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
rz(-2.4706609) q[3];
sx q[3];
rz(-1.4064515) q[3];
sx q[3];
rz(1.575904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
