OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.2494994) q[0];
sx q[0];
rz(-1.9568994) q[0];
sx q[0];
rz(1.1353593) q[0];
rz(-0.42179498) q[1];
sx q[1];
rz(3.9219895) q[1];
sx q[1];
rz(11.94413) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61257768) q[0];
sx q[0];
rz(-1.6941931) q[0];
sx q[0];
rz(-1.4234659) q[0];
rz(1.6517992) q[2];
sx q[2];
rz(-0.44854823) q[2];
sx q[2];
rz(-2.0392921) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92889228) q[1];
sx q[1];
rz(-1.5824741) q[1];
sx q[1];
rz(2.547193) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22316605) q[3];
sx q[3];
rz(-1.0112178) q[3];
sx q[3];
rz(2.1055438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3368095) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(-0.08610227) q[2];
rz(2.6297249) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(-0.69387236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43857273) q[0];
sx q[0];
rz(-0.95412552) q[0];
sx q[0];
rz(-0.91617209) q[0];
rz(0.29197261) q[1];
sx q[1];
rz(-0.59140721) q[1];
sx q[1];
rz(-0.77436647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2949897) q[0];
sx q[0];
rz(-1.2520391) q[0];
sx q[0];
rz(-1.2002719) q[0];
rz(-pi) q[1];
rz(-3.1355895) q[2];
sx q[2];
rz(-1.9161388) q[2];
sx q[2];
rz(0.83568043) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7978002) q[1];
sx q[1];
rz(-1.6722073) q[1];
sx q[1];
rz(3.0180458) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1470469) q[3];
sx q[3];
rz(-0.835383) q[3];
sx q[3];
rz(-1.6540838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3463717) q[2];
sx q[2];
rz(-0.56196153) q[2];
sx q[2];
rz(-2.5566761) q[2];
rz(-1.6631205) q[3];
sx q[3];
rz(-1.1774747) q[3];
sx q[3];
rz(-1.6519206) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65603489) q[0];
sx q[0];
rz(-1.0233044) q[0];
sx q[0];
rz(-0.61505944) q[0];
rz(-2.8133605) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(-1.0438017) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9639643) q[0];
sx q[0];
rz(-1.5657122) q[0];
sx q[0];
rz(3.0891434) q[0];
x q[1];
rz(1.5756025) q[2];
sx q[2];
rz(-1.9991181) q[2];
sx q[2];
rz(-2.7132963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9658326) q[1];
sx q[1];
rz(-1.1783667) q[1];
sx q[1];
rz(0.85939851) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2112593) q[3];
sx q[3];
rz(-1.1973527) q[3];
sx q[3];
rz(-1.5869753) q[3];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798645) q[0];
sx q[0];
rz(-2.3966615) q[0];
sx q[0];
rz(-0.47759011) q[0];
rz(-1.2848805) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(-2.6755948) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1168524) q[0];
sx q[0];
rz(-0.55944397) q[0];
sx q[0];
rz(-0.14450216) q[0];
x q[1];
rz(-2.6902507) q[2];
sx q[2];
rz(-1.9867431) q[2];
sx q[2];
rz(-2.0030264) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6324959) q[1];
sx q[1];
rz(-1.6057311) q[1];
sx q[1];
rz(2.3216802) q[1];
rz(-pi) q[2];
rz(-2.1804564) q[3];
sx q[3];
rz(-0.98126047) q[3];
sx q[3];
rz(1.2948376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3232702) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(1.6935879) q[2];
rz(1.7285534) q[3];
sx q[3];
rz(-1.6084684) q[3];
sx q[3];
rz(-1.1477227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84561658) q[0];
sx q[0];
rz(-0.71346658) q[0];
sx q[0];
rz(-0.071685858) q[0];
rz(1.4777615) q[1];
sx q[1];
rz(-1.7774589) q[1];
sx q[1];
rz(-0.62612265) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043163096) q[0];
sx q[0];
rz(-1.6610896) q[0];
sx q[0];
rz(-2.3066269) q[0];
rz(-3.0813498) q[2];
sx q[2];
rz(-2.2186096) q[2];
sx q[2];
rz(-0.94442764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71264919) q[1];
sx q[1];
rz(-0.82368851) q[1];
sx q[1];
rz(-0.23703833) q[1];
rz(-1.7528698) q[3];
sx q[3];
rz(-1.9183049) q[3];
sx q[3];
rz(-0.073542882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6757297) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(1.6443058) q[2];
rz(-0.40694445) q[3];
sx q[3];
rz(-2.1502083) q[3];
sx q[3];
rz(2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96754718) q[0];
sx q[0];
rz(-1.8438735) q[0];
sx q[0];
rz(-2.8795854) q[0];
rz(-1.920248) q[1];
sx q[1];
rz(-1.5120993) q[1];
sx q[1];
rz(-1.2006203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9438433) q[0];
sx q[0];
rz(-0.89559865) q[0];
sx q[0];
rz(2.4989216) q[0];
rz(-pi) q[1];
rz(2.2827705) q[2];
sx q[2];
rz(-1.1505652) q[2];
sx q[2];
rz(2.4855297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38528864) q[1];
sx q[1];
rz(-2.0722162) q[1];
sx q[1];
rz(-2.8629567) q[1];
x q[2];
rz(-1.6626705) q[3];
sx q[3];
rz(-1.535461) q[3];
sx q[3];
rz(1.0954344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0210586) q[2];
sx q[2];
rz(-2.6594682) q[2];
sx q[2];
rz(2.6684707) q[2];
rz(1.037723) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(1.13824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4859908) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(0.77014357) q[0];
rz(-0.4153525) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(1.279668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3949385) q[0];
sx q[0];
rz(-1.5459367) q[0];
sx q[0];
rz(-2.3034873) q[0];
rz(-1.2998313) q[2];
sx q[2];
rz(-0.44579664) q[2];
sx q[2];
rz(2.6482171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61381147) q[1];
sx q[1];
rz(-1.7281869) q[1];
sx q[1];
rz(1.792683) q[1];
rz(-1.5804251) q[3];
sx q[3];
rz(-2.0498247) q[3];
sx q[3];
rz(2.0183394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8659849) q[2];
sx q[2];
rz(-1.7476247) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(1.9870997) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6605717) q[0];
sx q[0];
rz(-0.72089973) q[0];
sx q[0];
rz(0.50203669) q[0];
rz(-1.0188811) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(0.92207164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4785504) q[0];
sx q[0];
rz(-1.6660569) q[0];
sx q[0];
rz(-0.24597286) q[0];
rz(-pi) q[1];
rz(0.33635535) q[2];
sx q[2];
rz(-0.65170641) q[2];
sx q[2];
rz(-2.327988) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0623297) q[1];
sx q[1];
rz(-0.070186071) q[1];
sx q[1];
rz(0.90026469) q[1];
x q[2];
rz(-1.4430674) q[3];
sx q[3];
rz(-1.5683577) q[3];
sx q[3];
rz(1.7170441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33793882) q[2];
sx q[2];
rz(-1.8622082) q[2];
sx q[2];
rz(-0.12019084) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-2.7836697) q[3];
sx q[3];
rz(3.1357989) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1015162) q[0];
sx q[0];
rz(-0.67779556) q[0];
sx q[0];
rz(-1.5090322) q[0];
rz(0.82935968) q[1];
sx q[1];
rz(-0.50193915) q[1];
sx q[1];
rz(-2.0109743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7735159) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(0.83004029) q[0];
x q[1];
rz(2.5700966) q[2];
sx q[2];
rz(-0.80249107) q[2];
sx q[2];
rz(2.9545138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.045550195) q[1];
sx q[1];
rz(-2.1438731) q[1];
sx q[1];
rz(2.6144774) q[1];
rz(-2.8895414) q[3];
sx q[3];
rz(-2.0743311) q[3];
sx q[3];
rz(1.9084712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.074038) q[2];
sx q[2];
rz(-2.270547) q[2];
sx q[2];
rz(-1.3646431) q[2];
rz(-3.126295) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(-1.3219705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.44052723) q[0];
sx q[0];
rz(-0.65182132) q[0];
sx q[0];
rz(-2.9891253) q[0];
rz(1.3854148) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(-0.32858953) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2300174) q[0];
sx q[0];
rz(-1.5571638) q[0];
sx q[0];
rz(1.5799119) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25224884) q[2];
sx q[2];
rz(-2.3636732) q[2];
sx q[2];
rz(-1.8053448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32544245) q[1];
sx q[1];
rz(-0.24005664) q[1];
sx q[1];
rz(-1.6961967) q[1];
rz(-2.4963107) q[3];
sx q[3];
rz(-1.8580164) q[3];
sx q[3];
rz(2.9336799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3370328) q[2];
sx q[2];
rz(-1.1351265) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(0.95450258) q[3];
sx q[3];
rz(-1.7007622) q[3];
sx q[3];
rz(-0.29630989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8148707) q[0];
sx q[0];
rz(-0.9878511) q[0];
sx q[0];
rz(-1.7350154) q[0];
rz(-1.3181435) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(0.043371928) q[2];
sx q[2];
rz(-1.6256754) q[2];
sx q[2];
rz(-0.91204249) q[2];
rz(0.26067692) q[3];
sx q[3];
rz(-0.68772975) q[3];
sx q[3];
rz(-2.9332193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
