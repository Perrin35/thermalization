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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6506158) q[0];
sx q[0];
rz(-2.9497006) q[0];
sx q[0];
rz(-2.2723115) q[0];
x q[1];
rz(3.102671) q[2];
sx q[2];
rz(-1.1238255) q[2];
sx q[2];
rz(1.0124504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4824145) q[1];
sx q[1];
rz(-2.5470922) q[1];
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
rz(-1.8047831) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(-0.08610227) q[2];
rz(0.51186776) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(-2.4477203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7030199) q[0];
sx q[0];
rz(-2.1874671) q[0];
sx q[0];
rz(2.2254206) q[0];
rz(0.29197261) q[1];
sx q[1];
rz(-0.59140721) q[1];
sx q[1];
rz(-0.77436647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2949897) q[0];
sx q[0];
rz(-1.2520391) q[0];
sx q[0];
rz(-1.2002719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9161445) q[2];
sx q[2];
rz(-1.5651476) q[2];
sx q[2];
rz(2.408509) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7978002) q[1];
sx q[1];
rz(-1.4693854) q[1];
sx q[1];
rz(-0.12354688) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5995042) q[3];
sx q[3];
rz(-0.89975587) q[3];
sx q[3];
rz(-2.4235291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3463717) q[2];
sx q[2];
rz(-2.5796311) q[2];
sx q[2];
rz(-0.58491659) q[2];
rz(-1.6631205) q[3];
sx q[3];
rz(-1.1774747) q[3];
sx q[3];
rz(-1.6519206) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65603489) q[0];
sx q[0];
rz(-1.0233044) q[0];
sx q[0];
rz(0.61505944) q[0];
rz(0.3282322) q[1];
sx q[1];
rz(-2.131772) q[1];
sx q[1];
rz(1.0438017) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1776284) q[0];
sx q[0];
rz(-1.5657122) q[0];
sx q[0];
rz(-0.052449277) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7132665) q[2];
sx q[2];
rz(-1.5664243) q[2];
sx q[2];
rz(1.1405038) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97745132) q[1];
sx q[1];
rz(-0.79557997) q[1];
sx q[1];
rz(1.005791) q[1];
x q[2];
rz(2.7451594) q[3];
sx q[3];
rz(-1.2370438) q[3];
sx q[3];
rz(-0.15245701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7528167) q[2];
sx q[2];
rz(-0.39377585) q[2];
sx q[2];
rz(0.30889312) q[2];
rz(0.4872407) q[3];
sx q[3];
rz(-1.4332708) q[3];
sx q[3];
rz(-0.88700956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76172817) q[0];
sx q[0];
rz(-0.74493113) q[0];
sx q[0];
rz(0.47759011) q[0];
rz(-1.8567122) q[1];
sx q[1];
rz(-1.711859) q[1];
sx q[1];
rz(0.46599785) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1168524) q[0];
sx q[0];
rz(-2.5821487) q[0];
sx q[0];
rz(-0.14450216) q[0];
rz(-pi) q[1];
rz(1.1144756) q[2];
sx q[2];
rz(-1.9812036) q[2];
sx q[2];
rz(2.9027651) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0424846) q[1];
sx q[1];
rz(-2.3900552) q[1];
sx q[1];
rz(1.5196176) q[1];
x q[2];
rz(-2.1804564) q[3];
sx q[3];
rz(-0.98126047) q[3];
sx q[3];
rz(-1.846755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-1.9938699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2959761) q[0];
sx q[0];
rz(-2.4281261) q[0];
sx q[0];
rz(3.0699068) q[0];
rz(-1.4777615) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(-0.62612265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043163096) q[0];
sx q[0];
rz(-1.6610896) q[0];
sx q[0];
rz(0.8349658) q[0];
rz(-1.4914054) q[2];
sx q[2];
rz(-2.491386) q[2];
sx q[2];
rz(-2.2967867) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0542293) q[1];
sx q[1];
rz(-2.3647671) q[1];
sx q[1];
rz(-1.3224949) q[1];
rz(-pi) q[2];
rz(1.7528698) q[3];
sx q[3];
rz(-1.9183049) q[3];
sx q[3];
rz(-3.0680498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46586299) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(-1.4972868) q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-0.99138433) q[3];
sx q[3];
rz(2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1740455) q[0];
sx q[0];
rz(-1.2977192) q[0];
sx q[0];
rz(-2.8795854) q[0];
rz(1.2213446) q[1];
sx q[1];
rz(-1.6294934) q[1];
sx q[1];
rz(-1.9409723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8190099) q[0];
sx q[0];
rz(-2.2457425) q[0];
sx q[0];
rz(2.2133166) q[0];
rz(-pi) q[1];
rz(-0.97091834) q[2];
sx q[2];
rz(-0.80764233) q[2];
sx q[2];
rz(-1.7852448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.756304) q[1];
sx q[1];
rz(-2.0722162) q[1];
sx q[1];
rz(2.8629567) q[1];
rz(-1.4789222) q[3];
sx q[3];
rz(-1.535461) q[3];
sx q[3];
rz(-1.0954344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0210586) q[2];
sx q[2];
rz(-2.6594682) q[2];
sx q[2];
rz(2.6684707) q[2];
rz(2.1038697) q[3];
sx q[3];
rz(-0.40268746) q[3];
sx q[3];
rz(1.13824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4859908) q[0];
sx q[0];
rz(-2.7779873) q[0];
sx q[0];
rz(-2.3714491) q[0];
rz(-0.4153525) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(1.279668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19821985) q[0];
sx q[0];
rz(-0.83838338) q[0];
sx q[0];
rz(-3.1081568) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0143731) q[2];
sx q[2];
rz(-1.1423649) q[2];
sx q[2];
rz(0.7920533) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99233124) q[1];
sx q[1];
rz(-1.3516973) q[1];
sx q[1];
rz(-2.9803139) q[1];
x q[2];
rz(-1.5804251) q[3];
sx q[3];
rz(-2.0498247) q[3];
sx q[3];
rz(2.0183394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27560774) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(-2.7371791) q[2];
rz(-1.154493) q[3];
sx q[3];
rz(-1.5279852) q[3];
sx q[3];
rz(-1.9563458) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6605717) q[0];
sx q[0];
rz(-0.72089973) q[0];
sx q[0];
rz(2.639556) q[0];
rz(-1.0188811) q[1];
sx q[1];
rz(-1.4733682) q[1];
sx q[1];
rz(2.219521) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2577137) q[0];
sx q[0];
rz(-1.3259616) q[0];
sx q[0];
rz(-1.6689945) q[0];
rz(2.8052373) q[2];
sx q[2];
rz(-2.4898862) q[2];
sx q[2];
rz(-2.327988) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3193911) q[1];
sx q[1];
rz(-1.5272045) q[1];
sx q[1];
rz(1.6258214) q[1];
x q[2];
rz(1.5899379) q[3];
sx q[3];
rz(-0.12775207) q[3];
sx q[3];
rz(-2.9763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8036538) q[2];
sx q[2];
rz(-1.8622082) q[2];
sx q[2];
rz(3.0214018) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-0.35792297) q[3];
sx q[3];
rz(0.005793747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1015162) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(-1.5090322) q[0];
rz(-2.312233) q[1];
sx q[1];
rz(-2.6396535) q[1];
sx q[1];
rz(2.0109743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7735159) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(-0.83004029) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0605325) q[2];
sx q[2];
rz(-2.2203373) q[2];
sx q[2];
rz(-0.559597) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.045550195) q[1];
sx q[1];
rz(-2.1438731) q[1];
sx q[1];
rz(-0.5271153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0535766) q[3];
sx q[3];
rz(-1.791009) q[3];
sx q[3];
rz(2.6802879) q[3];
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
rz(-0.015297628) q[3];
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
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7010654) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(0.1524674) q[0];
rz(1.3854148) q[1];
sx q[1];
rz(-2.0962174) q[1];
sx q[1];
rz(-2.8130031) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5009896) q[0];
sx q[0];
rz(-3.1251934) q[0];
sx q[0];
rz(-0.58935215) q[0];
rz(-0.25224884) q[2];
sx q[2];
rz(-2.3636732) q[2];
sx q[2];
rz(-1.3362479) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45450452) q[1];
sx q[1];
rz(-1.3326613) q[1];
sx q[1];
rz(3.1109875) q[1];
x q[2];
rz(-0.64528193) q[3];
sx q[3];
rz(-1.2835763) q[3];
sx q[3];
rz(2.9336799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80455989) q[2];
sx q[2];
rz(-2.0064662) q[2];
sx q[2];
rz(-1.0767153) q[2];
rz(-2.1870901) q[3];
sx q[3];
rz(-1.7007622) q[3];
sx q[3];
rz(2.8452828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3267219) q[0];
sx q[0];
rz(-0.9878511) q[0];
sx q[0];
rz(-1.7350154) q[0];
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
rz(1.7794505) q[3];
sx q[3];
rz(-2.2310774) q[3];
sx q[3];
rz(3.0175573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
