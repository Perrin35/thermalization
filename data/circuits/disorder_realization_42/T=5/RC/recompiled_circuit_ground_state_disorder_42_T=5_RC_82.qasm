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
rz(3.9219895) q[1];
sx q[1];
rz(11.94413) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.529015) q[0];
sx q[0];
rz(-1.4473995) q[0];
sx q[0];
rz(1.7181267) q[0];
rz(-pi) q[1];
x q[1];
rz(0.038921629) q[2];
sx q[2];
rz(-1.1238255) q[2];
sx q[2];
rz(-1.0124504) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4824145) q[1];
sx q[1];
rz(-0.59450048) q[1];
sx q[1];
rz(0.020850734) q[1];
rz(-pi) q[2];
rz(1.9104426) q[3];
sx q[3];
rz(-0.59800076) q[3];
sx q[3];
rz(-1.7015308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3368095) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(-3.0554904) q[2];
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
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.7030199) q[0];
sx q[0];
rz(-0.95412552) q[0];
sx q[0];
rz(2.2254206) q[0];
rz(2.84962) q[1];
sx q[1];
rz(-2.5501854) q[1];
sx q[1];
rz(2.3672262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385696) q[0];
sx q[0];
rz(-2.6576808) q[0];
sx q[0];
rz(0.83173521) q[0];
rz(-1.2254481) q[2];
sx q[2];
rz(-1.576445) q[2];
sx q[2];
rz(-2.408509) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.214433) q[1];
sx q[1];
rz(-1.4478874) q[1];
sx q[1];
rz(1.6729808) q[1];
rz(2.5995042) q[3];
sx q[3];
rz(-2.2418368) q[3];
sx q[3];
rz(-0.71806353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3463717) q[2];
sx q[2];
rz(-2.5796311) q[2];
sx q[2];
rz(0.58491659) q[2];
rz(1.6631205) q[3];
sx q[3];
rz(-1.9641179) q[3];
sx q[3];
rz(-1.6519206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4855578) q[0];
sx q[0];
rz(-1.0233044) q[0];
sx q[0];
rz(0.61505944) q[0];
rz(2.8133605) q[1];
sx q[1];
rz(-2.131772) q[1];
sx q[1];
rz(2.097791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8449677) q[0];
sx q[0];
rz(-3.0888978) q[0];
sx q[0];
rz(-0.096676306) q[0];
rz(-1.5659901) q[2];
sx q[2];
rz(-1.9991181) q[2];
sx q[2];
rz(-2.7132963) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4281299) q[1];
sx q[1];
rz(-0.9232115) q[1];
sx q[1];
rz(2.6414899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4099156) q[3];
sx q[3];
rz(-2.6291375) q[3];
sx q[3];
rz(2.3872914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7528167) q[2];
sx q[2];
rz(-0.39377585) q[2];
sx q[2];
rz(0.30889312) q[2];
rz(2.654352) q[3];
sx q[3];
rz(-1.4332708) q[3];
sx q[3];
rz(0.88700956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798645) q[0];
sx q[0];
rz(-0.74493113) q[0];
sx q[0];
rz(0.47759011) q[0];
rz(-1.8567122) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(-0.46599785) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14529252) q[0];
sx q[0];
rz(-1.0178653) q[0];
sx q[0];
rz(1.6607222) q[0];
rz(-pi) q[1];
rz(-1.1144756) q[2];
sx q[2];
rz(-1.9812036) q[2];
sx q[2];
rz(0.23882751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.029115) q[1];
sx q[1];
rz(-2.3211109) q[1];
sx q[1];
rz(0.047767834) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4336045) q[3];
sx q[3];
rz(-0.82105112) q[3];
sx q[3];
rz(2.193303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3232702) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(1.6935879) q[2];
rz(-1.4130392) q[3];
sx q[3];
rz(-1.6084684) q[3];
sx q[3];
rz(-1.1477227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2959761) q[0];
sx q[0];
rz(-2.4281261) q[0];
sx q[0];
rz(-0.071685858) q[0];
rz(-1.6638311) q[1];
sx q[1];
rz(-1.7774589) q[1];
sx q[1];
rz(2.51547) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6090995) q[0];
sx q[0];
rz(-0.83864826) q[0];
sx q[0];
rz(3.0200543) q[0];
rz(-pi) q[1];
rz(1.6501872) q[2];
sx q[2];
rz(-2.491386) q[2];
sx q[2];
rz(-2.2967867) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0542293) q[1];
sx q[1];
rz(-2.3647671) q[1];
sx q[1];
rz(-1.3224949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.463571) q[3];
sx q[3];
rz(-0.39059535) q[3];
sx q[3];
rz(0.56915802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6757297) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(-1.6443058) q[2];
rz(0.40694445) q[3];
sx q[3];
rz(-0.99138433) q[3];
sx q[3];
rz(2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.96754718) q[0];
sx q[0];
rz(-1.8438735) q[0];
sx q[0];
rz(0.26200727) q[0];
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
rz(2.8190099) q[0];
sx q[0];
rz(-0.8958502) q[0];
sx q[0];
rz(-0.92827601) q[0];
rz(-pi) q[1];
rz(-2.1706743) q[2];
sx q[2];
rz(-0.80764233) q[2];
sx q[2];
rz(-1.3563479) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3221601) q[1];
sx q[1];
rz(-1.3272078) q[1];
sx q[1];
rz(-2.0889682) q[1];
x q[2];
rz(-0.035484826) q[3];
sx q[3];
rz(-1.4789797) q[3];
sx q[3];
rz(0.47210708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1205341) q[2];
sx q[2];
rz(-0.48212442) q[2];
sx q[2];
rz(-2.6684707) q[2];
rz(2.1038697) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(-1.13824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4859908) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(-0.77014357) q[0];
rz(-0.4153525) q[1];
sx q[1];
rz(-1.5130006) q[1];
sx q[1];
rz(-1.279668) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9433728) q[0];
sx q[0];
rz(-0.83838338) q[0];
sx q[0];
rz(-3.1081568) q[0];
rz(-1.8417613) q[2];
sx q[2];
rz(-2.695796) q[2];
sx q[2];
rz(-0.4933756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7918167) q[1];
sx q[1];
rz(-2.8703051) q[1];
sx q[1];
rz(0.94601814) q[1];
rz(0.47904738) q[3];
sx q[3];
rz(-1.5622514) q[3];
sx q[3];
rz(-0.44310492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27560774) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(-1.9870997) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(1.1852468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.481021) q[0];
sx q[0];
rz(-0.72089973) q[0];
sx q[0];
rz(-0.50203669) q[0];
rz(1.0188811) q[1];
sx q[1];
rz(-1.4733682) q[1];
sx q[1];
rz(0.92207164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2577137) q[0];
sx q[0];
rz(-1.3259616) q[0];
sx q[0];
rz(-1.4725982) q[0];
rz(-1.3241295) q[2];
sx q[2];
rz(-2.1804159) q[2];
sx q[2];
rz(0.39928809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3905976) q[1];
sx q[1];
rz(-1.5158236) q[1];
sx q[1];
rz(-3.0979348) q[1];
x q[2];
rz(-1.6985252) q[3];
sx q[3];
rz(-1.5683577) q[3];
sx q[3];
rz(1.4245486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33793882) q[2];
sx q[2];
rz(-1.2793845) q[2];
sx q[2];
rz(3.0214018) q[2];
rz(3.1386612) q[3];
sx q[3];
rz(-2.7836697) q[3];
sx q[3];
rz(0.005793747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1015162) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(-1.6325604) q[0];
rz(-0.82935968) q[1];
sx q[1];
rz(-2.6396535) q[1];
sx q[1];
rz(1.1306184) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7735159) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(2.3115524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4254029) q[2];
sx q[2];
rz(-1.9702868) q[2];
sx q[2];
rz(-1.3376118) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8309511) q[1];
sx q[1];
rz(-2.0071908) q[1];
sx q[1];
rz(-0.92943014) q[1];
rz(1.9958903) q[3];
sx q[3];
rz(-2.5833874) q[3];
sx q[3];
rz(-2.3987215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.074038) q[2];
sx q[2];
rz(-2.270547) q[2];
sx q[2];
rz(1.3646431) q[2];
rz(-0.015297628) q[3];
sx q[3];
rz(-0.8316032) q[3];
sx q[3];
rz(-1.3219705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010654) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(2.9891253) q[0];
rz(-1.3854148) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(0.32858953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65909678) q[0];
sx q[0];
rz(-1.5616816) q[0];
sx q[0];
rz(3.1279596) q[0];
rz(-2.8893438) q[2];
sx q[2];
rz(-2.3636732) q[2];
sx q[2];
rz(-1.8053448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1235134) q[1];
sx q[1];
rz(-1.6005375) q[1];
sx q[1];
rz(1.8090387) q[1];
rz(-pi) q[2];
rz(0.64528193) q[3];
sx q[3];
rz(-1.8580164) q[3];
sx q[3];
rz(2.9336799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80455989) q[2];
sx q[2];
rz(-1.1351265) q[2];
sx q[2];
rz(-2.0648773) q[2];
rz(-2.1870901) q[3];
sx q[3];
rz(-1.7007622) q[3];
sx q[3];
rz(-0.29630989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(-2.238965) q[2];
sx q[2];
rz(-3.0716574) q[2];
sx q[2];
rz(2.8989094) q[2];
rz(-1.7794505) q[3];
sx q[3];
rz(-0.91051523) q[3];
sx q[3];
rz(-0.12403535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
