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
rz(-2.5193522) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4909768) q[0];
sx q[0];
rz(-0.19189206) q[0];
sx q[0];
rz(0.86928113) q[0];
x q[1];
rz(-1.4897935) q[2];
sx q[2];
rz(-2.6930444) q[2];
sx q[2];
rz(-1.1023006) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2127004) q[1];
sx q[1];
rz(-1.5591185) q[1];
sx q[1];
rz(-0.59439962) q[1];
x q[2];
rz(-2.9184266) q[3];
sx q[3];
rz(-2.1303749) q[3];
sx q[3];
rz(-2.1055438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8047831) q[2];
sx q[2];
rz(-1.4562891) q[2];
sx q[2];
rz(3.0554904) q[2];
rz(-2.6297249) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(0.69387236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43857273) q[0];
sx q[0];
rz(-0.95412552) q[0];
sx q[0];
rz(0.91617209) q[0];
rz(2.84962) q[1];
sx q[1];
rz(-0.59140721) q[1];
sx q[1];
rz(0.77436647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.403023) q[0];
sx q[0];
rz(-0.48391184) q[0];
sx q[0];
rz(2.3098574) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5541114) q[2];
sx q[2];
rz(-2.7962001) q[2];
sx q[2];
rz(2.2881803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3437924) q[1];
sx q[1];
rz(-1.6722073) q[1];
sx q[1];
rz(-3.0180458) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82335681) q[3];
sx q[3];
rz(-1.9867479) q[3];
sx q[3];
rz(0.49440629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.795221) q[2];
sx q[2];
rz(-0.56196153) q[2];
sx q[2];
rz(-2.5566761) q[2];
rz(-1.4784721) q[3];
sx q[3];
rz(-1.1774747) q[3];
sx q[3];
rz(1.6519206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4855578) q[0];
sx q[0];
rz(-2.1182883) q[0];
sx q[0];
rz(-2.5265332) q[0];
rz(-2.8133605) q[1];
sx q[1];
rz(-1.0098207) q[1];
sx q[1];
rz(2.097791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39343483) q[0];
sx q[0];
rz(-1.6232449) q[0];
sx q[0];
rz(-1.5657052) q[0];
rz(0.010525816) q[2];
sx q[2];
rz(-2.7132456) q[2];
sx q[2];
rz(2.7017252) q[2];
rz(pi/2) q[3];
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
rz(-0.39643328) q[3];
sx q[3];
rz(-1.2370438) q[3];
sx q[3];
rz(2.9891356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.388776) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(-0.30889312) q[2];
rz(0.4872407) q[3];
sx q[3];
rz(-1.7083218) q[3];
sx q[3];
rz(0.88700956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76172817) q[0];
sx q[0];
rz(-0.74493113) q[0];
sx q[0];
rz(-0.47759011) q[0];
rz(-1.8567122) q[1];
sx q[1];
rz(-1.711859) q[1];
sx q[1];
rz(-2.6755948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024740273) q[0];
sx q[0];
rz(-0.55944397) q[0];
sx q[0];
rz(2.9970905) q[0];
rz(-pi) q[1];
x q[1];
rz(2.027117) q[2];
sx q[2];
rz(-1.9812036) q[2];
sx q[2];
rz(0.23882751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0424846) q[1];
sx q[1];
rz(-2.3900552) q[1];
sx q[1];
rz(1.5196176) q[1];
rz(-pi) q[2];
rz(0.68434207) q[3];
sx q[3];
rz(-1.0747693) q[3];
sx q[3];
rz(-3.0471714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8183225) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(1.4480048) q[2];
rz(-1.4130392) q[3];
sx q[3];
rz(-1.5331242) q[3];
sx q[3];
rz(1.1477227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84561658) q[0];
sx q[0];
rz(-0.71346658) q[0];
sx q[0];
rz(0.071685858) q[0];
rz(-1.4777615) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(2.51547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0984296) q[0];
sx q[0];
rz(-1.4805031) q[0];
sx q[0];
rz(0.8349658) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.060242816) q[2];
sx q[2];
rz(-0.92298302) q[2];
sx q[2];
rz(-0.94442764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4461527) q[1];
sx q[1];
rz(-1.7439411) q[1];
sx q[1];
rz(-0.80953627) q[1];
rz(-pi) q[2];
rz(-1.7528698) q[3];
sx q[3];
rz(-1.9183049) q[3];
sx q[3];
rz(-0.073542882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46586299) q[2];
sx q[2];
rz(-2.6079123) q[2];
sx q[2];
rz(1.4972868) q[2];
rz(-0.40694445) q[3];
sx q[3];
rz(-2.1502083) q[3];
sx q[3];
rz(2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96754718) q[0];
sx q[0];
rz(-1.2977192) q[0];
sx q[0];
rz(-2.8795854) q[0];
rz(-1.2213446) q[1];
sx q[1];
rz(-1.6294934) q[1];
sx q[1];
rz(-1.2006203) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8190099) q[0];
sx q[0];
rz(-0.8958502) q[0];
sx q[0];
rz(-0.92827601) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97091834) q[2];
sx q[2];
rz(-0.80764233) q[2];
sx q[2];
rz(1.7852448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.756304) q[1];
sx q[1];
rz(-1.0693764) q[1];
sx q[1];
rz(-0.27863592) q[1];
x q[2];
rz(0.035484826) q[3];
sx q[3];
rz(-1.662613) q[3];
sx q[3];
rz(-2.6694856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0210586) q[2];
sx q[2];
rz(-0.48212442) q[2];
sx q[2];
rz(0.47312197) q[2];
rz(2.1038697) q[3];
sx q[3];
rz(-0.40268746) q[3];
sx q[3];
rz(-2.0033526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6556019) q[0];
sx q[0];
rz(-2.7779873) q[0];
sx q[0];
rz(-2.3714491) q[0];
rz(-0.4153525) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(-1.8619246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466542) q[0];
sx q[0];
rz(-1.5956559) q[0];
sx q[0];
rz(2.3034873) q[0];
x q[1];
rz(2.0023022) q[2];
sx q[2];
rz(-1.6864634) q[2];
sx q[2];
rz(2.3097599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1492614) q[1];
sx q[1];
rz(-1.7898954) q[1];
sx q[1];
rz(2.9803139) q[1];
rz(-1.5804251) q[3];
sx q[3];
rz(-2.0498247) q[3];
sx q[3];
rz(-1.1232532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8659849) q[2];
sx q[2];
rz(-1.393968) q[2];
sx q[2];
rz(-0.40441355) q[2];
rz(1.9870997) q[3];
sx q[3];
rz(-1.5279852) q[3];
sx q[3];
rz(-1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(2.639556) q[0];
rz(1.0188811) q[1];
sx q[1];
rz(-1.4733682) q[1];
sx q[1];
rz(-2.219521) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4785504) q[0];
sx q[0];
rz(-1.4755357) q[0];
sx q[0];
rz(2.8956198) q[0];
x q[1];
rz(2.8052373) q[2];
sx q[2];
rz(-2.4898862) q[2];
sx q[2];
rz(-2.327988) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82220158) q[1];
sx q[1];
rz(-1.5272045) q[1];
sx q[1];
rz(-1.5157713) q[1];
rz(-1.5516547) q[3];
sx q[3];
rz(-0.12775207) q[3];
sx q[3];
rz(-2.9763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33793882) q[2];
sx q[2];
rz(-1.2793845) q[2];
sx q[2];
rz(3.0214018) q[2];
rz(-3.1386612) q[3];
sx q[3];
rz(-0.35792297) q[3];
sx q[3];
rz(-3.1357989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040076479) q[0];
sx q[0];
rz(-2.4637971) q[0];
sx q[0];
rz(-1.5090322) q[0];
rz(-2.312233) q[1];
sx q[1];
rz(-0.50193915) q[1];
sx q[1];
rz(-2.0109743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36807674) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(0.83004029) q[0];
x q[1];
rz(-0.71618979) q[2];
sx q[2];
rz(-1.1713059) q[2];
sx q[2];
rz(-1.3376118) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.045550195) q[1];
sx q[1];
rz(-0.99771959) q[1];
sx q[1];
rz(-2.6144774) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9958903) q[3];
sx q[3];
rz(-0.55820528) q[3];
sx q[3];
rz(-0.74287117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.074038) q[2];
sx q[2];
rz(-0.87104565) q[2];
sx q[2];
rz(-1.3646431) q[2];
rz(0.015297628) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(1.8196222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7010654) q[0];
sx q[0];
rz(-0.65182132) q[0];
sx q[0];
rz(-2.9891253) q[0];
rz(-1.3854148) q[1];
sx q[1];
rz(-1.0453753) q[1];
sx q[1];
rz(0.32858953) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2300174) q[0];
sx q[0];
rz(-1.5844288) q[0];
sx q[0];
rz(-1.5616807) q[0];
rz(-pi) q[1];
rz(-2.8893438) q[2];
sx q[2];
rz(-0.77791947) q[2];
sx q[2];
rz(-1.3362479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0180792) q[1];
sx q[1];
rz(-1.5410551) q[1];
sx q[1];
rz(-1.3325539) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2166546) q[3];
sx q[3];
rz(-0.956007) q[3];
sx q[3];
rz(-1.572991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80455989) q[2];
sx q[2];
rz(-2.0064662) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(2.1870901) q[3];
sx q[3];
rz(-1.7007622) q[3];
sx q[3];
rz(0.29630989) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8148707) q[0];
sx q[0];
rz(-0.9878511) q[0];
sx q[0];
rz(-1.7350154) q[0];
rz(1.8234491) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(1.6257269) q[2];
sx q[2];
rz(-1.5274897) q[2];
sx q[2];
rz(-2.4804583) q[2];
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
