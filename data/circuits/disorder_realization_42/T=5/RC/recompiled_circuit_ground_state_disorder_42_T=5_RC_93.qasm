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
rz(4.3262859) q[0];
sx q[0];
rz(10.560137) q[0];
rz(-0.42179498) q[1];
sx q[1];
rz(3.9219895) q[1];
sx q[1];
rz(11.94413) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6506158) q[0];
sx q[0];
rz(-0.19189206) q[0];
sx q[0];
rz(0.86928113) q[0];
rz(0.038921629) q[2];
sx q[2];
rz(-2.0177671) q[2];
sx q[2];
rz(1.0124504) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2127004) q[1];
sx q[1];
rz(-1.5824741) q[1];
sx q[1];
rz(-0.59439962) q[1];
rz(-pi) q[2];
rz(-2.1417326) q[3];
sx q[3];
rz(-1.7594764) q[3];
sx q[3];
rz(2.7267371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8047831) q[2];
sx q[2];
rz(-1.6853036) q[2];
sx q[2];
rz(-3.0554904) q[2];
rz(0.51186776) q[3];
sx q[3];
rz(-1.0443001) q[3];
sx q[3];
rz(2.4477203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43857273) q[0];
sx q[0];
rz(-2.1874671) q[0];
sx q[0];
rz(2.2254206) q[0];
rz(-2.84962) q[1];
sx q[1];
rz(-0.59140721) q[1];
sx q[1];
rz(2.3672262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60304927) q[0];
sx q[0];
rz(-1.2197681) q[0];
sx q[0];
rz(-0.34026628) q[0];
rz(-pi) q[1];
rz(-1.5541114) q[2];
sx q[2];
rz(-2.7962001) q[2];
sx q[2];
rz(2.2881803) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3437924) q[1];
sx q[1];
rz(-1.4693854) q[1];
sx q[1];
rz(3.0180458) q[1];
rz(-pi) q[2];
rz(2.5995042) q[3];
sx q[3];
rz(-0.89975587) q[3];
sx q[3];
rz(0.71806353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3463717) q[2];
sx q[2];
rz(-0.56196153) q[2];
sx q[2];
rz(-2.5566761) q[2];
rz(1.4784721) q[3];
sx q[3];
rz(-1.9641179) q[3];
sx q[3];
rz(-1.4896721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4855578) q[0];
sx q[0];
rz(-2.1182883) q[0];
sx q[0];
rz(-0.61505944) q[0];
rz(-2.8133605) q[1];
sx q[1];
rz(-2.131772) q[1];
sx q[1];
rz(-2.097791) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8449677) q[0];
sx q[0];
rz(-3.0888978) q[0];
sx q[0];
rz(3.0449163) q[0];
x q[1];
rz(0.42832613) q[2];
sx q[2];
rz(-1.5664243) q[2];
sx q[2];
rz(-1.1405038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7134628) q[1];
sx q[1];
rz(-0.9232115) q[1];
sx q[1];
rz(2.6414899) q[1];
rz(-1.9303334) q[3];
sx q[3];
rz(-1.94424) q[3];
sx q[3];
rz(-1.5869753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7528167) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(2.8326995) q[2];
rz(0.4872407) q[3];
sx q[3];
rz(-1.7083218) q[3];
sx q[3];
rz(0.88700956) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798645) q[0];
sx q[0];
rz(-2.3966615) q[0];
sx q[0];
rz(-0.47759011) q[0];
rz(1.2848805) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(2.6755948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024740273) q[0];
sx q[0];
rz(-0.55944397) q[0];
sx q[0];
rz(0.14450216) q[0];
x q[1];
rz(0.79171574) q[2];
sx q[2];
rz(-2.5377065) q[2];
sx q[2];
rz(2.0144659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1124777) q[1];
sx q[1];
rz(-0.82048172) q[1];
sx q[1];
rz(3.0938248) q[1];
rz(2.4572506) q[3];
sx q[3];
rz(-2.0668233) q[3];
sx q[3];
rz(-3.0471714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3232702) q[2];
sx q[2];
rz(-0.44760901) q[2];
sx q[2];
rz(-1.4480048) q[2];
rz(-1.4130392) q[3];
sx q[3];
rz(-1.6084684) q[3];
sx q[3];
rz(1.9938699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2959761) q[0];
sx q[0];
rz(-0.71346658) q[0];
sx q[0];
rz(0.071685858) q[0];
rz(1.4777615) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(-2.51547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0984296) q[0];
sx q[0];
rz(-1.6610896) q[0];
sx q[0];
rz(-2.3066269) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2194835) q[2];
sx q[2];
rz(-1.6188237) q[2];
sx q[2];
rz(-0.66274984) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69543998) q[1];
sx q[1];
rz(-1.7439411) q[1];
sx q[1];
rz(-0.80953627) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3887229) q[3];
sx q[3];
rz(-1.2232878) q[3];
sx q[3];
rz(3.0680498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6757297) q[2];
sx q[2];
rz(-2.6079123) q[2];
sx q[2];
rz(-1.6443058) q[2];
rz(2.7346482) q[3];
sx q[3];
rz(-2.1502083) q[3];
sx q[3];
rz(-0.86010325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.1740455) q[0];
sx q[0];
rz(-1.2977192) q[0];
sx q[0];
rz(0.26200727) q[0];
rz(1.2213446) q[1];
sx q[1];
rz(-1.5120993) q[1];
sx q[1];
rz(1.9409723) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3308418) q[0];
sx q[0];
rz(-1.0839607) q[0];
sx q[0];
rz(2.3563516) q[0];
x q[1];
rz(2.2827705) q[2];
sx q[2];
rz(-1.9910275) q[2];
sx q[2];
rz(0.65606299) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8194325) q[1];
sx q[1];
rz(-1.3272078) q[1];
sx q[1];
rz(-1.0526245) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2030199) q[3];
sx q[3];
rz(-3.0431755) q[3];
sx q[3];
rz(-2.3000788) q[3];
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
rz(0.47312197) q[2];
rz(-1.037723) q[3];
sx q[3];
rz(-2.7389052) q[3];
sx q[3];
rz(-1.13824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6556019) q[0];
sx q[0];
rz(-0.36360535) q[0];
sx q[0];
rz(-0.77014357) q[0];
rz(-2.7262402) q[1];
sx q[1];
rz(-1.628592) q[1];
sx q[1];
rz(-1.279668) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19821985) q[0];
sx q[0];
rz(-2.3032093) q[0];
sx q[0];
rz(3.1081568) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2998313) q[2];
sx q[2];
rz(-2.695796) q[2];
sx q[2];
rz(0.4933756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99233124) q[1];
sx q[1];
rz(-1.7898954) q[1];
sx q[1];
rz(-2.9803139) q[1];
x q[2];
rz(-0.47904738) q[3];
sx q[3];
rz(-1.5793413) q[3];
sx q[3];
rz(-0.44310492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27560774) q[2];
sx q[2];
rz(-1.7476247) q[2];
sx q[2];
rz(-2.7371791) q[2];
rz(-1.9870997) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(-1.9563458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.481021) q[0];
sx q[0];
rz(-2.4206929) q[0];
sx q[0];
rz(2.639556) q[0];
rz(-1.0188811) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(0.92207164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8716836) q[0];
sx q[0];
rz(-2.8781664) q[0];
sx q[0];
rz(-2.7676537) q[0];
rz(-2.8052373) q[2];
sx q[2];
rz(-0.65170641) q[2];
sx q[2];
rz(0.8136047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82220158) q[1];
sx q[1];
rz(-1.5272045) q[1];
sx q[1];
rz(-1.6258214) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5899379) q[3];
sx q[3];
rz(-0.12775207) q[3];
sx q[3];
rz(-2.9763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33793882) q[2];
sx q[2];
rz(-1.8622082) q[2];
sx q[2];
rz(3.0214018) q[2];
rz(0.0029314824) q[3];
sx q[3];
rz(-2.7836697) q[3];
sx q[3];
rz(-0.005793747) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040076479) q[0];
sx q[0];
rz(-0.67779556) q[0];
sx q[0];
rz(-1.6325604) q[0];
rz(0.82935968) q[1];
sx q[1];
rz(-0.50193915) q[1];
sx q[1];
rz(1.1306184) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36807674) q[0];
sx q[0];
rz(-2.6720071) q[0];
sx q[0];
rz(2.3115524) q[0];
rz(-pi) q[1];
rz(1.0605325) q[2];
sx q[2];
rz(-0.9212554) q[2];
sx q[2];
rz(-0.559597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.045550195) q[1];
sx q[1];
rz(-2.1438731) q[1];
sx q[1];
rz(2.6144774) q[1];
x q[2];
rz(-2.088016) q[3];
sx q[3];
rz(-1.3505837) q[3];
sx q[3];
rz(2.6802879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.067554615) q[2];
sx q[2];
rz(-0.87104565) q[2];
sx q[2];
rz(1.7769495) q[2];
rz(3.126295) q[3];
sx q[3];
rz(-2.3099895) q[3];
sx q[3];
rz(-1.8196222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.44052723) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(2.9891253) q[0];
rz(1.7561779) q[1];
sx q[1];
rz(-2.0962174) q[1];
sx q[1];
rz(2.8130031) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6406031) q[0];
sx q[0];
rz(-0.016399212) q[0];
sx q[0];
rz(-2.5522405) q[0];
x q[1];
rz(-0.76184892) q[2];
sx q[2];
rz(-1.3947316) q[2];
sx q[2];
rz(-0.052964199) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0180792) q[1];
sx q[1];
rz(-1.6005375) q[1];
sx q[1];
rz(1.8090387) q[1];
x q[2];
rz(0.45654065) q[3];
sx q[3];
rz(-0.69788633) q[3];
sx q[3];
rz(2.1386351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3370328) q[2];
sx q[2];
rz(-1.1351265) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(-0.95450258) q[3];
sx q[3];
rz(-1.7007622) q[3];
sx q[3];
rz(-2.8452828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(1.8148707) q[0];
sx q[0];
rz(-2.1537415) q[0];
sx q[0];
rz(1.4065773) q[0];
rz(-1.3181435) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(0.90262765) q[2];
sx q[2];
rz(-3.0716574) q[2];
sx q[2];
rz(2.8989094) q[2];
rz(-1.3621422) q[3];
sx q[3];
rz(-2.2310774) q[3];
sx q[3];
rz(3.0175573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
