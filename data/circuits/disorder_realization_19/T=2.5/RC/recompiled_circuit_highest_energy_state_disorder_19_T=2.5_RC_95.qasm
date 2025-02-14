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
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(0.81508842) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(-1.7307245) q[1];
sx q[1];
rz(-1.0676395) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0995432) q[0];
sx q[0];
rz(-2.0021515) q[0];
sx q[0];
rz(3.1154446) q[0];
rz(-pi) q[1];
rz(1.9733623) q[2];
sx q[2];
rz(-2.0595008) q[2];
sx q[2];
rz(-0.89948621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5112665) q[1];
sx q[1];
rz(-2.8095803) q[1];
sx q[1];
rz(-1.8707596) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88495636) q[3];
sx q[3];
rz(-1.7772632) q[3];
sx q[3];
rz(-1.5790758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69433576) q[2];
sx q[2];
rz(-2.3981636) q[2];
sx q[2];
rz(1.2747964) q[2];
rz(-0.46191195) q[3];
sx q[3];
rz(-2.4670944) q[3];
sx q[3];
rz(1.1326724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79134113) q[0];
sx q[0];
rz(-2.8332062) q[0];
sx q[0];
rz(1.8885008) q[0];
rz(-2.99627) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(1.0911509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6515132) q[0];
sx q[0];
rz(-1.7848178) q[0];
sx q[0];
rz(2.0064615) q[0];
rz(-pi) q[1];
rz(1.2107466) q[2];
sx q[2];
rz(-2.0031679) q[2];
sx q[2];
rz(0.52456174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4636865) q[1];
sx q[1];
rz(-2.8645201) q[1];
sx q[1];
rz(-0.67986791) q[1];
rz(1.9698148) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-0.53129133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48065177) q[2];
sx q[2];
rz(-1.3903684) q[2];
sx q[2];
rz(-2.6118028) q[2];
rz(-0.79331136) q[3];
sx q[3];
rz(-1.5373693) q[3];
sx q[3];
rz(-2.5180499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48524258) q[0];
sx q[0];
rz(-2.1915477) q[0];
sx q[0];
rz(-0.52870885) q[0];
rz(-2.5953925) q[1];
sx q[1];
rz(-0.9587973) q[1];
sx q[1];
rz(-0.34034696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211283) q[0];
sx q[0];
rz(-1.2437151) q[0];
sx q[0];
rz(1.2820679) q[0];
x q[1];
rz(-1.4252001) q[2];
sx q[2];
rz(-1.8032089) q[2];
sx q[2];
rz(2.6111141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7527377) q[1];
sx q[1];
rz(-1.5675263) q[1];
sx q[1];
rz(-3.0470303) q[1];
rz(-1.163655) q[3];
sx q[3];
rz(-0.7837067) q[3];
sx q[3];
rz(0.78898417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(-0.42984143) q[2];
rz(-0.75628453) q[3];
sx q[3];
rz(-2.9679306) q[3];
sx q[3];
rz(0.82591301) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569358) q[0];
sx q[0];
rz(-1.6357559) q[0];
sx q[0];
rz(1.6480308) q[0];
rz(2.3736296) q[1];
sx q[1];
rz(-0.47052828) q[1];
sx q[1];
rz(0.54642645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9523729) q[0];
sx q[0];
rz(-1.5060987) q[0];
sx q[0];
rz(1.4995262) q[0];
rz(-pi) q[1];
rz(2.3713263) q[2];
sx q[2];
rz(-0.67521836) q[2];
sx q[2];
rz(0.16083052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.029303251) q[1];
sx q[1];
rz(-1.7186856) q[1];
sx q[1];
rz(1.4724971) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2222142) q[3];
sx q[3];
rz(-2.1438144) q[3];
sx q[3];
rz(2.1185377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7118608) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(0.30409733) q[2];
rz(1.7763304) q[3];
sx q[3];
rz(-1.2208166) q[3];
sx q[3];
rz(2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51233184) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(-3.1410134) q[0];
rz(2.0594788) q[1];
sx q[1];
rz(-0.52434701) q[1];
sx q[1];
rz(2.2023315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88193306) q[0];
sx q[0];
rz(-2.006001) q[0];
sx q[0];
rz(0.70910221) q[0];
x q[1];
rz(-2.0505191) q[2];
sx q[2];
rz(-0.66893286) q[2];
sx q[2];
rz(0.26593966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68785948) q[1];
sx q[1];
rz(-0.61111585) q[1];
sx q[1];
rz(1.9495717) q[1];
rz(1.2573832) q[3];
sx q[3];
rz(-1.7872918) q[3];
sx q[3];
rz(0.27319841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.038736343) q[2];
sx q[2];
rz(-1.2612217) q[2];
sx q[2];
rz(-0.2612513) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.4120925) q[3];
sx q[3];
rz(-0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29702) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(0.20508668) q[0];
rz(-2.0948441) q[1];
sx q[1];
rz(-1.7457242) q[1];
sx q[1];
rz(2.7632025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8240067) q[0];
sx q[0];
rz(-1.3841624) q[0];
sx q[0];
rz(1.9621137) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4055875) q[2];
sx q[2];
rz(-2.0370315) q[2];
sx q[2];
rz(-1.8777443) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4503895) q[1];
sx q[1];
rz(-0.65237633) q[1];
sx q[1];
rz(-2.1689586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69244416) q[3];
sx q[3];
rz(-1.9117578) q[3];
sx q[3];
rz(-2.8657262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64688524) q[2];
sx q[2];
rz(-1.1581706) q[2];
sx q[2];
rz(1.0542487) q[2];
rz(-0.65822893) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(-0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9996416) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(-2.5010338) q[0];
rz(2.7236252) q[1];
sx q[1];
rz(-1.9211946) q[1];
sx q[1];
rz(0.80498615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98950878) q[0];
sx q[0];
rz(-1.594127) q[0];
sx q[0];
rz(0.70744608) q[0];
rz(-1.9897132) q[2];
sx q[2];
rz(-0.81406128) q[2];
sx q[2];
rz(2.4474622) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0143483) q[1];
sx q[1];
rz(-1.5509203) q[1];
sx q[1];
rz(-1.5424506) q[1];
x q[2];
rz(1.2748313) q[3];
sx q[3];
rz(-2.1419542) q[3];
sx q[3];
rz(-2.7829602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(2.3804046) q[2];
rz(2.393764) q[3];
sx q[3];
rz(-1.3176094) q[3];
sx q[3];
rz(0.67659155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51650301) q[0];
sx q[0];
rz(-2.857132) q[0];
sx q[0];
rz(1.4768584) q[0];
rz(-0.45627108) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(-2.2705073) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81646279) q[0];
sx q[0];
rz(-1.5832381) q[0];
sx q[0];
rz(-2.5274656) q[0];
rz(-pi) q[1];
rz(1.7256868) q[2];
sx q[2];
rz(-1.911507) q[2];
sx q[2];
rz(2.129385) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91034094) q[1];
sx q[1];
rz(-1.3518855) q[1];
sx q[1];
rz(0.33556767) q[1];
rz(-pi) q[2];
rz(1.1444451) q[3];
sx q[3];
rz(-1.5160069) q[3];
sx q[3];
rz(-2.3940115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90290922) q[2];
sx q[2];
rz(-0.52672714) q[2];
sx q[2];
rz(0.61863679) q[2];
rz(-0.020261852) q[3];
sx q[3];
rz(-0.94347763) q[3];
sx q[3];
rz(-2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0778377) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(0.42386398) q[0];
rz(-2.9546812) q[1];
sx q[1];
rz(-0.82004768) q[1];
sx q[1];
rz(1.4580457) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20987147) q[0];
sx q[0];
rz(-0.11882028) q[0];
sx q[0];
rz(-3.0221536) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8700022) q[2];
sx q[2];
rz(-1.0729861) q[2];
sx q[2];
rz(-1.9099727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.600454) q[1];
sx q[1];
rz(-2.9935874) q[1];
sx q[1];
rz(-0.35548325) q[1];
x q[2];
rz(-2.7856876) q[3];
sx q[3];
rz(-1.6545466) q[3];
sx q[3];
rz(2.8474396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3905048) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(1.0941774) q[2];
rz(1.68082) q[3];
sx q[3];
rz(-1.3760309) q[3];
sx q[3];
rz(-1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8193034) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(1.5018916) q[0];
rz(-0.27885258) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(-1.0241114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49660027) q[0];
sx q[0];
rz(-1.4787714) q[0];
sx q[0];
rz(-2.651398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2959997) q[2];
sx q[2];
rz(-1.6195903) q[2];
sx q[2];
rz(0.21367197) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1361724) q[1];
sx q[1];
rz(-1.7917624) q[1];
sx q[1];
rz(-1.3653838) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82641469) q[3];
sx q[3];
rz(-1.3152988) q[3];
sx q[3];
rz(2.5662553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9378822) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(0.53696519) q[2];
rz(-1.2282061) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(2.8502407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0252329) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(-2.4346726) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(-0.42264414) q[2];
sx q[2];
rz(-0.44036897) q[2];
sx q[2];
rz(-1.8457495) q[2];
rz(-1.3523921) q[3];
sx q[3];
rz(-1.4256178) q[3];
sx q[3];
rz(-2.3322879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
