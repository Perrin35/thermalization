OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0903836) q[0];
sx q[0];
rz(3.9301498) q[0];
sx q[0];
rz(9.7333662) q[0];
rz(0.42449549) q[1];
sx q[1];
rz(4.2257809) q[1];
sx q[1];
rz(9.0492166) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449995) q[0];
sx q[0];
rz(-0.29175943) q[0];
sx q[0];
rz(-1.3048521) q[0];
rz(1.6029066) q[2];
sx q[2];
rz(-1.2689723) q[2];
sx q[2];
rz(0.15649199) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1218222) q[1];
sx q[1];
rz(-2.0149061) q[1];
sx q[1];
rz(-1.7239718) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4272903) q[3];
sx q[3];
rz(-1.5670027) q[3];
sx q[3];
rz(0.11303813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8033119) q[2];
sx q[2];
rz(-0.25140005) q[2];
sx q[2];
rz(-2.5982762) q[2];
rz(-1.733689) q[3];
sx q[3];
rz(-1.7270154) q[3];
sx q[3];
rz(0.89731115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498085) q[0];
sx q[0];
rz(-1.1487288) q[0];
sx q[0];
rz(-2.5445004) q[0];
rz(-1.0626556) q[1];
sx q[1];
rz(-2.5025044) q[1];
sx q[1];
rz(2.7426681) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2960348) q[0];
sx q[0];
rz(-0.71656967) q[0];
sx q[0];
rz(-0.96591732) q[0];
rz(-2.127803) q[2];
sx q[2];
rz(-1.5913871) q[2];
sx q[2];
rz(-2.2228732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.084735076) q[1];
sx q[1];
rz(-1.0001518) q[1];
sx q[1];
rz(-2.4407376) q[1];
x q[2];
rz(0.90090294) q[3];
sx q[3];
rz(-2.381061) q[3];
sx q[3];
rz(0.72872114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4661633) q[2];
sx q[2];
rz(-1.2389946) q[2];
sx q[2];
rz(2.0002401) q[2];
rz(0.91836786) q[3];
sx q[3];
rz(-2.4310302) q[3];
sx q[3];
rz(-0.58951283) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.454575) q[0];
sx q[0];
rz(-2.7641251) q[0];
sx q[0];
rz(0.17534176) q[0];
rz(-1.31458) q[1];
sx q[1];
rz(-1.250123) q[1];
sx q[1];
rz(-2.4741727) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1940031) q[0];
sx q[0];
rz(-2.0663446) q[0];
sx q[0];
rz(3.010514) q[0];
rz(-1.0834789) q[2];
sx q[2];
rz(-1.0333152) q[2];
sx q[2];
rz(1.4097241) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35386723) q[1];
sx q[1];
rz(-1.8866061) q[1];
sx q[1];
rz(1.5701957) q[1];
x q[2];
rz(0.36881558) q[3];
sx q[3];
rz(-0.64743667) q[3];
sx q[3];
rz(-1.2596697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0257001) q[2];
sx q[2];
rz(-2.1860055) q[2];
sx q[2];
rz(2.058775) q[2];
rz(-1.8358021) q[3];
sx q[3];
rz(-0.81597733) q[3];
sx q[3];
rz(-0.79735565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1378655) q[0];
sx q[0];
rz(-0.19345134) q[0];
sx q[0];
rz(0.61099148) q[0];
rz(-2.5773279) q[1];
sx q[1];
rz(-1.0287501) q[1];
sx q[1];
rz(0.70704031) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96933041) q[0];
sx q[0];
rz(-2.3670787) q[0];
sx q[0];
rz(-1.5277083) q[0];
rz(-pi) q[1];
rz(-1.2390751) q[2];
sx q[2];
rz(-1.1807071) q[2];
sx q[2];
rz(1.6029101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1752211) q[1];
sx q[1];
rz(-2.4680158) q[1];
sx q[1];
rz(2.3994956) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48816008) q[3];
sx q[3];
rz(-0.55503856) q[3];
sx q[3];
rz(-0.51923534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6881312) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(-1.0180265) q[2];
rz(-0.38797837) q[3];
sx q[3];
rz(-1.4464902) q[3];
sx q[3];
rz(-0.33897266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46868789) q[0];
sx q[0];
rz(-2.0141116) q[0];
sx q[0];
rz(-2.6197523) q[0];
rz(-0.28779596) q[1];
sx q[1];
rz(-0.58240533) q[1];
sx q[1];
rz(0.60307455) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0007684) q[0];
sx q[0];
rz(-0.10061564) q[0];
sx q[0];
rz(0.510143) q[0];
x q[1];
rz(-2.1289044) q[2];
sx q[2];
rz(-2.1753722) q[2];
sx q[2];
rz(0.1451491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9085616) q[1];
sx q[1];
rz(-2.7082293) q[1];
sx q[1];
rz(2.9410081) q[1];
rz(-pi) q[2];
rz(1.2279922) q[3];
sx q[3];
rz(-1.3206284) q[3];
sx q[3];
rz(-2.7968593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4294943) q[2];
sx q[2];
rz(-2.3488729) q[2];
sx q[2];
rz(-2.402044) q[2];
rz(-2.126501) q[3];
sx q[3];
rz(-0.52777094) q[3];
sx q[3];
rz(0.019006193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58448142) q[0];
sx q[0];
rz(-0.33753532) q[0];
sx q[0];
rz(2.2767516) q[0];
rz(-0.64019126) q[1];
sx q[1];
rz(-0.65018153) q[1];
sx q[1];
rz(-2.6472299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1304754) q[0];
sx q[0];
rz(-2.0164765) q[0];
sx q[0];
rz(0.22486873) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.041448822) q[2];
sx q[2];
rz(-2.1728656) q[2];
sx q[2];
rz(-0.96884251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4130654) q[1];
sx q[1];
rz(-1.5134437) q[1];
sx q[1];
rz(-1.5363218) q[1];
x q[2];
rz(-1.0406144) q[3];
sx q[3];
rz(-2.0145855) q[3];
sx q[3];
rz(0.63999301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.87259) q[2];
sx q[2];
rz(-2.5073017) q[2];
sx q[2];
rz(0.96160257) q[2];
rz(1.5087992) q[3];
sx q[3];
rz(-1.6274933) q[3];
sx q[3];
rz(3.0646724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85840571) q[0];
sx q[0];
rz(-0.34342331) q[0];
sx q[0];
rz(-2.0984233) q[0];
rz(-2.4647602) q[1];
sx q[1];
rz(-0.64462858) q[1];
sx q[1];
rz(-2.4136995) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71804172) q[0];
sx q[0];
rz(-1.4598993) q[0];
sx q[0];
rz(-1.9865722) q[0];
rz(-pi) q[1];
rz(-0.0038996242) q[2];
sx q[2];
rz(-2.3995288) q[2];
sx q[2];
rz(-0.33483792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3549865) q[1];
sx q[1];
rz(-0.97026541) q[1];
sx q[1];
rz(0.40734603) q[1];
rz(-pi) q[2];
rz(1.9123825) q[3];
sx q[3];
rz(-1.0083958) q[3];
sx q[3];
rz(-1.4319358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68208838) q[2];
sx q[2];
rz(-0.75152087) q[2];
sx q[2];
rz(-0.53306836) q[2];
rz(1.8367977) q[3];
sx q[3];
rz(-2.6663836) q[3];
sx q[3];
rz(-0.78002012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.099139) q[0];
sx q[0];
rz(-3.070153) q[0];
sx q[0];
rz(0.024600994) q[0];
rz(-0.050361659) q[1];
sx q[1];
rz(-2.5854526) q[1];
sx q[1];
rz(-0.76739001) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31130107) q[0];
sx q[0];
rz(-2.6008252) q[0];
sx q[0];
rz(-3.0917653) q[0];
rz(2.6219756) q[2];
sx q[2];
rz(-2.765968) q[2];
sx q[2];
rz(2.5779998) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.394262) q[1];
sx q[1];
rz(-2.39191) q[1];
sx q[1];
rz(-2.5452627) q[1];
rz(-1.1786091) q[3];
sx q[3];
rz(-1.0309891) q[3];
sx q[3];
rz(-1.3572049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6699827) q[2];
sx q[2];
rz(-0.47405425) q[2];
sx q[2];
rz(-1.6101884) q[2];
rz(-2.1286185) q[3];
sx q[3];
rz(-1.4644724) q[3];
sx q[3];
rz(0.60232919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7811964) q[0];
sx q[0];
rz(-0.66806) q[0];
sx q[0];
rz(0.49041954) q[0];
rz(-2.7410653) q[1];
sx q[1];
rz(-0.35919765) q[1];
sx q[1];
rz(1.5708956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2119031) q[0];
sx q[0];
rz(-2.4725998) q[0];
sx q[0];
rz(0.59697911) q[0];
x q[1];
rz(0.8428085) q[2];
sx q[2];
rz(-0.47177663) q[2];
sx q[2];
rz(-2.9695784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8388673) q[1];
sx q[1];
rz(-1.6714393) q[1];
sx q[1];
rz(-1.2194281) q[1];
x q[2];
rz(-0.83028806) q[3];
sx q[3];
rz(-0.92819302) q[3];
sx q[3];
rz(-2.682529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57764292) q[2];
sx q[2];
rz(-2.9190013) q[2];
sx q[2];
rz(-0.27601784) q[2];
rz(0.66463071) q[3];
sx q[3];
rz(-0.97739995) q[3];
sx q[3];
rz(2.9450534) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68263245) q[0];
sx q[0];
rz(-2.9660872) q[0];
sx q[0];
rz(1.578414) q[0];
rz(2.3873734) q[1];
sx q[1];
rz(-2.1140607) q[1];
sx q[1];
rz(0.064090699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9784713) q[0];
sx q[0];
rz(-0.45867094) q[0];
sx q[0];
rz(-1.6244829) q[0];
rz(-pi) q[1];
rz(-1.1214549) q[2];
sx q[2];
rz(-2.1394512) q[2];
sx q[2];
rz(1.519738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8386993) q[1];
sx q[1];
rz(-2.6461227) q[1];
sx q[1];
rz(1.7458781) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.552382) q[3];
sx q[3];
rz(-2.2452045) q[3];
sx q[3];
rz(2.8545344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2854707) q[2];
sx q[2];
rz(-2.8557114) q[2];
sx q[2];
rz(-0.90606436) q[2];
rz(1.9237349) q[3];
sx q[3];
rz(-1.2713615) q[3];
sx q[3];
rz(1.0539894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8080407) q[0];
sx q[0];
rz(-1.6186436) q[0];
sx q[0];
rz(2.9343395) q[0];
rz(-2.8605657) q[1];
sx q[1];
rz(-1.749975) q[1];
sx q[1];
rz(-0.71798807) q[1];
rz(2.9502921) q[2];
sx q[2];
rz(-0.40861599) q[2];
sx q[2];
rz(1.7392639) q[2];
rz(3.1140399) q[3];
sx q[3];
rz(-2.6644628) q[3];
sx q[3];
rz(0.62678496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
