OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.92575443) q[0];
sx q[0];
rz(2.821142) q[0];
sx q[0];
rz(9.5531643) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(2.554472) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7328661) q[0];
sx q[0];
rz(-0.48983296) q[0];
sx q[0];
rz(0.31087713) q[0];
x q[1];
rz(-1.112818) q[2];
sx q[2];
rz(-1.0931778) q[2];
sx q[2];
rz(-2.7736563) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11553227) q[1];
sx q[1];
rz(-0.75775131) q[1];
sx q[1];
rz(2.4819786) q[1];
rz(-pi) q[2];
rz(1.5641603) q[3];
sx q[3];
rz(-2.757098) q[3];
sx q[3];
rz(-0.16669434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2653653) q[2];
sx q[2];
rz(-0.23398016) q[2];
sx q[2];
rz(0.9210251) q[2];
rz(1.6169351) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(-2.9558712) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28144535) q[0];
sx q[0];
rz(-0.39491072) q[0];
sx q[0];
rz(-1.7922147) q[0];
rz(0.34740627) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.4030392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7424514) q[0];
sx q[0];
rz(-1.7871532) q[0];
sx q[0];
rz(0.59451806) q[0];
x q[1];
rz(0.58248708) q[2];
sx q[2];
rz(-0.91160027) q[2];
sx q[2];
rz(-0.90345736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35866061) q[1];
sx q[1];
rz(-0.73050776) q[1];
sx q[1];
rz(0.23083861) q[1];
x q[2];
rz(1.2469134) q[3];
sx q[3];
rz(-2.4270456) q[3];
sx q[3];
rz(-1.1199878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3226402) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(-2.9827706) q[2];
rz(-0.21765503) q[3];
sx q[3];
rz(-1.3889775) q[3];
sx q[3];
rz(1.7910819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36667103) q[0];
sx q[0];
rz(-1.3008302) q[0];
sx q[0];
rz(-1.2694673) q[0];
rz(0.41995755) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(1.6023844) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29572648) q[0];
sx q[0];
rz(-2.5877366) q[0];
sx q[0];
rz(-0.94828301) q[0];
rz(-pi) q[1];
rz(1.7544657) q[2];
sx q[2];
rz(-1.3546126) q[2];
sx q[2];
rz(-2.4930097) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6091385) q[1];
sx q[1];
rz(-1.517432) q[1];
sx q[1];
rz(1.9664001) q[1];
x q[2];
rz(-0.49323757) q[3];
sx q[3];
rz(-2.5158511) q[3];
sx q[3];
rz(1.431825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0594844) q[2];
sx q[2];
rz(-0.84222811) q[2];
sx q[2];
rz(0.19169894) q[2];
rz(2.5890403) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(2.1171169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7315652) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(0.60436526) q[0];
rz(1.8961689) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(-2.2818458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.788279) q[0];
sx q[0];
rz(-1.2333567) q[0];
sx q[0];
rz(2.2252625) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4687237) q[2];
sx q[2];
rz(-1.5849176) q[2];
sx q[2];
rz(-0.532224) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5936232) q[1];
sx q[1];
rz(-1.7523578) q[1];
sx q[1];
rz(-1.6214002) q[1];
x q[2];
rz(-0.28713496) q[3];
sx q[3];
rz(-0.74112149) q[3];
sx q[3];
rz(-0.59313074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2802281) q[2];
sx q[2];
rz(-2.1275529) q[2];
sx q[2];
rz(1.7029943) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-2.3137213) q[3];
sx q[3];
rz(2.0869702) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0629405) q[0];
sx q[0];
rz(-0.27314726) q[0];
sx q[0];
rz(2.8237421) q[0];
rz(-1.802313) q[1];
sx q[1];
rz(-2.7236718) q[1];
sx q[1];
rz(-0.97871614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.933799) q[0];
sx q[0];
rz(-1.9716671) q[0];
sx q[0];
rz(1.6036556) q[0];
rz(-1.6597346) q[2];
sx q[2];
rz(-1.3761889) q[2];
sx q[2];
rz(2.9165675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5086245) q[1];
sx q[1];
rz(-2.3787254) q[1];
sx q[1];
rz(-0.18277017) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.948165) q[3];
sx q[3];
rz(-1.3816091) q[3];
sx q[3];
rz(-0.54127579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19834441) q[2];
sx q[2];
rz(-2.1965616) q[2];
sx q[2];
rz(1.1172969) q[2];
rz(2.9605588) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(1.1972637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.70500526) q[0];
sx q[0];
rz(-2.8745108) q[0];
sx q[0];
rz(-1.3096814) q[0];
rz(-2.1197223) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(1.4885363) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2383136) q[0];
sx q[0];
rz(-0.92172232) q[0];
sx q[0];
rz(0.42239503) q[0];
x q[1];
rz(2.4632881) q[2];
sx q[2];
rz(-1.3875752) q[2];
sx q[2];
rz(-1.2228633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0206581) q[1];
sx q[1];
rz(-1.7132732) q[1];
sx q[1];
rz(2.6441036) q[1];
rz(-2.1316577) q[3];
sx q[3];
rz(-1.6672242) q[3];
sx q[3];
rz(-1.7280462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4788907) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(2.877511) q[2];
rz(-1.665202) q[3];
sx q[3];
rz(-0.68015209) q[3];
sx q[3];
rz(-0.41980729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61178094) q[0];
sx q[0];
rz(-1.5008858) q[0];
sx q[0];
rz(2.0773326) q[0];
rz(0.09919676) q[1];
sx q[1];
rz(-1.7925037) q[1];
sx q[1];
rz(1.681021) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6457466) q[0];
sx q[0];
rz(-2.6529713) q[0];
sx q[0];
rz(-0.13617985) q[0];
rz(-1.6978092) q[2];
sx q[2];
rz(-0.30108967) q[2];
sx q[2];
rz(1.8179229) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3702195) q[1];
sx q[1];
rz(-1.8533195) q[1];
sx q[1];
rz(-0.4346967) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6163382) q[3];
sx q[3];
rz(-1.958985) q[3];
sx q[3];
rz(-1.7565911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3508241) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(-0.59949818) q[2];
rz(2.2439469) q[3];
sx q[3];
rz(-1.41956) q[3];
sx q[3];
rz(-3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3277603) q[0];
sx q[0];
rz(-2.0812415) q[0];
sx q[0];
rz(1.7401485) q[0];
rz(3.0701045) q[1];
sx q[1];
rz(-1.46336) q[1];
sx q[1];
rz(-2.1404526) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545549) q[0];
sx q[0];
rz(-1.1251377) q[0];
sx q[0];
rz(-1.8513239) q[0];
x q[1];
rz(-1.9580574) q[2];
sx q[2];
rz(-0.94410482) q[2];
sx q[2];
rz(-1.837871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82187958) q[1];
sx q[1];
rz(-0.94187573) q[1];
sx q[1];
rz(-1.7340388) q[1];
x q[2];
rz(0.70434477) q[3];
sx q[3];
rz(-0.80432361) q[3];
sx q[3];
rz(-0.75596228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1772168) q[2];
sx q[2];
rz(-1.5528677) q[2];
sx q[2];
rz(2.4231518) q[2];
rz(0.046772379) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(-0.56979805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319594) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(0.69044789) q[0];
rz(0.18028232) q[1];
sx q[1];
rz(-2.0803662) q[1];
sx q[1];
rz(2.8188425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13234102) q[0];
sx q[0];
rz(-1.4173495) q[0];
sx q[0];
rz(1.4645534) q[0];
rz(-pi) q[1];
rz(-1.5813002) q[2];
sx q[2];
rz(-0.5521419) q[2];
sx q[2];
rz(-0.10607468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96092827) q[1];
sx q[1];
rz(-0.81392787) q[1];
sx q[1];
rz(2.3673173) q[1];
rz(0.68122562) q[3];
sx q[3];
rz(-2.109057) q[3];
sx q[3];
rz(-2.1665239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47101578) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(-1.0178817) q[2];
rz(-3.1380623) q[3];
sx q[3];
rz(-2.1719833) q[3];
sx q[3];
rz(1.2118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.3280535) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(-2.2156583) q[0];
rz(0.13735859) q[1];
sx q[1];
rz(-0.67247144) q[1];
sx q[1];
rz(0.59991178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5900676) q[0];
sx q[0];
rz(-2.5976564) q[0];
sx q[0];
rz(0.005225709) q[0];
x q[1];
rz(2.4358168) q[2];
sx q[2];
rz(-2.3786491) q[2];
sx q[2];
rz(-2.0313789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39420262) q[1];
sx q[1];
rz(-2.5521186) q[1];
sx q[1];
rz(-1.516045) q[1];
rz(-2.152485) q[3];
sx q[3];
rz(-2.5910561) q[3];
sx q[3];
rz(1.4359695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4124734) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(2.4998383) q[2];
rz(1.9852091) q[3];
sx q[3];
rz(-2.3803847) q[3];
sx q[3];
rz(-1.6183287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079035096) q[0];
sx q[0];
rz(-2.2618444) q[0];
sx q[0];
rz(-2.3656144) q[0];
rz(-1.4398126) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(-0.25664888) q[2];
sx q[2];
rz(-1.7636074) q[2];
sx q[2];
rz(-3.1008813) q[2];
rz(-0.84550459) q[3];
sx q[3];
rz(-0.68544023) q[3];
sx q[3];
rz(-1.1123085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
