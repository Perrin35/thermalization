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
rz(-0.32045066) q[0];
sx q[0];
rz(-0.12838636) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(2.554472) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7328661) q[0];
sx q[0];
rz(-2.6517597) q[0];
sx q[0];
rz(-2.8307155) q[0];
rz(2.618216) q[2];
sx q[2];
rz(-1.9743154) q[2];
sx q[2];
rz(1.425682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93361357) q[1];
sx q[1];
rz(-2.1449403) q[1];
sx q[1];
rz(1.0453392) q[1];
x q[2];
rz(-1.9552833) q[3];
sx q[3];
rz(-1.5732854) q[3];
sx q[3];
rz(1.3979504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2653653) q[2];
sx q[2];
rz(-0.23398016) q[2];
sx q[2];
rz(0.9210251) q[2];
rz(-1.5246576) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(0.18572148) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28144535) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(1.349378) q[0];
rz(2.7941864) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(-1.4030392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058075) q[0];
sx q[0];
rz(-2.5134235) q[0];
sx q[0];
rz(-0.37395333) q[0];
rz(0.82282339) q[2];
sx q[2];
rz(-1.1209271) q[2];
sx q[2];
rz(-2.8577386) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.782932) q[1];
sx q[1];
rz(-0.73050776) q[1];
sx q[1];
rz(-2.910754) q[1];
rz(-pi) q[2];
rz(2.8722309) q[3];
sx q[3];
rz(-2.2410763) q[3];
sx q[3];
rz(1.538185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8189524) q[2];
sx q[2];
rz(-1.8337269) q[2];
sx q[2];
rz(2.9827706) q[2];
rz(0.21765503) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(1.7910819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749216) q[0];
sx q[0];
rz(-1.3008302) q[0];
sx q[0];
rz(-1.2694673) q[0];
rz(-0.41995755) q[1];
sx q[1];
rz(-1.0008413) q[1];
sx q[1];
rz(-1.5392083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7363971) q[0];
sx q[0];
rz(-2.0123082) q[0];
sx q[0];
rz(2.7955152) q[0];
x q[1];
rz(-2.4478618) q[2];
sx q[2];
rz(-0.28273928) q[2];
sx q[2];
rz(0.065134243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53245413) q[1];
sx q[1];
rz(-1.517432) q[1];
sx q[1];
rz(1.9664001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49323757) q[3];
sx q[3];
rz(-0.62574157) q[3];
sx q[3];
rz(1.7097676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0594844) q[2];
sx q[2];
rz(-2.2993645) q[2];
sx q[2];
rz(2.9498937) q[2];
rz(0.55255237) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4100274) q[0];
sx q[0];
rz(-0.64842328) q[0];
sx q[0];
rz(-0.60436526) q[0];
rz(1.8961689) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(-2.2818458) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.788279) q[0];
sx q[0];
rz(-1.9082359) q[0];
sx q[0];
rz(2.2252625) q[0];
rz(-0.014195125) q[2];
sx q[2];
rz(-1.6728587) q[2];
sx q[2];
rz(1.0400187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27450505) q[1];
sx q[1];
rz(-2.9531859) q[1];
sx q[1];
rz(-0.26885689) q[1];
rz(-pi) q[2];
rz(-1.8243916) q[3];
sx q[3];
rz(-2.275064) q[3];
sx q[3];
rz(2.9292143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2802281) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(-1.4385983) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(-2.0869702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.0786521) q[0];
sx q[0];
rz(-0.27314726) q[0];
sx q[0];
rz(-2.8237421) q[0];
rz(-1.3392797) q[1];
sx q[1];
rz(-2.7236718) q[1];
sx q[1];
rz(-2.1628765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178378) q[0];
sx q[0];
rz(-2.7394501) q[0];
sx q[0];
rz(3.0642302) q[0];
rz(-pi) q[1];
rz(-1.481858) q[2];
sx q[2];
rz(-1.3761889) q[2];
sx q[2];
rz(0.22502514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3826344) q[1];
sx q[1];
rz(-2.3178604) q[1];
sx q[1];
rz(1.7428223) q[1];
x q[2];
rz(-1.3781015) q[3];
sx q[3];
rz(-1.7607302) q[3];
sx q[3];
rz(0.99270051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19834441) q[2];
sx q[2];
rz(-2.1965616) q[2];
sx q[2];
rz(2.0242958) q[2];
rz(-0.18103389) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(-1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4365874) q[0];
sx q[0];
rz(-0.26708189) q[0];
sx q[0];
rz(-1.8319112) q[0];
rz(-1.0218703) q[1];
sx q[1];
rz(-0.836687) q[1];
sx q[1];
rz(1.4885363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5426503) q[0];
sx q[0];
rz(-0.7573973) q[0];
sx q[0];
rz(-1.0754343) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4632881) q[2];
sx q[2];
rz(-1.3875752) q[2];
sx q[2];
rz(-1.9187294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0206581) q[1];
sx q[1];
rz(-1.7132732) q[1];
sx q[1];
rz(-2.6441036) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1316577) q[3];
sx q[3];
rz(-1.4743685) q[3];
sx q[3];
rz(-1.4135464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66270193) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(-2.877511) q[2];
rz(1.665202) q[3];
sx q[3];
rz(-0.68015209) q[3];
sx q[3];
rz(0.41980729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5298117) q[0];
sx q[0];
rz(-1.5008858) q[0];
sx q[0];
rz(-1.06426) q[0];
rz(-3.0423959) q[1];
sx q[1];
rz(-1.7925037) q[1];
sx q[1];
rz(1.681021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6457466) q[0];
sx q[0];
rz(-0.48862132) q[0];
sx q[0];
rz(-3.0054128) q[0];
x q[1];
rz(-1.869603) q[2];
sx q[2];
rz(-1.5332216) q[2];
sx q[2];
rz(0.36848289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3702195) q[1];
sx q[1];
rz(-1.8533195) q[1];
sx q[1];
rz(0.4346967) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1292633) q[3];
sx q[3];
rz(-1.0881502) q[3];
sx q[3];
rz(-3.1114374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3508241) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(-0.59949818) q[2];
rz(-0.89764578) q[3];
sx q[3];
rz(-1.7220327) q[3];
sx q[3];
rz(3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81383234) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(-1.7401485) q[0];
rz(-3.0701045) q[1];
sx q[1];
rz(-1.6782327) q[1];
sx q[1];
rz(-2.1404526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6342696) q[0];
sx q[0];
rz(-1.8232913) q[0];
sx q[0];
rz(0.4613614) q[0];
rz(-0.48076081) q[2];
sx q[2];
rz(-0.72274739) q[2];
sx q[2];
rz(-1.9113845) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82187958) q[1];
sx q[1];
rz(-2.1997169) q[1];
sx q[1];
rz(1.7340388) q[1];
x q[2];
rz(-0.66949943) q[3];
sx q[3];
rz(-1.0855128) q[3];
sx q[3];
rz(2.8593331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1772168) q[2];
sx q[2];
rz(-1.5528677) q[2];
sx q[2];
rz(2.4231518) q[2];
rz(3.0948203) q[3];
sx q[3];
rz(-2.434157) q[3];
sx q[3];
rz(-0.56979805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8319594) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(2.9613103) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(2.8188425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13234102) q[0];
sx q[0];
rz(-1.4173495) q[0];
sx q[0];
rz(-1.4645534) q[0];
rz(-pi) q[1];
rz(1.5602925) q[2];
sx q[2];
rz(-2.5894508) q[2];
sx q[2];
rz(-3.035518) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2219991) q[1];
sx q[1];
rz(-2.1173491) q[1];
sx q[1];
rz(-0.93354694) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3828245) q[3];
sx q[3];
rz(-0.840671) q[3];
sx q[3];
rz(1.1595977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6705769) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(1.0178817) q[2];
rz(3.1380623) q[3];
sx q[3];
rz(-2.1719833) q[3];
sx q[3];
rz(-1.2118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8135391) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(-0.92593431) q[0];
rz(0.13735859) q[1];
sx q[1];
rz(-0.67247144) q[1];
sx q[1];
rz(-2.5416809) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900676) q[0];
sx q[0];
rz(-0.54393629) q[0];
sx q[0];
rz(-3.1363669) q[0];
rz(2.4358168) q[2];
sx q[2];
rz(-2.3786491) q[2];
sx q[2];
rz(1.1102138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39420262) q[1];
sx q[1];
rz(-2.5521186) q[1];
sx q[1];
rz(-1.6255476) q[1];
rz(-pi) q[2];
rz(-1.0968911) q[3];
sx q[3];
rz(-1.2792493) q[3];
sx q[3];
rz(-0.64563676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7291193) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(2.4998383) q[2];
rz(-1.1563835) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(1.6183287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625576) q[0];
sx q[0];
rz(-2.2618444) q[0];
sx q[0];
rz(-2.3656144) q[0];
rz(-1.4398126) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(2.4859602) q[2];
sx q[2];
rz(-2.8218695) q[2];
sx q[2];
rz(-0.8994076) q[2];
rz(0.84550459) q[3];
sx q[3];
rz(-2.4561524) q[3];
sx q[3];
rz(2.0292841) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
