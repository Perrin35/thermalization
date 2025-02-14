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
rz(3.0132063) q[0];
rz(0.95500359) q[1];
sx q[1];
rz(-2.4817012) q[1];
sx q[1];
rz(0.58712062) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059537236) q[0];
sx q[0];
rz(-2.035241) q[0];
sx q[0];
rz(-1.4091253) q[0];
rz(2.4346515) q[2];
sx q[2];
rz(-2.4924008) q[2];
sx q[2];
rz(-2.6892218) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0260604) q[1];
sx q[1];
rz(-0.75775131) q[1];
sx q[1];
rz(0.65961403) q[1];
rz(-pi) q[2];
rz(3.1389075) q[3];
sx q[3];
rz(-1.955282) q[3];
sx q[3];
rz(0.17385305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8762274) q[2];
sx q[2];
rz(-0.23398016) q[2];
sx q[2];
rz(-2.2205676) q[2];
rz(1.5246576) q[3];
sx q[3];
rz(-1.8504668) q[3];
sx q[3];
rz(-2.9558712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8601473) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(1.349378) q[0];
rz(-0.34740627) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.7385534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31579298) q[0];
sx q[0];
rz(-2.1496338) q[0];
sx q[0];
rz(-1.3114503) q[0];
rz(-2.3187693) q[2];
sx q[2];
rz(-2.0206656) q[2];
sx q[2];
rz(-0.28385401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35866061) q[1];
sx q[1];
rz(-2.4110849) q[1];
sx q[1];
rz(0.23083861) q[1];
x q[2];
rz(-0.88256695) q[3];
sx q[3];
rz(-1.360713) q[3];
sx q[3];
rz(-2.9391367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8189524) q[2];
sx q[2];
rz(-1.8337269) q[2];
sx q[2];
rz(-2.9827706) q[2];
rz(-0.21765503) q[3];
sx q[3];
rz(-1.3889775) q[3];
sx q[3];
rz(-1.3505107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.7749216) q[0];
sx q[0];
rz(-1.3008302) q[0];
sx q[0];
rz(1.8721254) q[0];
rz(-2.7216351) q[1];
sx q[1];
rz(-1.0008413) q[1];
sx q[1];
rz(-1.6023844) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29572648) q[0];
sx q[0];
rz(-2.5877366) q[0];
sx q[0];
rz(0.94828301) q[0];
rz(-pi) q[1];
rz(0.2197651) q[2];
sx q[2];
rz(-1.3914492) q[2];
sx q[2];
rz(-2.1795535) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91130891) q[1];
sx q[1];
rz(-2.742594) q[1];
sx q[1];
rz(-1.7085275) q[1];
rz(-pi) q[2];
rz(0.49323757) q[3];
sx q[3];
rz(-0.62574157) q[3];
sx q[3];
rz(-1.7097676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.082108214) q[2];
sx q[2];
rz(-0.84222811) q[2];
sx q[2];
rz(-2.9498937) q[2];
rz(-0.55255237) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(2.1171169) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7315652) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(0.60436526) q[0];
rz(-1.2454237) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(0.85974685) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5167546) q[0];
sx q[0];
rz(-0.72480145) q[0];
sx q[0];
rz(-1.047931) q[0];
x q[1];
rz(-1.672869) q[2];
sx q[2];
rz(-1.5849176) q[2];
sx q[2];
rz(0.532224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0319718) q[1];
sx q[1];
rz(-1.6205677) q[1];
sx q[1];
rz(-0.18178908) q[1];
rz(-0.72034658) q[3];
sx q[3];
rz(-1.3784153) q[3];
sx q[3];
rz(1.1921574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2802281) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(1.4385983) q[2];
rz(1.8493308) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(-2.0869702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0629405) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(-2.8237421) q[0];
rz(1.3392797) q[1];
sx q[1];
rz(-0.41792089) q[1];
sx q[1];
rz(-2.1628765) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1237549) q[0];
sx q[0];
rz(-2.7394501) q[0];
sx q[0];
rz(3.0642302) q[0];
x q[1];
rz(-1.6597346) q[2];
sx q[2];
rz(-1.7654037) q[2];
sx q[2];
rz(0.22502514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7589583) q[1];
sx q[1];
rz(-0.82373229) q[1];
sx q[1];
rz(1.7428223) q[1];
rz(-2.948165) q[3];
sx q[3];
rz(-1.3816091) q[3];
sx q[3];
rz(2.6003169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9432482) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70500526) q[0];
sx q[0];
rz(-0.26708189) q[0];
sx q[0];
rz(1.3096814) q[0];
rz(2.1197223) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(1.6530564) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5989424) q[0];
sx q[0];
rz(-2.3841954) q[0];
sx q[0];
rz(2.0661584) q[0];
x q[1];
rz(-1.3371633) q[2];
sx q[2];
rz(-0.90590817) q[2];
sx q[2];
rz(0.49371749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52681953) q[1];
sx q[1];
rz(-2.062791) q[1];
sx q[1];
rz(1.4089877) q[1];
rz(-pi) q[2];
x q[2];
rz(1.390914) q[3];
sx q[3];
rz(-0.56821403) q[3];
sx q[3];
rz(0.0051604963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66270193) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(2.877511) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5298117) q[0];
sx q[0];
rz(-1.6407069) q[0];
sx q[0];
rz(-2.0773326) q[0];
rz(-3.0423959) q[1];
sx q[1];
rz(-1.7925037) q[1];
sx q[1];
rz(1.681021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997064) q[0];
sx q[0];
rz(-2.0545022) q[0];
sx q[0];
rz(1.4987491) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2719897) q[2];
sx q[2];
rz(-1.6083711) q[2];
sx q[2];
rz(-0.36848289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7713732) q[1];
sx q[1];
rz(-1.8533195) q[1];
sx q[1];
rz(-0.4346967) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0123294) q[3];
sx q[3];
rz(-1.0881502) q[3];
sx q[3];
rz(-3.1114374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7907685) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(-2.5420945) q[2];
rz(2.2439469) q[3];
sx q[3];
rz(-1.41956) q[3];
sx q[3];
rz(0.038014855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3277603) q[0];
sx q[0];
rz(-2.0812415) q[0];
sx q[0];
rz(-1.4014442) q[0];
rz(0.071488149) q[1];
sx q[1];
rz(-1.6782327) q[1];
sx q[1];
rz(1.00114) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18703774) q[0];
sx q[0];
rz(-1.1251377) q[0];
sx q[0];
rz(-1.8513239) q[0];
rz(-pi) q[1];
rz(1.1835353) q[2];
sx q[2];
rz(-0.94410482) q[2];
sx q[2];
rz(1.3037217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65232507) q[1];
sx q[1];
rz(-1.4389924) q[1];
sx q[1];
rz(2.5062923) q[1];
x q[2];
rz(2.4372479) q[3];
sx q[3];
rz(-2.337269) q[3];
sx q[3];
rz(-0.75596228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1772168) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(0.71844086) q[2];
rz(0.046772379) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30963323) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(-2.9613103) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(0.32275018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6868387) q[0];
sx q[0];
rz(-1.4658064) q[0];
sx q[0];
rz(0.15430321) q[0];
rz(1.5813002) q[2];
sx q[2];
rz(-2.5894508) q[2];
sx q[2];
rz(-0.10607468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1230939) q[1];
sx q[1];
rz(-1.0375751) q[1];
sx q[1];
rz(0.64792222) q[1];
rz(-0.75876816) q[3];
sx q[3];
rz(-2.3009217) q[3];
sx q[3];
rz(-1.981995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6705769) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(-2.1237109) q[2];
rz(-0.0035303591) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(1.2118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3280535) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(0.92593431) q[0];
rz(-0.13735859) q[1];
sx q[1];
rz(-0.67247144) q[1];
sx q[1];
rz(-0.59991178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5515251) q[0];
sx q[0];
rz(-0.54393629) q[0];
sx q[0];
rz(0.005225709) q[0];
rz(-pi) q[1];
rz(0.7057759) q[2];
sx q[2];
rz(-2.3786491) q[2];
sx q[2];
rz(-1.1102138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8132282) q[1];
sx q[1];
rz(-2.1592685) q[1];
sx q[1];
rz(0.036582734) q[1];
rz(-1.0968911) q[3];
sx q[3];
rz(-1.8623433) q[3];
sx q[3];
rz(-2.4959559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4124734) q[2];
sx q[2];
rz(-1.0189265) q[2];
sx q[2];
rz(0.64175433) q[2];
rz(1.1563835) q[3];
sx q[3];
rz(-0.76120794) q[3];
sx q[3];
rz(1.523264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625576) q[0];
sx q[0];
rz(-0.87974822) q[0];
sx q[0];
rz(0.77597822) q[0];
rz(1.4398126) q[1];
sx q[1];
rz(-1.4714614) q[1];
sx q[1];
rz(0.068838483) q[1];
rz(-2.4859602) q[2];
sx q[2];
rz(-0.31972319) q[2];
sx q[2];
rz(2.2421851) q[2];
rz(-2.2960881) q[3];
sx q[3];
rz(-2.4561524) q[3];
sx q[3];
rz(2.0292841) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
