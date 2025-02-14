OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2158382) q[0];
sx q[0];
rz(-2.821142) q[0];
sx q[0];
rz(-3.0132063) q[0];
rz(4.0965962) q[1];
sx q[1];
rz(8.7648865) q[1];
sx q[1];
rz(13.153491) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059537236) q[0];
sx q[0];
rz(-1.1063516) q[0];
sx q[0];
rz(-1.4091253) q[0];
rz(-pi) q[1];
rz(0.52337661) q[2];
sx q[2];
rz(-1.1672772) q[2];
sx q[2];
rz(-1.7159107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93361357) q[1];
sx q[1];
rz(-0.99665239) q[1];
sx q[1];
rz(1.0453392) q[1];
rz(-3.1389075) q[3];
sx q[3];
rz(-1.1863106) q[3];
sx q[3];
rz(0.17385305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8762274) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(2.2205676) q[2];
rz(-1.6169351) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(2.9558712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28144535) q[0];
sx q[0];
rz(-0.39491072) q[0];
sx q[0];
rz(1.7922147) q[0];
rz(2.7941864) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.7385534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3991412) q[0];
sx q[0];
rz(-1.7871532) q[0];
sx q[0];
rz(-2.5470746) q[0];
rz(-pi) q[1];
rz(-2.3187693) q[2];
sx q[2];
rz(-2.0206656) q[2];
sx q[2];
rz(-0.28385401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.053005347) q[1];
sx q[1];
rz(-0.8638051) q[1];
sx q[1];
rz(-1.3686351) q[1];
x q[2];
rz(2.2590257) q[3];
sx q[3];
rz(-1.360713) q[3];
sx q[3];
rz(-2.9391367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3226402) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(2.9827706) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36667103) q[0];
sx q[0];
rz(-1.8407624) q[0];
sx q[0];
rz(-1.8721254) q[0];
rz(2.7216351) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(1.5392083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8458662) q[0];
sx q[0];
rz(-2.5877366) q[0];
sx q[0];
rz(2.1933096) q[0];
rz(1.3871269) q[2];
sx q[2];
rz(-1.78698) q[2];
sx q[2];
rz(0.64858299) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.060614) q[1];
sx q[1];
rz(-1.1757869) q[1];
sx q[1];
rz(-3.0837713) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5747803) q[3];
sx q[3];
rz(-1.2897964) q[3];
sx q[3];
rz(-2.5917201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0594844) q[2];
sx q[2];
rz(-0.84222811) q[2];
sx q[2];
rz(2.9498937) q[2];
rz(-2.5890403) q[3];
sx q[3];
rz(-1.5250165) q[3];
sx q[3];
rz(2.1171169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4100274) q[0];
sx q[0];
rz(-0.64842328) q[0];
sx q[0];
rz(-2.5372274) q[0];
rz(-1.8961689) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(2.2818458) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3533136) q[0];
sx q[0];
rz(-1.2333567) q[0];
sx q[0];
rz(-2.2252625) q[0];
rz(-1.672869) q[2];
sx q[2];
rz(-1.5566751) q[2];
sx q[2];
rz(-0.532224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1096209) q[1];
sx q[1];
rz(-1.6205677) q[1];
sx q[1];
rz(0.18178908) q[1];
rz(-1.3172011) q[3];
sx q[3];
rz(-0.86652866) q[3];
sx q[3];
rz(-0.21237838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8613646) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(-1.4385983) q[2];
rz(-1.2922618) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(-2.0869702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0786521) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(2.8237421) q[0];
rz(1.3392797) q[1];
sx q[1];
rz(-2.7236718) q[1];
sx q[1];
rz(-0.97871614) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1237549) q[0];
sx q[0];
rz(-0.40214254) q[0];
sx q[0];
rz(-3.0642302) q[0];
rz(-pi) q[1];
rz(0.42338223) q[2];
sx q[2];
rz(-2.9278594) q[2];
sx q[2];
rz(-0.20704421) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7589583) q[1];
sx q[1];
rz(-2.3178604) q[1];
sx q[1];
rz(-1.7428223) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3781015) q[3];
sx q[3];
rz(-1.7607302) q[3];
sx q[3];
rz(-2.1488921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9432482) q[2];
sx q[2];
rz(-2.1965616) q[2];
sx q[2];
rz(-1.1172969) q[2];
rz(-0.18103389) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(1.1972637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(1.3096814) q[0];
rz(1.0218703) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(1.4885363) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90327901) q[0];
sx q[0];
rz(-0.92172232) q[0];
sx q[0];
rz(2.7191976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6783046) q[2];
sx q[2];
rz(-1.7540175) q[2];
sx q[2];
rz(1.2228633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9474839) q[1];
sx q[1];
rz(-0.51583898) q[1];
sx q[1];
rz(-0.2920002) q[1];
x q[2];
rz(-1.009935) q[3];
sx q[3];
rz(-1.4743685) q[3];
sx q[3];
rz(-1.7280462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4788907) q[2];
sx q[2];
rz(-1.3871437) q[2];
sx q[2];
rz(-2.877511) q[2];
rz(-1.665202) q[3];
sx q[3];
rz(-0.68015209) q[3];
sx q[3];
rz(-0.41980729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61178094) q[0];
sx q[0];
rz(-1.6407069) q[0];
sx q[0];
rz(2.0773326) q[0];
rz(-0.09919676) q[1];
sx q[1];
rz(-1.3490889) q[1];
sx q[1];
rz(-1.4605716) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.195358) q[0];
sx q[0];
rz(-1.6345662) q[0];
sx q[0];
rz(2.6568165) q[0];
rz(1.2719897) q[2];
sx q[2];
rz(-1.5332216) q[2];
sx q[2];
rz(-2.7731098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34014186) q[1];
sx q[1];
rz(-2.6281118) q[1];
sx q[1];
rz(-0.60350792) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1292633) q[3];
sx q[3];
rz(-2.0534424) q[3];
sx q[3];
rz(0.03015524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7907685) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(2.5420945) q[2];
rz(2.2439469) q[3];
sx q[3];
rz(-1.7220327) q[3];
sx q[3];
rz(3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81383234) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(-1.7401485) q[0];
rz(0.071488149) q[1];
sx q[1];
rz(-1.46336) q[1];
sx q[1];
rz(-1.00114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5073231) q[0];
sx q[0];
rz(-1.8232913) q[0];
sx q[0];
rz(-0.4613614) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4779448) q[2];
sx q[2];
rz(-1.8816684) q[2];
sx q[2];
rz(3.1093017) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3197131) q[1];
sx q[1];
rz(-2.1997169) q[1];
sx q[1];
rz(1.4075539) q[1];
rz(-0.66949943) q[3];
sx q[3];
rz(-2.0560798) q[3];
sx q[3];
rz(-2.8593331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9643758) q[2];
sx q[2];
rz(-1.588725) q[2];
sx q[2];
rz(2.4231518) q[2];
rz(-3.0948203) q[3];
sx q[3];
rz(-2.434157) q[3];
sx q[3];
rz(0.56979805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30963323) q[0];
sx q[0];
rz(-2.9871812) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(0.18028232) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(-2.8188425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0092516) q[0];
sx q[0];
rz(-1.7242431) q[0];
sx q[0];
rz(-1.4645534) q[0];
x q[1];
rz(1.0186791) q[2];
sx q[2];
rz(-1.5763057) q[2];
sx q[2];
rz(-1.4736648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2219991) q[1];
sx q[1];
rz(-1.0242435) q[1];
sx q[1];
rz(0.93354694) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75876816) q[3];
sx q[3];
rz(-2.3009217) q[3];
sx q[3];
rz(1.981995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47101578) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(2.1237109) q[2];
rz(0.0035303591) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(1.9297011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3280535) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(2.2156583) q[0];
rz(-0.13735859) q[1];
sx q[1];
rz(-2.4691212) q[1];
sx q[1];
rz(-2.5416809) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1178499) q[0];
sx q[0];
rz(-1.5735007) q[0];
sx q[0];
rz(-0.54393025) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5125347) q[2];
sx q[2];
rz(-2.0355844) q[2];
sx q[2];
rz(-2.128922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9194738) q[1];
sx q[1];
rz(-1.6012234) q[1];
sx q[1];
rz(-2.1595776) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32528721) q[3];
sx q[3];
rz(-2.023175) q[3];
sx q[3];
rz(-0.77879209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4124734) q[2];
sx q[2];
rz(-1.0189265) q[2];
sx q[2];
rz(-0.64175433) q[2];
rz(1.1563835) q[3];
sx q[3];
rz(-2.3803847) q[3];
sx q[3];
rz(-1.523264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625576) q[0];
sx q[0];
rz(-0.87974822) q[0];
sx q[0];
rz(0.77597822) q[0];
rz(1.7017801) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(-1.769968) q[2];
sx q[2];
rz(-1.3190075) q[2];
sx q[2];
rz(-1.5803303) q[2];
rz(-2.6445848) q[3];
sx q[3];
rz(-2.0642701) q[3];
sx q[3];
rz(-1.9652496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
