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
rz(0.95500359) q[1];
sx q[1];
rz(-2.4817012) q[1];
sx q[1];
rz(0.58712062) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40872657) q[0];
sx q[0];
rz(-0.48983296) q[0];
sx q[0];
rz(-0.31087713) q[0];
rz(-pi) q[1];
rz(2.4346515) q[2];
sx q[2];
rz(-2.4924008) q[2];
sx q[2];
rz(-2.6892218) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94226664) q[1];
sx q[1];
rz(-1.1360511) q[1];
sx q[1];
rz(-2.4995657) q[1];
rz(-pi) q[2];
rz(1.5774324) q[3];
sx q[3];
rz(-0.38449461) q[3];
sx q[3];
rz(2.9748983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8762274) q[2];
sx q[2];
rz(-2.9076125) q[2];
sx q[2];
rz(0.9210251) q[2];
rz(-1.5246576) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(-2.9558712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058075) q[0];
sx q[0];
rz(-0.62816915) q[0];
sx q[0];
rz(-0.37395333) q[0];
rz(-0.82282339) q[2];
sx q[2];
rz(-2.0206656) q[2];
sx q[2];
rz(0.28385401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.385437) q[1];
sx q[1];
rz(-1.7240579) q[1];
sx q[1];
rz(-2.4244244) q[1];
rz(-pi) q[2];
rz(-0.26936172) q[3];
sx q[3];
rz(-0.90051631) q[3];
sx q[3];
rz(-1.538185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3226402) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(2.9827706) q[2];
rz(-0.21765503) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(1.3505107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36667103) q[0];
sx q[0];
rz(-1.3008302) q[0];
sx q[0];
rz(1.2694673) q[0];
rz(0.41995755) q[1];
sx q[1];
rz(-1.0008413) q[1];
sx q[1];
rz(-1.6023844) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8231118) q[0];
sx q[0];
rz(-1.2590908) q[0];
sx q[0];
rz(-2.0363755) q[0];
rz(-2.9218276) q[2];
sx q[2];
rz(-1.3914492) q[2];
sx q[2];
rz(0.96203912) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.060614) q[1];
sx q[1];
rz(-1.9658057) q[1];
sx q[1];
rz(-3.0837713) q[1];
rz(-pi) q[2];
rz(2.6483551) q[3];
sx q[3];
rz(-2.5158511) q[3];
sx q[3];
rz(1.431825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.082108214) q[2];
sx q[2];
rz(-0.84222811) q[2];
sx q[2];
rz(0.19169894) q[2];
rz(-0.55255237) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(-1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.4100274) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(2.5372274) q[0];
rz(1.2454237) q[1];
sx q[1];
rz(-2.3912997) q[1];
sx q[1];
rz(-2.2818458) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62483803) q[0];
sx q[0];
rz(-0.72480145) q[0];
sx q[0];
rz(-2.0936616) q[0];
rz(1.672869) q[2];
sx q[2];
rz(-1.5566751) q[2];
sx q[2];
rz(-2.6093687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27450505) q[1];
sx q[1];
rz(-2.9531859) q[1];
sx q[1];
rz(-2.8727358) q[1];
rz(-pi) q[2];
rz(2.8544577) q[3];
sx q[3];
rz(-0.74112149) q[3];
sx q[3];
rz(-0.59313074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8613646) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(1.4385983) q[2];
rz(-1.2922618) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(1.0546225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0629405) q[0];
sx q[0];
rz(-0.27314726) q[0];
sx q[0];
rz(-2.8237421) q[0];
rz(-1.802313) q[1];
sx q[1];
rz(-0.41792089) q[1];
sx q[1];
rz(0.97871614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1237549) q[0];
sx q[0];
rz(-0.40214254) q[0];
sx q[0];
rz(-3.0642302) q[0];
rz(-pi) q[1];
rz(-2.7182104) q[2];
sx q[2];
rz(-2.9278594) q[2];
sx q[2];
rz(-0.20704421) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6329681) q[1];
sx q[1];
rz(-0.76286722) q[1];
sx q[1];
rz(0.18277017) q[1];
rz(-pi) q[2];
rz(2.3581402) q[3];
sx q[3];
rz(-2.8718567) q[3];
sx q[3];
rz(1.7945643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9432482) q[2];
sx q[2];
rz(-0.94503108) q[2];
sx q[2];
rz(2.0242958) q[2];
rz(2.9605588) q[3];
sx q[3];
rz(-1.9046611) q[3];
sx q[3];
rz(1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4365874) q[0];
sx q[0];
rz(-0.26708189) q[0];
sx q[0];
rz(-1.3096814) q[0];
rz(-1.0218703) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(1.6530564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7393417) q[0];
sx q[0];
rz(-1.9034804) q[0];
sx q[0];
rz(-2.2646623) q[0];
rz(2.8544442) q[2];
sx q[2];
rz(-0.69881546) q[2];
sx q[2];
rz(0.1255807) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.52681953) q[1];
sx q[1];
rz(-2.062791) q[1];
sx q[1];
rz(-1.4089877) q[1];
rz(1.7506787) q[3];
sx q[3];
rz(-2.5733786) q[3];
sx q[3];
rz(0.0051604963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4788907) q[2];
sx q[2];
rz(-1.3871437) q[2];
sx q[2];
rz(2.877511) q[2];
rz(-1.4763907) q[3];
sx q[3];
rz(-2.4614406) q[3];
sx q[3];
rz(-0.41980729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.61178094) q[0];
sx q[0];
rz(-1.5008858) q[0];
sx q[0];
rz(1.06426) q[0];
rz(0.09919676) q[1];
sx q[1];
rz(-1.7925037) q[1];
sx q[1];
rz(1.681021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6457466) q[0];
sx q[0];
rz(-0.48862132) q[0];
sx q[0];
rz(-0.13617985) q[0];
rz(-3.1022775) q[2];
sx q[2];
rz(-1.2722071) q[2];
sx q[2];
rz(-1.1907426) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7713732) q[1];
sx q[1];
rz(-1.2882731) q[1];
sx q[1];
rz(0.4346967) q[1];
rz(-pi) q[2];
rz(-2.6163382) q[3];
sx q[3];
rz(-1.958985) q[3];
sx q[3];
rz(1.3850016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3508241) q[2];
sx q[2];
rz(-1.1188353) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81383234) q[0];
sx q[0];
rz(-2.0812415) q[0];
sx q[0];
rz(-1.4014442) q[0];
rz(3.0701045) q[1];
sx q[1];
rz(-1.6782327) q[1];
sx q[1];
rz(-1.00114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6342696) q[0];
sx q[0];
rz(-1.3183013) q[0];
sx q[0];
rz(-0.4613614) q[0];
rz(-pi) q[1];
rz(-1.1835353) q[2];
sx q[2];
rz(-2.1974878) q[2];
sx q[2];
rz(-1.837871) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3197131) q[1];
sx q[1];
rz(-0.94187573) q[1];
sx q[1];
rz(1.4075539) q[1];
x q[2];
rz(-2.4720932) q[3];
sx q[3];
rz(-2.0560798) q[3];
sx q[3];
rz(-0.2822596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9643758) q[2];
sx q[2];
rz(-1.5528677) q[2];
sx q[2];
rz(-2.4231518) q[2];
rz(-0.046772379) q[3];
sx q[3];
rz(-0.7074357) q[3];
sx q[3];
rz(0.56979805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8319594) q[0];
sx q[0];
rz(-0.15441144) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(-0.18028232) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(-0.32275018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6868387) q[0];
sx q[0];
rz(-1.6757863) q[0];
sx q[0];
rz(0.15430321) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5602925) q[2];
sx q[2];
rz(-2.5894508) q[2];
sx q[2];
rz(0.10607468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1806644) q[1];
sx q[1];
rz(-0.81392787) q[1];
sx q[1];
rz(2.3673173) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75876816) q[3];
sx q[3];
rz(-2.3009217) q[3];
sx q[3];
rz(-1.981995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6705769) q[2];
sx q[2];
rz(-1.7936423) q[2];
sx q[2];
rz(1.0178817) q[2];
rz(-0.0035303591) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(-1.9297011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.3280535) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(0.92593431) q[0];
rz(-3.0042341) q[1];
sx q[1];
rz(-2.4691212) q[1];
sx q[1];
rz(2.5416809) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.545418) q[0];
sx q[0];
rz(-2.1147244) q[0];
sx q[0];
rz(-1.5676359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7057759) q[2];
sx q[2];
rz(-2.3786491) q[2];
sx q[2];
rz(-1.1102138) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9194738) q[1];
sx q[1];
rz(-1.6012234) q[1];
sx q[1];
rz(-0.98201507) q[1];
x q[2];
rz(0.98910768) q[3];
sx q[3];
rz(-0.55053655) q[3];
sx q[3];
rz(-1.4359695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4124734) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(0.64175433) q[2];
rz(-1.9852091) q[3];
sx q[3];
rz(-2.3803847) q[3];
sx q[3];
rz(-1.523264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079035096) q[0];
sx q[0];
rz(-0.87974822) q[0];
sx q[0];
rz(0.77597822) q[0];
rz(1.7017801) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(-2.8849438) q[2];
sx q[2];
rz(-1.3779852) q[2];
sx q[2];
rz(0.040711395) q[2];
rz(-2.1199119) q[3];
sx q[3];
rz(-2.0041448) q[3];
sx q[3];
rz(-0.14295391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
