OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(-0.50146377) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8034536) q[0];
sx q[0];
rz(-2.8087466) q[0];
sx q[0];
rz(1.8344318) q[0];
rz(-2.2519977) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(-1.7878469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7079561) q[1];
sx q[1];
rz(-2.4383713) q[1];
sx q[1];
rz(-1.0183079) q[1];
x q[2];
rz(2.6775042) q[3];
sx q[3];
rz(-2.1543192) q[3];
sx q[3];
rz(2.1634963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(2.9843176) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(0.10903407) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63576525) q[0];
sx q[0];
rz(-2.8553537) q[0];
sx q[0];
rz(-0.92606996) q[0];
rz(-pi) q[1];
rz(0.63983812) q[2];
sx q[2];
rz(-2.8176762) q[2];
sx q[2];
rz(1.6537635) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1034531) q[1];
sx q[1];
rz(-2.1557169) q[1];
sx q[1];
rz(1.000688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3974959) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-2.6599191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047086296) q[0];
sx q[0];
rz(-1.9425834) q[0];
sx q[0];
rz(-1.0466172) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1842314) q[2];
sx q[2];
rz(-1.8463677) q[2];
sx q[2];
rz(2.7797109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8183221) q[1];
sx q[1];
rz(-2.1774877) q[1];
sx q[1];
rz(-2.2491749) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5211469) q[3];
sx q[3];
rz(-2.1944322) q[3];
sx q[3];
rz(-0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(-2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.4612173) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9445779) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(1.5532956) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(-1.2447371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23403215) q[0];
sx q[0];
rz(-1.2313156) q[0];
sx q[0];
rz(-1.5725122) q[0];
rz(2.142698) q[2];
sx q[2];
rz(-2.9036387) q[2];
sx q[2];
rz(-1.1513125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4501805) q[1];
sx q[1];
rz(-2.1069063) q[1];
sx q[1];
rz(-0.24713534) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34282116) q[3];
sx q[3];
rz(-2.5431513) q[3];
sx q[3];
rz(-0.10859057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7446639) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(-0.98608195) q[0];
rz(-pi) q[1];
rz(2.224515) q[2];
sx q[2];
rz(-3.0055025) q[2];
sx q[2];
rz(1.2166785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4245783) q[1];
sx q[1];
rz(-1.4687521) q[1];
sx q[1];
rz(-1.2574408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2522069) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(1.1226908) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(-0.52156633) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-0.7985324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.417056) q[0];
sx q[0];
rz(-0.31053156) q[0];
sx q[0];
rz(-2.0470371) q[0];
rz(-2.6449634) q[2];
sx q[2];
rz(-1.2611258) q[2];
sx q[2];
rz(1.3336099) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16520271) q[1];
sx q[1];
rz(-0.573728) q[1];
sx q[1];
rz(-1.4742673) q[1];
rz(-pi) q[2];
rz(-0.60243209) q[3];
sx q[3];
rz(-0.89655399) q[3];
sx q[3];
rz(0.34365052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(0.19730332) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(0.46404776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-0.95373017) q[0];
sx q[0];
rz(-2.6130136) q[0];
x q[1];
rz(-3.1214141) q[2];
sx q[2];
rz(-1.8115461) q[2];
sx q[2];
rz(-1.9764331) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6205008) q[1];
sx q[1];
rz(-1.8720086) q[1];
sx q[1];
rz(-2.0786933) q[1];
rz(-pi) q[2];
rz(0.23630996) q[3];
sx q[3];
rz(-1.4910306) q[3];
sx q[3];
rz(0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.4415178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38948108) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(1.5556637) q[0];
rz(1.5324138) q[2];
sx q[2];
rz(-2.1830609) q[2];
sx q[2];
rz(-1.6910764) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2802906) q[1];
sx q[1];
rz(-2.960627) q[1];
sx q[1];
rz(0.7678395) q[1];
rz(3.0184047) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(1.5817225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(-0.22658919) q[2];
rz(2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9771377) q[0];
sx q[0];
rz(-1.7286321) q[0];
sx q[0];
rz(2.6437003) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95790205) q[2];
sx q[2];
rz(-2.3447678) q[2];
sx q[2];
rz(0.56607841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0341067) q[1];
sx q[1];
rz(-1.429538) q[1];
sx q[1];
rz(0.60955255) q[1];
x q[2];
rz(2.968623) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(3.1239307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(-1.7260889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6575359) q[0];
sx q[0];
rz(-0.41175479) q[0];
sx q[0];
rz(0.62396892) q[0];
rz(-pi) q[1];
rz(1.2070933) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(1.9050913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.43511697) q[1];
sx q[1];
rz(-0.26285989) q[1];
sx q[1];
rz(0.95107066) q[1];
x q[2];
rz(-2.9726082) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71904174) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(-0.66108274) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-2.8812257) q[2];
sx q[2];
rz(-1.3961062) q[2];
sx q[2];
rz(0.41873742) q[2];
rz(1.4205167) q[3];
sx q[3];
rz(-2.1223304) q[3];
sx q[3];
rz(0.35029678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
