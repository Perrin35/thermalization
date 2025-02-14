OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17570198) q[0];
sx q[0];
rz(5.0134563) q[0];
sx q[0];
rz(11.391517) q[0];
rz(2.0570316) q[1];
sx q[1];
rz(-1.6697872) q[1];
sx q[1];
rz(-0.56301277) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63516894) q[0];
sx q[0];
rz(-0.20316589) q[0];
sx q[0];
rz(-0.26655339) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1839574) q[2];
sx q[2];
rz(-2.4102825) q[2];
sx q[2];
rz(-0.72145185) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35511103) q[1];
sx q[1];
rz(-1.4110089) q[1];
sx q[1];
rz(0.33336877) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7607542) q[3];
sx q[3];
rz(-1.2361174) q[3];
sx q[3];
rz(-3.1370441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53434831) q[2];
sx q[2];
rz(-1.6515942) q[2];
sx q[2];
rz(-1.6687757) q[2];
rz(-1.8913174) q[3];
sx q[3];
rz(-1.5273124) q[3];
sx q[3];
rz(2.1691624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8984574) q[0];
sx q[0];
rz(-0.055483015) q[0];
sx q[0];
rz(-2.6836416) q[0];
rz(-0.10174879) q[1];
sx q[1];
rz(-1.9513756) q[1];
sx q[1];
rz(0.32624689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20083429) q[0];
sx q[0];
rz(-2.5837008) q[0];
sx q[0];
rz(-0.87403654) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63133755) q[2];
sx q[2];
rz(-2.7777618) q[2];
sx q[2];
rz(1.199786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7375398) q[1];
sx q[1];
rz(-2.0758325) q[1];
sx q[1];
rz(2.2801599) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2539034) q[3];
sx q[3];
rz(-2.2884048) q[3];
sx q[3];
rz(-0.16555351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.132167) q[2];
sx q[2];
rz(-2.1475809) q[2];
sx q[2];
rz(1.8780635) q[2];
rz(-0.020091232) q[3];
sx q[3];
rz(-2.3720436) q[3];
sx q[3];
rz(1.9115537) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1175213) q[0];
sx q[0];
rz(-3.0776403) q[0];
sx q[0];
rz(-2.5653895) q[0];
rz(1.6974712) q[1];
sx q[1];
rz(-1.8918119) q[1];
sx q[1];
rz(1.8796657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0281953) q[0];
sx q[0];
rz(-1.602766) q[0];
sx q[0];
rz(1.2616004) q[0];
rz(1.1088702) q[2];
sx q[2];
rz(-1.0903768) q[2];
sx q[2];
rz(-1.2269208) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8518119) q[1];
sx q[1];
rz(-1.3262311) q[1];
sx q[1];
rz(0.62965392) q[1];
rz(0.62278449) q[3];
sx q[3];
rz(-2.3238164) q[3];
sx q[3];
rz(-2.0501514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5835517) q[2];
sx q[2];
rz(-2.5463107) q[2];
sx q[2];
rz(2.9884647) q[2];
rz(-2.2554452) q[3];
sx q[3];
rz(-1.7821507) q[3];
sx q[3];
rz(1.9396293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2578289) q[0];
sx q[0];
rz(-1.9572636) q[0];
sx q[0];
rz(0.15876874) q[0];
rz(1.0348882) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(2.3854756) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0466111) q[0];
sx q[0];
rz(-1.621052) q[0];
sx q[0];
rz(-3.0836225) q[0];
rz(1.6768084) q[2];
sx q[2];
rz(-2.8869625) q[2];
sx q[2];
rz(-0.98843304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19263527) q[1];
sx q[1];
rz(-1.7757984) q[1];
sx q[1];
rz(0.24688772) q[1];
rz(-1.3642374) q[3];
sx q[3];
rz(-1.3404003) q[3];
sx q[3];
rz(-1.2588866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5990344) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(2.2912045) q[2];
rz(1.6709857) q[3];
sx q[3];
rz(-1.3266404) q[3];
sx q[3];
rz(2.3175122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4849512) q[0];
sx q[0];
rz(-1.4018207) q[0];
sx q[0];
rz(0.17855074) q[0];
rz(1.2483596) q[1];
sx q[1];
rz(-1.6105885) q[1];
sx q[1];
rz(2.607883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8781136) q[0];
sx q[0];
rz(-1.3878667) q[0];
sx q[0];
rz(-0.35533743) q[0];
x q[1];
rz(-2.1458964) q[2];
sx q[2];
rz(-0.786869) q[2];
sx q[2];
rz(2.5520419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0977749) q[1];
sx q[1];
rz(-2.2378451) q[1];
sx q[1];
rz(-0.066867877) q[1];
rz(-1.3946831) q[3];
sx q[3];
rz(-0.51100547) q[3];
sx q[3];
rz(-0.7738302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3137714) q[2];
sx q[2];
rz(-1.8905819) q[2];
sx q[2];
rz(2.0884464) q[2];
rz(-0.18038067) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(0.36261305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5617705) q[0];
sx q[0];
rz(-1.0012015) q[0];
sx q[0];
rz(-1.2286105) q[0];
rz(0.46663943) q[1];
sx q[1];
rz(-1.9231611) q[1];
sx q[1];
rz(3.012588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056696293) q[0];
sx q[0];
rz(-1.4375303) q[0];
sx q[0];
rz(-0.10198089) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3903046) q[2];
sx q[2];
rz(-0.73798385) q[2];
sx q[2];
rz(-2.5904417) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5552703) q[1];
sx q[1];
rz(-0.96809371) q[1];
sx q[1];
rz(-1.9307677) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0971076) q[3];
sx q[3];
rz(-0.29957459) q[3];
sx q[3];
rz(-2.4921474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1351607) q[2];
sx q[2];
rz(-2.5401523) q[2];
sx q[2];
rz(0.50430164) q[2];
rz(1.0935316) q[3];
sx q[3];
rz(-1.8176327) q[3];
sx q[3];
rz(0.95660153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.6495551) q[0];
sx q[0];
rz(-1.2832337) q[0];
sx q[0];
rz(1.5564224) q[0];
rz(-2.5632312) q[1];
sx q[1];
rz(-1.2588986) q[1];
sx q[1];
rz(-3.0392821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3123042) q[0];
sx q[0];
rz(-1.967038) q[0];
sx q[0];
rz(0.90950729) q[0];
rz(-1.8004216) q[2];
sx q[2];
rz(-1.9582677) q[2];
sx q[2];
rz(-1.4619399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8312586) q[1];
sx q[1];
rz(-1.4142904) q[1];
sx q[1];
rz(-0.308728) q[1];
rz(-2.3915406) q[3];
sx q[3];
rz(-1.0033718) q[3];
sx q[3];
rz(-2.9134083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9461225) q[2];
sx q[2];
rz(-1.8778233) q[2];
sx q[2];
rz(0.61402399) q[2];
rz(2.2139003) q[3];
sx q[3];
rz(-1.5427019) q[3];
sx q[3];
rz(2.4641002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880661) q[0];
sx q[0];
rz(-2.8480242) q[0];
sx q[0];
rz(-1.9422096) q[0];
rz(-1.1792432) q[1];
sx q[1];
rz(-2.113138) q[1];
sx q[1];
rz(2.1882122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74443734) q[0];
sx q[0];
rz(-2.1040475) q[0];
sx q[0];
rz(-0.6102194) q[0];
rz(-pi) q[1];
rz(-2.1087476) q[2];
sx q[2];
rz(-1.6643833) q[2];
sx q[2];
rz(0.46582281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.868321) q[1];
sx q[1];
rz(-1.6810604) q[1];
sx q[1];
rz(0.83302193) q[1];
rz(-2.6598831) q[3];
sx q[3];
rz(-1.6407971) q[3];
sx q[3];
rz(-0.21177975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74932468) q[2];
sx q[2];
rz(-3.0159123) q[2];
sx q[2];
rz(2.9607062) q[2];
rz(1.1130029) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(-0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0522633) q[0];
sx q[0];
rz(-0.90173975) q[0];
sx q[0];
rz(2.3774636) q[0];
rz(1.0230505) q[1];
sx q[1];
rz(-1.6108395) q[1];
sx q[1];
rz(-1.6463564) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8500124) q[0];
sx q[0];
rz(-0.33667013) q[0];
sx q[0];
rz(0.48179896) q[0];
x q[1];
rz(2.1713397) q[2];
sx q[2];
rz(-1.2926599) q[2];
sx q[2];
rz(0.32626611) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5519515) q[1];
sx q[1];
rz(-1.1855339) q[1];
sx q[1];
rz(-1.6478219) q[1];
rz(-0.71919433) q[3];
sx q[3];
rz(-0.56108863) q[3];
sx q[3];
rz(2.6230375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9493147) q[2];
sx q[2];
rz(-1.01769) q[2];
sx q[2];
rz(0.44753543) q[2];
rz(0.13628515) q[3];
sx q[3];
rz(-2.648573) q[3];
sx q[3];
rz(-0.59806699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78551453) q[0];
sx q[0];
rz(-1.0250174) q[0];
sx q[0];
rz(-2.6688975) q[0];
rz(-0.39974943) q[1];
sx q[1];
rz(-0.70416299) q[1];
sx q[1];
rz(2.3230816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8560573) q[0];
sx q[0];
rz(-1.7489079) q[0];
sx q[0];
rz(-2.8277603) q[0];
rz(0.50869663) q[2];
sx q[2];
rz(-2.3983068) q[2];
sx q[2];
rz(2.7270779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.012072969) q[1];
sx q[1];
rz(-2.8188779) q[1];
sx q[1];
rz(-1.5606776) q[1];
rz(-pi) q[2];
rz(2.6916041) q[3];
sx q[3];
rz(-2.7594824) q[3];
sx q[3];
rz(1.1799605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5588536) q[2];
sx q[2];
rz(-1.0820729) q[2];
sx q[2];
rz(1.2782798) q[2];
rz(1.9418955) q[3];
sx q[3];
rz(-1.7007098) q[3];
sx q[3];
rz(-0.27967683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3661135) q[0];
sx q[0];
rz(-1.4237325) q[0];
sx q[0];
rz(1.7970418) q[0];
rz(-1.4115502) q[1];
sx q[1];
rz(-0.81137864) q[1];
sx q[1];
rz(2.484533) q[1];
rz(-1.5700339) q[2];
sx q[2];
rz(-1.7117725) q[2];
sx q[2];
rz(2.6374809) q[2];
rz(-1.1468519) q[3];
sx q[3];
rz(-2.4320571) q[3];
sx q[3];
rz(0.71362389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
