OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(4.9217304) q[0];
sx q[0];
rz(11.187727) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987797) q[0];
sx q[0];
rz(-3.0516041) q[0];
sx q[0];
rz(0.16422693) q[0];
rz(-pi) q[1];
rz(1.7832463) q[2];
sx q[2];
rz(-2.6531086) q[2];
sx q[2];
rz(2.8993895) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7887468) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(-2.7098141) q[1];
x q[2];
rz(3.046002) q[3];
sx q[3];
rz(-0.88926892) q[3];
sx q[3];
rz(-0.77212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(-2.700581) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-2.5355693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5490897) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3620136) q[0];
sx q[0];
rz(-1.5148666) q[0];
sx q[0];
rz(1.5895784) q[0];
rz(-2.942191) q[2];
sx q[2];
rz(-1.6316895) q[2];
sx q[2];
rz(0.25564889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8801404) q[1];
sx q[1];
rz(-1.4457236) q[1];
sx q[1];
rz(0.90420453) q[1];
rz(-0.22963345) q[3];
sx q[3];
rz(-0.7080871) q[3];
sx q[3];
rz(-1.9830444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(0.31093591) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3804669) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(2.321373) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8504282) q[0];
sx q[0];
rz(-2.4463852) q[0];
sx q[0];
rz(-1.7466963) q[0];
rz(2.064346) q[2];
sx q[2];
rz(-1.764467) q[2];
sx q[2];
rz(0.94656241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0606544) q[1];
sx q[1];
rz(-1.2186236) q[1];
sx q[1];
rz(-2.0209795) q[1];
rz(1.1060171) q[3];
sx q[3];
rz(-1.6764063) q[3];
sx q[3];
rz(1.6184023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(-0.43831929) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.8211676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94538044) q[0];
sx q[0];
rz(-1.5604661) q[0];
sx q[0];
rz(-1.1611847) q[0];
rz(-pi) q[1];
rz(-0.45122066) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(0.95552432) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3418386) q[1];
sx q[1];
rz(-0.8952039) q[1];
sx q[1];
rz(1.231133) q[1];
rz(-0.96418013) q[3];
sx q[3];
rz(-2.2607431) q[3];
sx q[3];
rz(-2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-0.34238112) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(-0.86529055) q[0];
rz(-1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.3006166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39869719) q[0];
sx q[0];
rz(-1.5455751) q[0];
sx q[0];
rz(-2.7874649) q[0];
rz(-pi) q[1];
rz(-2.775035) q[2];
sx q[2];
rz(-1.5036811) q[2];
sx q[2];
rz(1.7442489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.865766) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(0.57002108) q[1];
rz(-pi) q[2];
rz(0.50118581) q[3];
sx q[3];
rz(-0.29705829) q[3];
sx q[3];
rz(0.81851573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(0.011750301) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5884429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9827305) q[0];
sx q[0];
rz(-1.2356865) q[0];
sx q[0];
rz(-1.9431252) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6278218) q[2];
sx q[2];
rz(-0.36643039) q[2];
sx q[2];
rz(1.0183522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4658818) q[1];
sx q[1];
rz(-1.6556182) q[1];
sx q[1];
rz(-2.431543) q[1];
rz(-2.2736069) q[3];
sx q[3];
rz(-1.6276974) q[3];
sx q[3];
rz(-1.9753319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(2.1645434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0316232) q[0];
sx q[0];
rz(-2.1620746) q[0];
sx q[0];
rz(-1.6615608) q[0];
rz(-0.00074978272) q[2];
sx q[2];
rz(-0.13204083) q[2];
sx q[2];
rz(-0.11242871) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.44927412) q[1];
sx q[1];
rz(-0.34808967) q[1];
sx q[1];
rz(1.7983789) q[1];
rz(-pi) q[2];
rz(-0.3427152) q[3];
sx q[3];
rz(-1.8748771) q[3];
sx q[3];
rz(-2.8976687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(-2.1288669) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-2.4429328) q[0];
rz(1.1220804) q[1];
sx q[1];
rz(-2.2955003) q[1];
sx q[1];
rz(-1.2493856) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69985336) q[0];
sx q[0];
rz(-0.2158567) q[0];
sx q[0];
rz(0.21780832) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84455372) q[2];
sx q[2];
rz(-0.65659467) q[2];
sx q[2];
rz(-0.086364634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74832143) q[1];
sx q[1];
rz(-1.5763348) q[1];
sx q[1];
rz(0.022151532) q[1];
rz(-pi) q[2];
rz(-1.1864248) q[3];
sx q[3];
rz(-2.4723408) q[3];
sx q[3];
rz(0.8347019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.9630986) q[2];
rz(-1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(-1.7199843) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(-0.59757772) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9874728) q[0];
sx q[0];
rz(-0.14053908) q[0];
sx q[0];
rz(-1.8494291) q[0];
x q[1];
rz(-2.3915646) q[2];
sx q[2];
rz(-2.4366597) q[2];
sx q[2];
rz(0.36545576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24039195) q[1];
sx q[1];
rz(-1.8933834) q[1];
sx q[1];
rz(2.8278973) q[1];
rz(-pi) q[2];
rz(-0.11573128) q[3];
sx q[3];
rz(-2.5525408) q[3];
sx q[3];
rz(2.2310886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(-0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(0.023660252) q[0];
rz(0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-0.67217174) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8199352) q[0];
sx q[0];
rz(-1.5760839) q[0];
sx q[0];
rz(1.5810285) q[0];
rz(-pi) q[1];
rz(2.6300738) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(-2.9715003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8476665) q[1];
sx q[1];
rz(-0.77984174) q[1];
sx q[1];
rz(0.49077175) q[1];
rz(0.32398128) q[3];
sx q[3];
rz(-2.6880662) q[3];
sx q[3];
rz(1.6739664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(-1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(0.72369408) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(1.9359246) q[2];
sx q[2];
rz(-2.7177313) q[2];
sx q[2];
rz(2.360366) q[2];
rz(-2.8430812) q[3];
sx q[3];
rz(-1.1805503) q[3];
sx q[3];
rz(-1.3211484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
