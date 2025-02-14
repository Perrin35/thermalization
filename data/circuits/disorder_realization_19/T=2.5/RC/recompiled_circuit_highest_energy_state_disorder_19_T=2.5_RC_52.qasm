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
rz(-1.018723) q[0];
sx q[0];
rz(-0.85917226) q[0];
sx q[0];
rz(2.3265042) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(-1.7307245) q[1];
sx q[1];
rz(-1.0676395) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1045221) q[0];
sx q[0];
rz(-2.7094954) q[0];
sx q[0];
rz(-1.62754) q[0];
x q[1];
rz(-2.5065866) q[2];
sx q[2];
rz(-2.5189812) q[2];
sx q[2];
rz(-2.9787763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9166482) q[1];
sx q[1];
rz(-1.6672581) q[1];
sx q[1];
rz(-1.8889844) q[1];
rz(-2.2566363) q[3];
sx q[3];
rz(-1.7772632) q[3];
sx q[3];
rz(-1.5790758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4472569) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(-1.2747964) q[2];
rz(-0.46191195) q[3];
sx q[3];
rz(-2.4670944) q[3];
sx q[3];
rz(1.1326724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79134113) q[0];
sx q[0];
rz(-2.8332062) q[0];
sx q[0];
rz(-1.8885008) q[0];
rz(0.14532267) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(1.0911509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4900794) q[0];
sx q[0];
rz(-1.3567748) q[0];
sx q[0];
rz(2.0064615) q[0];
rz(2.4895489) q[2];
sx q[2];
rz(-0.55527675) q[2];
sx q[2];
rz(2.9342143) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.021046358) q[1];
sx q[1];
rz(-1.7851549) q[1];
sx q[1];
rz(1.7477186) q[1];
x q[2];
rz(-2.6464822) q[3];
sx q[3];
rz(-1.9260599) q[3];
sx q[3];
rz(0.85384599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6609409) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(-2.6118028) q[2];
rz(2.3482813) q[3];
sx q[3];
rz(-1.5373693) q[3];
sx q[3];
rz(0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6563501) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(0.52870885) q[0];
rz(0.54620019) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(0.34034696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211283) q[0];
sx q[0];
rz(-1.2437151) q[0];
sx q[0];
rz(1.2820679) q[0];
rz(-pi) q[1];
rz(0.23480798) q[2];
sx q[2];
rz(-1.4291414) q[2];
sx q[2];
rz(1.0740785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9593411) q[1];
sx q[1];
rz(-1.6653582) q[1];
sx q[1];
rz(1.574081) q[1];
rz(-1.163655) q[3];
sx q[3];
rz(-0.7837067) q[3];
sx q[3];
rz(0.78898417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9693552) q[2];
sx q[2];
rz(-2.7293971) q[2];
sx q[2];
rz(-2.7117512) q[2];
rz(-0.75628453) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(1.6480308) q[0];
rz(0.76796302) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(0.54642645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38619216) q[0];
sx q[0];
rz(-1.4996756) q[0];
sx q[0];
rz(-3.0767308) q[0];
rz(-pi) q[1];
rz(-2.3713263) q[2];
sx q[2];
rz(-2.4663743) q[2];
sx q[2];
rz(-2.9807621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5560234) q[1];
sx q[1];
rz(-1.6680191) q[1];
sx q[1];
rz(0.14859622) q[1];
rz(-2.3872822) q[3];
sx q[3];
rz(-0.83900634) q[3];
sx q[3];
rz(3.0712002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7118608) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(2.8374953) q[2];
rz(1.7763304) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(-2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51233184) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(3.1410134) q[0];
rz(-1.0821139) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(-2.2023315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88193306) q[0];
sx q[0];
rz(-1.1355917) q[0];
sx q[0];
rz(2.4324904) q[0];
rz(-pi) q[1];
rz(2.0505191) q[2];
sx q[2];
rz(-2.4726598) q[2];
sx q[2];
rz(-2.875653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56783695) q[1];
sx q[1];
rz(-1.7845961) q[1];
sx q[1];
rz(0.9937728) q[1];
rz(-0.22721283) q[3];
sx q[3];
rz(-1.2649396) q[3];
sx q[3];
rz(-1.3671041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1028563) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(-2.8803414) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-0.84457266) q[0];
sx q[0];
rz(-1.427587) q[0];
sx q[0];
rz(-0.20508668) q[0];
rz(-1.0467485) q[1];
sx q[1];
rz(-1.7457242) q[1];
sx q[1];
rz(0.37839016) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16941026) q[0];
sx q[0];
rz(-0.43146389) q[0];
sx q[0];
rz(-1.1110825) q[0];
x q[1];
rz(0.31585898) q[2];
sx q[2];
rz(-0.49260413) q[2];
sx q[2];
rz(-1.5225747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7582963) q[1];
sx q[1];
rz(-1.9196916) q[1];
sx q[1];
rz(-2.1339244) q[1];
rz(-pi) q[2];
rz(0.50726733) q[3];
sx q[3];
rz(-2.3823822) q[3];
sx q[3];
rz(1.4634446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64688524) q[2];
sx q[2];
rz(-1.1581706) q[2];
sx q[2];
rz(2.0873439) q[2];
rz(2.4833637) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(-0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14195104) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(0.64055881) q[0];
rz(2.7236252) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(2.3366065) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55401257) q[0];
sx q[0];
rz(-0.70776429) q[0];
sx q[0];
rz(-3.1057024) q[0];
x q[1];
rz(2.7348324) q[2];
sx q[2];
rz(-2.2972339) q[2];
sx q[2];
rz(3.022829) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1678214) q[1];
sx q[1];
rz(-0.034618363) q[1];
sx q[1];
rz(-0.95914118) q[1];
rz(2.5500357) q[3];
sx q[3];
rz(-1.322896) q[3];
sx q[3];
rz(1.7660559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54887041) q[2];
sx q[2];
rz(-2.1465116) q[2];
sx q[2];
rz(2.3804046) q[2];
rz(2.393764) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(-0.67659155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51650301) q[0];
sx q[0];
rz(-0.28446063) q[0];
sx q[0];
rz(-1.6647343) q[0];
rz(-0.45627108) q[1];
sx q[1];
rz(-1.4197333) q[1];
sx q[1];
rz(-0.87108535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696156) q[0];
sx q[0];
rz(-2.5273558) q[0];
sx q[0];
rz(-0.02158879) q[0];
rz(-pi) q[1];
rz(0.34452166) q[2];
sx q[2];
rz(-1.4248751) q[2];
sx q[2];
rz(-2.5308756) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2312517) q[1];
sx q[1];
rz(-1.7897072) q[1];
sx q[1];
rz(0.33556767) q[1];
x q[2];
rz(1.7026448) q[3];
sx q[3];
rz(-0.42964298) q[3];
sx q[3];
rz(0.94320083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90290922) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(2.5229559) q[2];
rz(0.020261852) q[3];
sx q[3];
rz(-0.94347763) q[3];
sx q[3];
rz(2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063754931) q[0];
sx q[0];
rz(-1.5564593) q[0];
sx q[0];
rz(-0.42386398) q[0];
rz(0.18691143) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(1.6835469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9317212) q[0];
sx q[0];
rz(-3.0227724) q[0];
sx q[0];
rz(3.0221536) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0293268) q[2];
sx q[2];
rz(-2.5800309) q[2];
sx q[2];
rz(2.4379345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9002237) q[1];
sx q[1];
rz(-1.4321064) q[1];
sx q[1];
rz(-1.6226416) q[1];
rz(-1.4814754) q[3];
sx q[3];
rz(-1.2161939) q[3];
sx q[3];
rz(1.3077298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3905048) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(-2.0474153) q[2];
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
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.3222892) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(-1.5018916) q[0];
rz(0.27885258) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(1.0241114) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2447284) q[0];
sx q[0];
rz(-2.6435268) q[0];
sx q[0];
rz(-0.19356541) q[0];
x q[1];
rz(-3.0908998) q[2];
sx q[2];
rz(-1.8452574) q[2];
sx q[2];
rz(1.3433742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0054202) q[1];
sx q[1];
rz(-1.7917624) q[1];
sx q[1];
rz(-1.3653838) q[1];
rz(-pi) q[2];
rz(1.9387705) q[3];
sx q[3];
rz(-0.77897969) q[3];
sx q[3];
rz(-0.72768962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9378822) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(-0.53696519) q[2];
rz(-1.9133866) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(0.29135191) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1163597) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(-0.7069201) q[1];
sx q[1];
rz(-1.128935) q[1];
sx q[1];
rz(3.1405906) q[1];
rz(-1.7617284) q[2];
sx q[2];
rz(-1.9700865) q[2];
sx q[2];
rz(-1.3843591) q[2];
rz(0.97719396) q[3];
sx q[3];
rz(-2.8799812) q[3];
sx q[3];
rz(2.9577586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
