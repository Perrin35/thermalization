OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0030274) q[0];
sx q[0];
rz(-2.2632817) q[0];
sx q[0];
rz(0.83100975) q[0];
rz(-0.57549685) q[1];
sx q[1];
rz(3.8776445) q[1];
sx q[1];
rz(10.164645) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4736126) q[0];
sx q[0];
rz(-1.4987317) q[0];
sx q[0];
rz(-0.18911171) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14278463) q[2];
sx q[2];
rz(-2.1342284) q[2];
sx q[2];
rz(0.44445064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1581067) q[1];
sx q[1];
rz(-2.0465104) q[1];
sx q[1];
rz(2.9768894) q[1];
x q[2];
rz(1.3786198) q[3];
sx q[3];
rz(-1.9488244) q[3];
sx q[3];
rz(2.4227633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9326707) q[2];
sx q[2];
rz(-0.17503665) q[2];
sx q[2];
rz(2.8026061) q[2];
rz(2.637376) q[3];
sx q[3];
rz(-1.0176858) q[3];
sx q[3];
rz(-1.1241815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9826688) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(2.6320631) q[0];
rz(-1.6391899) q[1];
sx q[1];
rz(-0.27703151) q[1];
sx q[1];
rz(0.94430077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664638) q[0];
sx q[0];
rz(-0.24953609) q[0];
sx q[0];
rz(0.93231045) q[0];
rz(-pi) q[1];
rz(0.048205094) q[2];
sx q[2];
rz(-1.71016) q[2];
sx q[2];
rz(-0.73783079) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2067926) q[1];
sx q[1];
rz(-2.1046608) q[1];
sx q[1];
rz(-0.67558788) q[1];
rz(-pi) q[2];
rz(1.9924762) q[3];
sx q[3];
rz(-2.4184368) q[3];
sx q[3];
rz(-0.92484091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.79537359) q[2];
sx q[2];
rz(-1.5920762) q[2];
sx q[2];
rz(2.6056371) q[2];
rz(-2.0761944) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(-0.65142256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1123493) q[0];
sx q[0];
rz(-1.1493244) q[0];
sx q[0];
rz(-0.73597062) q[0];
rz(0.53572267) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(-0.89964286) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0871219) q[0];
sx q[0];
rz(-2.961714) q[0];
sx q[0];
rz(0.69152559) q[0];
x q[1];
rz(2.1754335) q[2];
sx q[2];
rz(-1.7636429) q[2];
sx q[2];
rz(1.9230587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6361743) q[1];
sx q[1];
rz(-0.50105011) q[1];
sx q[1];
rz(1.4062642) q[1];
rz(0.24921649) q[3];
sx q[3];
rz(-1.2050516) q[3];
sx q[3];
rz(-0.62376991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6003517) q[2];
sx q[2];
rz(-1.9466126) q[2];
sx q[2];
rz(0.80292732) q[2];
rz(2.3927355) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(-0.15933855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.7298537) q[0];
sx q[0];
rz(-3.0112322) q[0];
rz(-1.8761926) q[1];
sx q[1];
rz(-1.9644901) q[1];
sx q[1];
rz(0.1098384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43183655) q[0];
sx q[0];
rz(-1.1188112) q[0];
sx q[0];
rz(2.7664037) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5750138) q[2];
sx q[2];
rz(-1.8061122) q[2];
sx q[2];
rz(-1.7641774) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3567645) q[1];
sx q[1];
rz(-0.82883976) q[1];
sx q[1];
rz(2.1020562) q[1];
rz(2.5205344) q[3];
sx q[3];
rz(-0.73535669) q[3];
sx q[3];
rz(-1.3842954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.639223) q[2];
sx q[2];
rz(-1.9618192) q[2];
sx q[2];
rz(0.42312527) q[2];
rz(-2.1028178) q[3];
sx q[3];
rz(-0.3813425) q[3];
sx q[3];
rz(1.0264621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1410685) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(0.27960676) q[0];
rz(0.33310834) q[1];
sx q[1];
rz(-0.50222841) q[1];
sx q[1];
rz(0.28848973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56474876) q[0];
sx q[0];
rz(-1.1001407) q[0];
sx q[0];
rz(3.070757) q[0];
rz(1.9847068) q[2];
sx q[2];
rz(-2.8749646) q[2];
sx q[2];
rz(-1.8265343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0049952) q[1];
sx q[1];
rz(-1.4815287) q[1];
sx q[1];
rz(-3.0725293) q[1];
rz(-2.1937859) q[3];
sx q[3];
rz(-1.3419328) q[3];
sx q[3];
rz(-0.2589489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.592411) q[2];
sx q[2];
rz(-1.2191399) q[2];
sx q[2];
rz(-2.9782817) q[2];
rz(1.9773989) q[3];
sx q[3];
rz(-2.4653698) q[3];
sx q[3];
rz(-1.7827079) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0613681) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(-1.9153216) q[0];
rz(-1.3310883) q[1];
sx q[1];
rz(-2.2115579) q[1];
sx q[1];
rz(-0.61940449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7295655) q[0];
sx q[0];
rz(-2.4270227) q[0];
sx q[0];
rz(-2.3518635) q[0];
rz(0.2529958) q[2];
sx q[2];
rz(-1.7645451) q[2];
sx q[2];
rz(0.27744833) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0521026) q[1];
sx q[1];
rz(-1.5137496) q[1];
sx q[1];
rz(-1.4228348) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9483637) q[3];
sx q[3];
rz(-1.7488297) q[3];
sx q[3];
rz(-1.61249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9161393) q[2];
sx q[2];
rz(-1.658354) q[2];
sx q[2];
rz(0.43530604) q[2];
rz(0.12886038) q[3];
sx q[3];
rz(-1.83056) q[3];
sx q[3];
rz(-0.35663566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4211166) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(-0.55225736) q[0];
rz(-2.5241959) q[1];
sx q[1];
rz(-1.2115024) q[1];
sx q[1];
rz(-1.5026106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5867859) q[0];
sx q[0];
rz(-1.1830336) q[0];
sx q[0];
rz(-3.094502) q[0];
rz(1.065997) q[2];
sx q[2];
rz(-1.6994972) q[2];
sx q[2];
rz(1.7247049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0459111) q[1];
sx q[1];
rz(-2.0740182) q[1];
sx q[1];
rz(2.9614425) q[1];
x q[2];
rz(1.5911359) q[3];
sx q[3];
rz(-1.830066) q[3];
sx q[3];
rz(1.532497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5804533) q[2];
sx q[2];
rz(-2.9353607) q[2];
sx q[2];
rz(2.0533766) q[2];
rz(-3.1214118) q[3];
sx q[3];
rz(-1.8686649) q[3];
sx q[3];
rz(-1.5599686) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35767558) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(2.1141323) q[0];
rz(-2.0383535) q[1];
sx q[1];
rz(-1.5682861) q[1];
sx q[1];
rz(-1.9267513) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9068844) q[0];
sx q[0];
rz(-0.32433332) q[0];
sx q[0];
rz(-1.2417481) q[0];
x q[1];
rz(-2.0452477) q[2];
sx q[2];
rz(-1.9812968) q[2];
sx q[2];
rz(0.83555789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.852927) q[1];
sx q[1];
rz(-0.92354362) q[1];
sx q[1];
rz(-2.701328) q[1];
x q[2];
rz(-2.4818871) q[3];
sx q[3];
rz(-0.72835975) q[3];
sx q[3];
rz(-1.212478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2601629) q[2];
sx q[2];
rz(-1.2036999) q[2];
sx q[2];
rz(2.0737958) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-0.46025899) q[3];
sx q[3];
rz(2.0021745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282309) q[0];
sx q[0];
rz(-0.75465337) q[0];
sx q[0];
rz(-0.2555787) q[0];
rz(-2.3981587) q[1];
sx q[1];
rz(-1.5579222) q[1];
sx q[1];
rz(1.5240634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34786087) q[0];
sx q[0];
rz(-1.1208236) q[0];
sx q[0];
rz(1.0211358) q[0];
rz(-pi) q[1];
rz(-0.41478283) q[2];
sx q[2];
rz(-0.66294248) q[2];
sx q[2];
rz(0.65185968) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0164106) q[1];
sx q[1];
rz(-2.4577854) q[1];
sx q[1];
rz(-1.7829624) q[1];
rz(-pi) q[2];
rz(0.26542191) q[3];
sx q[3];
rz(-2.1947104) q[3];
sx q[3];
rz(-1.3359631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0301547) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(1.6072404) q[2];
rz(-1.0715019) q[3];
sx q[3];
rz(-1.5415618) q[3];
sx q[3];
rz(-1.967954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.8906422) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(2.6232134) q[0];
rz(-2.3333343) q[1];
sx q[1];
rz(-1.6812485) q[1];
sx q[1];
rz(1.4498651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0302561) q[0];
sx q[0];
rz(-1.5606631) q[0];
sx q[0];
rz(1.8596605) q[0];
x q[1];
rz(-2.6053455) q[2];
sx q[2];
rz(-1.9384346) q[2];
sx q[2];
rz(0.73208955) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5603578) q[1];
sx q[1];
rz(-2.5311845) q[1];
sx q[1];
rz(-0.44652744) q[1];
x q[2];
rz(0.57761044) q[3];
sx q[3];
rz(-2.0311714) q[3];
sx q[3];
rz(-2.1652997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.9064648) q[2];
sx q[2];
rz(-0.9355363) q[2];
rz(1.6701291) q[3];
sx q[3];
rz(-1.4751438) q[3];
sx q[3];
rz(-1.0940301) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4569296) q[0];
sx q[0];
rz(-2.5419432) q[0];
sx q[0];
rz(2.3210617) q[0];
rz(-0.44878557) q[1];
sx q[1];
rz(-1.1460591) q[1];
sx q[1];
rz(-1.9312327) q[1];
rz(1.6937428) q[2];
sx q[2];
rz(-2.7715383) q[2];
sx q[2];
rz(-0.72889974) q[2];
rz(-2.5815677) q[3];
sx q[3];
rz(-1.8631794) q[3];
sx q[3];
rz(1.2003492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
