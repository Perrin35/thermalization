OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3202138) q[0];
sx q[0];
rz(-0.76205215) q[0];
sx q[0];
rz(-2.2574545) q[0];
rz(-0.54028264) q[1];
sx q[1];
rz(-1.5958495) q[1];
sx q[1];
rz(-0.63313142) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5954689) q[0];
sx q[0];
rz(-1.4065388) q[0];
sx q[0];
rz(1.4558653) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9146637) q[2];
sx q[2];
rz(-0.49751505) q[2];
sx q[2];
rz(-2.016248) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31189141) q[1];
sx q[1];
rz(-0.30712488) q[1];
sx q[1];
rz(-2.7077371) q[1];
x q[2];
rz(-2.5217275) q[3];
sx q[3];
rz(-2.4976118) q[3];
sx q[3];
rz(-0.61686531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4151609) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(-2.4994182) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(-1.6665529) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.967763) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(-0.94671384) q[0];
rz(0.058454839) q[1];
sx q[1];
rz(-0.36403257) q[1];
sx q[1];
rz(2.2296925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6175323) q[0];
sx q[0];
rz(-2.723411) q[0];
sx q[0];
rz(1.4773744) q[0];
rz(-1.6591658) q[2];
sx q[2];
rz(-2.4659133) q[2];
sx q[2];
rz(0.32127221) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30787779) q[1];
sx q[1];
rz(-2.6243997) q[1];
sx q[1];
rz(-0.50719313) q[1];
rz(-pi) q[2];
rz(0.1095406) q[3];
sx q[3];
rz(-2.6975204) q[3];
sx q[3];
rz(2.9346497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.45137206) q[2];
sx q[2];
rz(-0.25190941) q[2];
sx q[2];
rz(-0.33736649) q[2];
rz(-2.1515324) q[3];
sx q[3];
rz(-1.8127352) q[3];
sx q[3];
rz(1.9199853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5628691) q[0];
sx q[0];
rz(-2.8674485) q[0];
sx q[0];
rz(1.8291871) q[0];
rz(-0.21526543) q[1];
sx q[1];
rz(-1.637641) q[1];
sx q[1];
rz(0.85830918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68033531) q[0];
sx q[0];
rz(-1.3051998) q[0];
sx q[0];
rz(3.0895825) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2300163) q[2];
sx q[2];
rz(-1.7834605) q[2];
sx q[2];
rz(-0.46175428) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2006113) q[1];
sx q[1];
rz(-1.1340144) q[1];
sx q[1];
rz(-1.753356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8859707) q[3];
sx q[3];
rz(-1.2838253) q[3];
sx q[3];
rz(-0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0158318) q[2];
sx q[2];
rz(-0.55700004) q[2];
sx q[2];
rz(-1.0090984) q[2];
rz(-2.5741746) q[3];
sx q[3];
rz(-1.5083454) q[3];
sx q[3];
rz(-1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0030982) q[0];
sx q[0];
rz(-2.4479285) q[0];
sx q[0];
rz(0.19288119) q[0];
rz(2.0774972) q[1];
sx q[1];
rz(-0.27609438) q[1];
sx q[1];
rz(-0.9562794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67502695) q[0];
sx q[0];
rz(-2.1888632) q[0];
sx q[0];
rz(1.3542372) q[0];
rz(-pi) q[1];
rz(0.85727875) q[2];
sx q[2];
rz(-2.3403717) q[2];
sx q[2];
rz(2.4080133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2697105) q[1];
sx q[1];
rz(-1.0607036) q[1];
sx q[1];
rz(-1.2126232) q[1];
x q[2];
rz(-2.8344735) q[3];
sx q[3];
rz(-0.92832091) q[3];
sx q[3];
rz(2.0406587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8320147) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(2.2008956) q[2];
rz(2.7637123) q[3];
sx q[3];
rz(-2.2418719) q[3];
sx q[3];
rz(2.5009724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85510319) q[0];
sx q[0];
rz(-1.6096977) q[0];
sx q[0];
rz(1.3252277) q[0];
rz(-0.27769604) q[1];
sx q[1];
rz(-0.9869906) q[1];
sx q[1];
rz(1.9353386) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0632478) q[0];
sx q[0];
rz(-1.4916911) q[0];
sx q[0];
rz(1.7633171) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5304753) q[2];
sx q[2];
rz(-1.9131994) q[2];
sx q[2];
rz(-0.26027203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5884456) q[1];
sx q[1];
rz(-0.35970769) q[1];
sx q[1];
rz(1.5313696) q[1];
rz(-2.8998781) q[3];
sx q[3];
rz(-1.4069342) q[3];
sx q[3];
rz(3.0047447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18579927) q[2];
sx q[2];
rz(-2.9904371) q[2];
sx q[2];
rz(-1.7929057) q[2];
rz(1.0079314) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0808829) q[0];
sx q[0];
rz(-3.1302629) q[0];
sx q[0];
rz(0.19705535) q[0];
rz(1.9673678) q[1];
sx q[1];
rz(-0.86001992) q[1];
sx q[1];
rz(0.42541447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933612) q[0];
sx q[0];
rz(-1.5299242) q[0];
sx q[0];
rz(-2.6806683) q[0];
rz(-1.7091016) q[2];
sx q[2];
rz(-0.30685654) q[2];
sx q[2];
rz(0.30447659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6007486) q[1];
sx q[1];
rz(-1.9862529) q[1];
sx q[1];
rz(1.284152) q[1];
x q[2];
rz(-0.80003512) q[3];
sx q[3];
rz(-0.87585319) q[3];
sx q[3];
rz(-1.8025229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60874385) q[2];
sx q[2];
rz(-3.0259735) q[2];
sx q[2];
rz(-2.4382584) q[2];
rz(-0.15130791) q[3];
sx q[3];
rz(-1.9688828) q[3];
sx q[3];
rz(-2.5942588) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587826) q[0];
sx q[0];
rz(-2.3334239) q[0];
sx q[0];
rz(0.95635828) q[0];
rz(-0.76902485) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(1.0395315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2295264) q[0];
sx q[0];
rz(-1.5943702) q[0];
sx q[0];
rz(-1.0957155) q[0];
x q[1];
rz(-2.4843487) q[2];
sx q[2];
rz(-2.1604156) q[2];
sx q[2];
rz(0.60962617) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36827792) q[1];
sx q[1];
rz(-1.8551143) q[1];
sx q[1];
rz(-1.8892688) q[1];
x q[2];
rz(-1.9045181) q[3];
sx q[3];
rz(-1.7705373) q[3];
sx q[3];
rz(-0.39272308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9397554) q[2];
sx q[2];
rz(-2.4079683) q[2];
sx q[2];
rz(1.6377829) q[2];
rz(0.87320295) q[3];
sx q[3];
rz(-0.88770294) q[3];
sx q[3];
rz(2.8100815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37209508) q[0];
sx q[0];
rz(-0.3898507) q[0];
sx q[0];
rz(-1.9918342) q[0];
rz(0.86504966) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(1.4335776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2385927) q[0];
sx q[0];
rz(-1.2051788) q[0];
sx q[0];
rz(-2.1793188) q[0];
rz(-pi) q[1];
rz(2.6636567) q[2];
sx q[2];
rz(-1.897287) q[2];
sx q[2];
rz(-0.70604347) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79093248) q[1];
sx q[1];
rz(-2.6106842) q[1];
sx q[1];
rz(-1.7450619) q[1];
rz(-2.0074437) q[3];
sx q[3];
rz(-0.45387196) q[3];
sx q[3];
rz(1.8271104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64234298) q[2];
sx q[2];
rz(-1.4433824) q[2];
sx q[2];
rz(2.3976105) q[2];
rz(-0.83114019) q[3];
sx q[3];
rz(-1.7404514) q[3];
sx q[3];
rz(-2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1019854) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(-2.4264477) q[0];
rz(-1.8483509) q[1];
sx q[1];
rz(-1.5512385) q[1];
sx q[1];
rz(2.901758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47166804) q[0];
sx q[0];
rz(-1.0201766) q[0];
sx q[0];
rz(2.1651833) q[0];
rz(-pi) q[1];
rz(-0.56279601) q[2];
sx q[2];
rz(-1.6220731) q[2];
sx q[2];
rz(-2.4382961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4308987) q[1];
sx q[1];
rz(-1.7393117) q[1];
sx q[1];
rz(0.88697834) q[1];
x q[2];
rz(-1.5375983) q[3];
sx q[3];
rz(-1.0799468) q[3];
sx q[3];
rz(2.6404078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0750601) q[2];
sx q[2];
rz(-1.5873453) q[2];
sx q[2];
rz(-1.6481605) q[2];
rz(-3.1028124) q[3];
sx q[3];
rz(-1.6438234) q[3];
sx q[3];
rz(-1.7356167) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99603242) q[0];
sx q[0];
rz(-1.4936438) q[0];
sx q[0];
rz(-2.5160568) q[0];
rz(-1.372555) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(0.58403429) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2997949) q[0];
sx q[0];
rz(-1.3158176) q[0];
sx q[0];
rz(0.38354446) q[0];
x q[1];
rz(2.4768684) q[2];
sx q[2];
rz(-1.692284) q[2];
sx q[2];
rz(0.45035502) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4702813) q[1];
sx q[1];
rz(-1.352942) q[1];
sx q[1];
rz(-2.4163428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32665334) q[3];
sx q[3];
rz(-0.71200899) q[3];
sx q[3];
rz(-2.6077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0293616) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(1.7870129) q[2];
rz(2.3692047) q[3];
sx q[3];
rz(-1.200501) q[3];
sx q[3];
rz(-1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1188485) q[0];
sx q[0];
rz(-0.72853715) q[0];
sx q[0];
rz(-1.5991612) q[0];
rz(2.6842885) q[1];
sx q[1];
rz(-2.2685469) q[1];
sx q[1];
rz(-0.1874478) q[1];
rz(-1.3832548) q[2];
sx q[2];
rz(-1.6482856) q[2];
sx q[2];
rz(2.0411574) q[2];
rz(-2.0304092) q[3];
sx q[3];
rz(-1.2818973) q[3];
sx q[3];
rz(-3.1218798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
