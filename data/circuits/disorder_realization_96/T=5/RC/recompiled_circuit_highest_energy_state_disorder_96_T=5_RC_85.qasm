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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(-2.9856292) q[0];
rz(-1.9744385) q[1];
sx q[1];
rz(-0.65216291) q[1];
sx q[1];
rz(-2.6551533) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4733361) q[0];
sx q[0];
rz(-1.0205246) q[0];
sx q[0];
rz(-2.2815198) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2335792) q[2];
sx q[2];
rz(-0.62866271) q[2];
sx q[2];
rz(0.44633055) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.26775441) q[1];
sx q[1];
rz(-2.7044786) q[1];
sx q[1];
rz(-0.74437352) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3144929) q[3];
sx q[3];
rz(-1.4129644) q[3];
sx q[3];
rz(0.99770297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8638986) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(2.7126183) q[2];
rz(-0.046836827) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(2.528791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2752537) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(0.66542768) q[0];
rz(-2.0003419) q[1];
sx q[1];
rz(-1.3803866) q[1];
sx q[1];
rz(-0.64249396) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25099558) q[0];
sx q[0];
rz(-2.3038376) q[0];
sx q[0];
rz(1.6392073) q[0];
x q[1];
rz(2.4464954) q[2];
sx q[2];
rz(-1.3386269) q[2];
sx q[2];
rz(2.1921981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68243945) q[1];
sx q[1];
rz(-0.87284589) q[1];
sx q[1];
rz(1.5345011) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1763686) q[3];
sx q[3];
rz(-2.1049728) q[3];
sx q[3];
rz(-1.4458223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3147754) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(2.5891506) q[2];
rz(-1.6636482) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(2.4655931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.8058572) q[0];
sx q[0];
rz(-2.4215846) q[0];
sx q[0];
rz(-0.82255256) q[0];
rz(-2.3303253) q[1];
sx q[1];
rz(-2.8370116) q[1];
sx q[1];
rz(-1.4623581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050348076) q[0];
sx q[0];
rz(-1.1093265) q[0];
sx q[0];
rz(1.1450923) q[0];
x q[1];
rz(1.64531) q[2];
sx q[2];
rz(-1.1587811) q[2];
sx q[2];
rz(-2.3262784) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87440675) q[1];
sx q[1];
rz(-1.864434) q[1];
sx q[1];
rz(-1.0974357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30841804) q[3];
sx q[3];
rz(-1.2756516) q[3];
sx q[3];
rz(0.67527387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86639577) q[2];
sx q[2];
rz(-2.3123645) q[2];
sx q[2];
rz(-0.59494507) q[2];
rz(-2.7368937) q[3];
sx q[3];
rz(-1.1070822) q[3];
sx q[3];
rz(1.8428724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73612708) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(-0.24293105) q[0];
rz(-2.2566707) q[1];
sx q[1];
rz(-2.4156069) q[1];
sx q[1];
rz(3.1387709) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0003687) q[0];
sx q[0];
rz(-1.1871115) q[0];
sx q[0];
rz(-2.3176718) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1921333) q[2];
sx q[2];
rz(-2.0120125) q[2];
sx q[2];
rz(0.37829933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5127848) q[1];
sx q[1];
rz(-1.6411726) q[1];
sx q[1];
rz(0.87035279) q[1];
rz(2.3798928) q[3];
sx q[3];
rz(-0.77328592) q[3];
sx q[3];
rz(-2.9536216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.60488492) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(-2.4082157) q[2];
rz(2.4865161) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(-0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68375278) q[0];
sx q[0];
rz(-0.34128749) q[0];
sx q[0];
rz(0.14232464) q[0];
rz(-1.7798452) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(0.49837643) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5866668) q[0];
sx q[0];
rz(-0.85882551) q[0];
sx q[0];
rz(-3.0601383) q[0];
x q[1];
rz(2.3827078) q[2];
sx q[2];
rz(-1.6079796) q[2];
sx q[2];
rz(0.55472022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8784762) q[1];
sx q[1];
rz(-0.85678393) q[1];
sx q[1];
rz(1.1787547) q[1];
x q[2];
rz(1.9047391) q[3];
sx q[3];
rz(-2.32226) q[3];
sx q[3];
rz(-2.5317095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0087697) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(0.69457561) q[2];
rz(-2.3092367) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(-0.25025234) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4271127) q[0];
sx q[0];
rz(-0.50943333) q[0];
sx q[0];
rz(-0.66202778) q[0];
rz(1.1854677) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(-1.1659291) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87905721) q[0];
sx q[0];
rz(-1.3692807) q[0];
sx q[0];
rz(-2.461947) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5334305) q[2];
sx q[2];
rz(-1.6972622) q[2];
sx q[2];
rz(-1.6236931) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37167376) q[1];
sx q[1];
rz(-3.0536302) q[1];
sx q[1];
rz(-2.1723353) q[1];
rz(-0.61975584) q[3];
sx q[3];
rz(-1.0527851) q[3];
sx q[3];
rz(0.56700804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7621496) q[2];
sx q[2];
rz(-1.9337485) q[2];
sx q[2];
rz(-2.4318648) q[2];
rz(-0.48192561) q[3];
sx q[3];
rz(-2.66633) q[3];
sx q[3];
rz(-0.019439241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703281) q[0];
sx q[0];
rz(-0.98099357) q[0];
sx q[0];
rz(-1.574466) q[0];
rz(-1.4944685) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(-2.2692197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.746856) q[0];
sx q[0];
rz(-2.7867705) q[0];
sx q[0];
rz(-2.0779993) q[0];
x q[1];
rz(-0.21193223) q[2];
sx q[2];
rz(-2.4080347) q[2];
sx q[2];
rz(-2.9447945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10441594) q[1];
sx q[1];
rz(-3.0311916) q[1];
sx q[1];
rz(-1.8502185) q[1];
rz(-pi) q[2];
rz(-1.5386343) q[3];
sx q[3];
rz(-0.65493203) q[3];
sx q[3];
rz(2.3155568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81277728) q[2];
sx q[2];
rz(-2.2810563) q[2];
sx q[2];
rz(1.1278197) q[2];
rz(0.40337107) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(-2.4283714) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9047852) q[0];
sx q[0];
rz(-1.1431563) q[0];
sx q[0];
rz(1.2971725) q[0];
rz(2.6203268) q[1];
sx q[1];
rz(-1.8788985) q[1];
sx q[1];
rz(-3.0788132) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15069467) q[0];
sx q[0];
rz(-2.6932062) q[0];
sx q[0];
rz(1.4175116) q[0];
rz(2.5812878) q[2];
sx q[2];
rz(-2.7231999) q[2];
sx q[2];
rz(-3.0458618) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52555841) q[1];
sx q[1];
rz(-1.8174606) q[1];
sx q[1];
rz(-2.9933369) q[1];
rz(-0.92446961) q[3];
sx q[3];
rz(-1.1271584) q[3];
sx q[3];
rz(0.2547338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89821833) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(-0.96837366) q[2];
rz(2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-2.8106522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53720713) q[0];
sx q[0];
rz(-0.66050118) q[0];
sx q[0];
rz(0.78395098) q[0];
rz(1.9369269) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(1.3752259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581552) q[0];
sx q[0];
rz(-2.7984507) q[0];
sx q[0];
rz(-2.457042) q[0];
rz(1.398574) q[2];
sx q[2];
rz(-2.9841514) q[2];
sx q[2];
rz(-0.96913183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7837401) q[1];
sx q[1];
rz(-0.87271089) q[1];
sx q[1];
rz(-0.74905209) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1777283) q[3];
sx q[3];
rz(-1.8788218) q[3];
sx q[3];
rz(-2.7975688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84460008) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(2.6958579) q[2];
rz(2.5388057) q[3];
sx q[3];
rz(-1.876839) q[3];
sx q[3];
rz(-0.84356892) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8730901) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(2.6173746) q[0];
rz(1.9408608) q[1];
sx q[1];
rz(-1.8914765) q[1];
sx q[1];
rz(-1.1842309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94468695) q[0];
sx q[0];
rz(-1.9160998) q[0];
sx q[0];
rz(1.7717591) q[0];
rz(-1.3648562) q[2];
sx q[2];
rz(-2.8194234) q[2];
sx q[2];
rz(1.3141156) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.47329119) q[1];
sx q[1];
rz(-1.4379825) q[1];
sx q[1];
rz(1.5831489) q[1];
rz(1.6053725) q[3];
sx q[3];
rz(-0.94217052) q[3];
sx q[3];
rz(0.93116597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21620096) q[2];
sx q[2];
rz(-2.6915458) q[2];
sx q[2];
rz(-1.0518543) q[2];
rz(2.9205186) q[3];
sx q[3];
rz(-1.3007921) q[3];
sx q[3];
rz(-0.93824798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5476407) q[0];
sx q[0];
rz(-1.6407536) q[0];
sx q[0];
rz(-3.092691) q[0];
rz(0.37638695) q[1];
sx q[1];
rz(-0.86224894) q[1];
sx q[1];
rz(-1.1663762) q[1];
rz(3.0920269) q[2];
sx q[2];
rz(-2.5709346) q[2];
sx q[2];
rz(1.554629) q[2];
rz(-2.6003414) q[3];
sx q[3];
rz(-0.79081906) q[3];
sx q[3];
rz(-2.9959903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
