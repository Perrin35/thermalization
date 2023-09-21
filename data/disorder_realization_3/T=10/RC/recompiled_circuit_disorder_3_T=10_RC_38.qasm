OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(-3.1103818) q[0];
sx q[0];
rz(-2.6565235) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.480455) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(0.35868355) q[0];
x q[1];
rz(1.4235052) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(1.5073843) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5536082) q[1];
sx q[1];
rz(-1.6349287) q[1];
sx q[1];
rz(-1.1448121) q[1];
rz(2.8768086) q[3];
sx q[3];
rz(-1.6695108) q[3];
sx q[3];
rz(-2.7717154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9238613) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-0.031575354) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(2.6020715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3300433) q[0];
sx q[0];
rz(-2.0199611) q[0];
sx q[0];
rz(-0.051785843) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9794481) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(1.0811999) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0961699) q[1];
sx q[1];
rz(-1.3816557) q[1];
sx q[1];
rz(-2.0778836) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5856421) q[3];
sx q[3];
rz(-2.2009938) q[3];
sx q[3];
rz(-0.48660183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(-2.8853436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51107823) q[0];
sx q[0];
rz(-2.179638) q[0];
sx q[0];
rz(1.8947253) q[0];
rz(-pi) q[1];
rz(2.696064) q[2];
sx q[2];
rz(-2.0085213) q[2];
sx q[2];
rz(0.23756269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33830723) q[1];
sx q[1];
rz(-0.76122621) q[1];
sx q[1];
rz(-1.085698) q[1];
x q[2];
rz(-1.6085298) q[3];
sx q[3];
rz(-1.747526) q[3];
sx q[3];
rz(2.7023774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.2333966) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-0.73192275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61011945) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(-0.5831232) q[0];
rz(-0.71728431) q[2];
sx q[2];
rz(-1.4826164) q[2];
sx q[2];
rz(-0.35536534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6432453) q[1];
sx q[1];
rz(-0.32532641) q[1];
sx q[1];
rz(2.494032) q[1];
x q[2];
rz(2.9210806) q[3];
sx q[3];
rz(-1.4482499) q[3];
sx q[3];
rz(-1.8054655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(0.79777065) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20194963) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(-1.4506838) q[0];
x q[1];
rz(0.30349489) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(1.7393302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3623912) q[1];
sx q[1];
rz(-0.63450846) q[1];
sx q[1];
rz(1.4125376) q[1];
rz(2.5161414) q[3];
sx q[3];
rz(-0.52895412) q[3];
sx q[3];
rz(-1.4048502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(1.9942412) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-3.0715122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018774059) q[0];
sx q[0];
rz(-2.9111324) q[0];
sx q[0];
rz(2.3019058) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55307936) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(1.4097708) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3190223) q[1];
sx q[1];
rz(-2.550107) q[1];
sx q[1];
rz(0.2925847) q[1];
x q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(2.7217334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-2.5202259) q[2];
rz(1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(0.85987464) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-0.0079356114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82848362) q[0];
sx q[0];
rz(-2.5223753) q[0];
sx q[0];
rz(-0.45066582) q[0];
x q[1];
rz(-1.0646348) q[2];
sx q[2];
rz(-2.3828265) q[2];
sx q[2];
rz(2.2366692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8994645) q[1];
sx q[1];
rz(-2.2367396) q[1];
sx q[1];
rz(-2.4813586) q[1];
x q[2];
rz(-1.1778529) q[3];
sx q[3];
rz(-1.279139) q[3];
sx q[3];
rz(-1.1100811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(0.41637862) q[2];
rz(-1.3683866) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72909268) q[0];
sx q[0];
rz(-1.3206498) q[0];
sx q[0];
rz(-1.581122) q[0];
x q[1];
rz(-2.6697568) q[2];
sx q[2];
rz(-0.89315692) q[2];
sx q[2];
rz(-3.133528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9598436) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(-0.77397857) q[1];
rz(1.3650465) q[3];
sx q[3];
rz(-1.660941) q[3];
sx q[3];
rz(-2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-1.0650744) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(0.41697821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96796658) q[0];
sx q[0];
rz(-2.2244503) q[0];
sx q[0];
rz(-0.047441479) q[0];
rz(-pi) q[1];
rz(-2.142799) q[2];
sx q[2];
rz(-0.98888328) q[2];
sx q[2];
rz(1.4243766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4638085) q[1];
sx q[1];
rz(-0.52075547) q[1];
sx q[1];
rz(1.2893454) q[1];
x q[2];
rz(-3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(-2.9337692) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.5243994) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(-1.8189925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8947637) q[0];
sx q[0];
rz(-1.449495) q[0];
sx q[0];
rz(-1.9508719) q[0];
x q[1];
rz(2.0595423) q[2];
sx q[2];
rz(-0.49968038) q[2];
sx q[2];
rz(2.3479455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0364914) q[1];
sx q[1];
rz(-1.0579915) q[1];
sx q[1];
rz(-2.3829616) q[1];
rz(1.7694468) q[3];
sx q[3];
rz(-0.66505265) q[3];
sx q[3];
rz(2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(0.017756391) q[2];
sx q[2];
rz(-2.6323071) q[2];
sx q[2];
rz(1.4321362) q[2];
rz(1.3397459) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
