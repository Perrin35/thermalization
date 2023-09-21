OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75738534) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(0.39461179) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73859282) q[2];
sx q[2];
rz(-1.9361708) q[2];
sx q[2];
rz(-1.5942758) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7100916) q[1];
sx q[1];
rz(-0.7845062) q[1];
sx q[1];
rz(1.1169408) q[1];
rz(-pi) q[2];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.1551757) q[3];
sx q[3];
rz(0.88413873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(0.1201771) q[2];
rz(1.1581356) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(2.2226298) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08081089) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(2.2297915) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(2.8149014) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90601901) q[0];
sx q[0];
rz(-2.2740907) q[0];
sx q[0];
rz(-2.2554382) q[0];
rz(1.8961043) q[2];
sx q[2];
rz(-0.46519687) q[2];
sx q[2];
rz(2.543769) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8596526) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(-2.2380026) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2585906) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(1.4955213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(1.3818285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-0.32989311) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1304504) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(1.2304473) q[0];
x q[1];
rz(-1.7705275) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(-1.0945601) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33144618) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(-1.0815094) q[1];
rz(-1.3789165) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.3624297) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9497313) q[0];
sx q[0];
rz(-2.9182069) q[0];
sx q[0];
rz(-2.1641157) q[0];
rz(-1.3894765) q[2];
sx q[2];
rz(-1.7315947) q[2];
sx q[2];
rz(0.60418512) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71388984) q[1];
sx q[1];
rz(-1.0298924) q[1];
sx q[1];
rz(1.1869903) q[1];
rz(-1.1449279) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(-1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797392) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97868279) q[0];
sx q[0];
rz(-1.1497578) q[0];
sx q[0];
rz(-1.0324423) q[0];
rz(-pi) q[1];
rz(0.69552341) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-0.34851375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1272808) q[1];
sx q[1];
rz(-0.31186549) q[1];
sx q[1];
rz(0.77906268) q[1];
rz(-pi) q[2];
rz(0.36851818) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.9011964) q[2];
rz(2.5455348) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(-2.9388359) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-0.99745497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8231682) q[0];
sx q[0];
rz(-0.71821763) q[0];
sx q[0];
rz(-1.7757925) q[0];
x q[1];
rz(-2.4100254) q[2];
sx q[2];
rz(-1.8533857) q[2];
sx q[2];
rz(0.76991316) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5207386) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(2.005307) q[1];
x q[2];
rz(0.84647471) q[3];
sx q[3];
rz(-1.9123565) q[3];
sx q[3];
rz(2.4404756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(2.5278032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261624) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(-1.2994231) q[0];
rz(-0.28366144) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(2.9279857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9871414) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-0.51388545) q[1];
x q[2];
rz(1.3404487) q[3];
sx q[3];
rz(-2.6655572) q[3];
sx q[3];
rz(-1.3766833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(-1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(0.4883782) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280633) q[0];
sx q[0];
rz(-0.54740471) q[0];
sx q[0];
rz(-0.071896032) q[0];
rz(-pi) q[1];
rz(2.7266399) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(-2.6331537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7998432) q[1];
sx q[1];
rz(-1.5768331) q[1];
sx q[1];
rz(2.7762665) q[1];
x q[2];
rz(-1.8833313) q[3];
sx q[3];
rz(-2.0274037) q[3];
sx q[3];
rz(-1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70665923) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(-2.6064176) q[2];
rz(-1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500279) q[0];
sx q[0];
rz(-0.92804337) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029322421) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(-0.14365833) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2885867) q[2];
sx q[2];
rz(-2.2249613) q[2];
sx q[2];
rz(-2.5329563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33854252) q[1];
sx q[1];
rz(-1.1330789) q[1];
sx q[1];
rz(0.16191698) q[1];
x q[2];
rz(1.4521452) q[3];
sx q[3];
rz(-1.3680397) q[3];
sx q[3];
rz(0.53393902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(-1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(1.013247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-2.9985715) q[0];
rz(2.7878694) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(1.313414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4229065) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(-0.42214091) q[1];
x q[2];
rz(-0.89684422) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(-0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6932678) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(-1.7761207) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41291819) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(0.076889597) q[2];
sx q[2];
rz(-2.0145609) q[2];
sx q[2];
rz(-1.2876074) q[2];
rz(-0.87293252) q[3];
sx q[3];
rz(-1.6193661) q[3];
sx q[3];
rz(-1.247874) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
