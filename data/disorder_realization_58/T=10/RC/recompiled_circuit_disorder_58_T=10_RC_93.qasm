OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(2.7352754) q[0];
sx q[0];
rz(8.6046594) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2031189) q[0];
sx q[0];
rz(-1.1955368) q[0];
sx q[0];
rz(1.24036) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0933502) q[2];
sx q[2];
rz(-0.89077836) q[2];
sx q[2];
rz(2.8505461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1066061) q[1];
sx q[1];
rz(-2.258746) q[1];
sx q[1];
rz(-0.41253849) q[1];
rz(-pi) q[2];
rz(1.1675695) q[3];
sx q[3];
rz(-1.9864169) q[3];
sx q[3];
rz(2.2574539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(1.9834571) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2355736) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(2.2554382) q[0];
rz(-pi) q[1];
rz(-2.0147444) q[2];
sx q[2];
rz(-1.4269281) q[2];
sx q[2];
rz(2.4614046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5888728) q[1];
sx q[1];
rz(-0.96541222) q[1];
sx q[1];
rz(2.6436716) q[1];
x q[2];
rz(-1.6750402) q[3];
sx q[3];
rz(-2.4511271) q[3];
sx q[3];
rz(-3.1363827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6313173) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(-0.32989311) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-2.0842016) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111423) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(-1.2304473) q[0];
rz(-0.51809394) q[2];
sx q[2];
rz(-2.7535541) q[2];
sx q[2];
rz(-0.54259091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33144618) q[1];
sx q[1];
rz(-0.76008893) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
rz(-0.4226513) q[3];
sx q[3];
rz(-1.7461516) q[3];
sx q[3];
rz(-0.56817504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7968984) q[0];
sx q[0];
rz(-1.3860774) q[0];
sx q[0];
rz(3.0152507) q[0];
x q[1];
rz(0.16343127) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-2.2043259) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.0488102) q[1];
sx q[1];
rz(-0.65199344) q[1];
sx q[1];
rz(0.55744967) q[1];
rz(1.1449279) q[3];
sx q[3];
rz(-2.5440359) q[3];
sx q[3];
rz(-1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39607221) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(-2.0139587) q[3];
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
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-2.9096471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9493128) q[0];
sx q[0];
rz(-0.67040196) q[0];
sx q[0];
rz(-2.288726) q[0];
rz(-pi) q[1];
rz(-1.7115713) q[2];
sx q[2];
rz(-2.2614334) q[2];
sx q[2];
rz(1.3123133) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8307454) q[1];
sx q[1];
rz(-1.7880882) q[1];
sx q[1];
rz(-2.9160935) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5512974) q[3];
sx q[3];
rz(-1.202342) q[3];
sx q[3];
rz(1.3080477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.0971138) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(2.1441377) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(-0.17605619) q[0];
x q[1];
rz(-0.41001292) q[2];
sx q[2];
rz(-2.3668681) q[2];
sx q[2];
rz(-2.0395525) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5207386) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(-2.005307) q[1];
rz(1.0783844) q[3];
sx q[3];
rz(-2.3541854) q[3];
sx q[3];
rz(1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.30615) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928771) q[0];
sx q[0];
rz(-2.2971417) q[0];
sx q[0];
rz(2.8918299) q[0];
rz(2.9250547) q[2];
sx q[2];
rz(-2.8515186) q[2];
sx q[2];
rz(-1.1494344) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9871414) q[1];
sx q[1];
rz(-3.0220991) q[1];
sx q[1];
rz(-2.6277072) q[1];
rz(-pi) q[2];
rz(0.11717637) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43341407) q[2];
sx q[2];
rz(-0.91579473) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(0.60837778) q[3];
sx q[3];
rz(-1.7374246) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1293837) q[0];
sx q[0];
rz(-2.1166271) q[0];
sx q[0];
rz(1.6145541) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41495277) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(-2.6331537) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9148548) q[1];
sx q[1];
rz(-1.2054772) q[1];
sx q[1];
rz(1.564333) q[1];
rz(-2.6650388) q[3];
sx q[3];
rz(-1.8503975) q[3];
sx q[3];
rz(2.9001146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-2.6064176) q[2];
rz(-1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(2.1323269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(0.92591441) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029322421) q[0];
sx q[0];
rz(-1.3591213) q[0];
sx q[0];
rz(0.14365833) q[0];
rz(-2.2885867) q[2];
sx q[2];
rz(-0.91663137) q[2];
sx q[2];
rz(-0.6086364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3013819) q[1];
sx q[1];
rz(-1.4242607) q[1];
sx q[1];
rz(1.128003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6191606) q[3];
sx q[3];
rz(-2.9070832) q[3];
sx q[3];
rz(1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-2.6055028) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(-1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(-1.013247) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.917406) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(-1.3361206) q[0];
rz(0.35372325) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(-1.8281787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71868616) q[1];
sx q[1];
rz(-1.6675073) q[1];
sx q[1];
rz(2.7194517) q[1];
rz(-pi) q[2];
rz(2.3771044) q[3];
sx q[3];
rz(-1.0478684) q[3];
sx q[3];
rz(-2.5564699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(2.3821793) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41291819) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-0.57327523) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(0.076889597) q[2];
sx q[2];
rz(-2.0145609) q[2];
sx q[2];
rz(-1.2876074) q[2];
rz(-1.4952954) q[3];
sx q[3];
rz(-2.4423238) q[3];
sx q[3];
rz(-2.8764976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
