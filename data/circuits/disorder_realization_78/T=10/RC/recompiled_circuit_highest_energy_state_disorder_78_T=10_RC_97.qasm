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
rz(1.3820833) q[0];
sx q[0];
rz(-0.74892646) q[0];
sx q[0];
rz(1.7480667) q[0];
rz(-2.3387609) q[1];
sx q[1];
rz(-2.3454911) q[1];
sx q[1];
rz(2.610745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40924588) q[0];
sx q[0];
rz(-2.0589925) q[0];
sx q[0];
rz(0.40803473) q[0];
rz(-pi) q[1];
rz(-0.16678236) q[2];
sx q[2];
rz(-1.5537401) q[2];
sx q[2];
rz(0.38558233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.894132) q[1];
sx q[1];
rz(-1.0399721) q[1];
sx q[1];
rz(-3.0546419) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10194709) q[3];
sx q[3];
rz(-2.0876214) q[3];
sx q[3];
rz(-2.068813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73534766) q[2];
sx q[2];
rz(-0.9333868) q[2];
sx q[2];
rz(-0.99785844) q[2];
rz(0.88456279) q[3];
sx q[3];
rz(-1.1793143) q[3];
sx q[3];
rz(0.71200371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3957921) q[0];
sx q[0];
rz(-2.6048248) q[0];
sx q[0];
rz(1.0774379) q[0];
rz(-0.61915818) q[1];
sx q[1];
rz(-1.2666603) q[1];
sx q[1];
rz(1.9177297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9400459) q[0];
sx q[0];
rz(-0.34149656) q[0];
sx q[0];
rz(1.0251371) q[0];
rz(-pi) q[1];
rz(-1.5101132) q[2];
sx q[2];
rz(-3.1113495) q[2];
sx q[2];
rz(-1.1497505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2457285) q[1];
sx q[1];
rz(-2.1460377) q[1];
sx q[1];
rz(2.0644581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7875062) q[3];
sx q[3];
rz(-1.0117784) q[3];
sx q[3];
rz(2.1632568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9549442) q[2];
sx q[2];
rz(-0.84543219) q[2];
sx q[2];
rz(-1.0364944) q[2];
rz(1.1896108) q[3];
sx q[3];
rz(-1.9018973) q[3];
sx q[3];
rz(-3.0724683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4619197) q[0];
sx q[0];
rz(-2.9000977) q[0];
sx q[0];
rz(-1.6756206) q[0];
rz(-1.8849323) q[1];
sx q[1];
rz(-2.4938221) q[1];
sx q[1];
rz(-2.5648436) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.758848) q[0];
sx q[0];
rz(-0.85475105) q[0];
sx q[0];
rz(-2.3780253) q[0];
rz(-1.5336214) q[2];
sx q[2];
rz(-2.0293183) q[2];
sx q[2];
rz(-0.22500817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27611342) q[1];
sx q[1];
rz(-2.7250368) q[1];
sx q[1];
rz(-1.6989442) q[1];
rz(3.0528487) q[3];
sx q[3];
rz(-1.7341494) q[3];
sx q[3];
rz(-3.0609012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3999148) q[2];
sx q[2];
rz(-0.33881131) q[2];
sx q[2];
rz(-0.76876387) q[2];
rz(0.1712884) q[3];
sx q[3];
rz(-1.7007217) q[3];
sx q[3];
rz(-2.7870074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6866368) q[0];
sx q[0];
rz(-2.2577715) q[0];
sx q[0];
rz(-0.53632847) q[0];
rz(0.80279154) q[1];
sx q[1];
rz(-1.2840459) q[1];
sx q[1];
rz(3.0435496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5216828) q[0];
sx q[0];
rz(-0.3371211) q[0];
sx q[0];
rz(-1.8127006) q[0];
rz(-pi) q[1];
rz(2.1613902) q[2];
sx q[2];
rz(-0.88104311) q[2];
sx q[2];
rz(-0.99926567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1817916) q[1];
sx q[1];
rz(-2.0804633) q[1];
sx q[1];
rz(2.9087429) q[1];
rz(-pi) q[2];
rz(-1.1418267) q[3];
sx q[3];
rz(-1.9460856) q[3];
sx q[3];
rz(-0.39794014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1596277) q[2];
sx q[2];
rz(-1.030913) q[2];
sx q[2];
rz(2.8221596) q[2];
rz(2.5751513) q[3];
sx q[3];
rz(-0.74685493) q[3];
sx q[3];
rz(-1.9802861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40120688) q[0];
sx q[0];
rz(-1.8219319) q[0];
sx q[0];
rz(0.55289406) q[0];
rz(-2.3623908) q[1];
sx q[1];
rz(-2.3190277) q[1];
sx q[1];
rz(2.353277) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9181128) q[0];
sx q[0];
rz(-1.8310552) q[0];
sx q[0];
rz(1.6421559) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9328038) q[2];
sx q[2];
rz(-1.5475071) q[2];
sx q[2];
rz(1.67729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8116469) q[1];
sx q[1];
rz(-1.1378985) q[1];
sx q[1];
rz(1.2654773) q[1];
rz(-pi) q[2];
rz(-2.8421229) q[3];
sx q[3];
rz(-1.8078929) q[3];
sx q[3];
rz(2.8194328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7439482) q[2];
sx q[2];
rz(-2.0136191) q[2];
sx q[2];
rz(-0.20565847) q[2];
rz(-1.3501984) q[3];
sx q[3];
rz(-2.4009027) q[3];
sx q[3];
rz(-3.1309483) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3288997) q[0];
sx q[0];
rz(-2.6943272) q[0];
sx q[0];
rz(-1.6269667) q[0];
rz(2.336592) q[1];
sx q[1];
rz(-1.6945508) q[1];
sx q[1];
rz(2.8256493) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8522911) q[0];
sx q[0];
rz(-1.2292432) q[0];
sx q[0];
rz(-0.41595398) q[0];
rz(-pi) q[1];
rz(-0.61152835) q[2];
sx q[2];
rz(-1.2622764) q[2];
sx q[2];
rz(0.53603822) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22623569) q[1];
sx q[1];
rz(-1.8457883) q[1];
sx q[1];
rz(-0.26696856) q[1];
rz(-pi) q[2];
rz(2.6942433) q[3];
sx q[3];
rz(-2.7998689) q[3];
sx q[3];
rz(-1.1319834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67779764) q[2];
sx q[2];
rz(-2.4424876) q[2];
sx q[2];
rz(0.15764906) q[2];
rz(0.53358233) q[3];
sx q[3];
rz(-1.6505417) q[3];
sx q[3];
rz(-1.6238448) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5428298) q[0];
sx q[0];
rz(-2.6451126) q[0];
sx q[0];
rz(2.1606523) q[0];
rz(0.44278231) q[1];
sx q[1];
rz(-1.1863703) q[1];
sx q[1];
rz(2.669899) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169051) q[0];
sx q[0];
rz(-2.0022831) q[0];
sx q[0];
rz(-2.8189075) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4002789) q[2];
sx q[2];
rz(-1.2954752) q[2];
sx q[2];
rz(1.6570651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.065475796) q[1];
sx q[1];
rz(-3.0055176) q[1];
sx q[1];
rz(3.0398366) q[1];
x q[2];
rz(0.34044644) q[3];
sx q[3];
rz(-1.7448815) q[3];
sx q[3];
rz(2.749884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20785759) q[2];
sx q[2];
rz(-2.2722878) q[2];
sx q[2];
rz(3.0233439) q[2];
rz(-2.5737428) q[3];
sx q[3];
rz(-1.7852781) q[3];
sx q[3];
rz(-2.3387199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589013) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(-2.4543104) q[0];
rz(-1.2673238) q[1];
sx q[1];
rz(-1.2143538) q[1];
sx q[1];
rz(-2.5158688) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1835577) q[0];
sx q[0];
rz(-1.5661217) q[0];
sx q[0];
rz(1.5497394) q[0];
x q[1];
rz(2.2082616) q[2];
sx q[2];
rz(-2.0431165) q[2];
sx q[2];
rz(2.460091) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-pi/4) q[1];
sx q[1];
rz(-0.57481261) q[1];
sx q[1];
rz(1.5936127) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3097624) q[3];
sx q[3];
rz(-2.1283669) q[3];
sx q[3];
rz(0.46771177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9132793) q[2];
sx q[2];
rz(-1.8931171) q[2];
sx q[2];
rz(-1.0949562) q[2];
rz(0.45371184) q[3];
sx q[3];
rz(-0.79115051) q[3];
sx q[3];
rz(2.9554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5943282) q[0];
sx q[0];
rz(-2.6481977) q[0];
sx q[0];
rz(0.067597978) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.2815579) q[1];
sx q[1];
rz(-0.13551113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9026176) q[0];
sx q[0];
rz(-2.1186566) q[0];
sx q[0];
rz(0.93642127) q[0];
x q[1];
rz(-0.74112792) q[2];
sx q[2];
rz(-1.4504045) q[2];
sx q[2];
rz(0.84726221) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1780336) q[1];
sx q[1];
rz(-1.8521223) q[1];
sx q[1];
rz(-0.60124361) q[1];
rz(-pi) q[2];
rz(1.816808) q[3];
sx q[3];
rz(-1.8455077) q[3];
sx q[3];
rz(1.6149855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4542666) q[2];
sx q[2];
rz(-0.25298515) q[2];
sx q[2];
rz(2.6542286) q[2];
rz(-0.38960114) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(-0.065936955) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6479263) q[0];
sx q[0];
rz(-1.0856029) q[0];
sx q[0];
rz(1.3847463) q[0];
rz(-2.0866277) q[1];
sx q[1];
rz(-2.2672548) q[1];
sx q[1];
rz(-2.8816282) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9956293) q[0];
sx q[0];
rz(-2.4267247) q[0];
sx q[0];
rz(-1.174182) q[0];
rz(2.2222338) q[2];
sx q[2];
rz(-0.16460379) q[2];
sx q[2];
rz(2.2155264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1321226) q[1];
sx q[1];
rz(-2.1984221) q[1];
sx q[1];
rz(0.043379003) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22180827) q[3];
sx q[3];
rz(-1.3768428) q[3];
sx q[3];
rz(1.8499595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42183033) q[2];
sx q[2];
rz(-0.43785849) q[2];
sx q[2];
rz(1.7913294) q[2];
rz(2.0566025) q[3];
sx q[3];
rz(-1.9271873) q[3];
sx q[3];
rz(-2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0709025) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(-1.6593973) q[1];
sx q[1];
rz(-1.2980325) q[1];
sx q[1];
rz(1.4549805) q[1];
rz(-0.58100159) q[2];
sx q[2];
rz(-1.1270564) q[2];
sx q[2];
rz(-1.8869274) q[2];
rz(-1.0559215) q[3];
sx q[3];
rz(-0.56031651) q[3];
sx q[3];
rz(0.65800695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
