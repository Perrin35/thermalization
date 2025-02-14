OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.48489269) q[0];
sx q[0];
rz(-3.0527053) q[0];
sx q[0];
rz(0.77699295) q[0];
rz(-2.4465893) q[1];
sx q[1];
rz(-0.13091317) q[1];
sx q[1];
rz(-0.34832365) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94960864) q[0];
sx q[0];
rz(-0.79553793) q[0];
sx q[0];
rz(2.863671) q[0];
x q[1];
rz(0.25599249) q[2];
sx q[2];
rz(-0.29689327) q[2];
sx q[2];
rz(-0.96742899) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49807032) q[1];
sx q[1];
rz(-0.76979945) q[1];
sx q[1];
rz(-0.78242015) q[1];
rz(-pi) q[2];
rz(-3.0248966) q[3];
sx q[3];
rz(-2.7435997) q[3];
sx q[3];
rz(1.3838878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4942887) q[2];
sx q[2];
rz(-2.5755197) q[2];
sx q[2];
rz(2.0658134) q[2];
rz(0.11476573) q[3];
sx q[3];
rz(-1.6598631) q[3];
sx q[3];
rz(2.4858294) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7455416) q[0];
sx q[0];
rz(-2.2069554) q[0];
sx q[0];
rz(2.4387687) q[0];
rz(-1.5463691) q[1];
sx q[1];
rz(-0.74235761) q[1];
sx q[1];
rz(-2.5667618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0718577) q[0];
sx q[0];
rz(-1.1848579) q[0];
sx q[0];
rz(1.2161847) q[0];
x q[1];
rz(0.6906688) q[2];
sx q[2];
rz(-2.4021742) q[2];
sx q[2];
rz(1.2659466) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65233764) q[1];
sx q[1];
rz(-1.4397419) q[1];
sx q[1];
rz(-0.62789692) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0466667) q[3];
sx q[3];
rz(-2.1364177) q[3];
sx q[3];
rz(-3.0622967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90291643) q[2];
sx q[2];
rz(-1.5006337) q[2];
sx q[2];
rz(-0.73671663) q[2];
rz(2.369407) q[3];
sx q[3];
rz(-2.9530647) q[3];
sx q[3];
rz(-2.7948715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9521028) q[0];
sx q[0];
rz(-2.0661856) q[0];
sx q[0];
rz(1.9814251) q[0];
rz(0.53933764) q[1];
sx q[1];
rz(-0.62349206) q[1];
sx q[1];
rz(0.070405237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6269913) q[0];
sx q[0];
rz(-0.87341269) q[0];
sx q[0];
rz(-0.91430362) q[0];
rz(-1.7504832) q[2];
sx q[2];
rz(-1.7298815) q[2];
sx q[2];
rz(0.87766872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.119795) q[1];
sx q[1];
rz(-2.7547495) q[1];
sx q[1];
rz(-0.17328823) q[1];
rz(2.4910624) q[3];
sx q[3];
rz(-1.1766947) q[3];
sx q[3];
rz(0.89571834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.77014852) q[2];
sx q[2];
rz(-1.4263209) q[2];
sx q[2];
rz(-0.16177978) q[2];
rz(-2.3700304) q[3];
sx q[3];
rz(-0.13650376) q[3];
sx q[3];
rz(-1.0426883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5914087) q[0];
sx q[0];
rz(-2.3593481) q[0];
sx q[0];
rz(2.8915306) q[0];
rz(2.1024044) q[1];
sx q[1];
rz(-2.3999374) q[1];
sx q[1];
rz(-0.80113062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.455832) q[0];
sx q[0];
rz(-1.555611) q[0];
sx q[0];
rz(2.9299829) q[0];
x q[1];
rz(-2.8829128) q[2];
sx q[2];
rz(-1.341267) q[2];
sx q[2];
rz(1.9088285) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.58248478) q[1];
sx q[1];
rz(-1.6402157) q[1];
sx q[1];
rz(-0.51377929) q[1];
x q[2];
rz(1.2053554) q[3];
sx q[3];
rz(-1.1574928) q[3];
sx q[3];
rz(-2.6843021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4074126) q[2];
sx q[2];
rz(-1.696618) q[2];
sx q[2];
rz(2.5511191) q[2];
rz(-0.34481314) q[3];
sx q[3];
rz(-2.7761288) q[3];
sx q[3];
rz(1.4709681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.19584) q[0];
sx q[0];
rz(-1.7029637) q[0];
sx q[0];
rz(1.6666743) q[0];
rz(1.5226115) q[1];
sx q[1];
rz(-1.5432065) q[1];
sx q[1];
rz(-2.7326857) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98981114) q[0];
sx q[0];
rz(-1.5623925) q[0];
sx q[0];
rz(3.0947797) q[0];
rz(2.2719745) q[2];
sx q[2];
rz(-1.2572351) q[2];
sx q[2];
rz(0.28570785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5572714) q[1];
sx q[1];
rz(-1.9964594) q[1];
sx q[1];
rz(-2.0549279) q[1];
x q[2];
rz(-2.4675995) q[3];
sx q[3];
rz(-1.4970652) q[3];
sx q[3];
rz(1.1300095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54395479) q[2];
sx q[2];
rz(-3.0089617) q[2];
sx q[2];
rz(-0.65072101) q[2];
rz(-1.4607653) q[3];
sx q[3];
rz(-1.7284349) q[3];
sx q[3];
rz(2.7248603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.8553218) q[0];
sx q[0];
rz(-1.0690419) q[0];
sx q[0];
rz(-1.2438114) q[0];
rz(2.6142201) q[1];
sx q[1];
rz(-1.6970044) q[1];
sx q[1];
rz(1.7605555) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46845159) q[0];
sx q[0];
rz(-1.8126789) q[0];
sx q[0];
rz(-1.558248) q[0];
x q[1];
rz(-0.63746286) q[2];
sx q[2];
rz(-0.12345498) q[2];
sx q[2];
rz(0.055094624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11070793) q[1];
sx q[1];
rz(-2.0763016) q[1];
sx q[1];
rz(2.5896789) q[1];
rz(-pi) q[2];
rz(-2.2562669) q[3];
sx q[3];
rz(-2.4866085) q[3];
sx q[3];
rz(-1.2745491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1018299) q[2];
sx q[2];
rz(-1.3527752) q[2];
sx q[2];
rz(-1.9977894) q[2];
rz(-0.77491289) q[3];
sx q[3];
rz(-1.9655922) q[3];
sx q[3];
rz(-0.2253783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5483122) q[0];
sx q[0];
rz(-2.6052124) q[0];
sx q[0];
rz(2.8450052) q[0];
rz(2.2071154) q[1];
sx q[1];
rz(-0.88948932) q[1];
sx q[1];
rz(2.8737822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4992384) q[0];
sx q[0];
rz(-1.7846414) q[0];
sx q[0];
rz(1.3539899) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8765062) q[2];
sx q[2];
rz(-0.62787442) q[2];
sx q[2];
rz(1.6123079) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6567475) q[1];
sx q[1];
rz(-2.2094927) q[1];
sx q[1];
rz(1.188422) q[1];
rz(-pi) q[2];
rz(2.9835864) q[3];
sx q[3];
rz(-1.4337162) q[3];
sx q[3];
rz(2.7003845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3725793) q[2];
sx q[2];
rz(-1.0472426) q[2];
sx q[2];
rz(-3.0229342) q[2];
rz(2.3246121) q[3];
sx q[3];
rz(-1.9575565) q[3];
sx q[3];
rz(2.584804) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939222) q[0];
sx q[0];
rz(-0.477328) q[0];
sx q[0];
rz(-1.073786) q[0];
rz(2.1444164) q[1];
sx q[1];
rz(-1.8938096) q[1];
sx q[1];
rz(2.1352077) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9341612) q[0];
sx q[0];
rz(-1.6994358) q[0];
sx q[0];
rz(3.1324803) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9053551) q[2];
sx q[2];
rz(-1.932497) q[2];
sx q[2];
rz(-1.6702056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83676618) q[1];
sx q[1];
rz(-1.8988892) q[1];
sx q[1];
rz(-1.3681218) q[1];
rz(-pi) q[2];
rz(-0.48000042) q[3];
sx q[3];
rz(-2.1566803) q[3];
sx q[3];
rz(2.4898477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5239079) q[2];
sx q[2];
rz(-1.074581) q[2];
sx q[2];
rz(-0.96134821) q[2];
rz(-2.5210467) q[3];
sx q[3];
rz(-2.4397662) q[3];
sx q[3];
rz(-2.3294241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470225) q[0];
sx q[0];
rz(-0.71806878) q[0];
sx q[0];
rz(-2.4156003) q[0];
rz(-0.30966169) q[1];
sx q[1];
rz(-2.1084712) q[1];
sx q[1];
rz(0.18212254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0326906) q[0];
sx q[0];
rz(-1.6918403) q[0];
sx q[0];
rz(0.75184476) q[0];
x q[1];
rz(-0.77960555) q[2];
sx q[2];
rz(-2.4692573) q[2];
sx q[2];
rz(2.5434932) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5096542) q[1];
sx q[1];
rz(-1.7052691) q[1];
sx q[1];
rz(-0.40841497) q[1];
rz(1.9791696) q[3];
sx q[3];
rz(-2.1675936) q[3];
sx q[3];
rz(-0.86270553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2187659) q[2];
sx q[2];
rz(-1.4309692) q[2];
sx q[2];
rz(-1.9975086) q[2];
rz(-0.73793441) q[3];
sx q[3];
rz(-0.72665557) q[3];
sx q[3];
rz(-1.1206333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69796193) q[0];
sx q[0];
rz(-1.8101036) q[0];
sx q[0];
rz(-2.4586082) q[0];
rz(1.5525275) q[1];
sx q[1];
rz(-0.88337675) q[1];
sx q[1];
rz(2.24486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67499298) q[0];
sx q[0];
rz(-1.3764125) q[0];
sx q[0];
rz(-0.43793508) q[0];
x q[1];
rz(-1.0227772) q[2];
sx q[2];
rz(-1.6923133) q[2];
sx q[2];
rz(0.97655988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0005714) q[1];
sx q[1];
rz(-1.440417) q[1];
sx q[1];
rz(1.6597553) q[1];
rz(-pi) q[2];
rz(0.50161485) q[3];
sx q[3];
rz(-1.9700431) q[3];
sx q[3];
rz(0.58540067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54301846) q[2];
sx q[2];
rz(-1.908058) q[2];
sx q[2];
rz(2.7198071) q[2];
rz(-2.5770523) q[3];
sx q[3];
rz(-0.83344236) q[3];
sx q[3];
rz(-2.25901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17557872) q[0];
sx q[0];
rz(-1.6015552) q[0];
sx q[0];
rz(1.8764499) q[0];
rz(1.235678) q[1];
sx q[1];
rz(-2.0912248) q[1];
sx q[1];
rz(0.32774027) q[1];
rz(0.36959481) q[2];
sx q[2];
rz(-1.7679749) q[2];
sx q[2];
rz(1.5883636) q[2];
rz(-0.81461904) q[3];
sx q[3];
rz(-0.64523166) q[3];
sx q[3];
rz(1.393691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
