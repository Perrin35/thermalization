OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(0.81737104) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(1.9411545) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5759597) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(-0.54310449) q[0];
rz(-1.3517411) q[2];
sx q[2];
rz(-2.5918505) q[2];
sx q[2];
rz(-1.654939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8564954) q[1];
sx q[1];
rz(-0.76438099) q[1];
sx q[1];
rz(-0.38715036) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2458385) q[3];
sx q[3];
rz(-2.339139) q[3];
sx q[3];
rz(1.5047531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(-1.9460829) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(2.6057459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8230096) q[0];
sx q[0];
rz(-2.5479925) q[0];
sx q[0];
rz(-1.3472605) q[0];
x q[1];
rz(1.3185805) q[2];
sx q[2];
rz(-2.2644342) q[2];
sx q[2];
rz(-2.0351621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1670926) q[1];
sx q[1];
rz(-1.8247316) q[1];
sx q[1];
rz(-2.7212935) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4278533) q[3];
sx q[3];
rz(-1.5385475) q[3];
sx q[3];
rz(1.8247719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0216996) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(2.9911175) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(0.025807468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3993527) q[0];
sx q[0];
rz(-0.25063801) q[0];
sx q[0];
rz(-1.4901194) q[0];
x q[1];
rz(1.8823207) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(0.66750079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64297134) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(-0.69795124) q[1];
rz(-pi) q[2];
rz(1.953655) q[3];
sx q[3];
rz(-1.448505) q[3];
sx q[3];
rz(2.0910636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0187443) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(2.5562111) q[2];
rz(-2.9600926) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.90081763) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(-0.17661072) q[0];
rz(-0.88090849) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49318424) q[0];
sx q[0];
rz(-2.3097561) q[0];
sx q[0];
rz(-2.1302845) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6279814) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(-1.1255217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9163461) q[1];
sx q[1];
rz(-0.9496453) q[1];
sx q[1];
rz(2.0398554) q[1];
rz(2.7809308) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(-3.0879471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(1.0446576) q[2];
rz(-0.70703834) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(2.7846591) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(1.6479187) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(3.0984745) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2152104) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(2.9125288) q[0];
rz(-pi) q[1];
rz(-2.8565065) q[2];
sx q[2];
rz(-2.0721772) q[2];
sx q[2];
rz(1.3146871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9777898) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(2.8603641) q[1];
rz(-pi) q[2];
rz(-0.14857265) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(-0.038392301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(-0.35219231) q[2];
rz(-2.5514065) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-2.0828784) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(2.7013742) q[0];
rz(-pi) q[1];
rz(-0.58745678) q[2];
sx q[2];
rz(-0.60411462) q[2];
sx q[2];
rz(0.47746745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2211654) q[1];
sx q[1];
rz(-1.3472918) q[1];
sx q[1];
rz(1.2616874) q[1];
rz(-pi) q[2];
rz(-0.22535725) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(-0.55707896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(-1.1550711) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-3.0999775) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.6301427) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(2.2999433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.323303) q[0];
sx q[0];
rz(-2.0879732) q[0];
sx q[0];
rz(-1.8316168) q[0];
rz(-1.3782578) q[2];
sx q[2];
rz(-0.95270573) q[2];
sx q[2];
rz(-2.3877909) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6517666) q[1];
sx q[1];
rz(-0.090432743) q[1];
sx q[1];
rz(2.138278) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82960415) q[3];
sx q[3];
rz(-1.4083574) q[3];
sx q[3];
rz(-0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(-0.049953071) q[2];
rz(-0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(-2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71224371) q[0];
sx q[0];
rz(-1.1375543) q[0];
sx q[0];
rz(-0.54746763) q[0];
x q[1];
rz(-1.8629486) q[2];
sx q[2];
rz(-2.7610965) q[2];
sx q[2];
rz(-2.4351062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.046125267) q[1];
sx q[1];
rz(-1.0091262) q[1];
sx q[1];
rz(-2.5820877) q[1];
rz(-pi) q[2];
rz(-0.91248625) q[3];
sx q[3];
rz(-0.45966002) q[3];
sx q[3];
rz(1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(2.0987089) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-2.1946857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3627975) q[0];
sx q[0];
rz(-2.6914094) q[0];
sx q[0];
rz(-2.3953526) q[0];
rz(-3.0445381) q[2];
sx q[2];
rz(-1.6104873) q[2];
sx q[2];
rz(-1.5511712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35518256) q[1];
sx q[1];
rz(-1.9043515) q[1];
sx q[1];
rz(1.4817609) q[1];
rz(-pi) q[2];
rz(-1.1831207) q[3];
sx q[3];
rz(-2.6689853) q[3];
sx q[3];
rz(1.8567059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(1.8927195) q[2];
rz(2.4272264) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-2.4972829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5376741) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(2.8331579) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8337433) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(0.76542379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58562467) q[1];
sx q[1];
rz(-1.854419) q[1];
sx q[1];
rz(2.7908299) q[1];
rz(-pi) q[2];
rz(0.44390042) q[3];
sx q[3];
rz(-1.1392986) q[3];
sx q[3];
rz(-0.75507009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(0.59990668) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(3.0145666) q[2];
sx q[2];
rz(-2.1913678) q[2];
sx q[2];
rz(0.36507228) q[2];
rz(-1.0944081) q[3];
sx q[3];
rz(-2.4088358) q[3];
sx q[3];
rz(0.81523304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
