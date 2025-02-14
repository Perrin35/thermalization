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
rz(-0.93936062) q[0];
sx q[0];
rz(-2.0225749) q[0];
sx q[0];
rz(0.033666704) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(1.508498) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38070883) q[0];
sx q[0];
rz(-3.0233726) q[0];
sx q[0];
rz(1.3143172) q[0];
x q[1];
rz(1.2766507) q[2];
sx q[2];
rz(-1.4987117) q[2];
sx q[2];
rz(-2.8800137) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.6623792) q[1];
sx q[1];
rz(-2.2594927) q[1];
sx q[1];
rz(-1.7583048) q[1];
x q[2];
rz(-0.90773031) q[3];
sx q[3];
rz(-2.2624348) q[3];
sx q[3];
rz(1.3724821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.17026751) q[2];
sx q[2];
rz(-2.2496932) q[2];
sx q[2];
rz(-2.3049221) q[2];
rz(0.80638805) q[3];
sx q[3];
rz(-0.92525768) q[3];
sx q[3];
rz(0.18536082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760913) q[0];
sx q[0];
rz(-2.0747023) q[0];
sx q[0];
rz(1.2170894) q[0];
rz(2.941046) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(-0.54380551) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018262176) q[0];
sx q[0];
rz(-1.2549434) q[0];
sx q[0];
rz(1.3297775) q[0];
x q[1];
rz(-1.5571652) q[2];
sx q[2];
rz(-2.1657155) q[2];
sx q[2];
rz(0.10228678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4098379) q[1];
sx q[1];
rz(-1.2827164) q[1];
sx q[1];
rz(-1.9695915) q[1];
rz(1.614566) q[3];
sx q[3];
rz(-2.4659116) q[3];
sx q[3];
rz(-2.0213493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4726938) q[2];
sx q[2];
rz(-2.3987179) q[2];
sx q[2];
rz(1.5354068) q[2];
rz(-1.6425447) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(-1.2523874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.0405618) q[0];
sx q[0];
rz(-2.5767548) q[0];
sx q[0];
rz(3.132013) q[0];
rz(-1.0293695) q[1];
sx q[1];
rz(-1.5734438) q[1];
sx q[1];
rz(-1.8052489) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1270008) q[0];
sx q[0];
rz(-2.3107502) q[0];
sx q[0];
rz(-0.71944351) q[0];
x q[1];
rz(2.1161386) q[2];
sx q[2];
rz(-0.57078123) q[2];
sx q[2];
rz(-2.2557543) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50130288) q[1];
sx q[1];
rz(-2.7510018) q[1];
sx q[1];
rz(1.0949446) q[1];
rz(-pi) q[2];
rz(1.1824781) q[3];
sx q[3];
rz(-0.42436436) q[3];
sx q[3];
rz(2.1774815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0336527) q[2];
sx q[2];
rz(-2.630271) q[2];
sx q[2];
rz(-1.8281724) q[2];
rz(3.1016453) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(-2.0136755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(3.0082598) q[0];
sx q[0];
rz(-1.8513716) q[0];
sx q[0];
rz(-2.7449352) q[0];
rz(0.89504129) q[1];
sx q[1];
rz(-1.4245234) q[1];
sx q[1];
rz(-0.76362124) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1930924) q[0];
sx q[0];
rz(-1.4903127) q[0];
sx q[0];
rz(-1.3250687) q[0];
rz(-pi) q[1];
rz(-0.12053804) q[2];
sx q[2];
rz(-1.5876895) q[2];
sx q[2];
rz(-2.8617045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8473007) q[1];
sx q[1];
rz(-2.9025893) q[1];
sx q[1];
rz(-1.4579178) q[1];
x q[2];
rz(1.5337297) q[3];
sx q[3];
rz(-1.8614113) q[3];
sx q[3];
rz(0.99729482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7274196) q[2];
sx q[2];
rz(-1.5199993) q[2];
sx q[2];
rz(1.1302036) q[2];
rz(2.0027509) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(1.1677008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30911699) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(-2.6600237) q[0];
rz(0.32749495) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(-2.17735) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835579) q[0];
sx q[0];
rz(-2.6917771) q[0];
sx q[0];
rz(-1.0490473) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98924182) q[2];
sx q[2];
rz(-1.4099401) q[2];
sx q[2];
rz(3.1015143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2468919) q[1];
sx q[1];
rz(-2.8671226) q[1];
sx q[1];
rz(1.1054705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6752969) q[3];
sx q[3];
rz(-0.43192902) q[3];
sx q[3];
rz(1.6598048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1714736) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(-1.4946651) q[2];
rz(-1.0129048) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(0.51500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2038912) q[0];
sx q[0];
rz(-2.1370115) q[0];
sx q[0];
rz(1.0585693) q[0];
rz(-2.2189498) q[1];
sx q[1];
rz(-2.1244996) q[1];
sx q[1];
rz(0.10291084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424868) q[0];
sx q[0];
rz(-2.1916336) q[0];
sx q[0];
rz(2.5607361) q[0];
rz(-pi) q[1];
rz(-0.40252588) q[2];
sx q[2];
rz(-2.5136097) q[2];
sx q[2];
rz(1.8875853) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6158689) q[1];
sx q[1];
rz(-1.3868887) q[1];
sx q[1];
rz(-2.1083819) q[1];
x q[2];
rz(3.135337) q[3];
sx q[3];
rz(-2.2872074) q[3];
sx q[3];
rz(2.7992918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.70204488) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(0.31583819) q[2];
rz(-0.99509197) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(-2.4163213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0072667) q[0];
sx q[0];
rz(-0.26679978) q[0];
sx q[0];
rz(2.5079978) q[0];
rz(2.7569356) q[1];
sx q[1];
rz(-1.9622012) q[1];
sx q[1];
rz(0.42919174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776183) q[0];
sx q[0];
rz(-2.0380479) q[0];
sx q[0];
rz(0.28470566) q[0];
rz(-3.09893) q[2];
sx q[2];
rz(-2.4555619) q[2];
sx q[2];
rz(-1.3801284) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7197912) q[1];
sx q[1];
rz(-2.4445577) q[1];
sx q[1];
rz(2.5747673) q[1];
x q[2];
rz(1.0121392) q[3];
sx q[3];
rz(-1.1691165) q[3];
sx q[3];
rz(2.5334849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.089243285) q[2];
sx q[2];
rz(-0.77755916) q[2];
sx q[2];
rz(2.8746129) q[2];
rz(-0.071768196) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(1.4192386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778022) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(1.4005533) q[0];
rz(1.1446704) q[1];
sx q[1];
rz(-1.0619699) q[1];
sx q[1];
rz(2.5544419) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3187744) q[0];
sx q[0];
rz(-2.1923034) q[0];
sx q[0];
rz(-1.0448635) q[0];
x q[1];
rz(2.7233974) q[2];
sx q[2];
rz(-1.4640239) q[2];
sx q[2];
rz(-2.622294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0441832) q[1];
sx q[1];
rz(-0.72260586) q[1];
sx q[1];
rz(-0.34383066) q[1];
rz(-pi) q[2];
rz(-0.95119344) q[3];
sx q[3];
rz(-1.8319329) q[3];
sx q[3];
rz(-1.2984683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2862386) q[2];
sx q[2];
rz(-1.9619532) q[2];
sx q[2];
rz(-0.41909763) q[2];
rz(1.6176443) q[3];
sx q[3];
rz(-2.2231299) q[3];
sx q[3];
rz(-1.6379448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986901) q[0];
sx q[0];
rz(-1.9200696) q[0];
sx q[0];
rz(-0.30938095) q[0];
rz(-1.7912553) q[1];
sx q[1];
rz(-2.5534936) q[1];
sx q[1];
rz(-2.5877171) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6375991) q[0];
sx q[0];
rz(-1.9862111) q[0];
sx q[0];
rz(-1.9204155) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.993169) q[2];
sx q[2];
rz(-1.191105) q[2];
sx q[2];
rz(2.6167999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9123037) q[1];
sx q[1];
rz(-2.4472651) q[1];
sx q[1];
rz(1.0780543) q[1];
rz(-pi) q[2];
rz(-0.70205261) q[3];
sx q[3];
rz(-1.3300196) q[3];
sx q[3];
rz(2.05621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89173335) q[2];
sx q[2];
rz(-1.6102108) q[2];
sx q[2];
rz(1.3472793) q[2];
rz(-2.3740718) q[3];
sx q[3];
rz(-1.1859505) q[3];
sx q[3];
rz(0.84848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42221853) q[0];
sx q[0];
rz(-0.30000559) q[0];
sx q[0];
rz(-1.2147709) q[0];
rz(2.1396554) q[1];
sx q[1];
rz(-1.5129713) q[1];
sx q[1];
rz(2.217206) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8633562) q[0];
sx q[0];
rz(-0.72310142) q[0];
sx q[0];
rz(-0.96295607) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1154249) q[2];
sx q[2];
rz(-1.5942425) q[2];
sx q[2];
rz(-3.0128827) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1781285) q[1];
sx q[1];
rz(-0.22451065) q[1];
sx q[1];
rz(0.1046532) q[1];
x q[2];
rz(-0.59238418) q[3];
sx q[3];
rz(-2.1796556) q[3];
sx q[3];
rz(-2.1453573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7944916) q[2];
sx q[2];
rz(-0.93239409) q[2];
sx q[2];
rz(-2.7744228) q[2];
rz(1.1651039) q[3];
sx q[3];
rz(-0.96787435) q[3];
sx q[3];
rz(-0.86618209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544871) q[0];
sx q[0];
rz(-2.4024873) q[0];
sx q[0];
rz(-3.0229229) q[0];
rz(-2.8754996) q[1];
sx q[1];
rz(-0.91231822) q[1];
sx q[1];
rz(-1.4516861) q[1];
rz(-1.8438392) q[2];
sx q[2];
rz(-1.8476877) q[2];
sx q[2];
rz(-2.3503204) q[2];
rz(-2.9333276) q[3];
sx q[3];
rz(-0.97934813) q[3];
sx q[3];
rz(-1.4772268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
