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
rz(-1.834637) q[0];
sx q[0];
rz(-0.87943465) q[0];
sx q[0];
rz(-0.49361324) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1373089) q[0];
sx q[0];
rz(-0.85544908) q[0];
sx q[0];
rz(-0.8887995) q[0];
rz(-pi) q[1];
rz(-1.2463208) q[2];
sx q[2];
rz(-2.4518584) q[2];
sx q[2];
rz(1.5943499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.263674) q[1];
sx q[1];
rz(-1.3527737) q[1];
sx q[1];
rz(0.71123567) q[1];
rz(-0.36840393) q[3];
sx q[3];
rz(-2.1737636) q[3];
sx q[3];
rz(1.6507698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1349858) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(-1.1084278) q[2];
rz(-1.3125575) q[3];
sx q[3];
rz(-0.24459845) q[3];
sx q[3];
rz(2.826622) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38158622) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(3.1257358) q[0];
rz(0.99041692) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(0.86318618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47967692) q[0];
sx q[0];
rz(-0.74115314) q[0];
sx q[0];
rz(-2.9984498) q[0];
rz(-pi) q[1];
rz(2.8065956) q[2];
sx q[2];
rz(-1.2588698) q[2];
sx q[2];
rz(0.27659135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3962209) q[1];
sx q[1];
rz(-2.6156151) q[1];
sx q[1];
rz(2.7891516) q[1];
x q[2];
rz(-0.14948577) q[3];
sx q[3];
rz(-2.2870697) q[3];
sx q[3];
rz(-3.1325185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70011675) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(-0.46879834) q[2];
rz(2.4539454) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1402635) q[0];
sx q[0];
rz(-2.2183473) q[0];
sx q[0];
rz(2.4053251) q[0];
rz(-0.34321347) q[1];
sx q[1];
rz(-0.88756573) q[1];
sx q[1];
rz(-2.9578178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4030754) q[0];
sx q[0];
rz(-0.74664298) q[0];
sx q[0];
rz(0.097642032) q[0];
rz(-pi) q[1];
rz(-2.0906581) q[2];
sx q[2];
rz(-2.1987872) q[2];
sx q[2];
rz(-1.7597212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0597093) q[1];
sx q[1];
rz(-2.3569111) q[1];
sx q[1];
rz(2.7126524) q[1];
rz(2.7815357) q[3];
sx q[3];
rz(-2.0814133) q[3];
sx q[3];
rz(-1.0086446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.78202128) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(-2.6256631) q[2];
rz(1.1082209) q[3];
sx q[3];
rz(-0.55086946) q[3];
sx q[3];
rz(1.9903323) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16294031) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(2.6287855) q[0];
rz(0.45979744) q[1];
sx q[1];
rz(-1.5922981) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054454) q[0];
sx q[0];
rz(-0.17373611) q[0];
sx q[0];
rz(-0.72189118) q[0];
rz(1.8945208) q[2];
sx q[2];
rz(-1.1737595) q[2];
sx q[2];
rz(0.70003245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7859753) q[1];
sx q[1];
rz(-1.9099059) q[1];
sx q[1];
rz(0.29735844) q[1];
rz(2.7598792) q[3];
sx q[3];
rz(-0.36892051) q[3];
sx q[3];
rz(2.5526508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1226471) q[2];
sx q[2];
rz(-0.15572369) q[2];
sx q[2];
rz(-2.8141008) q[2];
rz(-0.28199768) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9649488) q[0];
sx q[0];
rz(-0.38589859) q[0];
sx q[0];
rz(2.1642245) q[0];
rz(-2.7727959) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(1.4126973) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0127129) q[0];
sx q[0];
rz(-1.0815718) q[0];
sx q[0];
rz(0.21531824) q[0];
rz(-3.1277282) q[2];
sx q[2];
rz(-1.0540451) q[2];
sx q[2];
rz(0.93010974) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1407847) q[1];
sx q[1];
rz(-2.85726) q[1];
sx q[1];
rz(-1.6823014) q[1];
rz(-1.2431954) q[3];
sx q[3];
rz(-2.1246111) q[3];
sx q[3];
rz(-0.49529759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8480924) q[2];
sx q[2];
rz(-0.81587452) q[2];
sx q[2];
rz(-2.9917742) q[2];
rz(2.7430429) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(-0.15628763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77001101) q[0];
sx q[0];
rz(-0.12938975) q[0];
sx q[0];
rz(2.5883801) q[0];
rz(-2.5999056) q[1];
sx q[1];
rz(-2.0952416) q[1];
sx q[1];
rz(-2.6936626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49419379) q[0];
sx q[0];
rz(-1.9564391) q[0];
sx q[0];
rz(1.4146039) q[0];
x q[1];
rz(-2.3382285) q[2];
sx q[2];
rz(-2.5584872) q[2];
sx q[2];
rz(1.2606114) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3975911) q[1];
sx q[1];
rz(-0.41644704) q[1];
sx q[1];
rz(1.2214425) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16334857) q[3];
sx q[3];
rz(-2.8046097) q[3];
sx q[3];
rz(2.8548422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64122671) q[2];
sx q[2];
rz(-1.9772823) q[2];
sx q[2];
rz(0.79916239) q[2];
rz(-0.3847807) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(2.1414115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7810829) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(0.59180301) q[0];
rz(-2.7599755) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(-2.9891678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8417609) q[0];
sx q[0];
rz(-0.61474568) q[0];
sx q[0];
rz(-2.1257945) q[0];
x q[1];
rz(0.32976361) q[2];
sx q[2];
rz(-0.50055365) q[2];
sx q[2];
rz(2.9422174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3445216) q[1];
sx q[1];
rz(-2.7943369) q[1];
sx q[1];
rz(3.0879435) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4983879) q[3];
sx q[3];
rz(-0.99785757) q[3];
sx q[3];
rz(-0.3518897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3064208) q[2];
sx q[2];
rz(-2.1495337) q[2];
sx q[2];
rz(-1.4527808) q[2];
rz(-0.49550223) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.751916) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(2.9801242) q[0];
rz(0.031919315) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(1.8643103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93715765) q[0];
sx q[0];
rz(-1.3124518) q[0];
sx q[0];
rz(-0.72465557) q[0];
rz(-0.23060282) q[2];
sx q[2];
rz(-1.9829033) q[2];
sx q[2];
rz(-1.7526363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.016170597) q[1];
sx q[1];
rz(-1.4628264) q[1];
sx q[1];
rz(0.30534621) q[1];
rz(1.2767918) q[3];
sx q[3];
rz(-1.6075896) q[3];
sx q[3];
rz(2.5599673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6892467) q[2];
sx q[2];
rz(-2.2301058) q[2];
sx q[2];
rz(-0.39608836) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(-0.8766492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24092995) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(2.990429) q[0];
rz(-2.1647272) q[1];
sx q[1];
rz(-0.38115373) q[1];
sx q[1];
rz(2.3968598) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9827288) q[0];
sx q[0];
rz(-1.7773643) q[0];
sx q[0];
rz(1.2875272) q[0];
rz(-pi) q[1];
rz(-1.2538386) q[2];
sx q[2];
rz(-0.43192568) q[2];
sx q[2];
rz(0.077684075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3724461) q[1];
sx q[1];
rz(-2.219104) q[1];
sx q[1];
rz(1.5725971) q[1];
rz(-2.7895097) q[3];
sx q[3];
rz(-1.3439281) q[3];
sx q[3];
rz(2.9405103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5186844) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(2.1935479) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(-1.074056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0483765) q[0];
sx q[0];
rz(-2.8792448) q[0];
sx q[0];
rz(0.23396215) q[0];
rz(-1.9562862) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(-1.2497466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3435182) q[0];
sx q[0];
rz(-1.8600704) q[0];
sx q[0];
rz(-2.5630361) q[0];
rz(1.363867) q[2];
sx q[2];
rz(-0.96038681) q[2];
sx q[2];
rz(1.0640984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9712914) q[1];
sx q[1];
rz(-2.2066322) q[1];
sx q[1];
rz(0.88199349) q[1];
rz(-2.4602182) q[3];
sx q[3];
rz(-1.7077966) q[3];
sx q[3];
rz(-0.43230012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.027111) q[2];
sx q[2];
rz(-1.9805084) q[2];
sx q[2];
rz(1.4410045) q[2];
rz(-0.13127413) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-2.89768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(3.092939) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(-2.1009905) q[1];
sx q[1];
rz(-0.84747172) q[1];
sx q[1];
rz(-0.52451959) q[1];
rz(-2.8362989) q[2];
sx q[2];
rz(-2.1758781) q[2];
sx q[2];
rz(0.97347409) q[2];
rz(2.9374801) q[3];
sx q[3];
rz(-0.64328803) q[3];
sx q[3];
rz(-2.358123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
