OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(-1.8072577) q[1];
sx q[1];
rz(-0.72307888) q[1];
sx q[1];
rz(-1.877797) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9845147) q[0];
sx q[0];
rz(-2.4369168) q[0];
sx q[0];
rz(1.1873755) q[0];
rz(-0.42277067) q[2];
sx q[2];
rz(-1.0385993) q[2];
sx q[2];
rz(-1.339004) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8883724) q[1];
sx q[1];
rz(-1.651763) q[1];
sx q[1];
rz(-1.2451647) q[1];
rz(-pi) q[2];
rz(1.1321409) q[3];
sx q[3];
rz(-1.7925471) q[3];
sx q[3];
rz(2.7799299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15452142) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(1.7729574) q[2];
rz(-2.9030419) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(0.11428782) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7411165) q[0];
sx q[0];
rz(-2.0023161) q[0];
sx q[0];
rz(-1.1503295) q[0];
rz(0.087609619) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(1.5709343) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8920736) q[0];
sx q[0];
rz(-2.6457835) q[0];
sx q[0];
rz(0.22479381) q[0];
rz(-0.3377487) q[2];
sx q[2];
rz(-2.9523473) q[2];
sx q[2];
rz(0.35795975) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1024061) q[1];
sx q[1];
rz(-0.79823179) q[1];
sx q[1];
rz(-2.185553) q[1];
rz(-pi) q[2];
rz(2.6704795) q[3];
sx q[3];
rz(-2.2918252) q[3];
sx q[3];
rz(1.6278933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7972083) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(2.9551282) q[2];
rz(0.79948419) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(2.1605087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28476533) q[0];
sx q[0];
rz(-1.7600049) q[0];
sx q[0];
rz(-3.0322266) q[0];
rz(-0.60375396) q[1];
sx q[1];
rz(-1.3394638) q[1];
sx q[1];
rz(1.0341136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2952356) q[0];
sx q[0];
rz(-0.24938008) q[0];
sx q[0];
rz(2.1378822) q[0];
x q[1];
rz(1.6414406) q[2];
sx q[2];
rz(-1.6598668) q[2];
sx q[2];
rz(1.8235109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32419606) q[1];
sx q[1];
rz(-1.4651555) q[1];
sx q[1];
rz(-2.8608341) q[1];
rz(-1.4294741) q[3];
sx q[3];
rz(-2.3374448) q[3];
sx q[3];
rz(2.647561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(-0.80101454) q[3];
sx q[3];
rz(-2.2478734) q[3];
sx q[3];
rz(-3.0494704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48701778) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(-0.45904485) q[0];
rz(-1.1162988) q[1];
sx q[1];
rz(-0.79367677) q[1];
sx q[1];
rz(0.8265411) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61454489) q[0];
sx q[0];
rz(-1.2189208) q[0];
sx q[0];
rz(-1.1220758) q[0];
x q[1];
rz(-3.0439348) q[2];
sx q[2];
rz(-1.6689166) q[2];
sx q[2];
rz(-0.94933921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7956808) q[1];
sx q[1];
rz(-1.374578) q[1];
sx q[1];
rz(1.0190585) q[1];
rz(-pi) q[2];
rz(-2.8119254) q[3];
sx q[3];
rz(-2.7894944) q[3];
sx q[3];
rz(-2.6702777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35215846) q[2];
sx q[2];
rz(-2.9106079) q[2];
sx q[2];
rz(2.3918772) q[2];
rz(1.1226783) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(-2.1813755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4753251) q[0];
sx q[0];
rz(-2.1551977) q[0];
sx q[0];
rz(-3.0295897) q[0];
rz(-0.22459596) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(-2.3602643) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3488779) q[0];
sx q[0];
rz(-1.166806) q[0];
sx q[0];
rz(-3.1125665) q[0];
x q[1];
rz(2.8113643) q[2];
sx q[2];
rz(-1.4817186) q[2];
sx q[2];
rz(2.02998) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71415662) q[1];
sx q[1];
rz(-2.2793152) q[1];
sx q[1];
rz(2.5664213) q[1];
rz(-pi) q[2];
rz(2.9797793) q[3];
sx q[3];
rz(-2.1573632) q[3];
sx q[3];
rz(-0.88255054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3096699) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(-2.208948) q[2];
rz(-2.0617088) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(-0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90750736) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(1.750741) q[0];
rz(-2.8098409) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(-1.5914241) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0557779) q[0];
sx q[0];
rz(-1.5356488) q[0];
sx q[0];
rz(-2.2136392) q[0];
rz(-0.45959242) q[2];
sx q[2];
rz(-1.5377279) q[2];
sx q[2];
rz(3.1350373) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4188288) q[1];
sx q[1];
rz(-0.19954844) q[1];
sx q[1];
rz(-0.90252374) q[1];
rz(0.8121536) q[3];
sx q[3];
rz(-2.6366173) q[3];
sx q[3];
rz(0.92588378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8288237) q[2];
sx q[2];
rz(-0.88461107) q[2];
sx q[2];
rz(1.620232) q[2];
rz(-2.17365) q[3];
sx q[3];
rz(-1.5588375) q[3];
sx q[3];
rz(-2.2255285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62436002) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(2.970001) q[0];
rz(2.102237) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(1.4412057) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65214163) q[0];
sx q[0];
rz(-1.891937) q[0];
sx q[0];
rz(1.7016861) q[0];
rz(-2.6256297) q[2];
sx q[2];
rz(-1.668264) q[2];
sx q[2];
rz(-0.24701842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4525057) q[1];
sx q[1];
rz(-1.9259261) q[1];
sx q[1];
rz(1.0853115) q[1];
x q[2];
rz(-1.0073184) q[3];
sx q[3];
rz(-1.451056) q[3];
sx q[3];
rz(0.086825018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8695996) q[2];
sx q[2];
rz(-1.7337493) q[2];
sx q[2];
rz(2.5970411) q[2];
rz(-0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(-2.669529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8498103) q[0];
sx q[0];
rz(-2.6047459) q[0];
sx q[0];
rz(0.52919069) q[0];
rz(-2.5657907) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(2.0226488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75615785) q[0];
sx q[0];
rz(-1.5498112) q[0];
sx q[0];
rz(0.051647112) q[0];
x q[1];
rz(2.697102) q[2];
sx q[2];
rz(-0.40676446) q[2];
sx q[2];
rz(3.0009746) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24209472) q[1];
sx q[1];
rz(-0.78259727) q[1];
sx q[1];
rz(2.461754) q[1];
x q[2];
rz(1.0435853) q[3];
sx q[3];
rz(-1.7526748) q[3];
sx q[3];
rz(1.393569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.240856) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(0.71211234) q[2];
rz(-1.6093048) q[3];
sx q[3];
rz(-0.94523793) q[3];
sx q[3];
rz(-1.7942662) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79494548) q[0];
sx q[0];
rz(-1.1279339) q[0];
sx q[0];
rz(-1.0711063) q[0];
rz(1.2062997) q[1];
sx q[1];
rz(-2.1251528) q[1];
sx q[1];
rz(-1.6546904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8140355) q[0];
sx q[0];
rz(-0.90891713) q[0];
sx q[0];
rz(-1.3332086) q[0];
rz(-pi) q[1];
rz(0.73545154) q[2];
sx q[2];
rz(-2.0273773) q[2];
sx q[2];
rz(-2.9614855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2675954) q[1];
sx q[1];
rz(-1.2800299) q[1];
sx q[1];
rz(1.3895172) q[1];
rz(-2.8603691) q[3];
sx q[3];
rz(-1.1373925) q[3];
sx q[3];
rz(0.12896695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6335166) q[2];
sx q[2];
rz(-2.241892) q[2];
sx q[2];
rz(0.077795204) q[2];
rz(-2.9142694) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(1.0556861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530497) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(-2.7689834) q[0];
rz(0.44652069) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(2.2850697) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.646333) q[0];
sx q[0];
rz(-0.071852751) q[0];
sx q[0];
rz(-0.65793474) q[0];
rz(-pi) q[1];
rz(3.0356016) q[2];
sx q[2];
rz(-1.4376831) q[2];
sx q[2];
rz(0.0078545257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0622934) q[1];
sx q[1];
rz(-0.41626272) q[1];
sx q[1];
rz(1.2757311) q[1];
x q[2];
rz(2.4803745) q[3];
sx q[3];
rz(-0.91877979) q[3];
sx q[3];
rz(-0.87424247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38241688) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(2.5860533) q[2];
rz(2.7555452) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(2.4773795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1356708) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(-2.338943) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(-2.7886709) q[2];
sx q[2];
rz(-1.4581994) q[2];
sx q[2];
rz(-1.1117473) q[2];
rz(1.4925719) q[3];
sx q[3];
rz(-2.409392) q[3];
sx q[3];
rz(0.87421855) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
