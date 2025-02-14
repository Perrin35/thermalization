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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(-1.0519354) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(-0.63313484) q[1];
sx q[1];
rz(0.57300353) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7318037) q[0];
sx q[0];
rz(-0.62792009) q[0];
sx q[0];
rz(-1.6768006) q[0];
rz(0.51197211) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(1.3775502) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.68202735) q[1];
sx q[1];
rz(-2.6394156) q[1];
sx q[1];
rz(2.9952496) q[1];
rz(-pi) q[2];
rz(-1.6735733) q[3];
sx q[3];
rz(-2.2647018) q[3];
sx q[3];
rz(-2.4569907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8594592) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(-2.0852883) q[2];
rz(3.1193962) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2317155) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(-1.7040303) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24457045) q[0];
sx q[0];
rz(-1.8720227) q[0];
sx q[0];
rz(1.3006849) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.050927863) q[2];
sx q[2];
rz(-1.5136592) q[2];
sx q[2];
rz(-0.47540755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3900657) q[1];
sx q[1];
rz(-1.2428778) q[1];
sx q[1];
rz(-2.4310914) q[1];
rz(-pi) q[2];
rz(-0.35234612) q[3];
sx q[3];
rz(-2.5735613) q[3];
sx q[3];
rz(-1.8527135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.4613287) q[2];
sx q[2];
rz(2.6961668) q[2];
rz(-2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55260783) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(0.71480042) q[0];
rz(2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(-0.32726273) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48178534) q[0];
sx q[0];
rz(-2.5707173) q[0];
sx q[0];
rz(1.8828189) q[0];
rz(-1.0017603) q[2];
sx q[2];
rz(-1.941701) q[2];
sx q[2];
rz(-0.47874622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5048557) q[1];
sx q[1];
rz(-1.550192) q[1];
sx q[1];
rz(-1.4192753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1932262) q[3];
sx q[3];
rz(-2.7223058) q[3];
sx q[3];
rz(0.47800999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(1.6376015) q[2];
rz(-1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0729436) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-2.9856227) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(-1.2329996) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1240163) q[0];
sx q[0];
rz(-1.6917545) q[0];
sx q[0];
rz(1.8328387) q[0];
rz(-pi) q[1];
rz(-2.2188051) q[2];
sx q[2];
rz(-1.7402667) q[2];
sx q[2];
rz(0.10709912) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2501331) q[1];
sx q[1];
rz(-2.1910408) q[1];
sx q[1];
rz(-2.7930082) q[1];
rz(0.057179515) q[3];
sx q[3];
rz(-0.25186497) q[3];
sx q[3];
rz(-0.72309031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4431346) q[2];
sx q[2];
rz(-0.036018697) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.776942) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(-1.1530217) q[0];
rz(2.5979089) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(2.8797454) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9898997) q[0];
sx q[0];
rz(-1.6542871) q[0];
sx q[0];
rz(-0.033647353) q[0];
rz(-2.1543617) q[2];
sx q[2];
rz(-1.8529743) q[2];
sx q[2];
rz(-2.5099177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3186372) q[1];
sx q[1];
rz(-1.793211) q[1];
sx q[1];
rz(0.59808029) q[1];
rz(-0.20528593) q[3];
sx q[3];
rz(-0.79278273) q[3];
sx q[3];
rz(2.1403811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5220945) q[2];
sx q[2];
rz(-2.5358989) q[2];
sx q[2];
rz(-1.4052793) q[2];
rz(-2.3960579) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0010506823) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(2.7897799) q[0];
rz(-2.70128) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(-0.17131677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85799828) q[0];
sx q[0];
rz(-0.71526113) q[0];
sx q[0];
rz(0.37566988) q[0];
rz(1.8753042) q[2];
sx q[2];
rz(-2.300548) q[2];
sx q[2];
rz(-0.27308057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76207668) q[1];
sx q[1];
rz(-2.2691896) q[1];
sx q[1];
rz(1.5435436) q[1];
rz(-pi) q[2];
rz(-2.2619352) q[3];
sx q[3];
rz(-2.1591957) q[3];
sx q[3];
rz(-3.1373051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0773086) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(-0.94179955) q[2];
rz(-1.9866379) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(1.8010767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6963541) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(3.136694) q[0];
rz(2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-2.4536224) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5410091) q[0];
sx q[0];
rz(-2.2760411) q[0];
sx q[0];
rz(-1.3441104) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46779386) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(1.837567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6354949) q[1];
sx q[1];
rz(-0.48759547) q[1];
sx q[1];
rz(-2.0388842) q[1];
x q[2];
rz(0.87422411) q[3];
sx q[3];
rz(-2.1316574) q[3];
sx q[3];
rz(1.107533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61116162) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(-0.92998695) q[2];
rz(1.9200578) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(2.5482224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(0.68194836) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(-0.090106877) q[0];
rz(-1.7168761) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(-2.7065014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97532192) q[0];
sx q[0];
rz(-1.4251801) q[0];
sx q[0];
rz(2.6377489) q[0];
rz(-pi) q[1];
rz(1.7516364) q[2];
sx q[2];
rz(-0.9107843) q[2];
sx q[2];
rz(1.7969799) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96174091) q[1];
sx q[1];
rz(-1.6992178) q[1];
sx q[1];
rz(0.95203103) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7172312) q[3];
sx q[3];
rz(-1.5699982) q[3];
sx q[3];
rz(-2.3319946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6626176) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(2.5153861) q[2];
rz(-2.9446757) q[3];
sx q[3];
rz(-0.99072376) q[3];
sx q[3];
rz(1.3885434) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597252) q[0];
sx q[0];
rz(-1.4511061) q[0];
sx q[0];
rz(-3.0873121) q[0];
rz(2.718603) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(1.3351701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2847808) q[0];
sx q[0];
rz(-2.3216212) q[0];
sx q[0];
rz(0.7818082) q[0];
rz(-pi) q[1];
rz(2.1077431) q[2];
sx q[2];
rz(-0.62990153) q[2];
sx q[2];
rz(-0.065187188) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.685251) q[1];
sx q[1];
rz(-2.0944746) q[1];
sx q[1];
rz(-2.1728553) q[1];
x q[2];
rz(-1.1109675) q[3];
sx q[3];
rz(-2.9879192) q[3];
sx q[3];
rz(-0.75017649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52428952) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(2.7900901) q[2];
rz(-1.3737804) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(-2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3897301) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(0.75916284) q[0];
rz(2.0579386) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(-2.0738475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5006951) q[0];
sx q[0];
rz(-1.0984549) q[0];
sx q[0];
rz(-0.91820902) q[0];
x q[1];
rz(2.3876486) q[2];
sx q[2];
rz(-2.747376) q[2];
sx q[2];
rz(2.8682402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0800277) q[1];
sx q[1];
rz(-2.1127709) q[1];
sx q[1];
rz(-2.6656277) q[1];
rz(-2.2678947) q[3];
sx q[3];
rz(-1.8305959) q[3];
sx q[3];
rz(1.8099305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(-0.21027002) q[2];
rz(-2.353904) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(-0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44337153) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(1.5086077) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(0.40661033) q[2];
sx q[2];
rz(-0.92401531) q[2];
sx q[2];
rz(0.32377908) q[2];
rz(-1.8456712) q[3];
sx q[3];
rz(-1.4209005) q[3];
sx q[3];
rz(3.0358547) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
