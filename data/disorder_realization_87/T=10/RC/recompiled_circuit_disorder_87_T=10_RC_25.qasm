OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694323) q[0];
sx q[0];
rz(-1.1496815) q[0];
sx q[0];
rz(0.62299563) q[0];
rz(-pi) q[1];
x q[1];
rz(2.063077) q[2];
sx q[2];
rz(-2.1247851) q[2];
sx q[2];
rz(-0.89821494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.01552445) q[1];
sx q[1];
rz(-1.4926732) q[1];
sx q[1];
rz(-1.7541581) q[1];
rz(-pi) q[2];
rz(1.9107781) q[3];
sx q[3];
rz(-0.28453207) q[3];
sx q[3];
rz(-1.9576548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9822838) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(0.086159555) q[2];
rz(2.384095) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(1.3264867) q[0];
rz(1.2558698) q[1];
sx q[1];
rz(-1.565226) q[1];
sx q[1];
rz(-2.870141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4685681) q[0];
sx q[0];
rz(-2.4789171) q[0];
sx q[0];
rz(-2.2492692) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2618622) q[2];
sx q[2];
rz(-1.075282) q[2];
sx q[2];
rz(0.77906424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93310624) q[1];
sx q[1];
rz(-2.3098364) q[1];
sx q[1];
rz(-1.0528802) q[1];
rz(-0.35630393) q[3];
sx q[3];
rz(-2.2998527) q[3];
sx q[3];
rz(-0.19405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0469971) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24401027) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(-1.7890731) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(2.9325063) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3269743) q[0];
sx q[0];
rz(-1.0083535) q[0];
sx q[0];
rz(0.66280611) q[0];
rz(-pi) q[1];
rz(1.9833343) q[2];
sx q[2];
rz(-0.87450714) q[2];
sx q[2];
rz(-2.11889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0641388) q[1];
sx q[1];
rz(-1.1541379) q[1];
sx q[1];
rz(-1.0799079) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8385931) q[3];
sx q[3];
rz(-1.3649366) q[3];
sx q[3];
rz(0.96637615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(-2.1787815) q[0];
rz(0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(3.1406291) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92361802) q[0];
sx q[0];
rz(-1.0110564) q[0];
sx q[0];
rz(-2.2773507) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0041299303) q[2];
sx q[2];
rz(-3.0175856) q[2];
sx q[2];
rz(-0.73465243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8648659) q[1];
sx q[1];
rz(-0.58502561) q[1];
sx q[1];
rz(-1.3328972) q[1];
rz(0.56941454) q[3];
sx q[3];
rz(-1.4953519) q[3];
sx q[3];
rz(2.8358592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(-3.0299419) q[2];
rz(-0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(-2.8809663) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97317364) q[0];
sx q[0];
rz(-1.5289613) q[0];
sx q[0];
rz(-2.9996458) q[0];
x q[1];
rz(1.8243276) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(-1.4075116) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.063559859) q[1];
sx q[1];
rz(-1.5884001) q[1];
sx q[1];
rz(2.5672008) q[1];
rz(-pi) q[2];
rz(1.7138249) q[3];
sx q[3];
rz(-1.7149394) q[3];
sx q[3];
rz(3.0343057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9115209) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(1.0391327) q[1];
sx q[1];
rz(-1.8455448) q[1];
sx q[1];
rz(-1.7274436) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6468069) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(1.5216212) q[0];
x q[1];
rz(-0.44712375) q[2];
sx q[2];
rz(-0.56338718) q[2];
sx q[2];
rz(-0.901957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1097088) q[1];
sx q[1];
rz(-2.0307396) q[1];
sx q[1];
rz(2.1214532) q[1];
x q[2];
rz(1.7870951) q[3];
sx q[3];
rz(-1.4970475) q[3];
sx q[3];
rz(-2.4019965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1386537) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(-1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3843) q[0];
sx q[0];
rz(-1.5889656) q[0];
sx q[0];
rz(0.12167715) q[0];
rz(1.1514459) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(2.7391403) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84797317) q[0];
sx q[0];
rz(-1.8403887) q[0];
sx q[0];
rz(-2.2718391) q[0];
rz(3.0955663) q[2];
sx q[2];
rz(-1.5791025) q[2];
sx q[2];
rz(1.9866895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23558815) q[1];
sx q[1];
rz(-1.8530122) q[1];
sx q[1];
rz(1.7547592) q[1];
rz(-3.1137755) q[3];
sx q[3];
rz(-1.9189315) q[3];
sx q[3];
rz(0.58084014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1272614) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(-2.6640653) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(-0.91807085) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(-0.14006242) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(-2.1070811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0089793423) q[0];
sx q[0];
rz(-1.7404557) q[0];
sx q[0];
rz(-1.8011814) q[0];
rz(-1.6413692) q[2];
sx q[2];
rz(-1.6732054) q[2];
sx q[2];
rz(-0.18825738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67147493) q[1];
sx q[1];
rz(-2.5744994) q[1];
sx q[1];
rz(1.6398029) q[1];
rz(-0.23969527) q[3];
sx q[3];
rz(-2.7125159) q[3];
sx q[3];
rz(-2.808188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2404279) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(2.365716) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4416606) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(0.016816703) q[0];
rz(-3.1230714) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-0.7787849) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22395615) q[0];
sx q[0];
rz(-2.7178239) q[0];
sx q[0];
rz(-0.64026041) q[0];
rz(-pi) q[1];
rz(-1.0672827) q[2];
sx q[2];
rz(-0.52769606) q[2];
sx q[2];
rz(-0.50349456) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2494431) q[1];
sx q[1];
rz(-1.3378694) q[1];
sx q[1];
rz(1.1942785) q[1];
rz(0.064568297) q[3];
sx q[3];
rz(-1.7913006) q[3];
sx q[3];
rz(-1.8622423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4497711) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(-1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(-0.17627136) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(2.418628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13847362) q[0];
sx q[0];
rz(-1.4644633) q[0];
sx q[0];
rz(1.7220201) q[0];
x q[1];
rz(1.1935913) q[2];
sx q[2];
rz(-1.7282439) q[2];
sx q[2];
rz(0.27062705) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42230095) q[1];
sx q[1];
rz(-1.0675759) q[1];
sx q[1];
rz(-0.06312381) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78852699) q[3];
sx q[3];
rz(-1.2796113) q[3];
sx q[3];
rz(-2.9332719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(-0.39653683) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(2.6635086) q[2];
sx q[2];
rz(-0.95778428) q[2];
sx q[2];
rz(-0.20245353) q[2];
rz(0.59393926) q[3];
sx q[3];
rz(-2.754302) q[3];
sx q[3];
rz(2.1828628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];