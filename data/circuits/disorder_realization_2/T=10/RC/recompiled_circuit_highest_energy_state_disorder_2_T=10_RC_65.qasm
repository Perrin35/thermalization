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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
rz(-2.9479041) q[1];
sx q[1];
rz(-2.5084578) q[1];
sx q[1];
rz(2.5685891) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.409789) q[0];
sx q[0];
rz(-0.62792009) q[0];
sx q[0];
rz(1.464792) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51197211) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(-1.7640424) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51542789) q[1];
sx q[1];
rz(-2.0671131) q[1];
sx q[1];
rz(-1.4908916) q[1];
rz(-pi) q[2];
rz(-1.4680193) q[3];
sx q[3];
rz(-2.2647018) q[3];
sx q[3];
rz(2.4569907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28213349) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(-1.0563043) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(1.4200042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9098772) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(-2.7114482) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(1.7040303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50629726) q[0];
sx q[0];
rz(-0.40182913) q[0];
sx q[0];
rz(0.70962972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84354894) q[2];
sx q[2];
rz(-3.0650716) q[2];
sx q[2];
rz(2.8882972) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75152698) q[1];
sx q[1];
rz(-1.2428778) q[1];
sx q[1];
rz(0.71050127) q[1];
rz(-pi) q[2];
rz(2.6019179) q[3];
sx q[3];
rz(-1.3840578) q[3];
sx q[3];
rz(0.58247551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(-2.6961668) q[2];
rz(2.7323501) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(-0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889848) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(-0.71480042) q[0];
rz(2.0843166) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(0.32726273) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11571685) q[0];
sx q[0];
rz(-2.1109606) q[0];
sx q[0];
rz(2.9469423) q[0];
x q[1];
rz(-2.1398323) q[2];
sx q[2];
rz(-1.941701) q[2];
sx q[2];
rz(0.47874622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79993805) q[1];
sx q[1];
rz(-0.1529049) q[1];
sx q[1];
rz(1.4351109) q[1];
rz(-0.25423519) q[3];
sx q[3];
rz(-1.2336858) q[3];
sx q[3];
rz(1.997662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(1.5039911) q[2];
rz(1.2184294) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(-2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0729436) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(-2.9856227) q[0];
rz(-0.34863696) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.9085931) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6207244) q[0];
sx q[0];
rz(-1.8308795) q[0];
sx q[0];
rz(-0.12518945) q[0];
rz(-2.9301823) q[2];
sx q[2];
rz(-0.93358913) q[2];
sx q[2];
rz(1.5906972) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2501331) q[1];
sx q[1];
rz(-0.95055184) q[1];
sx q[1];
rz(-2.7930082) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8901222) q[3];
sx q[3];
rz(-1.5850388) q[3];
sx q[3];
rz(0.79232681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4431346) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(1.6301463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(-1.9885709) q[0];
rz(0.54368377) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(-0.26184729) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1516929) q[0];
sx q[0];
rz(-1.4873056) q[0];
sx q[0];
rz(-0.033647353) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8072281) q[2];
sx q[2];
rz(-1.0131256) q[2];
sx q[2];
rz(-1.1209436) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56526977) q[1];
sx q[1];
rz(-0.63332958) q[1];
sx q[1];
rz(0.38193462) q[1];
rz(-2.35942) q[3];
sx q[3];
rz(-1.7165136) q[3];
sx q[3];
rz(-2.4268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61949817) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(1.4052793) q[2];
rz(0.74553472) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(0.94329992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.140542) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(-2.7897799) q[0];
rz(-0.44031269) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(2.9702759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.139411) q[0];
sx q[0];
rz(-1.8137964) q[0];
sx q[0];
rz(0.67968207) q[0];
x q[1];
rz(-0.32342685) q[2];
sx q[2];
rz(-2.3617871) q[2];
sx q[2];
rz(0.16743539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.379516) q[1];
sx q[1];
rz(-0.87240309) q[1];
sx q[1];
rz(1.5435436) q[1];
rz(-pi) q[2];
rz(2.3790509) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(0.97555893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.064284023) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(0.94179955) q[2];
rz(-1.1549548) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(-1.340516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(-3.136694) q[0];
rz(-2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(2.4536224) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60058355) q[0];
sx q[0];
rz(-0.86555153) q[0];
sx q[0];
rz(1.3441104) q[0];
x q[1];
rz(0.46779386) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(1.3040257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35536218) q[1];
sx q[1];
rz(-1.3578102) q[1];
sx q[1];
rz(-2.0128241) q[1];
rz(-pi) q[2];
rz(0.68617188) q[3];
sx q[3];
rz(-2.1451575) q[3];
sx q[3];
rz(2.2597093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.61116162) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(0.92998695) q[2];
rz(-1.2215349) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596443) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(-3.0514858) q[0];
rz(1.4247165) q[1];
sx q[1];
rz(-2.5017891) q[1];
sx q[1];
rz(-0.43509126) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6259436) q[0];
sx q[0];
rz(-1.0727779) q[0];
sx q[0];
rz(-1.4048715) q[0];
rz(-pi) q[1];
rz(2.9138759) q[2];
sx q[2];
rz(-0.68074742) q[2];
sx q[2];
rz(-1.6344223) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1798517) q[1];
sx q[1];
rz(-1.4423749) q[1];
sx q[1];
rz(-2.1895616) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7172312) q[3];
sx q[3];
rz(-1.5699982) q[3];
sx q[3];
rz(0.80959807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6626176) q[2];
sx q[2];
rz(-1.0650029) q[2];
sx q[2];
rz(0.62620658) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(-3.0873121) q[0];
rz(0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(1.8064226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8322382) q[0];
sx q[0];
rz(-2.1119499) q[0];
sx q[0];
rz(2.4914129) q[0];
rz(-pi) q[1];
rz(-1.0338496) q[2];
sx q[2];
rz(-2.5116911) q[2];
sx q[2];
rz(0.065187188) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74354913) q[1];
sx q[1];
rz(-2.365592) q[1];
sx q[1];
rz(-2.3659124) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1109675) q[3];
sx q[3];
rz(-0.15367344) q[3];
sx q[3];
rz(-0.75017649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52428952) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(-1.3737804) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(-0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518625) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(-0.75916284) q[0];
rz(1.0836541) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(2.0738475) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6408975) q[0];
sx q[0];
rz(-2.0431378) q[0];
sx q[0];
rz(0.91820902) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75394404) q[2];
sx q[2];
rz(-0.39421668) q[2];
sx q[2];
rz(2.8682402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0800277) q[1];
sx q[1];
rz(-1.0288218) q[1];
sx q[1];
rz(-2.6656277) q[1];
rz(-2.8078733) q[3];
sx q[3];
rz(-2.2401143) q[3];
sx q[3];
rz(-3.1143318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7119673) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(-2.353904) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(-0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6982211) q[0];
sx q[0];
rz(-0.5589232) q[0];
sx q[0];
rz(1.9236175) q[0];
rz(1.632985) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(2.0532578) q[2];
sx q[2];
rz(-2.3934622) q[2];
sx q[2];
rz(-0.29665034) q[2];
rz(-1.2959215) q[3];
sx q[3];
rz(-1.7206921) q[3];
sx q[3];
rz(-0.10573798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
