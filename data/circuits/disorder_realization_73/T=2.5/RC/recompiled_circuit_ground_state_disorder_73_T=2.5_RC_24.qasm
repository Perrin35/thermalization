OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8483228) q[0];
sx q[0];
rz(-0.25265101) q[0];
sx q[0];
rz(1.9360315) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(-1.815058) q[1];
sx q[1];
rz(2.3958652) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099961258) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(-3.0801386) q[0];
rz(1.5775852) q[2];
sx q[2];
rz(-1.7951843) q[2];
sx q[2];
rz(0.16358384) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28994014) q[1];
sx q[1];
rz(-1.2356241) q[1];
sx q[1];
rz(-0.22214684) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0449355) q[3];
sx q[3];
rz(-1.385194) q[3];
sx q[3];
rz(1.4193648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2241609) q[2];
sx q[2];
rz(-1.5459583) q[2];
sx q[2];
rz(1.4707627) q[2];
rz(1.5077695) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(0.16954999) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(3.0068908) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8226263) q[0];
sx q[0];
rz(-1.9958226) q[0];
sx q[0];
rz(1.1898196) q[0];
x q[1];
rz(1.3943761) q[2];
sx q[2];
rz(-3.0718832) q[2];
sx q[2];
rz(1.7084054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.811482) q[1];
sx q[1];
rz(-0.54381424) q[1];
sx q[1];
rz(2.5303929) q[1];
rz(-pi) q[2];
rz(-1.991113) q[3];
sx q[3];
rz(-0.1945217) q[3];
sx q[3];
rz(1.9357301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1255101) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(0.15277319) q[2];
rz(-1.3656535) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.633054) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(0.48164865) q[0];
rz(0.18474361) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(-0.99536037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99670519) q[0];
sx q[0];
rz(-1.2249682) q[0];
sx q[0];
rz(0.24390999) q[0];
rz(-1.5103417) q[2];
sx q[2];
rz(-1.619307) q[2];
sx q[2];
rz(-1.8497576) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75835704) q[1];
sx q[1];
rz(-2.2192973) q[1];
sx q[1];
rz(1.9489669) q[1];
rz(-pi) q[2];
rz(-1.6158726) q[3];
sx q[3];
rz(-2.5451676) q[3];
sx q[3];
rz(-3.0395122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92622009) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(3.0769297) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30508405) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(2.3205561) q[0];
rz(0.07846421) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(-0.99748126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4560735) q[0];
sx q[0];
rz(-2.4302097) q[0];
sx q[0];
rz(-0.4943686) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036254739) q[2];
sx q[2];
rz(-1.5116351) q[2];
sx q[2];
rz(1.7922557) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23826829) q[1];
sx q[1];
rz(-1.0449755) q[1];
sx q[1];
rz(-0.30491288) q[1];
rz(-pi) q[2];
rz(0.11028744) q[3];
sx q[3];
rz(-2.6607041) q[3];
sx q[3];
rz(0.65910027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(1.9878261) q[2];
rz(0.18618259) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(-0.33128273) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368197) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(-1.1432884) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(-2.6236261) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85803131) q[0];
sx q[0];
rz(-2.1037344) q[0];
sx q[0];
rz(-2.3951247) q[0];
rz(-pi) q[1];
rz(-0.00092351726) q[2];
sx q[2];
rz(-1.5557655) q[2];
sx q[2];
rz(1.3236486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7593545) q[1];
sx q[1];
rz(-1.2975946) q[1];
sx q[1];
rz(1.2279991) q[1];
x q[2];
rz(0.2980026) q[3];
sx q[3];
rz(-1.4299222) q[3];
sx q[3];
rz(0.44671392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8225857) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(-2.8138568) q[2];
rz(1.94708) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(-0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1389393) q[0];
sx q[0];
rz(-0.26873538) q[0];
sx q[0];
rz(-2.5797381) q[0];
rz(-1.4909164) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(-3.0432826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033103) q[0];
sx q[0];
rz(-1.1489963) q[0];
sx q[0];
rz(-3.1027334) q[0];
rz(-pi) q[1];
rz(-1.570669) q[2];
sx q[2];
rz(-1.5714368) q[2];
sx q[2];
rz(-0.31631472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5206388) q[1];
sx q[1];
rz(-2.6700182) q[1];
sx q[1];
rz(3.0388799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0344072) q[3];
sx q[3];
rz(-1.1157728) q[3];
sx q[3];
rz(-0.95038271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(1.0657715) q[2];
rz(2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-0.081900224) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(-2.0211925) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(-0.33946005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0908112) q[0];
sx q[0];
rz(-3.0635186) q[0];
sx q[0];
rz(0.17103057) q[0];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5590057) q[2];
sx q[2];
rz(-1.5761216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.016638674) q[1];
sx q[1];
rz(-1.5724239) q[1];
sx q[1];
rz(1.485199) q[1];
rz(2.8651627) q[3];
sx q[3];
rz(-2.7904841) q[3];
sx q[3];
rz(-1.029226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7772943) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(-0.36336362) q[2];
rz(-1.0904788) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5217487) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(-0.92292619) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-0.62983477) q[1];
sx q[1];
rz(3.1373851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.370458) q[0];
sx q[0];
rz(-0.97386375) q[0];
sx q[0];
rz(0.93907066) q[0];
x q[1];
rz(0.0059295456) q[2];
sx q[2];
rz(-1.1599132) q[2];
sx q[2];
rz(3.1290999) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5397268) q[1];
sx q[1];
rz(-1.5063573) q[1];
sx q[1];
rz(1.0947202) q[1];
x q[2];
rz(-0.065807314) q[3];
sx q[3];
rz(-0.61344693) q[3];
sx q[3];
rz(2.891401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5795472) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(1.9407678) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(-0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414108) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(1.9218943) q[0];
rz(1.8118743) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(0.16389287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78744473) q[0];
sx q[0];
rz(-1.7423) q[0];
sx q[0];
rz(-1.9308912) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.741305) q[2];
sx q[2];
rz(-2.5185555) q[2];
sx q[2];
rz(-1.7965339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1023952) q[1];
sx q[1];
rz(-1.452924) q[1];
sx q[1];
rz(-3.0815691) q[1];
rz(-pi) q[2];
rz(-0.71299841) q[3];
sx q[3];
rz(-0.87942356) q[3];
sx q[3];
rz(2.612243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4280052) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-0.04341393) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(-0.77891946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.49122214) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(2.4468415) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(2.6659226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93906389) q[0];
sx q[0];
rz(-0.049204218) q[0];
sx q[0];
rz(-0.81728191) q[0];
rz(-0.86464244) q[2];
sx q[2];
rz(-1.5448031) q[2];
sx q[2];
rz(0.038250462) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1022288) q[1];
sx q[1];
rz(-2.1783531) q[1];
sx q[1];
rz(-1.5629014) q[1];
rz(-1.7559515) q[3];
sx q[3];
rz(-1.2655711) q[3];
sx q[3];
rz(-3.117962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67705578) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(-1.8439058) q[2];
rz(1.248598) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0155335) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(0.88232782) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(0.37226128) q[2];
sx q[2];
rz(-3.0426171) q[2];
sx q[2];
rz(-0.097886861) q[2];
rz(-0.23918693) q[3];
sx q[3];
rz(-1.5659955) q[3];
sx q[3];
rz(-1.5495054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
