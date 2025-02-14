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
rz(-2.1751997) q[0];
sx q[0];
rz(4.8973358) q[0];
sx q[0];
rz(10.020221) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(2.5505677) q[1];
sx q[1];
rz(12.443065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46245136) q[0];
sx q[0];
rz(-1.9872905) q[0];
sx q[0];
rz(-0.97521675) q[0];
rz(-0.78322905) q[2];
sx q[2];
rz(-0.95942623) q[2];
sx q[2];
rz(-2.2303344) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6001624) q[1];
sx q[1];
rz(-1.3568391) q[1];
sx q[1];
rz(0.096264953) q[1];
x q[2];
rz(-1.7233347) q[3];
sx q[3];
rz(-1.9328874) q[3];
sx q[3];
rz(0.83521508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82723242) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(3.0296791) q[2];
rz(1.3917475) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(-2.8351423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157442) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(-2.9916905) q[0];
rz(-2.8556178) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(2.012595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492831) q[0];
sx q[0];
rz(-2.8516927) q[0];
sx q[0];
rz(-1.4647746) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3311884) q[2];
sx q[2];
rz(-3.0802266) q[2];
sx q[2];
rz(-0.97739109) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6817878) q[1];
sx q[1];
rz(-2.7437609) q[1];
sx q[1];
rz(-0.43983207) q[1];
rz(-2.3150857) q[3];
sx q[3];
rz(-2.3621268) q[3];
sx q[3];
rz(0.41050467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67843208) q[2];
sx q[2];
rz(-0.9223991) q[2];
sx q[2];
rz(-1.9909667) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.4169644) q[3];
sx q[3];
rz(-0.020603389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(0.38159698) q[0];
rz(1.8474139) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(-0.19827422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37368942) q[0];
sx q[0];
rz(-1.7977909) q[0];
sx q[0];
rz(-2.9625708) q[0];
rz(0.90132331) q[2];
sx q[2];
rz(-1.633989) q[2];
sx q[2];
rz(-0.92869273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72858649) q[1];
sx q[1];
rz(-0.31794128) q[1];
sx q[1];
rz(0.7403265) q[1];
rz(3.1214633) q[3];
sx q[3];
rz(-1.745589) q[3];
sx q[3];
rz(-0.59390261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76564378) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(-2.622733) q[2];
rz(-1.8910003) q[3];
sx q[3];
rz(-2.0542681) q[3];
sx q[3];
rz(-2.5605104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64312235) q[0];
sx q[0];
rz(-2.9415218) q[0];
sx q[0];
rz(-2.6923687) q[0];
rz(1.6953702) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(-0.62072388) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3005488) q[0];
sx q[0];
rz(-0.67557205) q[0];
sx q[0];
rz(0.23209693) q[0];
rz(0.63346699) q[2];
sx q[2];
rz(-1.6704428) q[2];
sx q[2];
rz(-1.5758621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6171332) q[1];
sx q[1];
rz(-1.7726328) q[1];
sx q[1];
rz(-2.733888) q[1];
rz(0.053106282) q[3];
sx q[3];
rz(-2.1087007) q[3];
sx q[3];
rz(-2.2369179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0351403) q[2];
sx q[2];
rz(-1.8702714) q[2];
sx q[2];
rz(-0.16109666) q[2];
rz(-2.7437239) q[3];
sx q[3];
rz(-1.4784003) q[3];
sx q[3];
rz(-0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-0.75905269) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(-0.25704849) q[0];
rz(0.88047475) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(-1.8880728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.391025) q[0];
sx q[0];
rz(-1.5888968) q[0];
sx q[0];
rz(-1.6068293) q[0];
x q[1];
rz(-0.23373105) q[2];
sx q[2];
rz(-2.3178551) q[2];
sx q[2];
rz(-1.1540138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94379986) q[1];
sx q[1];
rz(-1.5256923) q[1];
sx q[1];
rz(0.79774858) q[1];
x q[2];
rz(2.7944466) q[3];
sx q[3];
rz(-1.5691363) q[3];
sx q[3];
rz(1.1356789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4132061) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(1.8750635) q[2];
rz(-0.62147102) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(2.6004041) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40015873) q[0];
sx q[0];
rz(-2.3024237) q[0];
sx q[0];
rz(-0.025064502) q[0];
rz(-1.2114245) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(-2.8866344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716041) q[0];
sx q[0];
rz(-1.5657664) q[0];
sx q[0];
rz(-3.0922573) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67398418) q[2];
sx q[2];
rz(-0.70491284) q[2];
sx q[2];
rz(0.92943905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7073696) q[1];
sx q[1];
rz(-2.1277782) q[1];
sx q[1];
rz(-0.93644489) q[1];
x q[2];
rz(0.99825777) q[3];
sx q[3];
rz(-1.3911595) q[3];
sx q[3];
rz(-1.4416173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24594626) q[2];
sx q[2];
rz(-1.0422372) q[2];
sx q[2];
rz(2.6825405) q[2];
rz(-1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46886214) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(-3.0555308) q[0];
rz(-0.91880265) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(-1.2109717) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54706193) q[0];
sx q[0];
rz(-1.8939202) q[0];
sx q[0];
rz(1.3573775) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2077256) q[2];
sx q[2];
rz(-2.4709765) q[2];
sx q[2];
rz(-0.071337168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0667292) q[1];
sx q[1];
rz(-1.7742429) q[1];
sx q[1];
rz(-1.2438846) q[1];
rz(-pi) q[2];
rz(-1.9350697) q[3];
sx q[3];
rz(-1.9134054) q[3];
sx q[3];
rz(-2.2628257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3226037) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(0.71713478) q[2];
rz(1.0273733) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(-2.7023442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9913919) q[0];
sx q[0];
rz(-0.99656492) q[0];
sx q[0];
rz(-0.014658654) q[0];
rz(2.390059) q[1];
sx q[1];
rz(-1.2208168) q[1];
sx q[1];
rz(1.470648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6122307) q[0];
sx q[0];
rz(-0.23461537) q[0];
sx q[0];
rz(-2.505872) q[0];
rz(-0.25663767) q[2];
sx q[2];
rz(-1.6483288) q[2];
sx q[2];
rz(2.3687349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9110221) q[1];
sx q[1];
rz(-1.9729905) q[1];
sx q[1];
rz(-0.12253094) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1378391) q[3];
sx q[3];
rz(-1.9158746) q[3];
sx q[3];
rz(1.0964637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26059255) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(-2.8680958) q[2];
rz(-1.1547487) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(-2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038789373) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(2.0661085) q[0];
rz(-3.0134046) q[1];
sx q[1];
rz(-2.2905541) q[1];
sx q[1];
rz(2.7959965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49508383) q[0];
sx q[0];
rz(-0.91143196) q[0];
sx q[0];
rz(2.93883) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1394452) q[2];
sx q[2];
rz(-1.8446184) q[2];
sx q[2];
rz(1.7632926) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84050485) q[1];
sx q[1];
rz(-1.6337758) q[1];
sx q[1];
rz(-2.5999864) q[1];
rz(-1.7365428) q[3];
sx q[3];
rz(-0.9791383) q[3];
sx q[3];
rz(-1.6935284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(1.5209939) q[2];
rz(1.3711551) q[3];
sx q[3];
rz(-2.3833279) q[3];
sx q[3];
rz(-0.41485205) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3138251) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-2.9058822) q[0];
rz(-0.57304263) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(1.7726353) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0018250759) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(2.1609309) q[0];
x q[1];
rz(2.4339381) q[2];
sx q[2];
rz(-2.3295412) q[2];
sx q[2];
rz(0.52834137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6073709) q[1];
sx q[1];
rz(-1.0080308) q[1];
sx q[1];
rz(0.20197473) q[1];
rz(-pi) q[2];
rz(-2.8970191) q[3];
sx q[3];
rz(-0.37062708) q[3];
sx q[3];
rz(-1.9300435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77091757) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(3.0885922) q[2];
rz(0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62405217) q[0];
sx q[0];
rz(-2.1255827) q[0];
sx q[0];
rz(-0.95638635) q[0];
rz(-1.3507631) q[1];
sx q[1];
rz(-2.7443934) q[1];
sx q[1];
rz(2.9846356) q[1];
rz(-2.9861957) q[2];
sx q[2];
rz(-2.4383895) q[2];
sx q[2];
rz(-0.40021605) q[2];
rz(-0.57784537) q[3];
sx q[3];
rz(-0.89580775) q[3];
sx q[3];
rz(2.5474594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
