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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(2.733736) q[0];
rz(3.1240533) q[1];
sx q[1];
rz(-1.2516302) q[1];
sx q[1];
rz(1.6012021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33690573) q[0];
sx q[0];
rz(-1.1338455) q[0];
sx q[0];
rz(-3.0880465) q[0];
x q[1];
rz(-1.7061739) q[2];
sx q[2];
rz(-1.1203566) q[2];
sx q[2];
rz(3.0036486) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1931128) q[1];
sx q[1];
rz(-2.6958637) q[1];
sx q[1];
rz(-2.6325001) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92272051) q[3];
sx q[3];
rz(-0.86258234) q[3];
sx q[3];
rz(-3.129209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5539598) q[2];
sx q[2];
rz(-1.0593869) q[2];
sx q[2];
rz(0.21767347) q[2];
rz(0.31146464) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(1.9757087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199453) q[0];
sx q[0];
rz(-0.2758652) q[0];
sx q[0];
rz(0.69197792) q[0];
rz(2.3852589) q[1];
sx q[1];
rz(-0.5223918) q[1];
sx q[1];
rz(1.5914894) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1082008) q[0];
sx q[0];
rz(-0.72275677) q[0];
sx q[0];
rz(0.18527822) q[0];
x q[1];
rz(1.3266356) q[2];
sx q[2];
rz(-2.7327644) q[2];
sx q[2];
rz(-0.93016184) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31526947) q[1];
sx q[1];
rz(-0.20139748) q[1];
sx q[1];
rz(-2.4884495) q[1];
rz(-pi) q[2];
rz(-0.49868985) q[3];
sx q[3];
rz(-2.0741921) q[3];
sx q[3];
rz(-1.4273912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99593607) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(0.23925979) q[2];
rz(-1.4907106) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(-1.9132805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.2759129) q[0];
sx q[0];
rz(-2.2353807) q[0];
sx q[0];
rz(-0.014634125) q[0];
rz(-0.89302653) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-2.9258974) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0365643) q[0];
sx q[0];
rz(-1.6220548) q[0];
sx q[0];
rz(3.1369484) q[0];
rz(-pi) q[1];
rz(-0.65962445) q[2];
sx q[2];
rz(-1.8462034) q[2];
sx q[2];
rz(0.026247488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5758504) q[1];
sx q[1];
rz(-1.8025959) q[1];
sx q[1];
rz(-1.4398458) q[1];
x q[2];
rz(-0.78554947) q[3];
sx q[3];
rz(-2.2758898) q[3];
sx q[3];
rz(2.7712703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0821685) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(2.1480985) q[2];
rz(-2.4332186) q[3];
sx q[3];
rz(-1.1185442) q[3];
sx q[3];
rz(0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.74743903) q[0];
sx q[0];
rz(-2.6469632) q[0];
sx q[0];
rz(-2.4476449) q[0];
rz(-0.28678647) q[1];
sx q[1];
rz(-1.6096121) q[1];
sx q[1];
rz(2.0538816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75380775) q[0];
sx q[0];
rz(-1.0669363) q[0];
sx q[0];
rz(0.14150454) q[0];
rz(1.0469646) q[2];
sx q[2];
rz(-2.0598542) q[2];
sx q[2];
rz(0.2917052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4162594) q[1];
sx q[1];
rz(-0.97914588) q[1];
sx q[1];
rz(-2.7736431) q[1];
rz(-0.71432738) q[3];
sx q[3];
rz(-1.5363524) q[3];
sx q[3];
rz(1.8689276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9168758) q[2];
sx q[2];
rz(-1.9852873) q[2];
sx q[2];
rz(-1.7737596) q[2];
rz(1.2380838) q[3];
sx q[3];
rz(-1.5690469) q[3];
sx q[3];
rz(2.9253173) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5424159) q[0];
sx q[0];
rz(-1.8734064) q[0];
sx q[0];
rz(2.5237778) q[0];
rz(-1.1082331) q[1];
sx q[1];
rz(-0.30312678) q[1];
sx q[1];
rz(0.88315001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9882433) q[0];
sx q[0];
rz(-2.3267965) q[0];
sx q[0];
rz(-1.5606461) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4083715) q[2];
sx q[2];
rz(-1.139763) q[2];
sx q[2];
rz(-1.169432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2715192) q[1];
sx q[1];
rz(-0.49527676) q[1];
sx q[1];
rz(1.7090767) q[1];
x q[2];
rz(-1.0155025) q[3];
sx q[3];
rz(-1.9262553) q[3];
sx q[3];
rz(-2.773284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21657476) q[2];
sx q[2];
rz(-2.8897132) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(1.2733634) q[3];
sx q[3];
rz(-1.7073771) q[3];
sx q[3];
rz(-0.98289615) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044416044) q[0];
sx q[0];
rz(-1.2460848) q[0];
sx q[0];
rz(-2.521305) q[0];
rz(-2.7445131) q[1];
sx q[1];
rz(-0.47080165) q[1];
sx q[1];
rz(-1.9420067) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.202855) q[0];
sx q[0];
rz(-1.5498383) q[0];
sx q[0];
rz(1.6771132) q[0];
rz(-1.7825837) q[2];
sx q[2];
rz(-1.7149936) q[2];
sx q[2];
rz(-1.5825001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5127677) q[1];
sx q[1];
rz(-2.5803452) q[1];
sx q[1];
rz(0.36833771) q[1];
rz(2.7971917) q[3];
sx q[3];
rz(-1.4592429) q[3];
sx q[3];
rz(1.9007746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8951796) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(-2.4708774) q[2];
rz(2.8504168) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(0.62180716) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42596844) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(0.17844644) q[0];
rz(-2.2715691) q[1];
sx q[1];
rz(-2.2585637) q[1];
sx q[1];
rz(-0.80284405) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74551302) q[0];
sx q[0];
rz(-0.95142309) q[0];
sx q[0];
rz(-1.6016763) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8151547) q[2];
sx q[2];
rz(-1.0006703) q[2];
sx q[2];
rz(-1.9889835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72910129) q[1];
sx q[1];
rz(-1.2454528) q[1];
sx q[1];
rz(0.71112432) q[1];
rz(0.25814806) q[3];
sx q[3];
rz(-2.1827563) q[3];
sx q[3];
rz(1.1560837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7650083) q[2];
sx q[2];
rz(-1.8439801) q[2];
sx q[2];
rz(0.89361781) q[2];
rz(1.5417967) q[3];
sx q[3];
rz(-1.6518281) q[3];
sx q[3];
rz(0.60683513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97745085) q[0];
sx q[0];
rz(-0.41005382) q[0];
sx q[0];
rz(0.2151016) q[0];
rz(2.2970301) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(2.3473158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71090224) q[0];
sx q[0];
rz(-1.5973418) q[0];
sx q[0];
rz(-1.558976) q[0];
x q[1];
rz(0.45659839) q[2];
sx q[2];
rz(-1.8482882) q[2];
sx q[2];
rz(0.014739787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2596742) q[1];
sx q[1];
rz(-1.6615903) q[1];
sx q[1];
rz(-1.0301678) q[1];
rz(-pi) q[2];
rz(2.8135962) q[3];
sx q[3];
rz(-0.80910002) q[3];
sx q[3];
rz(-2.3878333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8636785) q[2];
sx q[2];
rz(-1.9387551) q[2];
sx q[2];
rz(-1.4097144) q[2];
rz(0.75285161) q[3];
sx q[3];
rz(-1.1466305) q[3];
sx q[3];
rz(-2.1394155) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3442605) q[0];
sx q[0];
rz(-0.40818885) q[0];
sx q[0];
rz(-0.97287384) q[0];
rz(-1.6196039) q[1];
sx q[1];
rz(-1.7771959) q[1];
sx q[1];
rz(-0.63046986) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3469543) q[0];
sx q[0];
rz(-1.3890059) q[0];
sx q[0];
rz(1.7111045) q[0];
rz(-pi) q[1];
rz(-1.4138572) q[2];
sx q[2];
rz(-1.5782159) q[2];
sx q[2];
rz(2.3049624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7122102) q[1];
sx q[1];
rz(-0.80388596) q[1];
sx q[1];
rz(1.6401119) q[1];
rz(2.7961066) q[3];
sx q[3];
rz(-2.957323) q[3];
sx q[3];
rz(2.1994182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1260599) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(-1.3014334) q[2];
rz(-1.556501) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(-0.41527709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62366098) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(-2.6483722) q[0];
rz(1.2552235) q[1];
sx q[1];
rz(-1.247765) q[1];
sx q[1];
rz(-2.7769322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55281524) q[0];
sx q[0];
rz(-1.4328261) q[0];
sx q[0];
rz(-1.2347883) q[0];
rz(-pi) q[1];
rz(-2.1388059) q[2];
sx q[2];
rz(-0.73633654) q[2];
sx q[2];
rz(-1.9970837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2233878) q[1];
sx q[1];
rz(-2.3406696) q[1];
sx q[1];
rz(-3.0632698) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8871898) q[3];
sx q[3];
rz(-1.5324161) q[3];
sx q[3];
rz(2.5298713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6898592) q[2];
sx q[2];
rz(-2.1351337) q[2];
sx q[2];
rz(-2.149392) q[2];
rz(-2.774488) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(-2.4632857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.622396) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(0.92338152) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(1.3707163) q[2];
sx q[2];
rz(-1.4309819) q[2];
sx q[2];
rz(-0.33071721) q[2];
rz(-0.97040776) q[3];
sx q[3];
rz(-0.98894377) q[3];
sx q[3];
rz(0.43882216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
