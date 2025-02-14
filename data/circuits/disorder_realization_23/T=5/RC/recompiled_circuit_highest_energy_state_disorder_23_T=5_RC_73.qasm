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
rz(0.55573207) q[0];
sx q[0];
rz(-1.8620123) q[0];
sx q[0];
rz(-0.32787856) q[0];
rz(0.15283395) q[1];
sx q[1];
rz(-0.489355) q[1];
sx q[1];
rz(-2.1305003) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7578825) q[0];
sx q[0];
rz(-0.61311904) q[0];
sx q[0];
rz(0.28045537) q[0];
x q[1];
rz(1.7494124) q[2];
sx q[2];
rz(-1.9733323) q[2];
sx q[2];
rz(2.2312763) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0786405) q[1];
sx q[1];
rz(-2.302189) q[1];
sx q[1];
rz(-0.49830978) q[1];
x q[2];
rz(-2.070363) q[3];
sx q[3];
rz(-0.46854436) q[3];
sx q[3];
rz(0.51817453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27935394) q[2];
sx q[2];
rz(-0.78247672) q[2];
sx q[2];
rz(1.5581101) q[2];
rz(2.8067348) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(-2.8607232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25664499) q[0];
sx q[0];
rz(-0.56459752) q[0];
sx q[0];
rz(-0.33194342) q[0];
rz(2.7815107) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(0.28775451) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641787) q[0];
sx q[0];
rz(-2.8197643) q[0];
sx q[0];
rz(0.68049707) q[0];
x q[1];
rz(-0.73189484) q[2];
sx q[2];
rz(-1.8127155) q[2];
sx q[2];
rz(0.11920028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33720484) q[1];
sx q[1];
rz(-1.241695) q[1];
sx q[1];
rz(1.9371402) q[1];
rz(-0.61057456) q[3];
sx q[3];
rz(-0.27974162) q[3];
sx q[3];
rz(-2.4795462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.687279) q[2];
sx q[2];
rz(-1.5328898) q[2];
sx q[2];
rz(1.3909371) q[2];
rz(0.85919291) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(-0.19101645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.745382) q[0];
sx q[0];
rz(-1.5353545) q[0];
sx q[0];
rz(-2.9111653) q[0];
rz(-0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(-2.8527625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17265564) q[0];
sx q[0];
rz(-2.7285577) q[0];
sx q[0];
rz(-1.719219) q[0];
rz(-pi) q[1];
rz(2.6730113) q[2];
sx q[2];
rz(-1.5043253) q[2];
sx q[2];
rz(0.85487142) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17902148) q[1];
sx q[1];
rz(-1.5177392) q[1];
sx q[1];
rz(1.9808484) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2346441) q[3];
sx q[3];
rz(-0.28093279) q[3];
sx q[3];
rz(-0.30981608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1032224) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(0.040741097) q[2];
rz(3.0401958) q[3];
sx q[3];
rz(-1.3359759) q[3];
sx q[3];
rz(2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539702) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(2.907584) q[0];
rz(-2.9435844) q[1];
sx q[1];
rz(-1.0375236) q[1];
sx q[1];
rz(2.1580946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6925544) q[0];
sx q[0];
rz(-0.48547599) q[0];
sx q[0];
rz(1.4794769) q[0];
x q[1];
rz(-2.9806541) q[2];
sx q[2];
rz(-0.51921028) q[2];
sx q[2];
rz(-2.1759261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9069104) q[1];
sx q[1];
rz(-2.249472) q[1];
sx q[1];
rz(-2.3386416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4928905) q[3];
sx q[3];
rz(-0.39038218) q[3];
sx q[3];
rz(-1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6733751) q[2];
sx q[2];
rz(-2.8523291) q[2];
sx q[2];
rz(0.77486983) q[2];
rz(-1.8153927) q[3];
sx q[3];
rz(-1.9522791) q[3];
sx q[3];
rz(0.51138043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73959094) q[0];
sx q[0];
rz(-2.2690161) q[0];
sx q[0];
rz(-0.032489754) q[0];
rz(1.912311) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(0.083018735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1047528) q[0];
sx q[0];
rz(-2.6304465) q[0];
sx q[0];
rz(0.47170659) q[0];
x q[1];
rz(-1.4013702) q[2];
sx q[2];
rz(-1.8352785) q[2];
sx q[2];
rz(-1.6319909) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8578945) q[1];
sx q[1];
rz(-2.2304728) q[1];
sx q[1];
rz(-0.96925756) q[1];
rz(-pi) q[2];
x q[2];
rz(0.015542726) q[3];
sx q[3];
rz(-0.93399007) q[3];
sx q[3];
rz(-1.4525082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7897537) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(0.16072533) q[2];
rz(1.1927346) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6650894) q[0];
sx q[0];
rz(-2.022321) q[0];
sx q[0];
rz(0.29092586) q[0];
rz(-1.3833969) q[1];
sx q[1];
rz(-1.8427126) q[1];
sx q[1];
rz(2.6417522) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4668149) q[0];
sx q[0];
rz(-1.5272015) q[0];
sx q[0];
rz(1.6131667) q[0];
rz(-pi) q[1];
x q[1];
rz(2.585829) q[2];
sx q[2];
rz(-0.92257181) q[2];
sx q[2];
rz(-0.52582914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0840869) q[1];
sx q[1];
rz(-2.018714) q[1];
sx q[1];
rz(0.65051669) q[1];
rz(-2.1736702) q[3];
sx q[3];
rz(-2.335066) q[3];
sx q[3];
rz(-1.6953118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9342039) q[2];
sx q[2];
rz(-2.5298205) q[2];
sx q[2];
rz(0.70154166) q[2];
rz(2.1675341) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(-0.90132236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262254) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(2.4819964) q[0];
rz(1.9024128) q[1];
sx q[1];
rz(-1.0737123) q[1];
sx q[1];
rz(-1.0741796) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91759085) q[0];
sx q[0];
rz(-1.5573643) q[0];
sx q[0];
rz(-1.5852527) q[0];
x q[1];
rz(0.72788179) q[2];
sx q[2];
rz(-0.54833503) q[2];
sx q[2];
rz(2.9408583) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32666001) q[1];
sx q[1];
rz(-2.6505917) q[1];
sx q[1];
rz(-1.0990952) q[1];
rz(2.6857008) q[3];
sx q[3];
rz(-0.8178725) q[3];
sx q[3];
rz(-2.4867833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85019511) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(-0.58471739) q[2];
rz(1.0662474) q[3];
sx q[3];
rz(-1.9138347) q[3];
sx q[3];
rz(-2.7339981) q[3];
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
rz(-0.19121118) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(2.6608652) q[0];
rz(-1.2113384) q[1];
sx q[1];
rz(-1.6119266) q[1];
sx q[1];
rz(0.617625) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4803333) q[0];
sx q[0];
rz(-0.32513371) q[0];
sx q[0];
rz(-1.5576157) q[0];
rz(-pi) q[1];
rz(1.0477871) q[2];
sx q[2];
rz(-2.1126502) q[2];
sx q[2];
rz(2.7043992) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.94719515) q[1];
sx q[1];
rz(-2.7734967) q[1];
sx q[1];
rz(-0.93666623) q[1];
rz(-pi) q[2];
rz(2.9493773) q[3];
sx q[3];
rz(-1.6055067) q[3];
sx q[3];
rz(0.46958967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1787313) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(2.2507131) q[2];
rz(-1.0293055) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(-1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35128281) q[0];
sx q[0];
rz(-2.2064378) q[0];
sx q[0];
rz(1.1736322) q[0];
rz(-1.952518) q[1];
sx q[1];
rz(-2.3937841) q[1];
sx q[1];
rz(-0.48114166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089398459) q[0];
sx q[0];
rz(-1.7939095) q[0];
sx q[0];
rz(2.5714178) q[0];
rz(0.5598716) q[2];
sx q[2];
rz(-1.6563043) q[2];
sx q[2];
rz(-0.97031236) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9607842) q[1];
sx q[1];
rz(-2.167302) q[1];
sx q[1];
rz(-0.1779314) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88925006) q[3];
sx q[3];
rz(-0.64683952) q[3];
sx q[3];
rz(1.3726039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9159307) q[2];
sx q[2];
rz(-1.91232) q[2];
sx q[2];
rz(1.6602328) q[2];
rz(3.0952752) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(-1.9143298) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7939821) q[0];
sx q[0];
rz(-1.0337669) q[0];
sx q[0];
rz(3.0718497) q[0];
rz(2.8522885) q[1];
sx q[1];
rz(-1.8569088) q[1];
sx q[1];
rz(1.0522254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3892059) q[0];
sx q[0];
rz(-1.2121887) q[0];
sx q[0];
rz(-1.9176656) q[0];
rz(-pi) q[1];
rz(-2.2329726) q[2];
sx q[2];
rz(-1.7025456) q[2];
sx q[2];
rz(0.8140623) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47202808) q[1];
sx q[1];
rz(-1.2109562) q[1];
sx q[1];
rz(-1.9739975) q[1];
rz(-pi) q[2];
rz(0.42843786) q[3];
sx q[3];
rz(-1.9588545) q[3];
sx q[3];
rz(-0.65754277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5436486) q[2];
sx q[2];
rz(-1.9660549) q[2];
sx q[2];
rz(-0.0607461) q[2];
rz(-2.8525823) q[3];
sx q[3];
rz(-1.3649536) q[3];
sx q[3];
rz(1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.65171) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(-1.0501077) q[1];
sx q[1];
rz(-2.0790015) q[1];
sx q[1];
rz(1.9008295) q[1];
rz(1.1012668) q[2];
sx q[2];
rz(-1.4360089) q[2];
sx q[2];
rz(0.77965005) q[2];
rz(1.0084739) q[3];
sx q[3];
rz(-2.6431966) q[3];
sx q[3];
rz(0.24503844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
