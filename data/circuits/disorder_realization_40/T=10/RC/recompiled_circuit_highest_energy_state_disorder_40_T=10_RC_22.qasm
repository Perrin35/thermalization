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
rz(2.1514819) q[0];
sx q[0];
rz(2.4408542) q[0];
sx q[0];
rz(10.182575) q[0];
rz(-2.4611729) q[1];
sx q[1];
rz(-1.0935723) q[1];
sx q[1];
rz(-0.099686058) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282277) q[0];
sx q[0];
rz(-1.037957) q[0];
sx q[0];
rz(-1.9060978) q[0];
rz(2.6131479) q[2];
sx q[2];
rz(-0.64982729) q[2];
sx q[2];
rz(0.53534269) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.062391578) q[1];
sx q[1];
rz(-1.7686558) q[1];
sx q[1];
rz(-1.8366835) q[1];
rz(-pi) q[2];
rz(-1.2124196) q[3];
sx q[3];
rz(-1.626251) q[3];
sx q[3];
rz(0.90772861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0199355) q[2];
sx q[2];
rz(-2.0243553) q[2];
sx q[2];
rz(3.115444) q[2];
rz(-1.3610241) q[3];
sx q[3];
rz(-2.0676421) q[3];
sx q[3];
rz(2.0323544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015198135) q[0];
sx q[0];
rz(-2.4488738) q[0];
sx q[0];
rz(-1.3808845) q[0];
rz(2.8610002) q[1];
sx q[1];
rz(-2.3255489) q[1];
sx q[1];
rz(0.38995829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98509648) q[0];
sx q[0];
rz(-0.76129434) q[0];
sx q[0];
rz(3.0846315) q[0];
rz(-pi) q[1];
rz(2.2510914) q[2];
sx q[2];
rz(-1.599031) q[2];
sx q[2];
rz(-0.75301127) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7610723) q[1];
sx q[1];
rz(-1.7514844) q[1];
sx q[1];
rz(-1.5086345) q[1];
rz(0.51334629) q[3];
sx q[3];
rz(-1.5510552) q[3];
sx q[3];
rz(-1.6672368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0754764) q[2];
sx q[2];
rz(-2.7710997) q[2];
sx q[2];
rz(0.18744913) q[2];
rz(-0.97505331) q[3];
sx q[3];
rz(-1.8911898) q[3];
sx q[3];
rz(1.8072849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27075574) q[0];
sx q[0];
rz(-2.2360531) q[0];
sx q[0];
rz(1.4181197) q[0];
rz(-2.1899147) q[1];
sx q[1];
rz(-2.8680809) q[1];
sx q[1];
rz(-2.4804514) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1456297) q[0];
sx q[0];
rz(-1.055777) q[0];
sx q[0];
rz(1.6745695) q[0];
x q[1];
rz(1.9533526) q[2];
sx q[2];
rz(-0.91209665) q[2];
sx q[2];
rz(-0.11923583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60248261) q[1];
sx q[1];
rz(-1.6032673) q[1];
sx q[1];
rz(1.7131973) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6774353) q[3];
sx q[3];
rz(-2.1273945) q[3];
sx q[3];
rz(0.15579358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64267889) q[2];
sx q[2];
rz(-1.6940517) q[2];
sx q[2];
rz(-2.101669) q[2];
rz(-0.47636473) q[3];
sx q[3];
rz(-0.56291181) q[3];
sx q[3];
rz(-0.84144717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7845402) q[0];
sx q[0];
rz(-2.3226876) q[0];
sx q[0];
rz(-0.99405974) q[0];
rz(-1.1114063) q[1];
sx q[1];
rz(-1.2329654) q[1];
sx q[1];
rz(0.19817373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6918212) q[0];
sx q[0];
rz(-1.9110962) q[0];
sx q[0];
rz(1.4744076) q[0];
rz(-pi) q[1];
rz(2.5920766) q[2];
sx q[2];
rz(-1.1545815) q[2];
sx q[2];
rz(0.93577318) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4611283) q[1];
sx q[1];
rz(-1.8682419) q[1];
sx q[1];
rz(-2.4182085) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.010470114) q[3];
sx q[3];
rz(-1.0426077) q[3];
sx q[3];
rz(-1.8757602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.17890945) q[2];
sx q[2];
rz(-1.6390025) q[2];
sx q[2];
rz(-2.7965014) q[2];
rz(-2.1033449) q[3];
sx q[3];
rz(-2.0227261) q[3];
sx q[3];
rz(0.071917608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2800696) q[0];
sx q[0];
rz(-3.0298506) q[0];
sx q[0];
rz(-0.34568632) q[0];
rz(0.12738906) q[1];
sx q[1];
rz(-2.7368339) q[1];
sx q[1];
rz(1.5397127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72630779) q[0];
sx q[0];
rz(-1.7909701) q[0];
sx q[0];
rz(-2.9998695) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1131088) q[2];
sx q[2];
rz(-2.1808592) q[2];
sx q[2];
rz(-1.8230566) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7946164) q[1];
sx q[1];
rz(-1.8759751) q[1];
sx q[1];
rz(-2.6543314) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57347371) q[3];
sx q[3];
rz(-1.6454564) q[3];
sx q[3];
rz(2.4247313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4442399) q[2];
sx q[2];
rz(-1.8306066) q[2];
sx q[2];
rz(0.41651192) q[2];
rz(-0.087827772) q[3];
sx q[3];
rz(-0.4952687) q[3];
sx q[3];
rz(-0.4270359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86513615) q[0];
sx q[0];
rz(-0.75053954) q[0];
sx q[0];
rz(0.99391341) q[0];
rz(-1.3424501) q[1];
sx q[1];
rz(-1.7667222) q[1];
sx q[1];
rz(1.3329175) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11131824) q[0];
sx q[0];
rz(-1.9239178) q[0];
sx q[0];
rz(-1.4354491) q[0];
rz(-pi) q[1];
rz(1.7991884) q[2];
sx q[2];
rz(-1.1675541) q[2];
sx q[2];
rz(-2.5520812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7018334) q[1];
sx q[1];
rz(-1.6072131) q[1];
sx q[1];
rz(-1.4324051) q[1];
x q[2];
rz(-0.34324154) q[3];
sx q[3];
rz(-1.3784588) q[3];
sx q[3];
rz(2.1267935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6028613) q[2];
sx q[2];
rz(-1.2756602) q[2];
sx q[2];
rz(0.71962792) q[2];
rz(-2.1384278) q[3];
sx q[3];
rz(-1.4041308) q[3];
sx q[3];
rz(-2.7145332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3908983) q[0];
sx q[0];
rz(-0.4094795) q[0];
sx q[0];
rz(2.3636708) q[0];
rz(-1.3748417) q[1];
sx q[1];
rz(-0.96343103) q[1];
sx q[1];
rz(2.8533997) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2552528) q[0];
sx q[0];
rz(-1.1561516) q[0];
sx q[0];
rz(0.1840846) q[0];
rz(-pi) q[1];
rz(1.560236) q[2];
sx q[2];
rz(-1.6021072) q[2];
sx q[2];
rz(-2.8760394) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7544051) q[1];
sx q[1];
rz(-1.3164233) q[1];
sx q[1];
rz(0.53577975) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2533663) q[3];
sx q[3];
rz(-2.7446888) q[3];
sx q[3];
rz(-2.0407179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2253183) q[2];
sx q[2];
rz(-2.306873) q[2];
sx q[2];
rz(1.6342596) q[2];
rz(2.7507014) q[3];
sx q[3];
rz(-1.2619799) q[3];
sx q[3];
rz(-1.619537) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36562076) q[0];
sx q[0];
rz(-0.93637744) q[0];
sx q[0];
rz(2.9627664) q[0];
rz(1.7504182) q[1];
sx q[1];
rz(-1.4427253) q[1];
sx q[1];
rz(2.9475382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3349105) q[0];
sx q[0];
rz(-1.9437978) q[0];
sx q[0];
rz(-1.0399517) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4670635) q[2];
sx q[2];
rz(-0.87376587) q[2];
sx q[2];
rz(1.7473237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3755197) q[1];
sx q[1];
rz(-0.82036037) q[1];
sx q[1];
rz(-0.17800568) q[1];
x q[2];
rz(1.3369183) q[3];
sx q[3];
rz(-0.66576427) q[3];
sx q[3];
rz(3.105099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8630744) q[2];
sx q[2];
rz(-0.14894177) q[2];
sx q[2];
rz(-1.8982559) q[2];
rz(-0.30531034) q[3];
sx q[3];
rz(-1.6122626) q[3];
sx q[3];
rz(-0.39795157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0968001) q[0];
sx q[0];
rz(-0.069792062) q[0];
sx q[0];
rz(-0.10833649) q[0];
rz(2.4414252) q[1];
sx q[1];
rz(-1.533875) q[1];
sx q[1];
rz(2.7549423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2359206) q[0];
sx q[0];
rz(-2.2804944) q[0];
sx q[0];
rz(2.2527812) q[0];
rz(0.79603802) q[2];
sx q[2];
rz(-2.926119) q[2];
sx q[2];
rz(-0.79677671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0476847) q[1];
sx q[1];
rz(-2.3833774) q[1];
sx q[1];
rz(2.7205633) q[1];
x q[2];
rz(2.1709444) q[3];
sx q[3];
rz(-1.978756) q[3];
sx q[3];
rz(2.7768918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.75006968) q[2];
sx q[2];
rz(-1.1015026) q[2];
sx q[2];
rz(-0.92292419) q[2];
rz(-0.37724885) q[3];
sx q[3];
rz(-2.173285) q[3];
sx q[3];
rz(-1.7678461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39483505) q[0];
sx q[0];
rz(-0.45811284) q[0];
sx q[0];
rz(-0.75570345) q[0];
rz(0.79403263) q[1];
sx q[1];
rz(-0.67818063) q[1];
sx q[1];
rz(1.979535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0335498) q[0];
sx q[0];
rz(-2.6876515) q[0];
sx q[0];
rz(1.0832216) q[0];
rz(2.5799345) q[2];
sx q[2];
rz(-2.409081) q[2];
sx q[2];
rz(-1.9106302) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2311467) q[1];
sx q[1];
rz(-2.1953464) q[1];
sx q[1];
rz(-0.39448078) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1677577) q[3];
sx q[3];
rz(-1.3091901) q[3];
sx q[3];
rz(2.2261432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0572759) q[2];
sx q[2];
rz(-0.82547775) q[2];
sx q[2];
rz(-3.0153073) q[2];
rz(-0.59857357) q[3];
sx q[3];
rz(-1.5263661) q[3];
sx q[3];
rz(-1.9954782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79692688) q[0];
sx q[0];
rz(-1.4753337) q[0];
sx q[0];
rz(-1.4519806) q[0];
rz(-1.3852193) q[1];
sx q[1];
rz(-1.1744371) q[1];
sx q[1];
rz(-0.074443346) q[1];
rz(2.0695873) q[2];
sx q[2];
rz(-1.0236403) q[2];
sx q[2];
rz(-0.97579379) q[2];
rz(-1.810146) q[3];
sx q[3];
rz(-1.1905154) q[3];
sx q[3];
rz(2.1811335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
