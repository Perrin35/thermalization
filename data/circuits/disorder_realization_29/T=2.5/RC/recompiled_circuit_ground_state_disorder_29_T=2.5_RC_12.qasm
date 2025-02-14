OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.86201) q[0];
sx q[0];
rz(5.7481264) q[0];
sx q[0];
rz(9.6587107) q[0];
rz(-2.2995931) q[1];
sx q[1];
rz(-1.41398) q[1];
sx q[1];
rz(1.5185305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9568125) q[0];
sx q[0];
rz(-2.0818458) q[0];
sx q[0];
rz(1.2132116) q[0];
x q[1];
rz(-1.1359623) q[2];
sx q[2];
rz(-1.7644492) q[2];
sx q[2];
rz(0.47392148) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3564811) q[1];
sx q[1];
rz(-0.82725305) q[1];
sx q[1];
rz(-2.3712296) q[1];
rz(-pi) q[2];
rz(-0.47187658) q[3];
sx q[3];
rz(-0.96376824) q[3];
sx q[3];
rz(0.31479731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5412377) q[2];
sx q[2];
rz(-1.7011832) q[2];
sx q[2];
rz(2.3360628) q[2];
rz(-2.7496998) q[3];
sx q[3];
rz(-1.3561748) q[3];
sx q[3];
rz(2.384757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58981744) q[0];
sx q[0];
rz(-0.68142319) q[0];
sx q[0];
rz(2.7031194) q[0];
rz(1.8665727) q[1];
sx q[1];
rz(-1.4507111) q[1];
sx q[1];
rz(0.10800392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1923982) q[0];
sx q[0];
rz(-3.0279403) q[0];
sx q[0];
rz(-1.7920919) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82307016) q[2];
sx q[2];
rz(-2.0688754) q[2];
sx q[2];
rz(2.4790539) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0569563) q[1];
sx q[1];
rz(-0.26727391) q[1];
sx q[1];
rz(1.8932503) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2233954) q[3];
sx q[3];
rz(-0.51651556) q[3];
sx q[3];
rz(1.7361189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3652304) q[2];
sx q[2];
rz(-2.6680816) q[2];
sx q[2];
rz(-3.0253809) q[2];
rz(-0.41415563) q[3];
sx q[3];
rz(-1.073758) q[3];
sx q[3];
rz(-0.80348408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6666343) q[0];
sx q[0];
rz(-1.3131498) q[0];
sx q[0];
rz(-0.83589244) q[0];
rz(1.5180786) q[1];
sx q[1];
rz(-1.514879) q[1];
sx q[1];
rz(-1.2380884) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9758751) q[0];
sx q[0];
rz(-1.9414066) q[0];
sx q[0];
rz(2.1134645) q[0];
rz(-pi) q[1];
rz(-2.798271) q[2];
sx q[2];
rz(-1.3385217) q[2];
sx q[2];
rz(-0.49988036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7459384) q[1];
sx q[1];
rz(-1.6881583) q[1];
sx q[1];
rz(3.0836041) q[1];
x q[2];
rz(-1.8585303) q[3];
sx q[3];
rz(-1.8418752) q[3];
sx q[3];
rz(0.5036186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8572924) q[2];
sx q[2];
rz(-2.9967873) q[2];
sx q[2];
rz(0.71883744) q[2];
rz(-0.58088628) q[3];
sx q[3];
rz(-1.7234756) q[3];
sx q[3];
rz(-2.412793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179287) q[0];
sx q[0];
rz(-1.5939465) q[0];
sx q[0];
rz(-1.1267598) q[0];
rz(-1.843938) q[1];
sx q[1];
rz(-1.490386) q[1];
sx q[1];
rz(1.0430956) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9595056) q[0];
sx q[0];
rz(-2.8694177) q[0];
sx q[0];
rz(-1.4578117) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9773433) q[2];
sx q[2];
rz(-0.57465982) q[2];
sx q[2];
rz(-2.3324147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8759497) q[1];
sx q[1];
rz(-1.416872) q[1];
sx q[1];
rz(1.8791844) q[1];
x q[2];
rz(1.4851486) q[3];
sx q[3];
rz(-0.91721917) q[3];
sx q[3];
rz(-1.9854922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2780693) q[2];
sx q[2];
rz(-1.3845283) q[2];
sx q[2];
rz(2.1045904) q[2];
rz(2.1841124) q[3];
sx q[3];
rz(-2.0786395) q[3];
sx q[3];
rz(2.3069042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1400414) q[0];
sx q[0];
rz(-0.20714864) q[0];
sx q[0];
rz(3.0539883) q[0];
rz(0.43830782) q[1];
sx q[1];
rz(-2.0594845) q[1];
sx q[1];
rz(1.3137438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4632516) q[0];
sx q[0];
rz(-1.3366564) q[0];
sx q[0];
rz(2.0666422) q[0];
rz(1.472166) q[2];
sx q[2];
rz(-1.5336138) q[2];
sx q[2];
rz(1.94869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.61699762) q[1];
sx q[1];
rz(-1.4819078) q[1];
sx q[1];
rz(-1.3981546) q[1];
rz(-2.143712) q[3];
sx q[3];
rz(-2.8393203) q[3];
sx q[3];
rz(2.5469766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6614512) q[2];
sx q[2];
rz(-1.2378614) q[2];
sx q[2];
rz(2.5782149) q[2];
rz(-1.2207458) q[3];
sx q[3];
rz(-1.7505587) q[3];
sx q[3];
rz(0.63108546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0193943) q[0];
sx q[0];
rz(-2.4833184) q[0];
sx q[0];
rz(3.1296545) q[0];
rz(0.38966933) q[1];
sx q[1];
rz(-0.7020815) q[1];
sx q[1];
rz(0.24519244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.747904) q[0];
sx q[0];
rz(-2.1576799) q[0];
sx q[0];
rz(-0.52397195) q[0];
x q[1];
rz(0.71227534) q[2];
sx q[2];
rz(-0.33604188) q[2];
sx q[2];
rz(-3.0672376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6840382) q[1];
sx q[1];
rz(-2.0054617) q[1];
sx q[1];
rz(2.5358389) q[1];
rz(-2.526029) q[3];
sx q[3];
rz(-2.6608753) q[3];
sx q[3];
rz(1.2445039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7616854) q[2];
sx q[2];
rz(-1.3096389) q[2];
sx q[2];
rz(1.7725819) q[2];
rz(-2.6109429) q[3];
sx q[3];
rz(-0.68415087) q[3];
sx q[3];
rz(2.0556889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.948792) q[0];
sx q[0];
rz(-1.8185607) q[0];
sx q[0];
rz(-0.2970933) q[0];
rz(-1.5104712) q[1];
sx q[1];
rz(-0.62989569) q[1];
sx q[1];
rz(0.43103257) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1265701) q[0];
sx q[0];
rz(-1.3608785) q[0];
sx q[0];
rz(0.2951584) q[0];
x q[1];
rz(1.6234342) q[2];
sx q[2];
rz(-1.3402914) q[2];
sx q[2];
rz(2.8090257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7073133) q[1];
sx q[1];
rz(-0.61798862) q[1];
sx q[1];
rz(0.98843482) q[1];
rz(-2.0484974) q[3];
sx q[3];
rz(-0.58148958) q[3];
sx q[3];
rz(-3.1149816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3057574) q[2];
sx q[2];
rz(-2.367986) q[2];
sx q[2];
rz(2.8744899) q[2];
rz(-3.109572) q[3];
sx q[3];
rz(-1.1705541) q[3];
sx q[3];
rz(-0.10572461) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0778462) q[0];
sx q[0];
rz(-1.2910605) q[0];
sx q[0];
rz(-2.5015976) q[0];
rz(-2.2867639) q[1];
sx q[1];
rz(-1.6565485) q[1];
sx q[1];
rz(1.4124426) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8518191) q[0];
sx q[0];
rz(-0.69727325) q[0];
sx q[0];
rz(2.4757705) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0785901) q[2];
sx q[2];
rz(-1.5255431) q[2];
sx q[2];
rz(2.947888) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.392726) q[1];
sx q[1];
rz(-2.7905474) q[1];
sx q[1];
rz(-2.8944394) q[1];
rz(-pi) q[2];
rz(2.7667609) q[3];
sx q[3];
rz(-1.9697947) q[3];
sx q[3];
rz(0.77038902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.60014805) q[2];
sx q[2];
rz(-1.7057799) q[2];
sx q[2];
rz(-2.7195462) q[2];
rz(0.47354928) q[3];
sx q[3];
rz(-2.519042) q[3];
sx q[3];
rz(2.9464909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7710829) q[0];
sx q[0];
rz(-0.93733731) q[0];
sx q[0];
rz(1.3516082) q[0];
rz(2.8260258) q[1];
sx q[1];
rz(-1.6866997) q[1];
sx q[1];
rz(-1.1531166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0457927) q[0];
sx q[0];
rz(-1.519916) q[0];
sx q[0];
rz(-1.6232383) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0598203) q[2];
sx q[2];
rz(-2.6047629) q[2];
sx q[2];
rz(0.54807907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7866384) q[1];
sx q[1];
rz(-1.4435205) q[1];
sx q[1];
rz(-2.5580225) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38001506) q[3];
sx q[3];
rz(-1.9560062) q[3];
sx q[3];
rz(-2.4917701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45712581) q[2];
sx q[2];
rz(-1.9292597) q[2];
sx q[2];
rz(2.7678164) q[2];
rz(1.2260381) q[3];
sx q[3];
rz(-2.2235179) q[3];
sx q[3];
rz(-1.8019684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1821197) q[0];
sx q[0];
rz(-0.68030334) q[0];
sx q[0];
rz(-1.8208338) q[0];
rz(-0.93718115) q[1];
sx q[1];
rz(-1.3359741) q[1];
sx q[1];
rz(0.46930596) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8784638) q[0];
sx q[0];
rz(-1.6471905) q[0];
sx q[0];
rz(1.4532277) q[0];
rz(-pi) q[1];
rz(0.98913828) q[2];
sx q[2];
rz(-1.4114221) q[2];
sx q[2];
rz(-1.5053269) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7368245) q[1];
sx q[1];
rz(-1.4206845) q[1];
sx q[1];
rz(1.4895358) q[1];
rz(-2.4604581) q[3];
sx q[3];
rz(-0.38060846) q[3];
sx q[3];
rz(-1.8073624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4361973) q[2];
sx q[2];
rz(-2.2787978) q[2];
sx q[2];
rz(1.9311284) q[2];
rz(-0.56810275) q[3];
sx q[3];
rz(-1.3848687) q[3];
sx q[3];
rz(0.97584045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217011) q[0];
sx q[0];
rz(-2.2911063) q[0];
sx q[0];
rz(-1.6794857) q[0];
rz(2.2484491) q[1];
sx q[1];
rz(-1.3124663) q[1];
sx q[1];
rz(1.5193473) q[1];
rz(-0.22994269) q[2];
sx q[2];
rz(-1.2886921) q[2];
sx q[2];
rz(-3.1293426) q[2];
rz(1.0004956) q[3];
sx q[3];
rz(-1.9420997) q[3];
sx q[3];
rz(1.6122769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
