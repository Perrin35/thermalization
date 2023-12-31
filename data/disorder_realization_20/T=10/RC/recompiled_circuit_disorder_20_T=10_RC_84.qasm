OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(2.9291908) q[0];
rz(0.71495932) q[1];
sx q[1];
rz(3.9290805) q[1];
sx q[1];
rz(10.706283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4709269) q[0];
sx q[0];
rz(-0.39917329) q[0];
sx q[0];
rz(2.1098718) q[0];
rz(-pi) q[1];
rz(-1.6002866) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(-2.8147547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9077397) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(-1.185226) q[1];
rz(-1.1164066) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(2.4860399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3296559) q[0];
sx q[0];
rz(-2.5663079) q[0];
sx q[0];
rz(-0.29559691) q[0];
x q[1];
rz(2.4376051) q[2];
sx q[2];
rz(-0.96103243) q[2];
sx q[2];
rz(-3.0733382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21178791) q[1];
sx q[1];
rz(-2.1165407) q[1];
sx q[1];
rz(0.8916698) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21807166) q[3];
sx q[3];
rz(-3.0017188) q[3];
sx q[3];
rz(1.1034031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(0.29671159) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(-2.6112774) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.9525607) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2169164) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(2.6485373) q[0];
rz(1.5306773) q[2];
sx q[2];
rz(-1.8187858) q[2];
sx q[2];
rz(2.9330394) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6991068) q[1];
sx q[1];
rz(-0.92256303) q[1];
sx q[1];
rz(0.68251619) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0344814) q[3];
sx q[3];
rz(-0.92671219) q[3];
sx q[3];
rz(0.37582276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90551463) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(-0.45423147) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.8006515) q[0];
rz(-2.5355133) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.1846503) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9979447) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(-0.68302897) q[0];
rz(-pi) q[1];
rz(-2.2934154) q[2];
sx q[2];
rz(-0.64372534) q[2];
sx q[2];
rz(-1.4959469) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21936101) q[1];
sx q[1];
rz(-2.1582099) q[1];
sx q[1];
rz(-3.0012793) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70373669) q[3];
sx q[3];
rz(-2.6061213) q[3];
sx q[3];
rz(-0.4243917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(2.502029) q[2];
rz(-2.284164) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.5109167) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(-0.52039352) q[0];
rz(-3.0386472) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-0.72881126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33671185) q[0];
sx q[0];
rz(-1.6687972) q[0];
sx q[0];
rz(-2.9360807) q[0];
rz(-pi) q[1];
rz(-0.16291933) q[2];
sx q[2];
rz(-0.85328057) q[2];
sx q[2];
rz(-1.7999072) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.049206991) q[1];
sx q[1];
rz(-1.2485463) q[1];
sx q[1];
rz(-0.21578034) q[1];
rz(2.6008984) q[3];
sx q[3];
rz(-1.027642) q[3];
sx q[3];
rz(0.83609304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(-0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1028041) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(0.9206413) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52299243) q[0];
sx q[0];
rz(-0.97609659) q[0];
sx q[0];
rz(2.3408875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4913885) q[2];
sx q[2];
rz(-1.4646155) q[2];
sx q[2];
rz(-0.90027819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9787394) q[1];
sx q[1];
rz(-1.1161242) q[1];
sx q[1];
rz(-1.5233558) q[1];
x q[2];
rz(2.382155) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(1.9721608) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5116611) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(1.1849665) q[0];
rz(-pi) q[1];
rz(-2.7514236) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(0.90261501) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2028423) q[1];
sx q[1];
rz(-0.65333594) q[1];
sx q[1];
rz(0.57196879) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5285989) q[3];
sx q[3];
rz(-1.183126) q[3];
sx q[3];
rz(-1.8488415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.823267) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(-0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(0.56234223) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(2.6838578) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477752) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(1.093822) q[0];
rz(-2.2675603) q[2];
sx q[2];
rz(-2.0119785) q[2];
sx q[2];
rz(0.87768427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9131635) q[1];
sx q[1];
rz(-0.73691165) q[1];
sx q[1];
rz(1.5921028) q[1];
x q[2];
rz(0.3433414) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(0.2688558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(2.8590554) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(-2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.2618077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5489053) q[0];
sx q[0];
rz(-1.3015916) q[0];
sx q[0];
rz(0.44434987) q[0];
x q[1];
rz(2.8724573) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(0.035949635) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0314434) q[1];
sx q[1];
rz(-1.3182536) q[1];
sx q[1];
rz(0.7048216) q[1];
rz(-pi) q[2];
rz(-2.7997167) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(-3.0059909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70301473) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(1.6837439) q[2];
rz(0.66649377) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2465729) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.9816459) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(1.0704401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69660891) q[0];
sx q[0];
rz(-0.84060003) q[0];
sx q[0];
rz(-0.22459774) q[0];
rz(-pi) q[1];
rz(-3.1033377) q[2];
sx q[2];
rz(-0.80637156) q[2];
sx q[2];
rz(0.56626608) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1420319) q[1];
sx q[1];
rz(-2.3107717) q[1];
sx q[1];
rz(1.9097411) q[1];
x q[2];
rz(-1.3954193) q[3];
sx q[3];
rz(-1.8150107) q[3];
sx q[3];
rz(2.0004686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3728309) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(-1.6147511) q[2];
rz(-2.5620143) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(0.037739567) q[2];
sx q[2];
rz(-0.83913091) q[2];
sx q[2];
rz(-0.057374949) q[2];
rz(-2.3970849) q[3];
sx q[3];
rz(-2.4662938) q[3];
sx q[3];
rz(0.288356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
