OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7933529) q[0];
sx q[0];
rz(-1.5433595) q[0];
sx q[0];
rz(-1.6399075) q[0];
rz(-1.3950672) q[1];
sx q[1];
rz(2.7434064) q[1];
sx q[1];
rz(12.413496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9119551) q[0];
sx q[0];
rz(-2.9792333) q[0];
sx q[0];
rz(2.1549757) q[0];
x q[1];
rz(-2.7911573) q[2];
sx q[2];
rz(-1.7261243) q[2];
sx q[2];
rz(0.84463476) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7167048) q[1];
sx q[1];
rz(-0.198303) q[1];
sx q[1];
rz(1.2601869) q[1];
rz(-pi) q[2];
rz(2.128162) q[3];
sx q[3];
rz(-2.3881222) q[3];
sx q[3];
rz(-2.4581152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.07831002) q[2];
sx q[2];
rz(-1.4274884) q[2];
sx q[2];
rz(-0.59986344) q[2];
rz(2.6317224) q[3];
sx q[3];
rz(-0.32250914) q[3];
sx q[3];
rz(2.9738284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55709368) q[0];
sx q[0];
rz(-2.2826513) q[0];
sx q[0];
rz(2.9993045) q[0];
rz(1.2917057) q[1];
sx q[1];
rz(-2.6085491) q[1];
sx q[1];
rz(0.60089111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37764022) q[0];
sx q[0];
rz(-2.6598192) q[0];
sx q[0];
rz(1.0642306) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1942062) q[2];
sx q[2];
rz(-2.1731659) q[2];
sx q[2];
rz(-2.1927689) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6577158) q[1];
sx q[1];
rz(-0.88798385) q[1];
sx q[1];
rz(-1.8020991) q[1];
x q[2];
rz(-2.6272247) q[3];
sx q[3];
rz(-1.4050438) q[3];
sx q[3];
rz(-0.35781281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94704023) q[2];
sx q[2];
rz(-1.9779454) q[2];
sx q[2];
rz(-2.2770503) q[2];
rz(0.38327992) q[3];
sx q[3];
rz(-0.42208233) q[3];
sx q[3];
rz(-0.88467902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8466723) q[0];
sx q[0];
rz(-2.8546794) q[0];
sx q[0];
rz(0.82557803) q[0];
rz(0.83746743) q[1];
sx q[1];
rz(-1.0191963) q[1];
sx q[1];
rz(2.3501863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90368962) q[0];
sx q[0];
rz(-0.76399732) q[0];
sx q[0];
rz(1.4542411) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2242658) q[2];
sx q[2];
rz(-2.978108) q[2];
sx q[2];
rz(-1.2965073) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36230642) q[1];
sx q[1];
rz(-2.0842881) q[1];
sx q[1];
rz(-1.1720285) q[1];
x q[2];
rz(1.0867001) q[3];
sx q[3];
rz(-1.181385) q[3];
sx q[3];
rz(-2.2661346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5241663) q[2];
sx q[2];
rz(-0.41773057) q[2];
sx q[2];
rz(1.8903271) q[2];
rz(1.1795801) q[3];
sx q[3];
rz(-1.1359295) q[3];
sx q[3];
rz(0.65666667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11611045) q[0];
sx q[0];
rz(-1.2283093) q[0];
sx q[0];
rz(-0.81892282) q[0];
rz(0.081309155) q[1];
sx q[1];
rz(-2.5550877) q[1];
sx q[1];
rz(2.3611045) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7538502) q[0];
sx q[0];
rz(-1.6293793) q[0];
sx q[0];
rz(0.37466074) q[0];
rz(-pi) q[1];
rz(-1.444818) q[2];
sx q[2];
rz(-0.65444817) q[2];
sx q[2];
rz(-1.0140527) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7794521) q[1];
sx q[1];
rz(-1.929053) q[1];
sx q[1];
rz(1.254093) q[1];
rz(-pi) q[2];
rz(-1.9742161) q[3];
sx q[3];
rz(-1.2781004) q[3];
sx q[3];
rz(-0.044805275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2180194) q[2];
sx q[2];
rz(-1.4998481) q[2];
sx q[2];
rz(-1.8531331) q[2];
rz(1.6740547) q[3];
sx q[3];
rz(-2.5947184) q[3];
sx q[3];
rz(2.7047777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9913637) q[0];
sx q[0];
rz(-2.4910091) q[0];
sx q[0];
rz(2.9184166) q[0];
rz(-3.0617833) q[1];
sx q[1];
rz(-2.1640919) q[1];
sx q[1];
rz(-0.83597437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8337934) q[0];
sx q[0];
rz(-0.5029486) q[0];
sx q[0];
rz(-1.0566684) q[0];
rz(1.9718593) q[2];
sx q[2];
rz(-1.1163841) q[2];
sx q[2];
rz(-1.5961702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1389991) q[1];
sx q[1];
rz(-1.9789322) q[1];
sx q[1];
rz(2.1311174) q[1];
rz(-pi) q[2];
rz(-1.9014408) q[3];
sx q[3];
rz(-0.6370987) q[3];
sx q[3];
rz(1.482687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9072546) q[2];
sx q[2];
rz(-2.2865488) q[2];
sx q[2];
rz(-1.3912531) q[2];
rz(-1.2587345) q[3];
sx q[3];
rz(-1.3842868) q[3];
sx q[3];
rz(0.39906991) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2788972) q[0];
sx q[0];
rz(-0.47877043) q[0];
sx q[0];
rz(-0.71267772) q[0];
rz(1.124294) q[1];
sx q[1];
rz(-1.5713888) q[1];
sx q[1];
rz(-2.1052776) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9534822) q[0];
sx q[0];
rz(-1.6878016) q[0];
sx q[0];
rz(-2.5036158) q[0];
x q[1];
rz(1.8646303) q[2];
sx q[2];
rz(-1.7987408) q[2];
sx q[2];
rz(1.0874401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6129174) q[1];
sx q[1];
rz(-2.5096748) q[1];
sx q[1];
rz(1.9252434) q[1];
rz(1.7927945) q[3];
sx q[3];
rz(-2.432468) q[3];
sx q[3];
rz(0.17418451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0599646) q[2];
sx q[2];
rz(-2.5689503) q[2];
sx q[2];
rz(2.6247978) q[2];
rz(-1.146727) q[3];
sx q[3];
rz(-2.1490993) q[3];
sx q[3];
rz(-0.049588047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0741172) q[0];
sx q[0];
rz(-2.4570486) q[0];
sx q[0];
rz(-1.462498) q[0];
rz(-1.0985724) q[1];
sx q[1];
rz(-1.699828) q[1];
sx q[1];
rz(2.3625653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.532004) q[0];
sx q[0];
rz(-1.6250984) q[0];
sx q[0];
rz(1.5724584) q[0];
rz(-pi) q[1];
rz(-0.46709664) q[2];
sx q[2];
rz(-2.0990685) q[2];
sx q[2];
rz(-1.686765) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15493079) q[1];
sx q[1];
rz(-2.8720461) q[1];
sx q[1];
rz(0.69828548) q[1];
rz(-2.2813517) q[3];
sx q[3];
rz(-1.8651984) q[3];
sx q[3];
rz(1.8301726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3940008) q[2];
sx q[2];
rz(-2.438075) q[2];
sx q[2];
rz(-0.35959378) q[2];
rz(0.30073419) q[3];
sx q[3];
rz(-2.3263003) q[3];
sx q[3];
rz(0.99641478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71036285) q[0];
sx q[0];
rz(-1.2082986) q[0];
sx q[0];
rz(0.30074686) q[0];
rz(0.53384471) q[1];
sx q[1];
rz(-1.0736991) q[1];
sx q[1];
rz(3.0278382) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9917002) q[0];
sx q[0];
rz(-2.0684557) q[0];
sx q[0];
rz(-2.0927621) q[0];
rz(2.5328296) q[2];
sx q[2];
rz(-1.5396313) q[2];
sx q[2];
rz(-1.4733914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99234527) q[1];
sx q[1];
rz(-2.2015044) q[1];
sx q[1];
rz(-1.0058606) q[1];
x q[2];
rz(-1.9451009) q[3];
sx q[3];
rz(-1.4153719) q[3];
sx q[3];
rz(0.96159305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45812312) q[2];
sx q[2];
rz(-0.14294954) q[2];
sx q[2];
rz(0.45041034) q[2];
rz(2.2266375) q[3];
sx q[3];
rz(-0.92792845) q[3];
sx q[3];
rz(-1.3373059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8704855) q[0];
sx q[0];
rz(-0.73942375) q[0];
sx q[0];
rz(-2.7258605) q[0];
rz(0.43139002) q[1];
sx q[1];
rz(-0.58512551) q[1];
sx q[1];
rz(-2.364667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5901075) q[0];
sx q[0];
rz(-3.0868106) q[0];
sx q[0];
rz(1.807674) q[0];
x q[1];
rz(2.3007352) q[2];
sx q[2];
rz(-1.0094202) q[2];
sx q[2];
rz(-0.64602588) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5239656) q[1];
sx q[1];
rz(-2.1212148) q[1];
sx q[1];
rz(-1.074076) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.611972) q[3];
sx q[3];
rz(-0.52859113) q[3];
sx q[3];
rz(1.9749255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0100157) q[2];
sx q[2];
rz(-2.3516529) q[2];
sx q[2];
rz(1.5124403) q[2];
rz(-2.7029964) q[3];
sx q[3];
rz(-0.83493835) q[3];
sx q[3];
rz(-2.6252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877391) q[0];
sx q[0];
rz(-0.53936154) q[0];
sx q[0];
rz(-1.3158276) q[0];
rz(3.059803) q[1];
sx q[1];
rz(-1.1809228) q[1];
sx q[1];
rz(1.8705503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5499175) q[0];
sx q[0];
rz(-1.5999874) q[0];
sx q[0];
rz(-1.1468588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6432297) q[2];
sx q[2];
rz(-0.73144215) q[2];
sx q[2];
rz(0.7805025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76190799) q[1];
sx q[1];
rz(-0.81173249) q[1];
sx q[1];
rz(1.9736675) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2972352) q[3];
sx q[3];
rz(-3.0561746) q[3];
sx q[3];
rz(-2.2568995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79909331) q[2];
sx q[2];
rz(-1.7010331) q[2];
sx q[2];
rz(-1.052617) q[2];
rz(-1.8140225) q[3];
sx q[3];
rz(-1.9460287) q[3];
sx q[3];
rz(-2.6322406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256182) q[0];
sx q[0];
rz(-1.7617891) q[0];
sx q[0];
rz(-1.2431385) q[0];
rz(1.7983408) q[1];
sx q[1];
rz(-0.46011283) q[1];
sx q[1];
rz(-2.7759001) q[1];
rz(-2.0475564) q[2];
sx q[2];
rz(-1.7443716) q[2];
sx q[2];
rz(-2.3732408) q[2];
rz(0.41158493) q[3];
sx q[3];
rz(-1.5586018) q[3];
sx q[3];
rz(2.3476521) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
