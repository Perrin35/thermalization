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
rz(2.0496378) q[0];
sx q[0];
rz(-2.0593934) q[0];
sx q[0];
rz(1.9424633) q[0];
rz(0.24425976) q[1];
sx q[1];
rz(-1.3173988) q[1];
sx q[1];
rz(2.2720845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0713214) q[0];
sx q[0];
rz(-2.3247221) q[0];
sx q[0];
rz(2.8512297) q[0];
rz(0.86816048) q[2];
sx q[2];
rz(-2.5294315) q[2];
sx q[2];
rz(0.56589076) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8972733) q[1];
sx q[1];
rz(-2.2530375) q[1];
sx q[1];
rz(1.2644493) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.364089) q[3];
sx q[3];
rz(-2.1453259) q[3];
sx q[3];
rz(-0.053701775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1797552) q[2];
sx q[2];
rz(-1.9840252) q[2];
sx q[2];
rz(-1.301514) q[2];
rz(0.88422042) q[3];
sx q[3];
rz(-1.0853465) q[3];
sx q[3];
rz(-1.6940544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41534153) q[0];
sx q[0];
rz(-1.2301507) q[0];
sx q[0];
rz(0.05106654) q[0];
rz(-2.6452433) q[1];
sx q[1];
rz(-1.8615078) q[1];
sx q[1];
rz(-3.0275717) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43085264) q[0];
sx q[0];
rz(-0.56409696) q[0];
sx q[0];
rz(0.40379712) q[0];
rz(-pi) q[1];
rz(1.586726) q[2];
sx q[2];
rz(-0.6932879) q[2];
sx q[2];
rz(0.42809286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.025578827) q[1];
sx q[1];
rz(-2.4878628) q[1];
sx q[1];
rz(2.3313001) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3475655) q[3];
sx q[3];
rz(-0.79004254) q[3];
sx q[3];
rz(1.5727596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30782917) q[2];
sx q[2];
rz(-1.2276063) q[2];
sx q[2];
rz(2.7755136) q[2];
rz(2.6858373) q[3];
sx q[3];
rz(-2.3594806) q[3];
sx q[3];
rz(2.8428049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5251821) q[0];
sx q[0];
rz(-2.1822378) q[0];
sx q[0];
rz(0.33732238) q[0];
rz(1.6711383) q[1];
sx q[1];
rz(-0.93177876) q[1];
sx q[1];
rz(-0.09387389) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69031624) q[0];
sx q[0];
rz(-0.85448336) q[0];
sx q[0];
rz(-1.4716427) q[0];
x q[1];
rz(-0.44849918) q[2];
sx q[2];
rz(-1.8828159) q[2];
sx q[2];
rz(1.2879368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17901668) q[1];
sx q[1];
rz(-1.1454795) q[1];
sx q[1];
rz(-0.0089745402) q[1];
rz(-1.8444421) q[3];
sx q[3];
rz(-0.85262596) q[3];
sx q[3];
rz(-2.9349365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1845392) q[2];
sx q[2];
rz(-0.55593714) q[2];
sx q[2];
rz(2.8957193) q[2];
rz(0.18694123) q[3];
sx q[3];
rz(-1.4176466) q[3];
sx q[3];
rz(-0.42417446) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4251855) q[0];
sx q[0];
rz(-2.462429) q[0];
sx q[0];
rz(-2.9492522) q[0];
rz(-1.1771857) q[1];
sx q[1];
rz(-0.81383759) q[1];
sx q[1];
rz(0.40843931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3661097) q[0];
sx q[0];
rz(-2.977109) q[0];
sx q[0];
rz(2.7874095) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2871509) q[2];
sx q[2];
rz(-2.7506094) q[2];
sx q[2];
rz(-1.6765082) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0535897) q[1];
sx q[1];
rz(-1.1704967) q[1];
sx q[1];
rz(2.4884239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8242661) q[3];
sx q[3];
rz(-2.4400024) q[3];
sx q[3];
rz(1.6768368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8583777) q[2];
sx q[2];
rz(-2.242531) q[2];
sx q[2];
rz(2.9646207) q[2];
rz(-0.58961287) q[3];
sx q[3];
rz(-1.6556581) q[3];
sx q[3];
rz(2.1069215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.706834) q[0];
sx q[0];
rz(-2.28573) q[0];
sx q[0];
rz(2.6245497) q[0];
rz(-0.21273461) q[1];
sx q[1];
rz(-2.0566302) q[1];
sx q[1];
rz(2.3092666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9134923) q[0];
sx q[0];
rz(-1.6128165) q[0];
sx q[0];
rz(-1.5786922) q[0];
rz(0.8950142) q[2];
sx q[2];
rz(-0.38215206) q[2];
sx q[2];
rz(-2.7016751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4246108) q[1];
sx q[1];
rz(-1.9321793) q[1];
sx q[1];
rz(-0.34504621) q[1];
rz(-pi) q[2];
rz(2.66327) q[3];
sx q[3];
rz(-1.7167313) q[3];
sx q[3];
rz(-1.137382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0011255) q[2];
sx q[2];
rz(-0.8320063) q[2];
sx q[2];
rz(-2.7531085) q[2];
rz(1.2515986) q[3];
sx q[3];
rz(-1.3585217) q[3];
sx q[3];
rz(1.0640915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0851704) q[0];
sx q[0];
rz(-0.52260411) q[0];
sx q[0];
rz(1.9942459) q[0];
rz(-1.7583678) q[1];
sx q[1];
rz(-1.5478094) q[1];
sx q[1];
rz(0.17471084) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095681079) q[0];
sx q[0];
rz(-2.2457684) q[0];
sx q[0];
rz(1.573091) q[0];
rz(-pi) q[1];
rz(0.60947588) q[2];
sx q[2];
rz(-1.5646439) q[2];
sx q[2];
rz(0.61340082) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0905007) q[1];
sx q[1];
rz(-1.2686172) q[1];
sx q[1];
rz(-1.3626182) q[1];
rz(-pi) q[2];
rz(-2.9851341) q[3];
sx q[3];
rz(-2.9269013) q[3];
sx q[3];
rz(-0.36461339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2061578) q[2];
sx q[2];
rz(-1.9400699) q[2];
sx q[2];
rz(0.60714444) q[2];
rz(-2.855865) q[3];
sx q[3];
rz(-0.61107475) q[3];
sx q[3];
rz(0.40209517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6242591) q[0];
sx q[0];
rz(-1.3579955) q[0];
sx q[0];
rz(-0.72878033) q[0];
rz(-2.0023316) q[1];
sx q[1];
rz(-1.1092721) q[1];
sx q[1];
rz(-1.2713825) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2654289) q[0];
sx q[0];
rz(-0.012152925) q[0];
sx q[0];
rz(-1.7677714) q[0];
rz(-pi) q[1];
rz(-1.2954577) q[2];
sx q[2];
rz(-1.0022114) q[2];
sx q[2];
rz(0.24623539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9460341) q[1];
sx q[1];
rz(-2.5321153) q[1];
sx q[1];
rz(-0.35802623) q[1];
rz(-pi) q[2];
rz(-2.0470601) q[3];
sx q[3];
rz(-1.4753398) q[3];
sx q[3];
rz(2.3396479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83069688) q[2];
sx q[2];
rz(-1.7252012) q[2];
sx q[2];
rz(-2.8719416) q[2];
rz(1.0722748) q[3];
sx q[3];
rz(-2.8282073) q[3];
sx q[3];
rz(1.0543893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6894182) q[0];
sx q[0];
rz(-1.9502689) q[0];
sx q[0];
rz(2.7574975) q[0];
rz(2.0317888) q[1];
sx q[1];
rz(-0.88491076) q[1];
sx q[1];
rz(-1.5193411) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8728646) q[0];
sx q[0];
rz(-0.11188398) q[0];
sx q[0];
rz(3.0203793) q[0];
rz(-pi) q[1];
rz(0.89092358) q[2];
sx q[2];
rz(-1.884084) q[2];
sx q[2];
rz(3.0591722) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1311595) q[1];
sx q[1];
rz(-1.6747867) q[1];
sx q[1];
rz(1.7883744) q[1];
x q[2];
rz(-1.396809) q[3];
sx q[3];
rz(-1.0871776) q[3];
sx q[3];
rz(2.015851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9263837) q[2];
sx q[2];
rz(-1.0101725) q[2];
sx q[2];
rz(2.9877648) q[2];
rz(-2.202863) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(1.7083098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4140103) q[0];
sx q[0];
rz(-1.3879956) q[0];
sx q[0];
rz(-2.599732) q[0];
rz(1.2205203) q[1];
sx q[1];
rz(-2.4668756) q[1];
sx q[1];
rz(0.064124785) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7008931) q[0];
sx q[0];
rz(-1.6055853) q[0];
sx q[0];
rz(-0.056408806) q[0];
rz(0.86875963) q[2];
sx q[2];
rz(-1.6693309) q[2];
sx q[2];
rz(1.9256401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87919533) q[1];
sx q[1];
rz(-0.91934312) q[1];
sx q[1];
rz(1.5919973) q[1];
rz(-pi) q[2];
rz(2.7087791) q[3];
sx q[3];
rz(-2.6901931) q[3];
sx q[3];
rz(1.4572382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8378143) q[2];
sx q[2];
rz(-2.1393675) q[2];
sx q[2];
rz(3.0617867) q[2];
rz(2.5507353) q[3];
sx q[3];
rz(-2.8320524) q[3];
sx q[3];
rz(-0.045056067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2703242) q[0];
sx q[0];
rz(-0.41759434) q[0];
sx q[0];
rz(-0.25750345) q[0];
rz(-2.2389257) q[1];
sx q[1];
rz(-2.3011484) q[1];
sx q[1];
rz(-2.2534175) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42300007) q[0];
sx q[0];
rz(-0.42536727) q[0];
sx q[0];
rz(1.0235538) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4492017) q[2];
sx q[2];
rz(-2.2071725) q[2];
sx q[2];
rz(0.51511918) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1015472) q[1];
sx q[1];
rz(-2.2253621) q[1];
sx q[1];
rz(2.1957869) q[1];
x q[2];
rz(2.2769663) q[3];
sx q[3];
rz(-0.99643222) q[3];
sx q[3];
rz(-2.7780617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73457926) q[2];
sx q[2];
rz(-0.96757704) q[2];
sx q[2];
rz(-2.5261397) q[2];
rz(-0.19221273) q[3];
sx q[3];
rz(-2.2902316) q[3];
sx q[3];
rz(-1.6063469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2007582) q[0];
sx q[0];
rz(-1.2949018) q[0];
sx q[0];
rz(-1.9042263) q[0];
rz(-1.0279961) q[1];
sx q[1];
rz(-2.4528687) q[1];
sx q[1];
rz(-2.5797896) q[1];
rz(-1.4877612) q[2];
sx q[2];
rz(-2.5032264) q[2];
sx q[2];
rz(-1.6089321) q[2];
rz(1.9348223) q[3];
sx q[3];
rz(-1.0586998) q[3];
sx q[3];
rz(2.8431526) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
