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
rz(-2.2313843) q[0];
sx q[0];
rz(-2.3910523) q[0];
sx q[0];
rz(2.5810177) q[0];
rz(1.9380467) q[1];
sx q[1];
rz(-2.7653341) q[1];
sx q[1];
rz(-0.19317746) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28130075) q[0];
sx q[0];
rz(-2.0112462) q[0];
sx q[0];
rz(-2.4971636) q[0];
x q[1];
rz(3.019252) q[2];
sx q[2];
rz(-2.2723778) q[2];
sx q[2];
rz(2.0035494) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3456125) q[1];
sx q[1];
rz(-0.66540816) q[1];
sx q[1];
rz(-1.3434975) q[1];
rz(-pi) q[2];
rz(1.3985996) q[3];
sx q[3];
rz(-1.2463257) q[3];
sx q[3];
rz(2.8419122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91118497) q[2];
sx q[2];
rz(-2.513803) q[2];
sx q[2];
rz(2.5837303) q[2];
rz(-1.9134391) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(3.0465928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2896344) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(1.5648382) q[0];
rz(1.9948888) q[1];
sx q[1];
rz(-1.4253989) q[1];
sx q[1];
rz(2.1099405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0165591) q[0];
sx q[0];
rz(-1.3368589) q[0];
sx q[0];
rz(2.1159322) q[0];
rz(-pi) q[1];
rz(2.7286058) q[2];
sx q[2];
rz(-1.7888165) q[2];
sx q[2];
rz(-0.08139164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7647142) q[1];
sx q[1];
rz(-2.2441314) q[1];
sx q[1];
rz(0.44133472) q[1];
rz(-pi) q[2];
rz(0.50526039) q[3];
sx q[3];
rz(-1.7499515) q[3];
sx q[3];
rz(3.1324838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8138294) q[2];
sx q[2];
rz(-2.8507865) q[2];
sx q[2];
rz(-1.4168868) q[2];
rz(2.318577) q[3];
sx q[3];
rz(-1.6133285) q[3];
sx q[3];
rz(2.4970162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7773892) q[0];
sx q[0];
rz(-1.1193898) q[0];
sx q[0];
rz(-2.7893344) q[0];
rz(-2.7525355) q[1];
sx q[1];
rz(-1.9720826) q[1];
sx q[1];
rz(2.9152117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.987623) q[0];
sx q[0];
rz(-2.8806318) q[0];
sx q[0];
rz(1.8960192) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7938114) q[2];
sx q[2];
rz(-2.116733) q[2];
sx q[2];
rz(0.40218654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55646642) q[1];
sx q[1];
rz(-2.0279998) q[1];
sx q[1];
rz(-2.6207506) q[1];
x q[2];
rz(2.705615) q[3];
sx q[3];
rz(-1.6395901) q[3];
sx q[3];
rz(1.1747557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9809197) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(2.7799907) q[2];
rz(-2.8412039) q[3];
sx q[3];
rz(-1.3114248) q[3];
sx q[3];
rz(2.1786407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2196197) q[0];
sx q[0];
rz(-2.0491056) q[0];
sx q[0];
rz(-0.12119448) q[0];
rz(1.9425862) q[1];
sx q[1];
rz(-1.0795178) q[1];
sx q[1];
rz(1.2746864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.629144) q[0];
sx q[0];
rz(-1.3487909) q[0];
sx q[0];
rz(1.52284) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5869467) q[2];
sx q[2];
rz(-1.3653737) q[2];
sx q[2];
rz(-1.1753163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53654799) q[1];
sx q[1];
rz(-1.0031317) q[1];
sx q[1];
rz(-1.8325388) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.995756) q[3];
sx q[3];
rz(-1.9327628) q[3];
sx q[3];
rz(1.2545619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.092209665) q[2];
sx q[2];
rz(-1.4468687) q[2];
sx q[2];
rz(1.764074) q[2];
rz(-1.1285909) q[3];
sx q[3];
rz(-2.1994574) q[3];
sx q[3];
rz(-0.82730627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89151299) q[0];
sx q[0];
rz(-2.2703607) q[0];
sx q[0];
rz(-1.0846035) q[0];
rz(-2.4705823) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(-0.074140851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4808124) q[0];
sx q[0];
rz(-0.63711053) q[0];
sx q[0];
rz(-2.6694032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83168516) q[2];
sx q[2];
rz(-2.7234969) q[2];
sx q[2];
rz(-1.8036133) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74846327) q[1];
sx q[1];
rz(-0.65915758) q[1];
sx q[1];
rz(2.0214969) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0228593) q[3];
sx q[3];
rz(-0.30395711) q[3];
sx q[3];
rz(0.077465103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.025617754) q[2];
sx q[2];
rz(-2.3927549) q[2];
sx q[2];
rz(-2.1785114) q[2];
rz(1.3747619) q[3];
sx q[3];
rz(-2.0120967) q[3];
sx q[3];
rz(-0.53205427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9328203) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(1.2579086) q[0];
rz(0.84016291) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(-2.449583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0170572) q[0];
sx q[0];
rz(-1.3003674) q[0];
sx q[0];
rz(-1.9084683) q[0];
rz(-0.18628405) q[2];
sx q[2];
rz(-0.77464235) q[2];
sx q[2];
rz(-1.8587405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51541182) q[1];
sx q[1];
rz(-2.0686935) q[1];
sx q[1];
rz(0.84705686) q[1];
rz(2.1282084) q[3];
sx q[3];
rz(-1.0104826) q[3];
sx q[3];
rz(-2.2309365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.96364) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(0.5160416) q[2];
rz(-1.7959203) q[3];
sx q[3];
rz(-2.7005152) q[3];
sx q[3];
rz(-1.0015944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8733785) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(2.1116665) q[0];
rz(0.65451199) q[1];
sx q[1];
rz(-1.3348568) q[1];
sx q[1];
rz(-2.9676504) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27979461) q[0];
sx q[0];
rz(-1.0423941) q[0];
sx q[0];
rz(2.2938919) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3340764) q[2];
sx q[2];
rz(-1.1828198) q[2];
sx q[2];
rz(-0.18460759) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6831931) q[1];
sx q[1];
rz(-1.8432861) q[1];
sx q[1];
rz(2.2098455) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6967322) q[3];
sx q[3];
rz(-2.5813817) q[3];
sx q[3];
rz(-1.0154533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30764636) q[2];
sx q[2];
rz(-2.6070194) q[2];
sx q[2];
rz(1.8024811) q[2];
rz(0.80859679) q[3];
sx q[3];
rz(-1.9924889) q[3];
sx q[3];
rz(-1.8554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22063743) q[0];
sx q[0];
rz(-2.0957102) q[0];
sx q[0];
rz(1.165423) q[0];
rz(-0.19365817) q[1];
sx q[1];
rz(-2.3083189) q[1];
sx q[1];
rz(1.9906893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26955596) q[0];
sx q[0];
rz(-1.8241166) q[0];
sx q[0];
rz(-2.4756858) q[0];
x q[1];
rz(-2.3370158) q[2];
sx q[2];
rz(-1.6940306) q[2];
sx q[2];
rz(0.096913902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.833646) q[1];
sx q[1];
rz(-0.6816136) q[1];
sx q[1];
rz(2.2824085) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70827534) q[3];
sx q[3];
rz(-2.357956) q[3];
sx q[3];
rz(-2.2275138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2151486) q[2];
sx q[2];
rz(-1.0215267) q[2];
sx q[2];
rz(1.273217) q[2];
rz(2.6269954) q[3];
sx q[3];
rz(-2.8223346) q[3];
sx q[3];
rz(-0.74015051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1505245) q[0];
sx q[0];
rz(-3.0758698) q[0];
sx q[0];
rz(-0.42732987) q[0];
rz(-1.7006251) q[1];
sx q[1];
rz(-1.4375552) q[1];
sx q[1];
rz(2.2429121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60481614) q[0];
sx q[0];
rz(-1.1425848) q[0];
sx q[0];
rz(2.8003119) q[0];
x q[1];
rz(-2.6701791) q[2];
sx q[2];
rz(-2.1080394) q[2];
sx q[2];
rz(1.3211847) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8391621) q[1];
sx q[1];
rz(-1.7453472) q[1];
sx q[1];
rz(-1.3329092) q[1];
rz(0.69435223) q[3];
sx q[3];
rz(-2.4789841) q[3];
sx q[3];
rz(0.41678762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.906189) q[2];
sx q[2];
rz(-1.0455422) q[2];
sx q[2];
rz(1.7380627) q[2];
rz(-2.4710726) q[3];
sx q[3];
rz(-1.5454005) q[3];
sx q[3];
rz(1.8286573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32187605) q[0];
sx q[0];
rz(-0.96915594) q[0];
sx q[0];
rz(0.72625351) q[0];
rz(0.67604524) q[1];
sx q[1];
rz(-1.7773726) q[1];
sx q[1];
rz(-0.77686754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51190864) q[0];
sx q[0];
rz(-2.4998695) q[0];
sx q[0];
rz(1.4944906) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29669478) q[2];
sx q[2];
rz(-2.9888973) q[2];
sx q[2];
rz(-2.8489528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.59275) q[1];
sx q[1];
rz(-1.934086) q[1];
sx q[1];
rz(-2.8123358) q[1];
x q[2];
rz(-2.131072) q[3];
sx q[3];
rz(-2.0891857) q[3];
sx q[3];
rz(0.7553525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0016510222) q[2];
sx q[2];
rz(-0.17639128) q[2];
sx q[2];
rz(1.6528992) q[2];
rz(2.0148924) q[3];
sx q[3];
rz(-1.1688346) q[3];
sx q[3];
rz(-2.6039629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61459944) q[0];
sx q[0];
rz(-2.499883) q[0];
sx q[0];
rz(-1.5747621) q[0];
rz(-2.3552786) q[1];
sx q[1];
rz(-2.96824) q[1];
sx q[1];
rz(2.7460964) q[1];
rz(-0.6720856) q[2];
sx q[2];
rz(-1.0955878) q[2];
sx q[2];
rz(0.53447117) q[2];
rz(-0.24675225) q[3];
sx q[3];
rz(-2.3507531) q[3];
sx q[3];
rz(1.7274461) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
