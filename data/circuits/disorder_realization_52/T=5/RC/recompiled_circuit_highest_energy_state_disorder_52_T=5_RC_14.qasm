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
rz(1.7088543) q[0];
sx q[0];
rz(6.611293) q[0];
sx q[0];
rz(10.264504) q[0];
rz(1.4108763) q[1];
sx q[1];
rz(-1.5411935) q[1];
sx q[1];
rz(0.76959258) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1043848) q[0];
sx q[0];
rz(-1.6055371) q[0];
sx q[0];
rz(-1.6503235) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.316759) q[2];
sx q[2];
rz(-1.5519594) q[2];
sx q[2];
rz(-2.5530262) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.283764) q[1];
sx q[1];
rz(-2.789838) q[1];
sx q[1];
rz(1.5689538) q[1];
rz(-0.64842865) q[3];
sx q[3];
rz(-0.75639137) q[3];
sx q[3];
rz(-1.3677011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8771693) q[2];
sx q[2];
rz(-1.8921655) q[2];
sx q[2];
rz(-1.1733615) q[2];
rz(-0.42826432) q[3];
sx q[3];
rz(-2.5731125) q[3];
sx q[3];
rz(-0.65304023) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1332755) q[0];
sx q[0];
rz(-0.44047099) q[0];
sx q[0];
rz(1.0622729) q[0];
rz(1.5509037) q[1];
sx q[1];
rz(-0.56113243) q[1];
sx q[1];
rz(1.0947469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0355201) q[0];
sx q[0];
rz(-2.8815334) q[0];
sx q[0];
rz(-1.3706657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44610141) q[2];
sx q[2];
rz(-0.94233905) q[2];
sx q[2];
rz(0.43804178) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1909263) q[1];
sx q[1];
rz(-0.48787531) q[1];
sx q[1];
rz(2.4937954) q[1];
rz(-pi) q[2];
rz(2.9238104) q[3];
sx q[3];
rz(-2.3541303) q[3];
sx q[3];
rz(-1.9290719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67794472) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(-2.2410683) q[2];
rz(-2.1053704) q[3];
sx q[3];
rz(-3.0039054) q[3];
sx q[3];
rz(0.010995939) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66505945) q[0];
sx q[0];
rz(-2.0158975) q[0];
sx q[0];
rz(1.6735459) q[0];
rz(-0.92848575) q[1];
sx q[1];
rz(-1.0037582) q[1];
sx q[1];
rz(-1.7195864) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4741459) q[0];
sx q[0];
rz(-2.8774539) q[0];
sx q[0];
rz(0.35281065) q[0];
rz(-pi) q[1];
rz(2.3238238) q[2];
sx q[2];
rz(-1.0330794) q[2];
sx q[2];
rz(1.0420711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0735237) q[1];
sx q[1];
rz(-1.0367107) q[1];
sx q[1];
rz(1.1492386) q[1];
rz(-2.5676651) q[3];
sx q[3];
rz(-2.6477154) q[3];
sx q[3];
rz(-0.18791325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15098393) q[2];
sx q[2];
rz(-1.432212) q[2];
sx q[2];
rz(-0.81986156) q[2];
rz(3.0153583) q[3];
sx q[3];
rz(-1.9641967) q[3];
sx q[3];
rz(-1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4594629) q[0];
sx q[0];
rz(-1.7403025) q[0];
sx q[0];
rz(1.5392186) q[0];
rz(-1.0914717) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.9713255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71890408) q[0];
sx q[0];
rz(-2.6391374) q[0];
sx q[0];
rz(-0.86645856) q[0];
rz(-pi) q[1];
x q[1];
rz(1.865388) q[2];
sx q[2];
rz(-0.74199332) q[2];
sx q[2];
rz(2.2592253) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2728426) q[1];
sx q[1];
rz(-3.0210872) q[1];
sx q[1];
rz(2.4122448) q[1];
x q[2];
rz(-1.7863196) q[3];
sx q[3];
rz(-1.2053262) q[3];
sx q[3];
rz(-0.48367031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9200661) q[2];
sx q[2];
rz(-2.2130794) q[2];
sx q[2];
rz(0.083219223) q[2];
rz(0.65256882) q[3];
sx q[3];
rz(-3.1196399) q[3];
sx q[3];
rz(0.86161247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49587747) q[0];
sx q[0];
rz(-0.28384122) q[0];
sx q[0];
rz(0.19677095) q[0];
rz(1.1605877) q[1];
sx q[1];
rz(-2.6996758) q[1];
sx q[1];
rz(-1.388185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8854839) q[0];
sx q[0];
rz(-1.5692595) q[0];
sx q[0];
rz(2.9716757) q[0];
rz(-2.9391727) q[2];
sx q[2];
rz(-1.7192063) q[2];
sx q[2];
rz(2.9925516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95294556) q[1];
sx q[1];
rz(-1.7428696) q[1];
sx q[1];
rz(0.91305542) q[1];
rz(-pi) q[2];
rz(-0.36404534) q[3];
sx q[3];
rz(-0.87556404) q[3];
sx q[3];
rz(-2.1882868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5768726) q[2];
sx q[2];
rz(-1.2083961) q[2];
sx q[2];
rz(-1.2811071) q[2];
rz(2.7992904) q[3];
sx q[3];
rz(-0.050693158) q[3];
sx q[3];
rz(-2.2127693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.36011919) q[0];
sx q[0];
rz(-1.0973278) q[0];
sx q[0];
rz(2.1088364) q[0];
rz(2.4646941) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(2.3597609) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9572409) q[0];
sx q[0];
rz(-2.4515984) q[0];
sx q[0];
rz(2.5646567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.135649) q[2];
sx q[2];
rz(-3.0810591) q[2];
sx q[2];
rz(-0.50598991) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1429174) q[1];
sx q[1];
rz(-1.8774722) q[1];
sx q[1];
rz(-0.59345133) q[1];
x q[2];
rz(2.8036814) q[3];
sx q[3];
rz(-1.6114144) q[3];
sx q[3];
rz(-0.13536832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8402164) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(0.8737348) q[2];
rz(2.3524763) q[3];
sx q[3];
rz(-1.5566166) q[3];
sx q[3];
rz(0.48621392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18272045) q[0];
sx q[0];
rz(-3.0954269) q[0];
sx q[0];
rz(3.0507372) q[0];
rz(0.60984045) q[1];
sx q[1];
rz(-1.5515386) q[1];
sx q[1];
rz(-0.59744936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2029363) q[0];
sx q[0];
rz(-0.19492004) q[0];
sx q[0];
rz(-0.7233497) q[0];
rz(-pi) q[1];
rz(2.2944974) q[2];
sx q[2];
rz(-0.48111808) q[2];
sx q[2];
rz(1.6498977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3143464) q[1];
sx q[1];
rz(-2.0318077) q[1];
sx q[1];
rz(-0.48996144) q[1];
rz(-pi) q[2];
rz(0.037887033) q[3];
sx q[3];
rz(-1.1514613) q[3];
sx q[3];
rz(1.9333206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37683836) q[2];
sx q[2];
rz(-0.54543269) q[2];
sx q[2];
rz(0.1405912) q[2];
rz(1.8600672) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(-2.6144821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22290467) q[0];
sx q[0];
rz(-1.0162901) q[0];
sx q[0];
rz(-1.5284398) q[0];
rz(-0.81360045) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(-0.80748564) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4190237) q[0];
sx q[0];
rz(-1.3815855) q[0];
sx q[0];
rz(-1.7721227) q[0];
x q[1];
rz(-0.1296223) q[2];
sx q[2];
rz(-2.0106914) q[2];
sx q[2];
rz(0.18378809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2827816) q[1];
sx q[1];
rz(-2.3859897) q[1];
sx q[1];
rz(-2.1367461) q[1];
rz(0.92370478) q[3];
sx q[3];
rz(-2.4974303) q[3];
sx q[3];
rz(1.3180863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7878824) q[2];
sx q[2];
rz(-1.7870125) q[2];
sx q[2];
rz(2.9969969) q[2];
rz(-2.0374129) q[3];
sx q[3];
rz(-2.927533) q[3];
sx q[3];
rz(-2.1000699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60177326) q[0];
sx q[0];
rz(-2.0105392) q[0];
sx q[0];
rz(-2.5850776) q[0];
rz(-0.76983184) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(0.46129033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2568478) q[0];
sx q[0];
rz(-1.0007273) q[0];
sx q[0];
rz(-3.0850201) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3084505) q[2];
sx q[2];
rz(-1.034948) q[2];
sx q[2];
rz(-1.3539202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76595054) q[1];
sx q[1];
rz(-1.8643171) q[1];
sx q[1];
rz(-2.3883345) q[1];
rz(-1.1848524) q[3];
sx q[3];
rz(-1.0331153) q[3];
sx q[3];
rz(-0.73962921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43872675) q[2];
sx q[2];
rz(-2.8361969) q[2];
sx q[2];
rz(-2.3766282) q[2];
rz(1.9276098) q[3];
sx q[3];
rz(-0.58911222) q[3];
sx q[3];
rz(0.37863076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7162489) q[0];
sx q[0];
rz(-3.0916164) q[0];
sx q[0];
rz(-2.7783527) q[0];
rz(-1.4092457) q[1];
sx q[1];
rz(-0.98594085) q[1];
sx q[1];
rz(-0.22663103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5674681) q[0];
sx q[0];
rz(-0.13524817) q[0];
sx q[0];
rz(3.1194206) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4754074) q[2];
sx q[2];
rz(-0.86951423) q[2];
sx q[2];
rz(-1.9124857) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5780826) q[1];
sx q[1];
rz(-0.77111608) q[1];
sx q[1];
rz(-1.095849) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8037386) q[3];
sx q[3];
rz(-1.5808592) q[3];
sx q[3];
rz(-2.4110766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54214415) q[2];
sx q[2];
rz(-2.6207974) q[2];
sx q[2];
rz(-2.4929872) q[2];
rz(-0.058622807) q[3];
sx q[3];
rz(-0.62871814) q[3];
sx q[3];
rz(-2.4962943) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7764353) q[0];
sx q[0];
rz(-1.9226274) q[0];
sx q[0];
rz(2.6489039) q[0];
rz(1.0538712) q[1];
sx q[1];
rz(-2.5851879) q[1];
sx q[1];
rz(-2.7256706) q[1];
rz(0.024439288) q[2];
sx q[2];
rz(-1.8069488) q[2];
sx q[2];
rz(-1.9843742) q[2];
rz(1.966393) q[3];
sx q[3];
rz(-0.7351055) q[3];
sx q[3];
rz(3.0529589) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
