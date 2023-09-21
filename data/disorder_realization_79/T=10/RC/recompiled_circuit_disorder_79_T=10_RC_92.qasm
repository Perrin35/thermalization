OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8389991) q[0];
sx q[0];
rz(-1.5035045) q[0];
sx q[0];
rz(1.5124613) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2201266) q[2];
sx q[2];
rz(-1.7979243) q[2];
sx q[2];
rz(2.7820058) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9477387) q[1];
sx q[1];
rz(-0.96809909) q[1];
sx q[1];
rz(1.0268474) q[1];
rz(-0.37791032) q[3];
sx q[3];
rz(-2.1034735) q[3];
sx q[3];
rz(2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(-0.57693276) q[1];
sx q[1];
rz(-0.88795841) q[1];
sx q[1];
rz(1.4651441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11008308) q[0];
sx q[0];
rz(-0.84658909) q[0];
sx q[0];
rz(1.5741041) q[0];
x q[1];
rz(-0.25603489) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(-1.408996) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2966753) q[1];
sx q[1];
rz(-1.0122609) q[1];
sx q[1];
rz(-2.1293473) q[1];
rz(-2.39605) q[3];
sx q[3];
rz(-1.3127483) q[3];
sx q[3];
rz(0.52449709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.96317545) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(2.3410472) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24414177) q[0];
sx q[0];
rz(-1.4276917) q[0];
sx q[0];
rz(1.98154) q[0];
x q[1];
rz(-2.3163047) q[2];
sx q[2];
rz(-0.85916677) q[2];
sx q[2];
rz(-0.7691783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2785981) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(1.5652656) q[1];
rz(-pi) q[2];
x q[2];
rz(1.21739) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1151428) q[0];
sx q[0];
rz(-1.0466252) q[0];
sx q[0];
rz(-1.3036222) q[0];
x q[1];
rz(-0.24809804) q[2];
sx q[2];
rz(-1.2144226) q[2];
sx q[2];
rz(-2.6756289) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3738721) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(2.9994681) q[1];
x q[2];
rz(1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(-0.15743263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.408668) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(0.98168215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0504793) q[0];
sx q[0];
rz(-2.1149153) q[0];
sx q[0];
rz(-1.0247466) q[0];
rz(0.9972516) q[2];
sx q[2];
rz(-2.855636) q[2];
sx q[2];
rz(-0.85948247) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2213858) q[1];
sx q[1];
rz(-1.1077987) q[1];
sx q[1];
rz(-1.7358001) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4622299) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(1.1443646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.813252) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-2.7691675) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(-2.9763124) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80711354) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(1.6001742) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68768244) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(-2.8962367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0300671) q[1];
sx q[1];
rz(-2.1038052) q[1];
sx q[1];
rz(-1.0982606) q[1];
rz(1.1092471) q[3];
sx q[3];
rz(-2.647532) q[3];
sx q[3];
rz(-2.9500614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(0.73648891) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7613477) q[0];
sx q[0];
rz(-1.5878829) q[0];
sx q[0];
rz(-1.6197617) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1963504) q[2];
sx q[2];
rz(-0.40996273) q[2];
sx q[2];
rz(1.9117102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.10662096) q[1];
sx q[1];
rz(-1.391341) q[1];
sx q[1];
rz(-0.76645318) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.001708) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(-1.0678837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(-2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078995973) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661134) q[0];
sx q[0];
rz(-2.4371494) q[0];
sx q[0];
rz(-0.2091614) q[0];
rz(-1.2577031) q[2];
sx q[2];
rz(-2.3790092) q[2];
sx q[2];
rz(-1.8956172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39107716) q[1];
sx q[1];
rz(-2.8120496) q[1];
sx q[1];
rz(-2.4604172) q[1];
rz(-0.17955762) q[3];
sx q[3];
rz(-0.87053821) q[3];
sx q[3];
rz(1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(2.4618861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98201671) q[0];
sx q[0];
rz(-1.7734818) q[0];
sx q[0];
rz(-0.52635877) q[0];
rz(-1.3557415) q[2];
sx q[2];
rz(-2.4348223) q[2];
sx q[2];
rz(0.099345318) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.745984) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(2.4376274) q[1];
x q[2];
rz(2.1867832) q[3];
sx q[3];
rz(-0.98815742) q[3];
sx q[3];
rz(2.5481176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(0.91040197) q[2];
rz(0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(0.7146548) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.1766599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991853) q[0];
sx q[0];
rz(-2.9710794) q[0];
sx q[0];
rz(1.8605581) q[0];
rz(-2.5942957) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(-1.0123569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1365876) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(0.088076061) q[1];
rz(1.482974) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(-2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(0.73666112) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(-1.5626004) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];