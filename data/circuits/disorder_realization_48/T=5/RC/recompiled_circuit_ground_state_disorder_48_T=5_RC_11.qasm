OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(-2.7645219) q[1];
sx q[1];
rz(-1.3426251) q[1];
sx q[1];
rz(-1.0599729) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19883991) q[0];
sx q[0];
rz(-1.6918385) q[0];
sx q[0];
rz(3.0225656) q[0];
x q[1];
rz(2.4490812) q[2];
sx q[2];
rz(-1.5438809) q[2];
sx q[2];
rz(-1.7153502) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6457451) q[1];
sx q[1];
rz(-1.8105257) q[1];
sx q[1];
rz(1.354753) q[1];
x q[2];
rz(-0.59935244) q[3];
sx q[3];
rz(-0.17440344) q[3];
sx q[3];
rz(2.6389183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9609191) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(-1.712435) q[2];
rz(-0.6116496) q[3];
sx q[3];
rz(-1.6983756) q[3];
sx q[3];
rz(0.071368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83577689) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(-2.3961156) q[0];
rz(1.2019134) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(-2.1048996) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03464493) q[0];
sx q[0];
rz(-1.4781524) q[0];
sx q[0];
rz(1.9361467) q[0];
x q[1];
rz(-1.5971423) q[2];
sx q[2];
rz(-2.3374632) q[2];
sx q[2];
rz(-0.94601099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0889609) q[1];
sx q[1];
rz(-2.3792069) q[1];
sx q[1];
rz(-0.41684581) q[1];
rz(-2.8193982) q[3];
sx q[3];
rz(-2.5327842) q[3];
sx q[3];
rz(0.52458602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1408954) q[2];
sx q[2];
rz(-1.2709728) q[2];
sx q[2];
rz(0.97266436) q[2];
rz(-0.85727143) q[3];
sx q[3];
rz(-2.075115) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5001517) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(-0.80297536) q[0];
rz(-0.650644) q[1];
sx q[1];
rz(-0.79901189) q[1];
sx q[1];
rz(-2.9237936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7729491) q[0];
sx q[0];
rz(-1.1945219) q[0];
sx q[0];
rz(2.2479288) q[0];
rz(-pi) q[1];
rz(-0.62816633) q[2];
sx q[2];
rz(-1.7587314) q[2];
sx q[2];
rz(0.70580259) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4068702) q[1];
sx q[1];
rz(-2.0592505) q[1];
sx q[1];
rz(-3.0367756) q[1];
x q[2];
rz(-0.59813114) q[3];
sx q[3];
rz(-2.7418828) q[3];
sx q[3];
rz(2.3593115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95152068) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(-1.5104843) q[2];
rz(-1.2269646) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(1.0043043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25201061) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(-2.4776283) q[0];
rz(0.12531677) q[1];
sx q[1];
rz(-1.7071416) q[1];
sx q[1];
rz(-1.0391611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9665899) q[0];
sx q[0];
rz(-3.1313854) q[0];
sx q[0];
rz(0.27423476) q[0];
rz(1.9661994) q[2];
sx q[2];
rz(-1.8186453) q[2];
sx q[2];
rz(-1.6958267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3323556) q[1];
sx q[1];
rz(-2.0505822) q[1];
sx q[1];
rz(-0.23608853) q[1];
x q[2];
rz(1.4757206) q[3];
sx q[3];
rz(-2.0119442) q[3];
sx q[3];
rz(0.26218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30109721) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(0.077979716) q[2];
rz(-2.0237427) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(-2.1054721) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2074821) q[0];
sx q[0];
rz(-2.9353862) q[0];
sx q[0];
rz(0.41314405) q[0];
rz(1.235599) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(1.8280169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395555) q[0];
sx q[0];
rz(-1.7540723) q[0];
sx q[0];
rz(2.2048143) q[0];
rz(1.1059472) q[2];
sx q[2];
rz(-1.8951956) q[2];
sx q[2];
rz(1.731763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47258618) q[1];
sx q[1];
rz(-0.60429497) q[1];
sx q[1];
rz(-0.86493203) q[1];
rz(-0.78727874) q[3];
sx q[3];
rz(-1.2378927) q[3];
sx q[3];
rz(-1.5008139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1818587) q[2];
sx q[2];
rz(-1.9071969) q[2];
sx q[2];
rz(-0.44446298) q[2];
rz(-1.2005165) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(3.0357231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926517) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(3.0643903) q[0];
rz(-2.1791747) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(-3.107792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12594906) q[0];
sx q[0];
rz(-1.9140049) q[0];
sx q[0];
rz(2.6604466) q[0];
rz(-0.54767139) q[2];
sx q[2];
rz(-1.593809) q[2];
sx q[2];
rz(0.097824899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2359743) q[1];
sx q[1];
rz(-1.2575358) q[1];
sx q[1];
rz(2.226023) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.46557) q[3];
sx q[3];
rz(-2.3447373) q[3];
sx q[3];
rz(-1.6355255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17322156) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(-2.8823493) q[2];
rz(-0.94414583) q[3];
sx q[3];
rz(-2.9278432) q[3];
sx q[3];
rz(-1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(-2.5282705) q[0];
rz(-0.90336409) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(2.9936252) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44339339) q[0];
sx q[0];
rz(-0.86981378) q[0];
sx q[0];
rz(-2.2156391) q[0];
rz(-pi) q[1];
rz(-2.7025928) q[2];
sx q[2];
rz(-1.5996965) q[2];
sx q[2];
rz(-2.7759107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1175673) q[1];
sx q[1];
rz(-2.2541775) q[1];
sx q[1];
rz(-2.1272117) q[1];
rz(-2.1359753) q[3];
sx q[3];
rz(-1.5907173) q[3];
sx q[3];
rz(-1.7651778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2030187) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(-0.99986783) q[2];
rz(-1.8262919) q[3];
sx q[3];
rz(-2.857693) q[3];
sx q[3];
rz(-0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95241791) q[0];
sx q[0];
rz(-0.84680951) q[0];
sx q[0];
rz(2.1536105) q[0];
rz(-1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1647427) q[0];
sx q[0];
rz(-1.4535507) q[0];
sx q[0];
rz(1.4848723) q[0];
x q[1];
rz(2.0664178) q[2];
sx q[2];
rz(-1.3775423) q[2];
sx q[2];
rz(0.19710625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33262256) q[1];
sx q[1];
rz(-1.5676985) q[1];
sx q[1];
rz(-1.078152) q[1];
rz(-pi) q[2];
rz(2.2211391) q[3];
sx q[3];
rz(-1.4152539) q[3];
sx q[3];
rz(0.33442861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24171955) q[2];
sx q[2];
rz(-2.6417929) q[2];
sx q[2];
rz(1.0618173) q[2];
rz(0.1344943) q[3];
sx q[3];
rz(-0.89956346) q[3];
sx q[3];
rz(-3.0883446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6580842) q[0];
sx q[0];
rz(-2.4299419) q[0];
sx q[0];
rz(0.2051951) q[0];
rz(-2.9070053) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(0.1964143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4171643) q[0];
sx q[0];
rz(-2.1311893) q[0];
sx q[0];
rz(0.36720328) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5242879) q[2];
sx q[2];
rz(-1.7066188) q[2];
sx q[2];
rz(-0.4834396) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3141296) q[1];
sx q[1];
rz(-1.2478831) q[1];
sx q[1];
rz(0.18675487) q[1];
rz(2.6077723) q[3];
sx q[3];
rz(-0.8554503) q[3];
sx q[3];
rz(1.2792339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74497574) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.9632001) q[2];
rz(-2.7922503) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0013393764) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(2.4440675) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(-2.2360905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6385429) q[0];
sx q[0];
rz(-2.7075504) q[0];
sx q[0];
rz(-3.1014355) q[0];
rz(-pi) q[1];
rz(-0.72737965) q[2];
sx q[2];
rz(-0.47107163) q[2];
sx q[2];
rz(-1.7076275) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.057396279) q[1];
sx q[1];
rz(-1.2603972) q[1];
sx q[1];
rz(-1.6204024) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8373413) q[3];
sx q[3];
rz(-2.6681134) q[3];
sx q[3];
rz(2.480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2424348) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(0.71851292) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-0.1040641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.318442) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(1.5009343) q[1];
sx q[1];
rz(-1.2702912) q[1];
sx q[1];
rz(2.695695) q[1];
rz(-2.9983042) q[2];
sx q[2];
rz(-1.454151) q[2];
sx q[2];
rz(3.1232338) q[2];
rz(0.80397687) q[3];
sx q[3];
rz(-1.0729651) q[3];
sx q[3];
rz(3.1333522) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
