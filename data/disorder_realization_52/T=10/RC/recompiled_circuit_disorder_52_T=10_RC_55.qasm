OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(6.8318879) q[0];
sx q[0];
rz(5.3988342) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8774672) q[0];
sx q[0];
rz(-1.2116417) q[0];
sx q[0];
rz(2.9496664) q[0];
x q[1];
rz(0.0083562334) q[2];
sx q[2];
rz(-0.5651606) q[2];
sx q[2];
rz(-1.9985808) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.002418) q[1];
sx q[1];
rz(-0.25484172) q[1];
sx q[1];
rz(-2.8789218) q[1];
rz(-1.1611657) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(-0.15256552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(-1.9338699) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2475964) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(-0.22878376) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996961) q[0];
sx q[0];
rz(-1.6640088) q[0];
sx q[0];
rz(1.9641563) q[0];
rz(2.8480808) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(-0.40700618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.05939535) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(2.5897964) q[1];
rz(-pi) q[2];
rz(-0.47827999) q[3];
sx q[3];
rz(-0.96049958) q[3];
sx q[3];
rz(1.573521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(-0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-2.8799768) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3221489) q[0];
sx q[0];
rz(-1.561164) q[0];
sx q[0];
rz(-1.5887512) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7058271) q[2];
sx q[2];
rz(-0.74880744) q[2];
sx q[2];
rz(2.7111862) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6690327) q[1];
sx q[1];
rz(-1.6918039) q[1];
sx q[1];
rz(2.3585412) q[1];
x q[2];
rz(0.055734169) q[3];
sx q[3];
rz(-2.627943) q[3];
sx q[3];
rz(-1.0518215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6802784) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(-0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.5699566) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.8932231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6682537) q[0];
sx q[0];
rz(-3.0248397) q[0];
sx q[0];
rz(2.1658685) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3893045) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(-2.4455621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8716988) q[1];
sx q[1];
rz(-1.4596241) q[1];
sx q[1];
rz(2.486869) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72327153) q[3];
sx q[3];
rz(-1.372882) q[3];
sx q[3];
rz(-0.42641446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(2.3045585) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(2.2391589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31061253) q[0];
sx q[0];
rz(-1.5251625) q[0];
sx q[0];
rz(-2.5020585) q[0];
rz(-pi) q[1];
rz(0.4544223) q[2];
sx q[2];
rz(-1.7917969) q[2];
sx q[2];
rz(2.1511252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24689281) q[1];
sx q[1];
rz(-1.0122074) q[1];
sx q[1];
rz(-0.65319368) q[1];
x q[2];
rz(-1.8875214) q[3];
sx q[3];
rz(-1.5681019) q[3];
sx q[3];
rz(2.5731784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(-1.442391) q[3];
sx q[3];
rz(-1.8701575) q[3];
sx q[3];
rz(-1.6825914) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8114132) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(-2.6384171) q[0];
rz(1.6852089) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(2.9398289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148894) q[0];
sx q[0];
rz(-3.0198583) q[0];
sx q[0];
rz(2.142971) q[0];
rz(-pi) q[1];
rz(2.0189507) q[2];
sx q[2];
rz(-2.0070224) q[2];
sx q[2];
rz(-3.089038) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4469874) q[1];
sx q[1];
rz(-1.9896549) q[1];
sx q[1];
rz(0.5254196) q[1];
rz(-pi) q[2];
rz(2.6334409) q[3];
sx q[3];
rz(-2.3414632) q[3];
sx q[3];
rz(1.2451764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21489828) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(0.26724896) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-0.94271359) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38521117) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(-2.5368607) q[0];
rz(2.1812181) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(0.20197091) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49589866) q[1];
sx q[1];
rz(-1.8897448) q[1];
sx q[1];
rz(-1.1369399) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25992486) q[3];
sx q[3];
rz(-1.0726895) q[3];
sx q[3];
rz(-2.2181524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(-2.102397) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246178) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-2.8714645) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(2.8576635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547202) q[0];
sx q[0];
rz(-2.3754639) q[0];
sx q[0];
rz(-0.92673577) q[0];
x q[1];
rz(-1.4242886) q[2];
sx q[2];
rz(-2.1657145) q[2];
sx q[2];
rz(2.811424) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64360196) q[1];
sx q[1];
rz(-1.0712578) q[1];
sx q[1];
rz(-1.040578) q[1];
rz(-pi) q[2];
x q[2];
rz(0.038589434) q[3];
sx q[3];
rz(-0.66017294) q[3];
sx q[3];
rz(2.431228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(1.2109057) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.4607666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47855908) q[0];
sx q[0];
rz(-1.2872818) q[0];
sx q[0];
rz(1.5525596) q[0];
rz(-1.3912348) q[2];
sx q[2];
rz(-1.1968687) q[2];
sx q[2];
rz(-2.0342846) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31635346) q[1];
sx q[1];
rz(-1.5484527) q[1];
sx q[1];
rz(-0.52492001) q[1];
rz(-1.1838412) q[3];
sx q[3];
rz(-0.82434067) q[3];
sx q[3];
rz(-1.6798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4298657) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(-2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56363737) q[0];
sx q[0];
rz(-1.6407688) q[0];
sx q[0];
rz(1.3689343) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1516018) q[2];
sx q[2];
rz(-2.2238646) q[2];
sx q[2];
rz(0.68819118) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9951574) q[1];
sx q[1];
rz(-0.53968118) q[1];
sx q[1];
rz(-1.9567009) q[1];
rz(-pi) q[2];
rz(-0.59488876) q[3];
sx q[3];
rz(-0.66685646) q[3];
sx q[3];
rz(0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(1.1364737) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(1.0271172) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(0.74577352) q[2];
sx q[2];
rz(-0.33200982) q[2];
sx q[2];
rz(0.40340323) q[2];
rz(1.1626701) q[3];
sx q[3];
rz(-2.6682731) q[3];
sx q[3];
rz(-0.70262739) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
