OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.409531) q[0];
sx q[0];
rz(-1.3652029) q[0];
sx q[0];
rz(-2.1172297) q[0];
rz(-2.536474) q[1];
sx q[1];
rz(-2.6095698) q[1];
sx q[1];
rz(-1.1693118) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697203) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(-2.0277434) q[0];
x q[1];
rz(2.3697882) q[2];
sx q[2];
rz(-0.82320854) q[2];
sx q[2];
rz(-0.31847218) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6572666) q[1];
sx q[1];
rz(-2.6063759) q[1];
sx q[1];
rz(0.80429299) q[1];
rz(-pi) q[2];
rz(1.4952881) q[3];
sx q[3];
rz(-1.8823874) q[3];
sx q[3];
rz(2.989245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(1.3226091) q[2];
rz(0.30098513) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(-1.3809563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-2.6665376) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-0.95502949) q[1];
sx q[1];
rz(2.1038726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.422721) q[0];
sx q[0];
rz(-1.0178716) q[0];
sx q[0];
rz(-1.3119112) q[0];
rz(-0.53595397) q[2];
sx q[2];
rz(-2.0336656) q[2];
sx q[2];
rz(-3.0218389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71643752) q[1];
sx q[1];
rz(-1.5967224) q[1];
sx q[1];
rz(1.9772711) q[1];
rz(-pi) q[2];
rz(2.328697) q[3];
sx q[3];
rz(-1.1162236) q[3];
sx q[3];
rz(1.2113435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4804046) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-3.0569055) q[2];
rz(2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111074) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-0.57317615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37704913) q[0];
sx q[0];
rz(-1.1408313) q[0];
sx q[0];
rz(0.069805108) q[0];
rz(-pi) q[1];
rz(2.3689752) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(-0.85341838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1307615) q[1];
sx q[1];
rz(-1.2847932) q[1];
sx q[1];
rz(-2.7421013) q[1];
rz(-pi) q[2];
rz(0.8543386) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(-2.4454988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(2.9476681) q[2];
rz(-0.097269416) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(-0.55066806) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(2.7788924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.761128) q[0];
sx q[0];
rz(-1.8109545) q[0];
sx q[0];
rz(-2.2855177) q[0];
rz(-pi) q[1];
rz(3.1231819) q[2];
sx q[2];
rz(-1.5393179) q[2];
sx q[2];
rz(-1.5359495) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4421694) q[1];
sx q[1];
rz(-1.4675092) q[1];
sx q[1];
rz(2.5073754) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0691931) q[3];
sx q[3];
rz(-0.94745938) q[3];
sx q[3];
rz(-2.6423531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4576733) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-3.1385699) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(0.017379934) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(1.4594706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97809726) q[0];
sx q[0];
rz(-0.28143829) q[0];
sx q[0];
rz(2.1241758) q[0];
x q[1];
rz(0.8874138) q[2];
sx q[2];
rz(-1.9878584) q[2];
sx q[2];
rz(0.29287698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0046878) q[1];
sx q[1];
rz(-1.7675752) q[1];
sx q[1];
rz(2.911527) q[1];
rz(-pi) q[2];
rz(-2.0586117) q[3];
sx q[3];
rz(-1.2851464) q[3];
sx q[3];
rz(-2.3398427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.83546272) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-0.84189502) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-1.1152277) q[3];
sx q[3];
rz(-1.5766597) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(0.58419624) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-0.13866436) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16620557) q[0];
sx q[0];
rz(-1.5741325) q[0];
sx q[0];
rz(-0.020394527) q[0];
rz(1.4939098) q[2];
sx q[2];
rz(-1.3704408) q[2];
sx q[2];
rz(2.4424057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2885292) q[1];
sx q[1];
rz(-1.001734) q[1];
sx q[1];
rz(0.89785518) q[1];
x q[2];
rz(-2.3994282) q[3];
sx q[3];
rz(-1.8932749) q[3];
sx q[3];
rz(-0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(-0.095120393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(0.53652525) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-3.0032872) q[1];
sx q[1];
rz(2.5172863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1034575) q[0];
sx q[0];
rz(-2.3893444) q[0];
sx q[0];
rz(-2.8318846) q[0];
rz(0.46632669) q[2];
sx q[2];
rz(-1.000324) q[2];
sx q[2];
rz(-1.8898659) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1075322) q[1];
sx q[1];
rz(-2.1851087) q[1];
sx q[1];
rz(-0.88453102) q[1];
rz(-pi) q[2];
rz(-0.47847139) q[3];
sx q[3];
rz(-2.1395184) q[3];
sx q[3];
rz(-2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5252934) q[2];
sx q[2];
rz(-1.1181744) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.265825) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.6775999) q[1];
sx q[1];
rz(-2.174607) q[1];
sx q[1];
rz(0.12891842) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3710204) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(-1.8380941) q[0];
rz(-pi) q[1];
rz(0.89959896) q[2];
sx q[2];
rz(-2.0313782) q[2];
sx q[2];
rz(-2.7923498) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32313777) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(-2.6998991) q[1];
rz(-pi) q[2];
x q[2];
rz(1.586732) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(-2.9908138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.504618) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-0.23877731) q[3];
sx q[3];
rz(-1.5978565) q[3];
sx q[3];
rz(-2.074266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(-1.2497485) q[0];
rz(-3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.3508266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46247813) q[0];
sx q[0];
rz(-0.20972855) q[0];
sx q[0];
rz(-2.1049343) q[0];
x q[1];
rz(2.5647854) q[2];
sx q[2];
rz(-1.9153321) q[2];
sx q[2];
rz(0.21201269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1180601) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(-0.94783028) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25165598) q[3];
sx q[3];
rz(-2.3725384) q[3];
sx q[3];
rz(0.8151527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(-1.8814686) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(1.982622) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262912) q[0];
sx q[0];
rz(-1.6537332) q[0];
sx q[0];
rz(0.21324883) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-1.001819) q[1];
sx q[1];
rz(2.5949123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71781681) q[0];
sx q[0];
rz(-2.2795838) q[0];
sx q[0];
rz(0.29341673) q[0];
rz(-pi) q[1];
rz(0.47537739) q[2];
sx q[2];
rz(-1.7065085) q[2];
sx q[2];
rz(-3.140608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.039779546) q[1];
sx q[1];
rz(-1.2817849) q[1];
sx q[1];
rz(-1.5323557) q[1];
x q[2];
rz(-1.7095079) q[3];
sx q[3];
rz(-1.3472392) q[3];
sx q[3];
rz(-2.9885938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(-3.1372916) q[2];
rz(-2.1440078) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988408) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(2.6976363) q[1];
sx q[1];
rz(-0.28354357) q[1];
sx q[1];
rz(1.2734738) q[1];
rz(-1.8757204) q[2];
sx q[2];
rz(-0.049449895) q[2];
sx q[2];
rz(0.54686875) q[2];
rz(2.8125433) q[3];
sx q[3];
rz(-2.0206738) q[3];
sx q[3];
rz(-2.1856367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
