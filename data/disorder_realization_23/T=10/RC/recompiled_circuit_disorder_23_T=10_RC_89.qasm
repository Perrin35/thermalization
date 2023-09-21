OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(-1.707466) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(2.0531274) q[0];
x q[1];
rz(1.7069874) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(-1.1510804) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0093706) q[1];
sx q[1];
rz(-1.8144061) q[1];
sx q[1];
rz(2.0396712) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1316142) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-0.83797541) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(-0.078991927) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74719602) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(-1.2778506) q[0];
rz(-2.9648119) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(0.4321672) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5837001) q[0];
sx q[0];
rz(-1.0520792) q[0];
sx q[0];
rz(1.8899263) q[0];
rz(-pi) q[1];
rz(1.1439267) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(-1.2183684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.00694) q[1];
sx q[1];
rz(-1.4538527) q[1];
sx q[1];
rz(-2.0140531) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0057955) q[3];
sx q[3];
rz(-1.5917935) q[3];
sx q[3];
rz(2.3137623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5923578) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040722672) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(2.7242463) q[0];
rz(-1.4886645) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(-0.506385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15784141) q[0];
sx q[0];
rz(-0.39641532) q[0];
sx q[0];
rz(1.0979963) q[0];
rz(-0.73451368) q[2];
sx q[2];
rz(-1.7506415) q[2];
sx q[2];
rz(2.432446) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80458927) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(-0.6188436) q[1];
rz(-1.9403463) q[3];
sx q[3];
rz(-1.4003716) q[3];
sx q[3];
rz(0.91526645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4425519) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(-0.58602035) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-0.83918321) q[0];
rz(0.025578586) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4913113) q[0];
sx q[0];
rz(-0.45931739) q[0];
sx q[0];
rz(-1.7772872) q[0];
rz(-pi) q[1];
rz(0.80231248) q[2];
sx q[2];
rz(-1.017184) q[2];
sx q[2];
rz(-0.72788903) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3241987) q[1];
sx q[1];
rz(-2.3014268) q[1];
sx q[1];
rz(2.6170931) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2876835) q[3];
sx q[3];
rz(-1.9861756) q[3];
sx q[3];
rz(2.0302041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.5754159) q[0];
sx q[0];
rz(-1.7594936) q[0];
sx q[0];
rz(-2.8856522) q[0];
rz(2.6804965) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81239031) q[0];
sx q[0];
rz(-2.985552) q[0];
sx q[0];
rz(-0.72326707) q[0];
x q[1];
rz(1.4270093) q[2];
sx q[2];
rz(-0.41848768) q[2];
sx q[2];
rz(-2.8029122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45436817) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(1.441799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7759336) q[3];
sx q[3];
rz(-2.1846111) q[3];
sx q[3];
rz(2.1706276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53133416) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(1.1791139) q[0];
rz(-2.9367661) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(1.0669605) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.897192) q[0];
sx q[0];
rz(-3.1163437) q[0];
sx q[0];
rz(-0.61141725) q[0];
x q[1];
rz(-2.7166769) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(-1.0645107) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9251717) q[1];
sx q[1];
rz(-0.74837084) q[1];
sx q[1];
rz(-1.5009576) q[1];
x q[2];
rz(2.9818929) q[3];
sx q[3];
rz(-0.80280639) q[3];
sx q[3];
rz(-1.476895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-0.55994326) q[2];
rz(-2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2840435) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-2.9395318) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(2.2834159) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9265384) q[0];
sx q[0];
rz(-0.19556043) q[0];
sx q[0];
rz(2.2901448) q[0];
x q[1];
rz(1.1602976) q[2];
sx q[2];
rz(-1.6121284) q[2];
sx q[2];
rz(-1.2683887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3837636) q[1];
sx q[1];
rz(-1.541829) q[1];
sx q[1];
rz(-1.5860735) q[1];
x q[2];
rz(1.9167561) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(-0.35017761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28785607) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.3254335) q[2];
rz(-0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637852) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(-0.97887865) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(-1.3495061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50586787) q[0];
sx q[0];
rz(-0.74698193) q[0];
sx q[0];
rz(-2.5542459) q[0];
rz(-3.1153203) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(0.19933137) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8229891) q[1];
sx q[1];
rz(-0.18853304) q[1];
sx q[1];
rz(1.8382501) q[1];
x q[2];
rz(1.8324864) q[3];
sx q[3];
rz(-2.1013386) q[3];
sx q[3];
rz(1.4025276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(1.7780001) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(-1.6040241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3106829) q[0];
sx q[0];
rz(-1.791782) q[0];
sx q[0];
rz(-1.2587147) q[0];
rz(-0.41861694) q[2];
sx q[2];
rz(-1.6659684) q[2];
sx q[2];
rz(1.0375432) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9201587) q[1];
sx q[1];
rz(-0.71343525) q[1];
sx q[1];
rz(0.17381298) q[1];
rz(2.780464) q[3];
sx q[3];
rz(-1.1361406) q[3];
sx q[3];
rz(-1.1772732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-0.76134479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.3392357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(0.95473081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2749598) q[0];
sx q[0];
rz(-1.2506335) q[0];
sx q[0];
rz(0.76036705) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54951349) q[2];
sx q[2];
rz(-1.2064484) q[2];
sx q[2];
rz(0.074631045) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4090875) q[1];
sx q[1];
rz(-1.115683) q[1];
sx q[1];
rz(1.485977) q[1];
rz(-pi) q[2];
rz(-2.2050489) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(1.6983502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0163991) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(-2.184536) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64086296) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-0.8846994) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(-0.62467081) q[2];
sx q[2];
rz(-0.93745898) q[2];
sx q[2];
rz(-2.6425101) q[2];
rz(-1.9839722) q[3];
sx q[3];
rz(-2.6631841) q[3];
sx q[3];
rz(-0.14315179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];