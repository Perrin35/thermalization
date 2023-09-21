OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0697203) q[0];
sx q[0];
rz(-0.9099996) q[0];
sx q[0];
rz(1.1138492) q[0];
rz(2.4835303) q[2];
sx q[2];
rz(-1.0339289) q[2];
sx q[2];
rz(1.8368349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77820233) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(-1.9744639) q[1];
rz(-1.4952881) q[3];
sx q[3];
rz(-1.2592053) q[3];
sx q[3];
rz(2.989245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-2.3031394) q[2];
sx q[2];
rz(-1.8189836) q[2];
rz(-2.8406075) q[3];
sx q[3];
rz(-2.5299278) q[3];
sx q[3];
rz(1.7606364) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.009636119) q[0];
sx q[0];
rz(-2.8490503) q[0];
sx q[0];
rz(0.47505501) q[0];
rz(1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-2.1038726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1514725) q[0];
sx q[0];
rz(-1.7904141) q[0];
sx q[0];
rz(-0.56818509) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3677164) q[2];
sx q[2];
rz(-2.4485588) q[2];
sx q[2];
rz(-2.3351923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71643752) q[1];
sx q[1];
rz(-1.5967224) q[1];
sx q[1];
rz(-1.9772711) q[1];
x q[2];
rz(-0.59229895) q[3];
sx q[3];
rz(-2.2364738) q[3];
sx q[3];
rz(-2.388282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-0.084687106) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(1.9975196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304853) q[0];
sx q[0];
rz(-1.100891) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(-0.39069191) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(2.5684165) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.976982) q[0];
sx q[0];
rz(-1.6342388) q[0];
sx q[0];
rz(-1.1399067) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8163221) q[2];
sx q[2];
rz(-0.79954445) q[2];
sx q[2];
rz(0.48649597) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1307615) q[1];
sx q[1];
rz(-1.8567994) q[1];
sx q[1];
rz(-2.7421013) q[1];
rz(-pi) q[2];
rz(0.40804789) q[3];
sx q[3];
rz(-1.9983665) q[3];
sx q[3];
rz(-1.5023295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-2.9476681) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.2852185) q[3];
sx q[3];
rz(2.9320419) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28213421) q[0];
sx q[0];
rz(-2.5971446) q[0];
sx q[0];
rz(2.5909246) q[0];
rz(-1.1286873) q[1];
sx q[1];
rz(-2.0813324) q[1];
sx q[1];
rz(-0.36270025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.761128) q[0];
sx q[0];
rz(-1.3306381) q[0];
sx q[0];
rz(-0.856075) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0998459) q[2];
sx q[2];
rz(-3.1051271) q[2];
sx q[2];
rz(-1.0066102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94707045) q[1];
sx q[1];
rz(-0.94049373) q[1];
sx q[1];
rz(-1.6987726) q[1];
rz(-pi) q[2];
rz(-1.4705212) q[3];
sx q[3];
rz(-2.5146211) q[3];
sx q[3];
rz(0.37563045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-2.0726911) q[2];
sx q[2];
rz(-0.0030227946) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-2.7987517) q[3];
sx q[3];
rz(-2.3390884) q[3];
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
rz(-0.18810774) q[0];
sx q[0];
rz(-0.67512023) q[0];
sx q[0];
rz(-0.014199646) q[0];
rz(-3.1242127) q[1];
sx q[1];
rz(-0.94795579) q[1];
sx q[1];
rz(1.682122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1634954) q[0];
sx q[0];
rz(-0.28143829) q[0];
sx q[0];
rz(1.0174169) q[0];
x q[1];
rz(-0.95894496) q[2];
sx q[2];
rz(-2.358837) q[2];
sx q[2];
rz(0.81629717) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1294714) q[1];
sx q[1];
rz(-2.8399889) q[1];
sx q[1];
rz(-0.71838897) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1305389) q[3];
sx q[3];
rz(-2.5821745) q[3];
sx q[3];
rz(1.2572446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83546272) q[2];
sx q[2];
rz(-2.7042522) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275948) q[0];
sx q[0];
rz(-0.70677775) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(1.2305413) q[1];
sx q[1];
rz(-1.1122333) q[1];
sx q[1];
rz(-0.13866436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5667144) q[0];
sx q[0];
rz(-3.1209271) q[0];
sx q[0];
rz(-2.979435) q[0];
rz(0.3615985) q[2];
sx q[2];
rz(-0.21441678) q[2];
sx q[2];
rz(-0.32985652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2885292) q[1];
sx q[1];
rz(-2.1398586) q[1];
sx q[1];
rz(2.2437375) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74216446) q[3];
sx q[3];
rz(-1.2483178) q[3];
sx q[3];
rz(-0.025067586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(-1.3343875) q[2];
rz(-1.1602317) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60550624) q[0];
sx q[0];
rz(-2.0083997) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(-2.5560608) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(0.62430635) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666703) q[0];
sx q[0];
rz(-0.86219388) q[0];
sx q[0];
rz(-1.8486345) q[0];
rz(-pi) q[1];
rz(2.1820071) q[2];
sx q[2];
rz(-0.72003905) q[2];
sx q[2];
rz(-0.50146539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.149951) q[1];
sx q[1];
rz(-2.2551564) q[1];
sx q[1];
rz(2.409694) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94654406) q[3];
sx q[3];
rz(-0.72580273) q[3];
sx q[3];
rz(3.1012227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5252934) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(2.7590511) q[2];
rz(0.031490695) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(-3.0338045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87576762) q[0];
sx q[0];
rz(-2.8540397) q[0];
sx q[0];
rz(-3.0016622) q[0];
rz(1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(-0.12891842) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77057225) q[0];
sx q[0];
rz(-2.7311374) q[0];
sx q[0];
rz(1.3034986) q[0];
x q[1];
rz(0.89959896) q[2];
sx q[2];
rz(-1.1102144) q[2];
sx q[2];
rz(2.7923498) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8184549) q[1];
sx q[1];
rz(-1.2982681) q[1];
sx q[1];
rz(2.6998991) q[1];
rz(1.586732) q[3];
sx q[3];
rz(-0.83105479) q[3];
sx q[3];
rz(0.15077886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.63697469) q[2];
sx q[2];
rz(-1.0417754) q[2];
sx q[2];
rz(0.62409419) q[2];
rz(2.9028153) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(2.074266) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(3.1255787) q[1];
sx q[1];
rz(-0.7557973) q[1];
sx q[1];
rz(1.790766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.006376) q[0];
sx q[0];
rz(-1.3906286) q[0];
sx q[0];
rz(0.10794497) q[0];
rz(-1.1662912) q[2];
sx q[2];
rz(-2.1098237) q[2];
sx q[2];
rz(-1.5750969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0235325) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(2.1937624) q[1];
x q[2];
rz(-2.8899367) q[3];
sx q[3];
rz(-0.76905426) q[3];
sx q[3];
rz(0.8151527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.089036971) q[2];
sx q[2];
rz(-0.75110835) q[2];
sx q[2];
rz(-3.0272711) q[2];
rz(1.2601241) q[3];
sx q[3];
rz(-1.0269287) q[3];
sx q[3];
rz(-1.1589706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.91530144) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(0.21324883) q[0];
rz(2.7217216) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(-2.5949123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8583546) q[0];
sx q[0];
rz(-2.3843117) q[0];
sx q[0];
rz(1.8961294) q[0];
x q[1];
rz(-0.28995138) q[2];
sx q[2];
rz(-0.49294127) q[2];
sx q[2];
rz(1.314756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9676799) q[1];
sx q[1];
rz(-0.2914857) q[1];
sx q[1];
rz(3.013054) q[1];
x q[2];
rz(-2.5952026) q[3];
sx q[3];
rz(-0.26248172) q[3];
sx q[3];
rz(-2.4266092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0032349) q[2];
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
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99988408) q[0];
sx q[0];
rz(-1.1078436) q[0];
sx q[0];
rz(-2.1583337) q[0];
rz(-2.6976363) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(1.6179686) q[2];
sx q[2];
rz(-1.5856367) q[2];
sx q[2];
rz(-0.71935364) q[2];
rz(0.32904939) q[3];
sx q[3];
rz(-1.1209189) q[3];
sx q[3];
rz(0.95595595) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
