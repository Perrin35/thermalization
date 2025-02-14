OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8483228) q[0];
sx q[0];
rz(-0.25265101) q[0];
sx q[0];
rz(-1.2055612) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(-1.815058) q[1];
sx q[1];
rz(2.3958652) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416314) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(0.061454031) q[0];
rz(2.9171997) q[2];
sx q[2];
rz(-1.5641777) q[2];
sx q[2];
rz(1.4057019) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28994014) q[1];
sx q[1];
rz(-1.2356241) q[1];
sx q[1];
rz(-2.9194458) q[1];
rz(1.0966572) q[3];
sx q[3];
rz(-1.7563987) q[3];
sx q[3];
rz(-1.4193648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91743177) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(-1.6708299) q[2];
rz(1.6338232) q[3];
sx q[3];
rz(-3.1247415) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(-1.5731328) q[0];
rz(-0.16954999) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(3.0068908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7260774) q[0];
sx q[0];
rz(-1.2252136) q[0];
sx q[0];
rz(0.4536566) q[0];
x q[1];
rz(1.3943761) q[2];
sx q[2];
rz(-3.0718832) q[2];
sx q[2];
rz(-1.4331872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0162279) q[1];
sx q[1];
rz(-2.008359) q[1];
sx q[1];
rz(1.2368278) q[1];
x q[2];
rz(1.1504796) q[3];
sx q[3];
rz(-2.947071) q[3];
sx q[3];
rz(-1.9357301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0160825) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(-2.9888195) q[2];
rz(-1.7759391) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50853866) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(-0.48164865) q[0];
rz(2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(0.99536037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1448875) q[0];
sx q[0];
rz(-1.2249682) q[0];
sx q[0];
rz(-2.8976827) q[0];
rz(1.6312509) q[2];
sx q[2];
rz(-1.5222856) q[2];
sx q[2];
rz(-1.2918351) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0479483) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(2.4574952) q[1];
rz(-pi) q[2];
rz(2.1667492) q[3];
sx q[3];
rz(-1.5454834) q[3];
sx q[3];
rz(1.635575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92622009) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(0.064662956) q[2];
rz(-1.0489382) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(-1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8365086) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-2.3205561) q[0];
rz(3.0631284) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(2.1441114) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68551915) q[0];
sx q[0];
rz(-0.711383) q[0];
sx q[0];
rz(2.6472241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5115963) q[2];
sx q[2];
rz(-1.534605) q[2];
sx q[2];
rz(-0.22360392) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4892024) q[1];
sx q[1];
rz(-1.308131) q[1];
sx q[1];
rz(-2.1173544) q[1];
x q[2];
rz(2.663199) q[3];
sx q[3];
rz(-1.5198623) q[3];
sx q[3];
rz(-0.8138322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(1.9878261) q[2];
rz(0.18618259) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(-0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368197) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(-1.1432884) q[0];
rz(-2.000287) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(0.51796651) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85803131) q[0];
sx q[0];
rz(-1.0378583) q[0];
sx q[0];
rz(2.3951247) q[0];
x q[1];
rz(-1.5094366) q[2];
sx q[2];
rz(-3.1265335) q[2];
sx q[2];
rz(-1.3850152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8570366) q[1];
sx q[1];
rz(-1.9003881) q[1];
sx q[1];
rz(-0.28917851) q[1];
rz(-pi) q[2];
rz(-2.6916396) q[3];
sx q[3];
rz(-2.8128689) q[3];
sx q[3];
rz(-2.4462819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8225857) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(-0.32773584) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0026534) q[0];
sx q[0];
rz(-0.26873538) q[0];
sx q[0];
rz(2.5797381) q[0];
rz(1.4909164) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(-3.0432826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54360169) q[0];
sx q[0];
rz(-0.4234792) q[0];
sx q[0];
rz(1.6571568) q[0];
rz(-0.00064050015) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(-1.887111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6358306) q[1];
sx q[1];
rz(-2.039685) q[1];
sx q[1];
rz(1.6230349) q[1];
rz(1.1071854) q[3];
sx q[3];
rz(-1.1157728) q[3];
sx q[3];
rz(-2.1912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.085122846) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(-2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(-1.2679509) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999076) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(-2.0211925) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(0.33946005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0507815) q[0];
sx q[0];
rz(-0.078074038) q[0];
sx q[0];
rz(0.17103057) q[0];
rz(-pi) q[1];
rz(-3.1192084) q[2];
sx q[2];
rz(-2.5867992) q[2];
sx q[2];
rz(3.1172397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5875747) q[1];
sx q[1];
rz(-1.4851991) q[1];
sx q[1];
rz(-3.1399591) q[1];
x q[2];
rz(1.4711597) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(-1.3226313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(2.778229) q[2];
rz(2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(-1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5217487) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(-2.2186665) q[0];
rz(-1.482831) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-0.0042075687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7711346) q[0];
sx q[0];
rz(-2.1677289) q[0];
sx q[0];
rz(0.93907066) q[0];
rz(-pi) q[1];
rz(1.9816859) q[2];
sx q[2];
rz(-1.5762323) q[2];
sx q[2];
rz(-1.5856575) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15531047) q[1];
sx q[1];
rz(-2.661507) q[1];
sx q[1];
rz(1.4309149) q[1];
x q[2];
rz(3.0757853) q[3];
sx q[3];
rz(-0.61344693) q[3];
sx q[3];
rz(-0.25019161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.2008249) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414108) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(-1.2196983) q[0];
rz(-1.8118743) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78744473) q[0];
sx q[0];
rz(-1.7423) q[0];
sx q[0];
rz(-1.9308912) q[0];
x q[1];
rz(0.40028769) q[2];
sx q[2];
rz(-0.62303715) q[2];
sx q[2];
rz(1.7965339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6170609) q[1];
sx q[1];
rz(-1.6304029) q[1];
sx q[1];
rz(1.4527133) q[1];
rz(-pi) q[2];
rz(-2.2396062) q[3];
sx q[3];
rz(-0.94873442) q[3];
sx q[3];
rz(-0.40611503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4280052) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(0.62291992) q[2];
rz(0.04341393) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(2.3626732) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6503705) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(0.69475118) q[1];
sx q[1];
rz(-0.34725747) q[1];
sx q[1];
rz(2.6659226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2025288) q[0];
sx q[0];
rz(-3.0923884) q[0];
sx q[0];
rz(-2.3243107) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6108405) q[2];
sx q[2];
rz(-2.4350428) q[2];
sx q[2];
rz(1.6395115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.039363843) q[1];
sx q[1];
rz(-2.1783531) q[1];
sx q[1];
rz(-1.5786912) q[1];
rz(-pi) q[2];
rz(0.31020152) q[3];
sx q[3];
rz(-1.7473012) q[3];
sx q[3];
rz(1.4909397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67705578) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(-1.8439058) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(-2.2592648) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(1.6068983) q[2];
sx q[2];
rz(-1.4786199) q[2];
sx q[2];
rz(0.27603966) q[2];
rz(-3.121331) q[3];
sx q[3];
rz(-2.9023585) q[3];
sx q[3];
rz(0.0016061847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
