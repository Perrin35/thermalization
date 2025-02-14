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
rz(1.9360315) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(-1.815058) q[1];
sx q[1];
rz(2.3958652) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099961258) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(-3.0801386) q[0];
x q[1];
rz(-3.1118561) q[2];
sx q[2];
rz(-2.9171037) q[2];
sx q[2];
rz(2.9475074) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9348976) q[1];
sx q[1];
rz(-1.3612011) q[1];
sx q[1];
rz(-1.2278201) q[1];
rz(-pi) q[2];
rz(-2.0449355) q[3];
sx q[3];
rz(-1.7563987) q[3];
sx q[3];
rz(-1.4193648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2241609) q[2];
sx q[2];
rz(-1.5459583) q[2];
sx q[2];
rz(1.6708299) q[2];
rz(1.6338232) q[3];
sx q[3];
rz(-3.1247415) q[3];
sx q[3];
rz(-0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(-1.5731328) q[0];
rz(2.9720427) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(3.0068908) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8226263) q[0];
sx q[0];
rz(-1.14577) q[0];
sx q[0];
rz(1.1898196) q[0];
rz(-pi) q[1];
rz(1.3943761) q[2];
sx q[2];
rz(-0.069709502) q[2];
sx q[2];
rz(1.4331872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4410494) q[1];
sx q[1];
rz(-1.8722539) q[1];
sx q[1];
rz(0.45977199) q[1];
rz(-1.7487583) q[3];
sx q[3];
rz(-1.6497532) q[3];
sx q[3];
rz(-0.7782026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(-0.15277319) q[2];
rz(-1.7759391) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(-0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50853866) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(0.48164865) q[0];
rz(0.18474361) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(0.99536037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48993123) q[0];
sx q[0];
rz(-1.7999987) q[0];
sx q[0];
rz(-1.9263173) q[0];
rz(-0.89389385) q[2];
sx q[2];
rz(-0.077493103) q[2];
sx q[2];
rz(-2.187196) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0936443) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(2.4574952) q[1];
rz(0.97484346) q[3];
sx q[3];
rz(-1.5961093) q[3];
sx q[3];
rz(1.635575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2153726) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(0.064662956) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(-1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30508405) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(0.82103658) q[0];
rz(0.07846421) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(0.99748126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0745463) q[0];
sx q[0];
rz(-0.95854488) q[0];
sx q[0];
rz(-1.9590098) q[0];
rz(-pi) q[1];
rz(-3.1053379) q[2];
sx q[2];
rz(-1.6299575) q[2];
sx q[2];
rz(1.7922557) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6523903) q[1];
sx q[1];
rz(-1.8334616) q[1];
sx q[1];
rz(1.0242382) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47839368) q[3];
sx q[3];
rz(-1.5198623) q[3];
sx q[3];
rz(0.8138322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(1.1537665) q[2];
rz(-0.18618259) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004772923) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(1.9983043) q[0];
rz(-2.000287) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(-2.6236261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1521027) q[0];
sx q[0];
rz(-2.1955262) q[0];
sx q[0];
rz(-2.2476907) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.632156) q[2];
sx q[2];
rz(-3.1265335) q[2];
sx q[2];
rz(1.3850152) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3822381) q[1];
sx q[1];
rz(-1.2975946) q[1];
sx q[1];
rz(1.9135936) q[1];
rz(2.84359) q[3];
sx q[3];
rz(-1.4299222) q[3];
sx q[3];
rz(-0.44671392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(2.8138568) q[2];
rz(-1.94708) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(2.5797381) q[0];
rz(-1.4909164) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(-3.0432826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.597991) q[0];
sx q[0];
rz(-0.4234792) q[0];
sx q[0];
rz(-1.4844358) q[0];
x q[1];
rz(-0.00064050015) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(-1.887111) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6209539) q[1];
sx q[1];
rz(-0.47157447) q[1];
sx q[1];
rz(3.0388799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6410511) q[3];
sx q[3];
rz(-1.1574452) q[3];
sx q[3];
rz(-2.7374637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0564698) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(2.0758212) q[2];
rz(2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(-1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.6837233) q[0];
rz(-1.1204002) q[1];
sx q[1];
rz(-0.1381865) q[1];
sx q[1];
rz(0.33946005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4910879) q[0];
sx q[0];
rz(-1.5840713) q[0];
sx q[0];
rz(0.076939452) q[0];
rz(-pi) q[1];
x q[1];
rz(0.022384251) q[2];
sx q[2];
rz(-2.5867992) q[2];
sx q[2];
rz(-0.024352976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5875747) q[1];
sx q[1];
rz(-1.4851991) q[1];
sx q[1];
rz(3.1399591) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27642997) q[3];
sx q[3];
rz(-2.7904841) q[3];
sx q[3];
rz(-2.1123667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(-2.778229) q[2];
rz(1.0904788) q[3];
sx q[3];
rz(-0.57991475) q[3];
sx q[3];
rz(-1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61984396) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(0.92292619) q[0];
rz(1.6587616) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(3.1373851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7711346) q[0];
sx q[0];
rz(-0.97386375) q[0];
sx q[0];
rz(2.202522) q[0];
rz(-pi) q[1];
rz(1.5844052) q[2];
sx q[2];
rz(-2.7306692) q[2];
sx q[2];
rz(-0.027337242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6018659) q[1];
sx q[1];
rz(-1.5063573) q[1];
sx q[1];
rz(-1.0947202) q[1];
rz(-pi) q[2];
rz(0.61242731) q[3];
sx q[3];
rz(-1.53293) q[3];
sx q[3];
rz(-1.374439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5795472) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(1.9407678) q[2];
rz(-1.7480525) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(-2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414108) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(1.2196983) q[0];
rz(-1.8118743) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9326646) q[0];
sx q[0];
rz(-0.39723662) q[0];
sx q[0];
rz(1.1139289) q[0];
x q[1];
rz(-0.40028769) q[2];
sx q[2];
rz(-2.5185555) q[2];
sx q[2];
rz(-1.3450587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1023952) q[1];
sx q[1];
rz(-1.6886687) q[1];
sx q[1];
rz(-3.0815691) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4285942) q[3];
sx q[3];
rz(-2.2621691) q[3];
sx q[3];
rz(-0.52934968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71358744) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(0.62291992) q[2];
rz(3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(-2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(-0.69475118) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(-0.4756701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9566475) q[0];
sx q[0];
rz(-1.534919) q[0];
sx q[0];
rz(-0.033680276) q[0];
rz(-0.034157201) q[2];
sx q[2];
rz(-0.86493051) q[2];
sx q[2];
rz(1.5868843) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.039363843) q[1];
sx q[1];
rz(-2.1783531) q[1];
sx q[1];
rz(1.5786912) q[1];
rz(-pi) q[2];
rz(1.7559515) q[3];
sx q[3];
rz(-1.8760215) q[3];
sx q[3];
rz(-3.117962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4645369) q[2];
sx q[2];
rz(-0.060364351) q[2];
sx q[2];
rz(1.2976868) q[2];
rz(-1.8929947) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(-0.36778522) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.5346943) q[2];
sx q[2];
rz(-1.4786199) q[2];
sx q[2];
rz(0.27603966) q[2];
rz(1.5757379) q[3];
sx q[3];
rz(-1.3316122) q[3];
sx q[3];
rz(0.022461654) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
