OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(-1.7655656) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5698619) q[0];
sx q[0];
rz(-2.1751582) q[0];
sx q[0];
rz(0.5423003) q[0];
rz(-2.4260169) q[2];
sx q[2];
rz(-1.1872429) q[2];
sx q[2];
rz(-0.54563145) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4151167) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(0.16023689) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81935482) q[3];
sx q[3];
rz(-1.6091533) q[3];
sx q[3];
rz(1.0591782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.084289) q[2];
sx q[2];
rz(-2.7412582) q[2];
sx q[2];
rz(0.98891813) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-2.3108216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(0.99769366) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-0.72431272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9162826) q[0];
sx q[0];
rz(-2.3819469) q[0];
sx q[0];
rz(2.0347974) q[0];
x q[1];
rz(2.9421147) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(0.85180887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95485479) q[1];
sx q[1];
rz(-1.2192982) q[1];
sx q[1];
rz(-0.21912205) q[1];
x q[2];
rz(-0.23726666) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(-2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(-3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.3954337) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.2190855) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.483706) q[0];
sx q[0];
rz(-1.3077967) q[0];
sx q[0];
rz(-1.2067243) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8457882) q[2];
sx q[2];
rz(-0.89106262) q[2];
sx q[2];
rz(-2.5909397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(2.0134258) q[1];
x q[2];
rz(0.015467042) q[3];
sx q[3];
rz(-1.3028212) q[3];
sx q[3];
rz(-1.0166575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(-1.6960309) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(0.82675654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2430192) q[0];
sx q[0];
rz(-1.9802226) q[0];
sx q[0];
rz(-0.55222521) q[0];
x q[1];
rz(-2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(-1.3354288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7552232) q[1];
sx q[1];
rz(-1.9519023) q[1];
sx q[1];
rz(-1.2632881) q[1];
x q[2];
rz(-1.8825674) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7171705) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(-1.0958584) q[0];
rz(-1.1826913) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(-2.2270122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73357108) q[1];
sx q[1];
rz(-1.5703778) q[1];
sx q[1];
rz(-1.8838521) q[1];
rz(-pi) q[2];
rz(-2.9167446) q[3];
sx q[3];
rz(-2.7307011) q[3];
sx q[3];
rz(-1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(-0.53945333) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-0.8845194) q[3];
sx q[3];
rz(-0.45421281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25813112) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.7745811) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54129823) q[0];
sx q[0];
rz(-1.6389264) q[0];
sx q[0];
rz(3.1085204) q[0];
rz(-pi) q[1];
rz(0.96732803) q[2];
sx q[2];
rz(-2.2852995) q[2];
sx q[2];
rz(-2.8514903) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4821266) q[1];
sx q[1];
rz(-0.30059338) q[1];
sx q[1];
rz(-1.2947004) q[1];
x q[2];
rz(-2.8664687) q[3];
sx q[3];
rz(-0.36558662) q[3];
sx q[3];
rz(-0.58055731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8453025) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(0.40346754) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(0.8738628) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0872333) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(1.5646294) q[0];
x q[1];
rz(2.7295693) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(-0.7574946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9782941) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(-0.75978029) q[1];
rz(-2.6074334) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(-1.1155038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(-2.5308385) q[2];
rz(2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(0.92774123) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24191813) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(-2.8444667) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(2.4954605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532115) q[0];
sx q[0];
rz(-1.6717981) q[0];
sx q[0];
rz(-1.3636916) q[0];
rz(-pi) q[1];
rz(1.4023151) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(-0.70665765) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7568126) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(-0.53201075) q[1];
rz(2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(0.83412795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.064676553) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(2.440195) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(1.1340244) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.425697) q[0];
sx q[0];
rz(-1.513039) q[0];
sx q[0];
rz(-0.9066559) q[0];
x q[1];
rz(1.8234532) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(2.2968963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(2.2309169) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7897723) q[3];
sx q[3];
rz(-2.4253546) q[3];
sx q[3];
rz(0.66974528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0124399) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(0.41839504) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(-0.72559124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-2.4627731) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(3.0864339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96345761) q[0];
sx q[0];
rz(-1.194251) q[0];
sx q[0];
rz(1.6403273) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32785428) q[2];
sx q[2];
rz(-2.556986) q[2];
sx q[2];
rz(2.2962928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(2.39141) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1049214) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(2.7367221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(-2.8293777) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(-0.24003868) q[3];
sx q[3];
rz(-1.6143027) q[3];
sx q[3];
rz(-1.3054813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
