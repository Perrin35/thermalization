OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2603962) q[0];
sx q[0];
rz(-1.2328147) q[0];
sx q[0];
rz(0.32786274) q[0];
rz(-0.39363632) q[1];
sx q[1];
rz(-1.5235695) q[1];
sx q[1];
rz(-1.5394428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0547787) q[0];
sx q[0];
rz(-0.20113763) q[0];
sx q[0];
rz(-3.1166409) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27199409) q[2];
sx q[2];
rz(-2.1025206) q[2];
sx q[2];
rz(-3.0716999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2214927) q[1];
sx q[1];
rz(-1.3650238) q[1];
sx q[1];
rz(0.39452067) q[1];
rz(-0.39253916) q[3];
sx q[3];
rz(-1.919636) q[3];
sx q[3];
rz(-2.4667645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1800805) q[2];
sx q[2];
rz(-2.1891258) q[2];
sx q[2];
rz(2.3251779) q[2];
rz(-2.4444729) q[3];
sx q[3];
rz(-2.7643047) q[3];
sx q[3];
rz(2.2009946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78696752) q[0];
sx q[0];
rz(-2.0836199) q[0];
sx q[0];
rz(-0.70607287) q[0];
rz(-2.970703) q[1];
sx q[1];
rz(-0.51097521) q[1];
sx q[1];
rz(0.99542803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1814897) q[0];
sx q[0];
rz(-1.8892291) q[0];
sx q[0];
rz(-0.35950999) q[0];
rz(-pi) q[1];
rz(2.1336251) q[2];
sx q[2];
rz(-2.3555104) q[2];
sx q[2];
rz(2.2843334) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.625841) q[1];
sx q[1];
rz(-1.4179572) q[1];
sx q[1];
rz(-1.3550853) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0964504) q[3];
sx q[3];
rz(-0.63796872) q[3];
sx q[3];
rz(0.54813069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55676111) q[2];
sx q[2];
rz(-2.2274667) q[2];
sx q[2];
rz(2.8625028) q[2];
rz(0.68785214) q[3];
sx q[3];
rz(-2.9180241) q[3];
sx q[3];
rz(-1.0031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85182178) q[0];
sx q[0];
rz(-2.2230447) q[0];
sx q[0];
rz(-2.7626792) q[0];
rz(1.1907578) q[1];
sx q[1];
rz(-0.36668101) q[1];
sx q[1];
rz(-0.77622882) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2874173) q[0];
sx q[0];
rz(-2.6834134) q[0];
sx q[0];
rz(-1.2121546) q[0];
rz(1.4188781) q[2];
sx q[2];
rz(-1.3891274) q[2];
sx q[2];
rz(-0.40975299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6783733) q[1];
sx q[1];
rz(-1.8253981) q[1];
sx q[1];
rz(-0.71779743) q[1];
rz(3.0954563) q[3];
sx q[3];
rz(-2.0946119) q[3];
sx q[3];
rz(-2.5906117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0188521) q[2];
sx q[2];
rz(-1.2383702) q[2];
sx q[2];
rz(-2.5416601) q[2];
rz(-1.9020724) q[3];
sx q[3];
rz(-1.7299088) q[3];
sx q[3];
rz(0.6111353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3032853) q[0];
sx q[0];
rz(-2.4628283) q[0];
sx q[0];
rz(1.432206) q[0];
rz(2.2749061) q[1];
sx q[1];
rz(-1.7984093) q[1];
sx q[1];
rz(2.0492882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537395) q[0];
sx q[0];
rz(-2.1668375) q[0];
sx q[0];
rz(0.88576646) q[0];
rz(-pi) q[1];
rz(-2.1020911) q[2];
sx q[2];
rz(-2.4091125) q[2];
sx q[2];
rz(-1.8341482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0835613) q[1];
sx q[1];
rz(-1.5688854) q[1];
sx q[1];
rz(2.6772736) q[1];
x q[2];
rz(-2.3208612) q[3];
sx q[3];
rz(-2.1607275) q[3];
sx q[3];
rz(1.8128519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0852574) q[2];
sx q[2];
rz(-1.9838355) q[2];
sx q[2];
rz(1.2752656) q[2];
rz(-1.7665675) q[3];
sx q[3];
rz(-1.2396038) q[3];
sx q[3];
rz(-2.1084771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035868693) q[0];
sx q[0];
rz(-2.2585456) q[0];
sx q[0];
rz(1.7778273) q[0];
rz(-0.56963244) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(-0.5035351) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.106905) q[0];
sx q[0];
rz(-1.4464753) q[0];
sx q[0];
rz(0.19363815) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84837368) q[2];
sx q[2];
rz(-0.97992291) q[2];
sx q[2];
rz(-0.17308519) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40696661) q[1];
sx q[1];
rz(-1.5583724) q[1];
sx q[1];
rz(1.6520651) q[1];
rz(-pi) q[2];
rz(2.8315169) q[3];
sx q[3];
rz(-2.3469124) q[3];
sx q[3];
rz(0.39439699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.86294952) q[2];
sx q[2];
rz(-0.38662616) q[2];
sx q[2];
rz(1.53842) q[2];
rz(-0.1362416) q[3];
sx q[3];
rz(-1.7377661) q[3];
sx q[3];
rz(-1.9041825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.49029008) q[0];
sx q[0];
rz(-0.61838111) q[0];
sx q[0];
rz(2.8001617) q[0];
rz(2.4104207) q[1];
sx q[1];
rz(-0.57594222) q[1];
sx q[1];
rz(1.0715356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8971469) q[0];
sx q[0];
rz(-1.9715152) q[0];
sx q[0];
rz(-0.44735502) q[0];
rz(2.5265556) q[2];
sx q[2];
rz(-1.5205548) q[2];
sx q[2];
rz(-0.096252354) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2372969) q[1];
sx q[1];
rz(-2.0091506) q[1];
sx q[1];
rz(-2.3890952) q[1];
rz(-pi) q[2];
rz(3.1370779) q[3];
sx q[3];
rz(-2.0751953) q[3];
sx q[3];
rz(-2.2405091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7858481) q[2];
sx q[2];
rz(-1.637849) q[2];
sx q[2];
rz(-2.6453633) q[2];
rz(-1.5931386) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(1.0065856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16794311) q[0];
sx q[0];
rz(-1.4730467) q[0];
sx q[0];
rz(-3.1267401) q[0];
rz(-1.7282093) q[1];
sx q[1];
rz(-1.1057066) q[1];
sx q[1];
rz(0.76494876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0891554) q[0];
sx q[0];
rz(-0.90298072) q[0];
sx q[0];
rz(2.3994156) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4582389) q[2];
sx q[2];
rz(-1.2938479) q[2];
sx q[2];
rz(2.5954557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97100869) q[1];
sx q[1];
rz(-2.7448344) q[1];
sx q[1];
rz(3.1301857) q[1];
x q[2];
rz(-1.1469923) q[3];
sx q[3];
rz(-1.7693345) q[3];
sx q[3];
rz(1.0120704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0689653) q[2];
sx q[2];
rz(-1.5494538) q[2];
sx q[2];
rz(0.7507945) q[2];
rz(-0.99610656) q[3];
sx q[3];
rz(-2.3361497) q[3];
sx q[3];
rz(2.6944845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873416) q[0];
sx q[0];
rz(-0.60289201) q[0];
sx q[0];
rz(2.531429) q[0];
rz(-2.8918686) q[1];
sx q[1];
rz(-1.4796673) q[1];
sx q[1];
rz(-3.0009559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661903) q[0];
sx q[0];
rz(-1.1151214) q[0];
sx q[0];
rz(2.3123017) q[0];
rz(-pi) q[1];
rz(2.6360077) q[2];
sx q[2];
rz(-1.4806432) q[2];
sx q[2];
rz(-2.8204927) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7761155) q[1];
sx q[1];
rz(-1.8473313) q[1];
sx q[1];
rz(1.6460258) q[1];
rz(2.3811627) q[3];
sx q[3];
rz(-1.3730341) q[3];
sx q[3];
rz(-0.20748479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7481302) q[2];
sx q[2];
rz(-2.738214) q[2];
sx q[2];
rz(2.051579) q[2];
rz(-1.9213093) q[3];
sx q[3];
rz(-2.0956495) q[3];
sx q[3];
rz(0.88855851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1423993) q[0];
sx q[0];
rz(-2.0274473) q[0];
sx q[0];
rz(-0.92495579) q[0];
rz(-1.6389182) q[1];
sx q[1];
rz(-2.131856) q[1];
sx q[1];
rz(-2.9138873) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0696237) q[0];
sx q[0];
rz(-1.6771131) q[0];
sx q[0];
rz(-2.1031455) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57982071) q[2];
sx q[2];
rz(-1.3353773) q[2];
sx q[2];
rz(-1.5832242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0947837) q[1];
sx q[1];
rz(-1.9247479) q[1];
sx q[1];
rz(-1.7750398) q[1];
rz(2.5277522) q[3];
sx q[3];
rz(-1.7305529) q[3];
sx q[3];
rz(-2.0789119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1091252) q[2];
sx q[2];
rz(-2.511907) q[2];
sx q[2];
rz(0.49051782) q[2];
rz(1.0444752) q[3];
sx q[3];
rz(-2.677768) q[3];
sx q[3];
rz(3.131955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.57546416) q[0];
sx q[0];
rz(-2.6818891) q[0];
sx q[0];
rz(1.8756961) q[0];
rz(1.5220386) q[1];
sx q[1];
rz(-1.8892037) q[1];
sx q[1];
rz(0.44580805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73484008) q[0];
sx q[0];
rz(-1.2452083) q[0];
sx q[0];
rz(2.1181137) q[0];
rz(-pi) q[1];
rz(1.5238026) q[2];
sx q[2];
rz(-1.0596152) q[2];
sx q[2];
rz(-1.5509449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5490732) q[1];
sx q[1];
rz(-0.62738327) q[1];
sx q[1];
rz(0.23457867) q[1];
x q[2];
rz(0.021546797) q[3];
sx q[3];
rz(-1.2974707) q[3];
sx q[3];
rz(0.97238982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.012099115) q[2];
sx q[2];
rz(-0.98126498) q[2];
sx q[2];
rz(1.0274461) q[2];
rz(-1.8217314) q[3];
sx q[3];
rz(-2.8362995) q[3];
sx q[3];
rz(-0.52904883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32578596) q[0];
sx q[0];
rz(-1.3675714) q[0];
sx q[0];
rz(2.0969781) q[0];
rz(0.927399) q[1];
sx q[1];
rz(-1.2694042) q[1];
sx q[1];
rz(-0.70613695) q[1];
rz(-1.2583707) q[2];
sx q[2];
rz(-1.8835405) q[2];
sx q[2];
rz(0.082688511) q[2];
rz(-2.5221323) q[3];
sx q[3];
rz(-1.3347698) q[3];
sx q[3];
rz(2.6556849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
