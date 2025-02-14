OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0889283) q[0];
sx q[0];
rz(5.1533617) q[0];
sx q[0];
rz(12.552153) q[0];
rz(0.05834236) q[1];
sx q[1];
rz(0.74906936) q[1];
sx q[1];
rz(11.963117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.280312) q[0];
sx q[0];
rz(-2.2076108) q[0];
sx q[0];
rz(1.1293108) q[0];
rz(-pi) q[1];
rz(2.2255664) q[2];
sx q[2];
rz(-1.5226622) q[2];
sx q[2];
rz(0.88716884) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4976884) q[1];
sx q[1];
rz(-0.7086646) q[1];
sx q[1];
rz(-0.24073118) q[1];
x q[2];
rz(-0.18666191) q[3];
sx q[3];
rz(-0.87282729) q[3];
sx q[3];
rz(-1.9912793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3732036) q[2];
sx q[2];
rz(-1.5378121) q[2];
sx q[2];
rz(1.0962037) q[2];
rz(2.0632035) q[3];
sx q[3];
rz(-0.53467852) q[3];
sx q[3];
rz(0.53511846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32644367) q[0];
sx q[0];
rz(-1.0468227) q[0];
sx q[0];
rz(-2.6640418) q[0];
rz(1.011147) q[1];
sx q[1];
rz(-2.5944581) q[1];
sx q[1];
rz(1.0844885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732373) q[0];
sx q[0];
rz(-0.38738444) q[0];
sx q[0];
rz(1.8854467) q[0];
x q[1];
rz(1.162503) q[2];
sx q[2];
rz(-1.5973685) q[2];
sx q[2];
rz(1.0607189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66531679) q[1];
sx q[1];
rz(-1.3937104) q[1];
sx q[1];
rz(0.55426532) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73954496) q[3];
sx q[3];
rz(-1.7266183) q[3];
sx q[3];
rz(0.23529822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61887211) q[2];
sx q[2];
rz(-0.42417696) q[2];
sx q[2];
rz(2.0514533) q[2];
rz(2.5777396) q[3];
sx q[3];
rz(-0.68958759) q[3];
sx q[3];
rz(-3.1088945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3801333) q[0];
sx q[0];
rz(-0.021012336) q[0];
sx q[0];
rz(-1.6047961) q[0];
rz(1.8236632) q[1];
sx q[1];
rz(-1.0937966) q[1];
sx q[1];
rz(0.39753786) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1434162) q[0];
sx q[0];
rz(-2.2097124) q[0];
sx q[0];
rz(-2.8692416) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2961646) q[2];
sx q[2];
rz(-2.1284747) q[2];
sx q[2];
rz(-2.9147343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2469625) q[1];
sx q[1];
rz(-1.8254465) q[1];
sx q[1];
rz(-0.92856426) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0769172) q[3];
sx q[3];
rz(-2.0398413) q[3];
sx q[3];
rz(2.7739547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2354108) q[2];
sx q[2];
rz(-1.9637039) q[2];
sx q[2];
rz(0.57514352) q[2];
rz(-0.67342526) q[3];
sx q[3];
rz(-1.6005102) q[3];
sx q[3];
rz(-1.5448236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.6509318) q[0];
sx q[0];
rz(-1.1216102) q[0];
sx q[0];
rz(-0.90890539) q[0];
rz(3.124369) q[1];
sx q[1];
rz(-0.52809087) q[1];
sx q[1];
rz(-3.1294894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41866583) q[0];
sx q[0];
rz(-1.8932492) q[0];
sx q[0];
rz(-0.61005436) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4922089) q[2];
sx q[2];
rz(-1.1282136) q[2];
sx q[2];
rz(-1.336652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3693272) q[1];
sx q[1];
rz(-0.64913926) q[1];
sx q[1];
rz(1.9499798) q[1];
rz(-pi) q[2];
rz(-0.12676012) q[3];
sx q[3];
rz(-1.6157627) q[3];
sx q[3];
rz(-2.8823406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.072731344) q[2];
sx q[2];
rz(-0.27421811) q[2];
sx q[2];
rz(-2.7552628) q[2];
rz(-0.38935152) q[3];
sx q[3];
rz(-1.6081622) q[3];
sx q[3];
rz(1.0079591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7858081) q[0];
sx q[0];
rz(-2.4171827) q[0];
sx q[0];
rz(-0.32387787) q[0];
rz(0.38662275) q[1];
sx q[1];
rz(-1.3893678) q[1];
sx q[1];
rz(-1.5547543) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51832685) q[0];
sx q[0];
rz(-0.73438209) q[0];
sx q[0];
rz(0.096303864) q[0];
rz(-2.7681211) q[2];
sx q[2];
rz(-1.0658776) q[2];
sx q[2];
rz(2.2511456) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91134478) q[1];
sx q[1];
rz(-2.2007588) q[1];
sx q[1];
rz(-0.98736422) q[1];
x q[2];
rz(1.0500141) q[3];
sx q[3];
rz(-0.84283295) q[3];
sx q[3];
rz(1.4340374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3891478) q[2];
sx q[2];
rz(-2.2989595) q[2];
sx q[2];
rz(-0.064099224) q[2];
rz(-2.5373503) q[3];
sx q[3];
rz(-1.3925545) q[3];
sx q[3];
rz(1.909168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9391249) q[0];
sx q[0];
rz(-0.61102837) q[0];
sx q[0];
rz(-0.87919277) q[0];
rz(-0.35036707) q[1];
sx q[1];
rz(-1.8354974) q[1];
sx q[1];
rz(-2.9894357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6317196) q[0];
sx q[0];
rz(-1.3001967) q[0];
sx q[0];
rz(2.4218049) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9486676) q[2];
sx q[2];
rz(-1.0954276) q[2];
sx q[2];
rz(2.1383291) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85206807) q[1];
sx q[1];
rz(-1.3686603) q[1];
sx q[1];
rz(-2.7006671) q[1];
x q[2];
rz(2.426126) q[3];
sx q[3];
rz(-0.56583929) q[3];
sx q[3];
rz(2.6239606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1019613) q[2];
sx q[2];
rz(-2.8122718) q[2];
sx q[2];
rz(-2.9015818) q[2];
rz(2.4890684) q[3];
sx q[3];
rz(-1.1381166) q[3];
sx q[3];
rz(1.8664546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08808247) q[0];
sx q[0];
rz(-1.9048012) q[0];
sx q[0];
rz(-0.85813338) q[0];
rz(-1.3872604) q[1];
sx q[1];
rz(-2.6871082) q[1];
sx q[1];
rz(2.8087356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3368206) q[0];
sx q[0];
rz(-0.66786486) q[0];
sx q[0];
rz(2.0940381) q[0];
rz(-pi) q[1];
rz(-1.6295241) q[2];
sx q[2];
rz(-1.8203805) q[2];
sx q[2];
rz(-1.6317612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7499125) q[1];
sx q[1];
rz(-1.9894162) q[1];
sx q[1];
rz(-1.7961572) q[1];
rz(2.694533) q[3];
sx q[3];
rz(-2.4363228) q[3];
sx q[3];
rz(-1.3652319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7994999) q[2];
sx q[2];
rz(-0.19793333) q[2];
sx q[2];
rz(1.0510772) q[2];
rz(1.6445232) q[3];
sx q[3];
rz(-1.9688508) q[3];
sx q[3];
rz(-2.3433949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77338162) q[0];
sx q[0];
rz(-2.1147275) q[0];
sx q[0];
rz(-0.89212242) q[0];
rz(0.46961531) q[1];
sx q[1];
rz(-0.62925595) q[1];
sx q[1];
rz(-2.3687252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.597201) q[0];
sx q[0];
rz(-0.8071476) q[0];
sx q[0];
rz(2.3788484) q[0];
rz(-0.30420423) q[2];
sx q[2];
rz(-2.4149564) q[2];
sx q[2];
rz(-1.3192504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6438021) q[1];
sx q[1];
rz(-1.1082768) q[1];
sx q[1];
rz(1.2116572) q[1];
x q[2];
rz(0.77247932) q[3];
sx q[3];
rz(-2.8101343) q[3];
sx q[3];
rz(0.14642388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5272687) q[2];
sx q[2];
rz(-1.9169151) q[2];
sx q[2];
rz(0.69250715) q[2];
rz(-0.45037371) q[3];
sx q[3];
rz(-1.9362484) q[3];
sx q[3];
rz(-1.6374121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5699128) q[0];
sx q[0];
rz(-0.79140651) q[0];
sx q[0];
rz(0.94863844) q[0];
rz(-2.0750849) q[1];
sx q[1];
rz(-0.98943168) q[1];
sx q[1];
rz(-1.271064) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61425754) q[0];
sx q[0];
rz(-1.7760881) q[0];
sx q[0];
rz(0.094747748) q[0];
rz(-pi) q[1];
rz(-1.5367212) q[2];
sx q[2];
rz(-1.5928895) q[2];
sx q[2];
rz(1.3349229) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3349325) q[1];
sx q[1];
rz(-2.1927823) q[1];
sx q[1];
rz(0.16652624) q[1];
x q[2];
rz(-0.27008121) q[3];
sx q[3];
rz(-2.684786) q[3];
sx q[3];
rz(2.1169259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54242414) q[2];
sx q[2];
rz(-1.8343265) q[2];
sx q[2];
rz(-3.001281) q[2];
rz(-3.1091651) q[3];
sx q[3];
rz(-2.2166538) q[3];
sx q[3];
rz(-1.2062581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5133544) q[0];
sx q[0];
rz(-1.4435377) q[0];
sx q[0];
rz(-2.3640609) q[0];
rz(-2.2041722) q[1];
sx q[1];
rz(-1.9098858) q[1];
sx q[1];
rz(2.6984528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454058) q[0];
sx q[0];
rz(-1.4898132) q[0];
sx q[0];
rz(0.085025351) q[0];
rz(-pi) q[1];
rz(0.97340092) q[2];
sx q[2];
rz(-1.8466766) q[2];
sx q[2];
rz(-1.2593002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3819246) q[1];
sx q[1];
rz(-0.72205359) q[1];
sx q[1];
rz(2.6284435) q[1];
rz(-1.0983766) q[3];
sx q[3];
rz(-0.34675035) q[3];
sx q[3];
rz(-0.43671331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0605269) q[2];
sx q[2];
rz(-2.492283) q[2];
sx q[2];
rz(-3.013986) q[2];
rz(-2.8429032) q[3];
sx q[3];
rz(-1.9689711) q[3];
sx q[3];
rz(-2.7711788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61113197) q[0];
sx q[0];
rz(-1.2584542) q[0];
sx q[0];
rz(-2.8897814) q[0];
rz(-0.61027377) q[1];
sx q[1];
rz(-0.88519575) q[1];
sx q[1];
rz(-2.0559678) q[1];
rz(-0.064432003) q[2];
sx q[2];
rz(-0.85307912) q[2];
sx q[2];
rz(-1.0882594) q[2];
rz(1.1518703) q[3];
sx q[3];
rz(-0.48696951) q[3];
sx q[3];
rz(1.3305668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
