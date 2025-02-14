OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3353424) q[0];
sx q[0];
rz(3.5080533) q[0];
sx q[0];
rz(8.9656497) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(-1.4067283) q[1];
sx q[1];
rz(-0.23174098) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329862) q[0];
sx q[0];
rz(-2.7911515) q[0];
sx q[0];
rz(0.20163433) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.628078) q[2];
sx q[2];
rz(-2.1165032) q[2];
sx q[2];
rz(-1.6451665) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6218839) q[1];
sx q[1];
rz(-0.95920282) q[1];
sx q[1];
rz(-2.0382463) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5292589) q[3];
sx q[3];
rz(-1.9088863) q[3];
sx q[3];
rz(1.6281782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.090652466) q[2];
sx q[2];
rz(-2.6068164) q[2];
sx q[2];
rz(-1.862662) q[2];
rz(-0.88979641) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5559674) q[0];
sx q[0];
rz(-0.90921679) q[0];
sx q[0];
rz(-2.2157748) q[0];
rz(-1.0649118) q[1];
sx q[1];
rz(-1.5644667) q[1];
sx q[1];
rz(-0.085478641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8259698) q[0];
sx q[0];
rz(-1.1809826) q[0];
sx q[0];
rz(1.5672383) q[0];
rz(2.7932554) q[2];
sx q[2];
rz(-0.98789633) q[2];
sx q[2];
rz(-2.9027241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6798415) q[1];
sx q[1];
rz(-1.0437168) q[1];
sx q[1];
rz(-0.26694571) q[1];
x q[2];
rz(-2.45473) q[3];
sx q[3];
rz(-1.9429038) q[3];
sx q[3];
rz(-0.61249477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3257137) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(-2.7549287) q[2];
rz(1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(0.2200505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81095186) q[0];
sx q[0];
rz(-2.7702259) q[0];
sx q[0];
rz(-0.49705848) q[0];
rz(1.0423543) q[1];
sx q[1];
rz(-2.1940239) q[1];
sx q[1];
rz(0.29464468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9884856) q[0];
sx q[0];
rz(-1.9026347) q[0];
sx q[0];
rz(1.6436623) q[0];
rz(2.5025236) q[2];
sx q[2];
rz(-1.7099755) q[2];
sx q[2];
rz(0.39060171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7859438) q[1];
sx q[1];
rz(-0.79141599) q[1];
sx q[1];
rz(-0.9407465) q[1];
x q[2];
rz(2.9975843) q[3];
sx q[3];
rz(-2.1668808) q[3];
sx q[3];
rz(-0.86897396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7872539) q[2];
sx q[2];
rz(-0.86091176) q[2];
sx q[2];
rz(-1.5617237) q[2];
rz(-1.076237) q[3];
sx q[3];
rz(-1.8393686) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(2.9122747) q[0];
rz(-1.5062821) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(0.062072676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44767013) q[0];
sx q[0];
rz(-1.0167443) q[0];
sx q[0];
rz(2.5497132) q[0];
rz(-0.96458413) q[2];
sx q[2];
rz(-2.6754489) q[2];
sx q[2];
rz(-3.0678444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8530555) q[1];
sx q[1];
rz(-1.848135) q[1];
sx q[1];
rz(-0.71373516) q[1];
x q[2];
rz(0.2009521) q[3];
sx q[3];
rz(-1.6621328) q[3];
sx q[3];
rz(0.82543711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.092992358) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(1.830706) q[2];
rz(3.1033031) q[3];
sx q[3];
rz(-1.7891276) q[3];
sx q[3];
rz(-2.766975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467317) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(-2.2485961) q[0];
rz(2.112174) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(-0.095887862) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5487089) q[0];
sx q[0];
rz(-1.0619231) q[0];
sx q[0];
rz(2.6611317) q[0];
rz(-pi) q[1];
rz(-0.070438373) q[2];
sx q[2];
rz(-1.6595006) q[2];
sx q[2];
rz(-2.8005637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.694987) q[1];
sx q[1];
rz(-2.1220653) q[1];
sx q[1];
rz(-1.1405844) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6259817) q[3];
sx q[3];
rz(-1.2429951) q[3];
sx q[3];
rz(-2.9189381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9868077) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(-0.42547697) q[2];
rz(-2.7191539) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285889) q[0];
sx q[0];
rz(-2.509403) q[0];
sx q[0];
rz(-3.0244306) q[0];
rz(1.4940184) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(2.07043) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2881831) q[0];
sx q[0];
rz(-2.5700535) q[0];
sx q[0];
rz(-0.54067143) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6980459) q[2];
sx q[2];
rz(-0.58813349) q[2];
sx q[2];
rz(-3.0556222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6630604) q[1];
sx q[1];
rz(-2.2473865) q[1];
sx q[1];
rz(0.29037906) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4070687) q[3];
sx q[3];
rz(-1.8843972) q[3];
sx q[3];
rz(-1.1160451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5693207) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(-0.039483698) q[2];
rz(0.076400541) q[3];
sx q[3];
rz(-1.1698086) q[3];
sx q[3];
rz(2.6881645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.4385248) q[0];
sx q[0];
rz(-0.81145966) q[0];
sx q[0];
rz(2.1912498) q[0];
rz(1.1514661) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(-1.3075525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3111585) q[0];
sx q[0];
rz(-1.3518999) q[0];
sx q[0];
rz(-1.2375184) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9588863) q[2];
sx q[2];
rz(-0.38190834) q[2];
sx q[2];
rz(-2.4289102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2491409) q[1];
sx q[1];
rz(-2.0607407) q[1];
sx q[1];
rz(2.9016414) q[1];
x q[2];
rz(0.96486196) q[3];
sx q[3];
rz(-1.6814702) q[3];
sx q[3];
rz(2.0697274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3211956) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(-2.5879228) q[2];
rz(0.6959483) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.11947908) q[0];
sx q[0];
rz(-1.6895634) q[0];
sx q[0];
rz(-2.8346862) q[0];
rz(1.8652929) q[1];
sx q[1];
rz(-0.46638322) q[1];
sx q[1];
rz(-1.9006405) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0690445) q[0];
sx q[0];
rz(-0.47397754) q[0];
sx q[0];
rz(-1.2233673) q[0];
rz(-pi) q[1];
rz(2.5503301) q[2];
sx q[2];
rz(-2.4545672) q[2];
sx q[2];
rz(-2.9369773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7684264) q[1];
sx q[1];
rz(-1.2618999) q[1];
sx q[1];
rz(-1.7684446) q[1];
rz(1.1697024) q[3];
sx q[3];
rz(-2.1425284) q[3];
sx q[3];
rz(-0.16756646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3500195) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(-1.8512858) q[2];
rz(2.4604515) q[3];
sx q[3];
rz(-0.9001503) q[3];
sx q[3];
rz(1.6397938) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660368) q[0];
sx q[0];
rz(-1.0365726) q[0];
sx q[0];
rz(-1.1736897) q[0];
rz(-0.58700079) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(0.079708286) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147238) q[0];
sx q[0];
rz(-2.7023315) q[0];
sx q[0];
rz(-1.0160116) q[0];
x q[1];
rz(-0.22144145) q[2];
sx q[2];
rz(-1.6802854) q[2];
sx q[2];
rz(1.8851282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4707798) q[1];
sx q[1];
rz(-0.97320405) q[1];
sx q[1];
rz(2.2714991) q[1];
rz(-pi) q[2];
rz(0.2853197) q[3];
sx q[3];
rz(-0.82672182) q[3];
sx q[3];
rz(-0.43275012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6012663) q[2];
sx q[2];
rz(-1.9129632) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(2.9529115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(0.48225668) q[0];
sx q[0];
rz(-1.5784669) q[0];
sx q[0];
rz(-1.851086) q[0];
rz(0.036529649) q[1];
sx q[1];
rz(-1.1643658) q[1];
sx q[1];
rz(-2.0972924) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1636476) q[0];
sx q[0];
rz(-0.53102101) q[0];
sx q[0];
rz(0.95300302) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2453111) q[2];
sx q[2];
rz(-2.7283629) q[2];
sx q[2];
rz(0.21597029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3022233) q[1];
sx q[1];
rz(-1.7249134) q[1];
sx q[1];
rz(0.83854143) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32907069) q[3];
sx q[3];
rz(-1.6794723) q[3];
sx q[3];
rz(1.4377126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9145987) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(1.8168137) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(2.2911086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.666438) q[0];
sx q[0];
rz(-1.403724) q[0];
sx q[0];
rz(-2.4028461) q[0];
rz(1.4229763) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(-1.1463317) q[2];
sx q[2];
rz(-0.88256114) q[2];
sx q[2];
rz(0.12086856) q[2];
rz(2.041009) q[3];
sx q[3];
rz(-2.1775424) q[3];
sx q[3];
rz(1.2829124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
