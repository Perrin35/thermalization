OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0377169) q[0];
sx q[0];
rz(-1.2020943) q[0];
sx q[0];
rz(-1.9934959) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(-1.3936477) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6205665) q[0];
sx q[0];
rz(-1.5728083) q[0];
sx q[0];
rz(2.8319671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.355028) q[2];
sx q[2];
rz(-0.85430356) q[2];
sx q[2];
rz(0.63209817) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4204971) q[1];
sx q[1];
rz(-1.8407397) q[1];
sx q[1];
rz(0.47081486) q[1];
rz(-pi) q[2];
rz(-2.1288793) q[3];
sx q[3];
rz(-1.3289684) q[3];
sx q[3];
rz(-1.1954607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48774886) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(-0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-3.0088186) q[0];
rz(1.4615387) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(0.24138385) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3791083) q[0];
sx q[0];
rz(-2.603841) q[0];
sx q[0];
rz(-2.3178029) q[0];
rz(-0.74626211) q[2];
sx q[2];
rz(-0.66821874) q[2];
sx q[2];
rz(1.8909188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.452255) q[1];
sx q[1];
rz(-0.66879767) q[1];
sx q[1];
rz(1.4237088) q[1];
x q[2];
rz(1.4199735) q[3];
sx q[3];
rz(-1.5966468) q[3];
sx q[3];
rz(-0.54102708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-1.1068809) q[2];
rz(1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7805507) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(-0.74044359) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(-2.0887451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42441472) q[0];
sx q[0];
rz(-2.4632235) q[0];
sx q[0];
rz(2.5413187) q[0];
rz(-2.8430976) q[2];
sx q[2];
rz(-2.2224991) q[2];
sx q[2];
rz(0.088108206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4098674) q[1];
sx q[1];
rz(-1.6907398) q[1];
sx q[1];
rz(2.9883595) q[1];
rz(-pi) q[2];
rz(3.1316109) q[3];
sx q[3];
rz(-2.0008893) q[3];
sx q[3];
rz(0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25049245) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(2.5946674) q[2];
rz(-2.8524103) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96994394) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(1.3522211) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073233152) q[0];
sx q[0];
rz(-1.7422361) q[0];
sx q[0];
rz(-1.839848) q[0];
rz(-2.7912555) q[2];
sx q[2];
rz(-1.0812034) q[2];
sx q[2];
rz(-1.1794773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(3.1349036) q[1];
x q[2];
rz(-1.719791) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(-1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(2.5881055) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(2.4374403) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(0.59590894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0740972) q[0];
sx q[0];
rz(-1.4795545) q[0];
sx q[0];
rz(-1.3542961) q[0];
rz(-pi) q[1];
rz(0.67363588) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(3.0631531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51922819) q[1];
sx q[1];
rz(-1.16751) q[1];
sx q[1];
rz(-1.7720023) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7549873) q[3];
sx q[3];
rz(-1.9362238) q[3];
sx q[3];
rz(-2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(-0.55111432) q[2];
rz(-0.20714949) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(2.1898848) q[0];
rz(-0.47479182) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808522) q[0];
sx q[0];
rz(-1.7096925) q[0];
sx q[0];
rz(1.172658) q[0];
rz(-pi) q[1];
rz(2.8304808) q[2];
sx q[2];
rz(-1.7179486) q[2];
sx q[2];
rz(-1.0263024) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5438248) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(-2.3901229) q[1];
rz(-pi) q[2];
rz(1.4548364) q[3];
sx q[3];
rz(-2.6280858) q[3];
sx q[3];
rz(-0.49125571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(0.93377101) q[2];
rz(1.6823403) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(-2.4767955) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-2.738293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032265183) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(-1.1778675) q[0];
x q[1];
rz(2.7935739) q[2];
sx q[2];
rz(-2.4325779) q[2];
sx q[2];
rz(1.6955171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94782103) q[1];
sx q[1];
rz(-1.3606451) q[1];
sx q[1];
rz(-3.139688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7014916) q[3];
sx q[3];
rz(-0.19154597) q[3];
sx q[3];
rz(2.6768315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91281259) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(2.1408391) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1338761) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(0.49466053) q[0];
rz(1.6330632) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(0.60044926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8884044) q[0];
sx q[0];
rz(-2.0300403) q[0];
sx q[0];
rz(-2.0183536) q[0];
x q[1];
rz(1.0457439) q[2];
sx q[2];
rz(-1.5935491) q[2];
sx q[2];
rz(-0.75726985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2272051) q[1];
sx q[1];
rz(-1.3890508) q[1];
sx q[1];
rz(1.6919796) q[1];
x q[2];
rz(0.30386691) q[3];
sx q[3];
rz(-1.3715203) q[3];
sx q[3];
rz(0.58135939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0042469) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329926) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(2.7240662) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5696213) q[0];
sx q[0];
rz(-0.53335359) q[0];
sx q[0];
rz(-1.4825975) q[0];
rz(-1.3067901) q[2];
sx q[2];
rz(-1.4679113) q[2];
sx q[2];
rz(-1.4357944) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36996499) q[1];
sx q[1];
rz(-1.4039478) q[1];
sx q[1];
rz(0.19373993) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0017654) q[3];
sx q[3];
rz(-1.7472072) q[3];
sx q[3];
rz(2.924502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(2.647906) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(-0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(2.0102274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2162227) q[0];
sx q[0];
rz(-2.4664481) q[0];
sx q[0];
rz(2.569909) q[0];
rz(-pi) q[1];
rz(-2.4234424) q[2];
sx q[2];
rz(-0.8469204) q[2];
sx q[2];
rz(1.3494327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86913925) q[1];
sx q[1];
rz(-1.9483882) q[1];
sx q[1];
rz(1.3294405) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0919308) q[3];
sx q[3];
rz(-2.1225956) q[3];
sx q[3];
rz(2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3161105) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(2.4342009) q[2];
rz(2.3317544) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.390776) q[2];
sx q[2];
rz(-2.3859947) q[2];
sx q[2];
rz(1.044556) q[2];
rz(2.2255696) q[3];
sx q[3];
rz(-2.3229204) q[3];
sx q[3];
rz(-1.7182072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
