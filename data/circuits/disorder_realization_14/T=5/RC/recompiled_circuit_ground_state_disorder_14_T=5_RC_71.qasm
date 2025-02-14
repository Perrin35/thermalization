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
rz(1.908778) q[0];
sx q[0];
rz(9.0969152) q[0];
rz(2.7479563) q[1];
sx q[1];
rz(-1.6180232) q[1];
sx q[1];
rz(-1.6021498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0293139) q[0];
sx q[0];
rz(-1.7718705) q[0];
sx q[0];
rz(-1.5657094) q[0];
rz(-pi) q[1];
rz(1.0225564) q[2];
sx q[2];
rz(-1.3371144) q[2];
sx q[2];
rz(1.3604239) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5757619) q[1];
sx q[1];
rz(-1.9565493) q[1];
sx q[1];
rz(-1.7931531) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76018192) q[3];
sx q[3];
rz(-2.6225445) q[3];
sx q[3];
rz(0.20582919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9615122) q[2];
sx q[2];
rz(-2.1891258) q[2];
sx q[2];
rz(-2.3251779) q[2];
rz(2.4444729) q[3];
sx q[3];
rz(-2.7643047) q[3];
sx q[3];
rz(-2.2009946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78696752) q[0];
sx q[0];
rz(-2.0836199) q[0];
sx q[0];
rz(-0.70607287) q[0];
rz(2.970703) q[1];
sx q[1];
rz(-2.6306174) q[1];
sx q[1];
rz(0.99542803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083849523) q[0];
sx q[0];
rz(-2.665976) q[0];
sx q[0];
rz(2.3887079) q[0];
rz(2.6508743) q[2];
sx q[2];
rz(-2.2123537) q[2];
sx q[2];
rz(3.0131755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4793206) q[1];
sx q[1];
rz(-0.26367981) q[1];
sx q[1];
rz(-0.94698833) q[1];
x q[2];
rz(-0.045142249) q[3];
sx q[3];
rz(-0.63796872) q[3];
sx q[3];
rz(-0.54813069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5848315) q[2];
sx q[2];
rz(-2.2274667) q[2];
sx q[2];
rz(0.27908984) q[2];
rz(2.4537405) q[3];
sx q[3];
rz(-0.22356859) q[3];
sx q[3];
rz(2.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2897709) q[0];
sx q[0];
rz(-2.2230447) q[0];
sx q[0];
rz(-2.7626792) q[0];
rz(1.9508349) q[1];
sx q[1];
rz(-0.36668101) q[1];
sx q[1];
rz(0.77622882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8541753) q[0];
sx q[0];
rz(-2.6834134) q[0];
sx q[0];
rz(-1.929438) q[0];
x q[1];
rz(-1.7227145) q[2];
sx q[2];
rz(-1.3891274) q[2];
sx q[2];
rz(2.7318397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6783733) q[1];
sx q[1];
rz(-1.3161945) q[1];
sx q[1];
rz(0.71779743) q[1];
rz(-pi) q[2];
rz(-1.6504694) q[3];
sx q[3];
rz(-2.6159379) q[3];
sx q[3];
rz(2.6826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0188521) q[2];
sx q[2];
rz(-1.2383702) q[2];
sx q[2];
rz(2.5416601) q[2];
rz(1.9020724) q[3];
sx q[3];
rz(-1.7299088) q[3];
sx q[3];
rz(2.5304573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83830738) q[0];
sx q[0];
rz(-2.4628283) q[0];
sx q[0];
rz(1.7093866) q[0];
rz(-2.2749061) q[1];
sx q[1];
rz(-1.7984093) q[1];
sx q[1];
rz(-2.0492882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019072559) q[0];
sx q[0];
rz(-2.2664223) q[0];
sx q[0];
rz(0.75059661) q[0];
rz(-pi) q[1];
rz(-2.714027) q[2];
sx q[2];
rz(-0.95631546) q[2];
sx q[2];
rz(-0.63860047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48627794) q[1];
sx q[1];
rz(-1.1064782) q[1];
sx q[1];
rz(-1.5686591) q[1];
x q[2];
rz(2.3471429) q[3];
sx q[3];
rz(-0.91728079) q[3];
sx q[3];
rz(0.78032035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0563353) q[2];
sx q[2];
rz(-1.1577572) q[2];
sx q[2];
rz(1.2752656) q[2];
rz(-1.7665675) q[3];
sx q[3];
rz(-1.9019889) q[3];
sx q[3];
rz(2.1084771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035868693) q[0];
sx q[0];
rz(-2.2585456) q[0];
sx q[0];
rz(-1.7778273) q[0];
rz(-0.56963244) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(-0.5035351) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0277594) q[0];
sx q[0];
rz(-0.22969023) q[0];
sx q[0];
rz(-0.57595797) q[0];
rz(-pi) q[1];
rz(2.3634144) q[2];
sx q[2];
rz(-0.89820895) q[2];
sx q[2];
rz(-1.1802117) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9787748) q[1];
sx q[1];
rz(-1.4895338) q[1];
sx q[1];
rz(-3.1291276) q[1];
x q[2];
rz(-0.77025099) q[3];
sx q[3];
rz(-1.7903084) q[3];
sx q[3];
rz(-1.3971922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86294952) q[2];
sx q[2];
rz(-0.38662616) q[2];
sx q[2];
rz(-1.53842) q[2];
rz(-0.1362416) q[3];
sx q[3];
rz(-1.4038266) q[3];
sx q[3];
rz(-1.2374102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49029008) q[0];
sx q[0];
rz(-2.5232115) q[0];
sx q[0];
rz(-0.34143099) q[0];
rz(0.731172) q[1];
sx q[1];
rz(-0.57594222) q[1];
sx q[1];
rz(-1.0715356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6302231) q[0];
sx q[0];
rz(-1.1611188) q[0];
sx q[0];
rz(2.0100586) q[0];
x q[1];
rz(-2.5265556) q[2];
sx q[2];
rz(-1.5205548) q[2];
sx q[2];
rz(-3.0453403) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2372969) q[1];
sx q[1];
rz(-2.0091506) q[1];
sx q[1];
rz(2.3890952) q[1];
x q[2];
rz(-1.066393) q[3];
sx q[3];
rz(-1.5747488) q[3];
sx q[3];
rz(2.469698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3557446) q[2];
sx q[2];
rz(-1.5037437) q[2];
sx q[2];
rz(-0.49622932) q[2];
rz(-1.5931386) q[3];
sx q[3];
rz(-1.691247) q[3];
sx q[3];
rz(-2.1350071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9736495) q[0];
sx q[0];
rz(-1.668546) q[0];
sx q[0];
rz(-0.014852511) q[0];
rz(-1.4133833) q[1];
sx q[1];
rz(-1.1057066) q[1];
sx q[1];
rz(2.3766439) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1434484) q[0];
sx q[0];
rz(-1.0113749) q[0];
sx q[0];
rz(-0.75152121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4582389) q[2];
sx q[2];
rz(-1.8477447) q[2];
sx q[2];
rz(-0.54613699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.98337631) q[1];
sx q[1];
rz(-1.1740654) q[1];
sx q[1];
rz(-1.5660172) q[1];
rz(-pi) q[2];
rz(2.9243605) q[3];
sx q[3];
rz(-1.9857556) q[3];
sx q[3];
rz(-0.46997786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0726274) q[2];
sx q[2];
rz(-1.5494538) q[2];
sx q[2];
rz(-0.7507945) q[2];
rz(2.1454861) q[3];
sx q[3];
rz(-2.3361497) q[3];
sx q[3];
rz(-0.44710818) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65425101) q[0];
sx q[0];
rz(-2.5387006) q[0];
sx q[0];
rz(0.61016369) q[0];
rz(0.24972406) q[1];
sx q[1];
rz(-1.6619253) q[1];
sx q[1];
rz(3.0009559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2785169) q[0];
sx q[0];
rz(-0.91916537) q[0];
sx q[0];
rz(2.5550575) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9570661) q[2];
sx q[2];
rz(-0.51287309) q[2];
sx q[2];
rz(-1.41092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.096123213) q[1];
sx q[1];
rz(-0.28633196) q[1];
sx q[1];
rz(-2.8827122) q[1];
rz(1.8406156) q[3];
sx q[3];
rz(-2.3128445) q[3];
sx q[3];
rz(-1.5480812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3934624) q[2];
sx q[2];
rz(-2.738214) q[2];
sx q[2];
rz(-2.051579) q[2];
rz(-1.2202834) q[3];
sx q[3];
rz(-2.0956495) q[3];
sx q[3];
rz(2.2530341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9991934) q[0];
sx q[0];
rz(-1.1141454) q[0];
sx q[0];
rz(-0.92495579) q[0];
rz(-1.6389182) q[1];
sx q[1];
rz(-1.0097367) q[1];
sx q[1];
rz(-0.22770539) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3229436) q[0];
sx q[0];
rz(-0.54185757) q[0];
sx q[0];
rz(1.3635554) q[0];
rz(-pi) q[1];
rz(1.8500344) q[2];
sx q[2];
rz(-1.0089356) q[2];
sx q[2];
rz(0.13915881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.046809) q[1];
sx q[1];
rz(-1.2168447) q[1];
sx q[1];
rz(-1.3665529) q[1];
x q[2];
rz(-2.5277522) q[3];
sx q[3];
rz(-1.7305529) q[3];
sx q[3];
rz(2.0789119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0324675) q[2];
sx q[2];
rz(-0.6296857) q[2];
sx q[2];
rz(2.6510748) q[2];
rz(2.0971175) q[3];
sx q[3];
rz(-0.46382469) q[3];
sx q[3];
rz(-0.0096376816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5661285) q[0];
sx q[0];
rz(-2.6818891) q[0];
sx q[0];
rz(-1.2658966) q[0];
rz(1.5220386) q[1];
sx q[1];
rz(-1.252389) q[1];
sx q[1];
rz(-0.44580805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8222919) q[0];
sx q[0];
rz(-2.5133488) q[0];
sx q[0];
rz(-0.99530812) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0580373) q[2];
sx q[2];
rz(-2.6284468) q[2];
sx q[2];
rz(1.4551077) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3109772) q[1];
sx q[1];
rz(-1.4339245) q[1];
sx q[1];
rz(-0.61424436) q[1];
rz(-pi) q[2];
rz(-3.1200459) q[3];
sx q[3];
rz(-1.8441219) q[3];
sx q[3];
rz(-0.97238982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1294935) q[2];
sx q[2];
rz(-2.1603277) q[2];
sx q[2];
rz(1.0274461) q[2];
rz(-1.3198613) q[3];
sx q[3];
rz(-0.3052932) q[3];
sx q[3];
rz(-0.52904883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32578596) q[0];
sx q[0];
rz(-1.3675714) q[0];
sx q[0];
rz(2.0969781) q[0];
rz(-2.2141937) q[1];
sx q[1];
rz(-1.2694042) q[1];
sx q[1];
rz(-0.70613695) q[1];
rz(0.32756424) q[2];
sx q[2];
rz(-1.2739999) q[2];
sx q[2];
rz(-1.3890541) q[2];
rz(1.8580243) q[3];
sx q[3];
rz(-2.1706222) q[3];
sx q[3];
rz(-2.2219347) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
