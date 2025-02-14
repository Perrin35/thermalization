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
rz(2.7479563) q[1];
sx q[1];
rz(-1.6180232) q[1];
sx q[1];
rz(-1.6021498) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868139) q[0];
sx q[0];
rz(-0.20113763) q[0];
sx q[0];
rz(-0.024951713) q[0];
x q[1];
rz(-1.9992158) q[2];
sx q[2];
rz(-2.5503469) q[2];
sx q[2];
rz(-2.5687886) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92009993) q[1];
sx q[1];
rz(-1.3650238) q[1];
sx q[1];
rz(-0.39452067) q[1];
rz(-pi) q[2];
rz(-2.3814107) q[3];
sx q[3];
rz(-2.6225445) q[3];
sx q[3];
rz(0.20582919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1800805) q[2];
sx q[2];
rz(-0.95246685) q[2];
sx q[2];
rz(2.3251779) q[2];
rz(2.4444729) q[3];
sx q[3];
rz(-2.7643047) q[3];
sx q[3];
rz(0.94059801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78696752) q[0];
sx q[0];
rz(-1.0579728) q[0];
sx q[0];
rz(0.70607287) q[0];
rz(-0.1708897) q[1];
sx q[1];
rz(-0.51097521) q[1];
sx q[1];
rz(2.1461646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1814897) q[0];
sx q[0];
rz(-1.8892291) q[0];
sx q[0];
rz(2.7820827) q[0];
rz(-pi) q[1];
rz(-2.1336251) q[2];
sx q[2];
rz(-0.78608222) q[2];
sx q[2];
rz(-0.85725923) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4793206) q[1];
sx q[1];
rz(-2.8779128) q[1];
sx q[1];
rz(-0.94698833) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.045142249) q[3];
sx q[3];
rz(-0.63796872) q[3];
sx q[3];
rz(-0.54813069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55676111) q[2];
sx q[2];
rz(-0.91412592) q[2];
sx q[2];
rz(-0.27908984) q[2];
rz(2.4537405) q[3];
sx q[3];
rz(-0.22356859) q[3];
sx q[3];
rz(-1.0031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85182178) q[0];
sx q[0];
rz(-2.2230447) q[0];
sx q[0];
rz(2.7626792) q[0];
rz(1.9508349) q[1];
sx q[1];
rz(-2.7749116) q[1];
sx q[1];
rz(2.3653638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8541753) q[0];
sx q[0];
rz(-0.45817927) q[0];
sx q[0];
rz(-1.2121546) q[0];
rz(-pi) q[1];
rz(0.68910901) q[2];
sx q[2];
rz(-2.9053134) q[2];
sx q[2];
rz(2.8483727) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8175428) q[1];
sx q[1];
rz(-0.88081283) q[1];
sx q[1];
rz(1.903456) q[1];
x q[2];
rz(-2.0950731) q[3];
sx q[3];
rz(-1.5308497) q[3];
sx q[3];
rz(2.0986882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0188521) q[2];
sx q[2];
rz(-1.2383702) q[2];
sx q[2];
rz(0.59993258) q[2];
rz(-1.9020724) q[3];
sx q[3];
rz(-1.4116838) q[3];
sx q[3];
rz(-0.6111353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3032853) q[0];
sx q[0];
rz(-2.4628283) q[0];
sx q[0];
rz(-1.432206) q[0];
rz(-0.86668658) q[1];
sx q[1];
rz(-1.3431834) q[1];
sx q[1];
rz(-2.0492882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0129654) q[0];
sx q[0];
rz(-2.1218637) q[0];
sx q[0];
rz(-2.4222213) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42756568) q[2];
sx q[2];
rz(-0.95631546) q[2];
sx q[2];
rz(0.63860047) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6553147) q[1];
sx q[1];
rz(-2.0351145) q[1];
sx q[1];
rz(-1.5729335) q[1];
rz(-pi) q[2];
rz(-2.400551) q[3];
sx q[3];
rz(-0.96858409) q[3];
sx q[3];
rz(-2.9054537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0852574) q[2];
sx q[2];
rz(-1.1577572) q[2];
sx q[2];
rz(1.866327) q[2];
rz(1.7665675) q[3];
sx q[3];
rz(-1.9019889) q[3];
sx q[3];
rz(-2.1084771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.105724) q[0];
sx q[0];
rz(-0.88304702) q[0];
sx q[0];
rz(1.7778273) q[0];
rz(-2.5719602) q[1];
sx q[1];
rz(-0.49915794) q[1];
sx q[1];
rz(-2.6380576) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1138332) q[0];
sx q[0];
rz(-2.9119024) q[0];
sx q[0];
rz(0.57595797) q[0];
rz(-pi) q[1];
rz(0.72959186) q[2];
sx q[2];
rz(-0.9895095) q[2];
sx q[2];
rz(-0.94129291) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40696661) q[1];
sx q[1];
rz(-1.5583724) q[1];
sx q[1];
rz(1.4895275) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31007575) q[3];
sx q[3];
rz(-2.3469124) q[3];
sx q[3];
rz(-0.39439699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2786431) q[2];
sx q[2];
rz(-2.7549665) q[2];
sx q[2];
rz(1.53842) q[2];
rz(0.1362416) q[3];
sx q[3];
rz(-1.4038266) q[3];
sx q[3];
rz(1.2374102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513026) q[0];
sx q[0];
rz(-2.5232115) q[0];
sx q[0];
rz(-2.8001617) q[0];
rz(-0.731172) q[1];
sx q[1];
rz(-2.5656504) q[1];
sx q[1];
rz(2.0700571) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6302231) q[0];
sx q[0];
rz(-1.9804738) q[0];
sx q[0];
rz(2.0100586) q[0];
rz(-pi) q[1];
rz(-3.0546636) q[2];
sx q[2];
rz(-0.61682103) q[2];
sx q[2];
rz(-1.4035743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9042957) q[1];
sx q[1];
rz(-1.132442) q[1];
sx q[1];
rz(-0.75249747) q[1];
x q[2];
rz(1.066393) q[3];
sx q[3];
rz(-1.5747488) q[3];
sx q[3];
rz(-2.469698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3557446) q[2];
sx q[2];
rz(-1.5037437) q[2];
sx q[2];
rz(-2.6453633) q[2];
rz(1.5931386) q[3];
sx q[3];
rz(-1.4503456) q[3];
sx q[3];
rz(-2.1350071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16794311) q[0];
sx q[0];
rz(-1.4730467) q[0];
sx q[0];
rz(0.014852511) q[0];
rz(-1.7282093) q[1];
sx q[1];
rz(-2.035886) q[1];
sx q[1];
rz(-0.76494876) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1434484) q[0];
sx q[0];
rz(-2.1302178) q[0];
sx q[0];
rz(2.3900714) q[0];
x q[1];
rz(0.68335376) q[2];
sx q[2];
rz(-1.2938479) q[2];
sx q[2];
rz(0.54613699) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5892667) q[1];
sx q[1];
rz(-1.5663884) q[1];
sx q[1];
rz(-2.7448576) q[1];
rz(-1.9946003) q[3];
sx q[3];
rz(-1.3722581) q[3];
sx q[3];
rz(1.0120704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0726274) q[2];
sx q[2];
rz(-1.5921389) q[2];
sx q[2];
rz(2.3907982) q[2];
rz(2.1454861) q[3];
sx q[3];
rz(-2.3361497) q[3];
sx q[3];
rz(2.6944845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873416) q[0];
sx q[0];
rz(-2.5387006) q[0];
sx q[0];
rz(0.61016369) q[0];
rz(0.24972406) q[1];
sx q[1];
rz(-1.4796673) q[1];
sx q[1];
rz(-3.0009559) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2785169) q[0];
sx q[0];
rz(-2.2224273) q[0];
sx q[0];
rz(-0.58653513) q[0];
rz(-1.6737559) q[2];
sx q[2];
rz(-2.0741346) q[2];
sx q[2];
rz(-1.1998985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.096123213) q[1];
sx q[1];
rz(-2.8552607) q[1];
sx q[1];
rz(2.8827122) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28293149) q[3];
sx q[3];
rz(-2.3608876) q[3];
sx q[3];
rz(-1.9819575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7481302) q[2];
sx q[2];
rz(-0.40337864) q[2];
sx q[2];
rz(2.051579) q[2];
rz(-1.9213093) q[3];
sx q[3];
rz(-1.0459432) q[3];
sx q[3];
rz(-0.88855851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1423993) q[0];
sx q[0];
rz(-1.1141454) q[0];
sx q[0];
rz(-0.92495579) q[0];
rz(-1.6389182) q[1];
sx q[1];
rz(-1.0097367) q[1];
sx q[1];
rz(2.9138873) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3229436) q[0];
sx q[0];
rz(-0.54185757) q[0];
sx q[0];
rz(-1.7780373) q[0];
rz(-pi) q[1];
rz(-2.7289189) q[2];
sx q[2];
rz(-2.5209171) q[2];
sx q[2];
rz(-0.35457573) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.633435) q[1];
sx q[1];
rz(-2.7351008) q[1];
sx q[1];
rz(0.50199957) q[1];
x q[2];
rz(-1.7654159) q[3];
sx q[3];
rz(-2.1756919) q[3];
sx q[3];
rz(-2.5218487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0324675) q[2];
sx q[2];
rz(-2.511907) q[2];
sx q[2];
rz(-2.6510748) q[2];
rz(-1.0444752) q[3];
sx q[3];
rz(-0.46382469) q[3];
sx q[3];
rz(-0.0096376816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5661285) q[0];
sx q[0];
rz(-2.6818891) q[0];
sx q[0];
rz(1.2658966) q[0];
rz(1.6195541) q[1];
sx q[1];
rz(-1.252389) q[1];
sx q[1];
rz(0.44580805) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73484008) q[0];
sx q[0];
rz(-1.8963843) q[0];
sx q[0];
rz(2.1181137) q[0];
rz(0.083555315) q[2];
sx q[2];
rz(-0.5131459) q[2];
sx q[2];
rz(1.6864849) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3109772) q[1];
sx q[1];
rz(-1.4339245) q[1];
sx q[1];
rz(-2.5273483) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.021546797) q[3];
sx q[3];
rz(-1.2974707) q[3];
sx q[3];
rz(2.1692028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1294935) q[2];
sx q[2];
rz(-2.1603277) q[2];
sx q[2];
rz(-2.1141466) q[2];
rz(1.8217314) q[3];
sx q[3];
rz(-2.8362995) q[3];
sx q[3];
rz(0.52904883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.2583707) q[2];
sx q[2];
rz(-1.2580522) q[2];
sx q[2];
rz(-3.0589041) q[2];
rz(-0.61946034) q[3];
sx q[3];
rz(-1.8068228) q[3];
sx q[3];
rz(-0.48590776) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
