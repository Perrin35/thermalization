OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(1.9934959) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(-2.2009067) q[1];
sx q[1];
rz(1.3936477) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049126547) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(1.5729088) q[0];
rz(-pi) q[1];
rz(-0.78656466) q[2];
sx q[2];
rz(-0.85430356) q[2];
sx q[2];
rz(2.5094945) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72109556) q[1];
sx q[1];
rz(-1.8407397) q[1];
sx q[1];
rz(2.6707778) q[1];
rz(2.0066891) q[3];
sx q[3];
rz(-2.5385058) q[3];
sx q[3];
rz(-3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(3.0207108) q[2];
rz(-0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(3.0088186) q[0];
rz(-1.6800539) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-0.24138385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.1855159) q[0];
sx q[0];
rz(-0.38498621) q[0];
rz(-pi) q[1];
rz(1.0788467) q[2];
sx q[2];
rz(-1.098512) q[2];
sx q[2];
rz(0.3837331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0028210359) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(0.9072733) q[1];
rz(-pi) q[2];
rz(1.7216191) q[3];
sx q[3];
rz(-1.5966468) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(-1.0528475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2967865) q[0];
sx q[0];
rz(-1.0264945) q[0];
sx q[0];
rz(1.143572) q[0];
x q[1];
rz(-2.244433) q[2];
sx q[2];
rz(-1.8067915) q[2];
sx q[2];
rz(1.8434075) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3215072) q[1];
sx q[1];
rz(-2.9472889) q[1];
sx q[1];
rz(-0.66837515) q[1];
rz(-pi) q[2];
rz(1.5490407) q[3];
sx q[3];
rz(-0.43020159) q[3];
sx q[3];
rz(0.83963001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25049245) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-2.5946674) q[2];
rz(-2.8524103) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(-3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1716487) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.3522211) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073233152) q[0];
sx q[0];
rz(-1.3993565) q[0];
sx q[0];
rz(-1.839848) q[0];
rz(-pi) q[1];
rz(0.9985853) q[2];
sx q[2];
rz(-0.59362715) q[2];
sx q[2];
rz(-1.8400536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.786495) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(3.1349036) q[1];
x q[2];
rz(-1.4218016) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(-1.6001584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2410879) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(-2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6699162) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(-2.4374403) q[0];
rz(-2.0856805) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-2.5456837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0376315) q[0];
sx q[0];
rz(-0.2346633) q[0];
sx q[0];
rz(-1.9734567) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17275177) q[2];
sx q[2];
rz(-2.4606332) q[2];
sx q[2];
rz(-1.7839884) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.039918598) q[1];
sx q[1];
rz(-2.6933751) q[1];
sx q[1];
rz(-2.7035473) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3866053) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(2.3900677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(2.5904783) q[2];
rz(-2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15774396) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(2.8463083) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8897032) q[0];
sx q[0];
rz(-1.9648874) q[0];
sx q[0];
rz(-2.9910812) q[0];
rz(1.7252543) q[2];
sx q[2];
rz(-1.2631577) q[2];
sx q[2];
rz(2.6442106) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5438248) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(-2.3901229) q[1];
rz(-pi) q[2];
rz(1.6867562) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(-0.49125571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(-1.4592524) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(2.899535) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(2.738293) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032265183) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(1.9637252) q[0];
x q[1];
rz(-2.7935739) q[2];
sx q[2];
rz(-2.4325779) q[2];
sx q[2];
rz(1.4460756) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1937716) q[1];
sx q[1];
rz(-1.3606451) q[1];
sx q[1];
rz(0.0019046849) q[1];
rz(0.025267406) q[3];
sx q[3];
rz(-1.7606887) q[3];
sx q[3];
rz(2.5437298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(0.51030695) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(2.5411434) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8884044) q[0];
sx q[0];
rz(-1.1115523) q[0];
sx q[0];
rz(-2.0183536) q[0];
rz(-1.0457439) q[2];
sx q[2];
rz(-1.5935491) q[2];
sx q[2];
rz(-2.3843228) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3214896) q[1];
sx q[1];
rz(-0.21807018) q[1];
sx q[1];
rz(0.58184187) q[1];
rz(-pi) q[2];
rz(1.7793711) q[3];
sx q[3];
rz(-1.2731291) q[3];
sx q[3];
rz(1.0514333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(3.0252769) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.6631888) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(2.7240662) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0644181) q[0];
sx q[0];
rz(-1.6155956) q[0];
sx q[0];
rz(2.1024465) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0350424) q[2];
sx q[2];
rz(-1.3082192) q[2];
sx q[2];
rz(-3.0343461) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7716277) q[1];
sx q[1];
rz(-1.4039478) q[1];
sx q[1];
rz(2.9478527) q[1];
rz(-pi) q[2];
rz(1.1674676) q[3];
sx q[3];
rz(-2.6780193) q[3];
sx q[3];
rz(1.7183255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4801165) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(0.49368668) q[2];
rz(-0.72475973) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17393728) q[0];
sx q[0];
rz(-0.65615654) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(-1.1313653) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5269276) q[0];
sx q[0];
rz(-2.1242495) q[0];
sx q[0];
rz(1.1620031) q[0];
rz(-pi) q[1];
rz(-2.2106902) q[2];
sx q[2];
rz(-0.97133342) q[2];
sx q[2];
rz(0.42720366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8611421) q[1];
sx q[1];
rz(-2.6965953) q[1];
sx q[1];
rz(0.54235561) q[1];
rz(-pi) q[2];
rz(-1.0919308) q[3];
sx q[3];
rz(-1.0189971) q[3];
sx q[3];
rz(-2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0512882) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.5384197) q[2];
sx q[2];
rz(-2.0576253) q[2];
sx q[2];
rz(-1.1228592) q[2];
rz(-0.91602305) q[3];
sx q[3];
rz(-2.3229204) q[3];
sx q[3];
rz(-1.7182072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];