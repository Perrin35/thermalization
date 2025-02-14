OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8213788) q[0];
sx q[0];
rz(-2.3795405) q[0];
sx q[0];
rz(2.2574545) q[0];
rz(2.60131) q[1];
sx q[1];
rz(-1.5457431) q[1];
sx q[1];
rz(-2.5084612) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54612371) q[0];
sx q[0];
rz(-1.7350539) q[0];
sx q[0];
rz(1.6857273) q[0];
rz(-pi) q[1];
rz(0.48670003) q[2];
sx q[2];
rz(-1.4632157) q[2];
sx q[2];
rz(0.24525596) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3772845) q[1];
sx q[1];
rz(-1.2929243) q[1];
sx q[1];
rz(-1.4382526) q[1];
x q[2];
rz(-1.9820561) q[3];
sx q[3];
rz(-1.0602128) q[3];
sx q[3];
rz(3.0298281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4151609) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(-0.66698709) q[2];
rz(-0.64217448) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(-1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17382962) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(-0.94671384) q[0];
rz(3.0831378) q[1];
sx q[1];
rz(-0.36403257) q[1];
sx q[1];
rz(0.91190016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13214853) q[0];
sx q[0];
rz(-1.6086888) q[0];
sx q[0];
rz(-1.1542341) q[0];
rz(3.0709708) q[2];
sx q[2];
rz(-2.2433519) q[2];
sx q[2];
rz(2.9333851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2649041) q[1];
sx q[1];
rz(-2.0177245) q[1];
sx q[1];
rz(-1.3012215) q[1];
rz(1.6227609) q[3];
sx q[3];
rz(-2.012019) q[3];
sx q[3];
rz(-3.0558464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45137206) q[2];
sx q[2];
rz(-0.25190941) q[2];
sx q[2];
rz(2.8042262) q[2];
rz(-0.99006027) q[3];
sx q[3];
rz(-1.8127352) q[3];
sx q[3];
rz(-1.9199853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5628691) q[0];
sx q[0];
rz(-2.8674485) q[0];
sx q[0];
rz(1.8291871) q[0];
rz(-0.21526543) q[1];
sx q[1];
rz(-1.637641) q[1];
sx q[1];
rz(-2.2832835) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2374683) q[0];
sx q[0];
rz(-1.6209813) q[0];
sx q[0];
rz(-1.8367358) q[0];
rz(-pi) q[1];
rz(-0.22521727) q[2];
sx q[2];
rz(-1.2379939) q[2];
sx q[2];
rz(1.0343347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2006113) q[1];
sx q[1];
rz(-1.1340144) q[1];
sx q[1];
rz(1.753356) q[1];
x q[2];
rz(-0.86238548) q[3];
sx q[3];
rz(-0.38194712) q[3];
sx q[3];
rz(-2.2464584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0158318) q[2];
sx q[2];
rz(-2.5845926) q[2];
sx q[2];
rz(2.1324943) q[2];
rz(2.5741746) q[3];
sx q[3];
rz(-1.5083454) q[3];
sx q[3];
rz(-1.6217888) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384945) q[0];
sx q[0];
rz(-2.4479285) q[0];
sx q[0];
rz(-0.19288119) q[0];
rz(-2.0774972) q[1];
sx q[1];
rz(-0.27609438) q[1];
sx q[1];
rz(-2.1853133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.372626) q[0];
sx q[0];
rz(-1.3947693) q[0];
sx q[0];
rz(2.5123216) q[0];
x q[1];
rz(-0.85727875) q[2];
sx q[2];
rz(-2.3403717) q[2];
sx q[2];
rz(0.73357936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6156441) q[1];
sx q[1];
rz(-0.61405776) q[1];
sx q[1];
rz(0.55974069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2364081) q[3];
sx q[3];
rz(-1.815261) q[3];
sx q[3];
rz(-2.8595231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3095779) q[2];
sx q[2];
rz(-2.7509406) q[2];
sx q[2];
rz(2.2008956) q[2];
rz(2.7637123) q[3];
sx q[3];
rz(-0.89972073) q[3];
sx q[3];
rz(0.64062029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85510319) q[0];
sx q[0];
rz(-1.531895) q[0];
sx q[0];
rz(1.3252277) q[0];
rz(-0.27769604) q[1];
sx q[1];
rz(-0.9869906) q[1];
sx q[1];
rz(1.9353386) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078344854) q[0];
sx q[0];
rz(-1.4916911) q[0];
sx q[0];
rz(1.7633171) q[0];
rz(-pi) q[1];
rz(0.61111738) q[2];
sx q[2];
rz(-1.9131994) q[2];
sx q[2];
rz(0.26027203) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5531471) q[1];
sx q[1];
rz(-0.35970769) q[1];
sx q[1];
rz(1.5313696) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24171452) q[3];
sx q[3];
rz(-1.7346584) q[3];
sx q[3];
rz(3.0047447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9557934) q[2];
sx q[2];
rz(-2.9904371) q[2];
sx q[2];
rz(-1.7929057) q[2];
rz(-2.1336613) q[3];
sx q[3];
rz(-1.734684) q[3];
sx q[3];
rz(0.036788363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0607097) q[0];
sx q[0];
rz(-3.1302629) q[0];
sx q[0];
rz(-0.19705535) q[0];
rz(-1.9673678) q[1];
sx q[1];
rz(-2.2815727) q[1];
sx q[1];
rz(0.42541447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0393164) q[0];
sx q[0];
rz(-2.031306) q[0];
sx q[0];
rz(1.5251681) q[0];
rz(-pi) q[1];
rz(0.043656761) q[2];
sx q[2];
rz(-1.8746286) q[2];
sx q[2];
rz(2.6921261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96996355) q[1];
sx q[1];
rz(-0.49996129) q[1];
sx q[1];
rz(0.5699711) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2813536) q[3];
sx q[3];
rz(-1.0061534) q[3];
sx q[3];
rz(2.3535239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60874385) q[2];
sx q[2];
rz(-0.11561919) q[2];
sx q[2];
rz(0.70333424) q[2];
rz(-0.15130791) q[3];
sx q[3];
rz(-1.9688828) q[3];
sx q[3];
rz(-2.5942588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828101) q[0];
sx q[0];
rz(-0.8081688) q[0];
sx q[0];
rz(-0.95635828) q[0];
rz(-0.76902485) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(1.0395315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7881986) q[0];
sx q[0];
rz(-2.0457342) q[0];
sx q[0];
rz(0.026508412) q[0];
rz(-pi) q[1];
rz(-2.4843487) q[2];
sx q[2];
rz(-2.1604156) q[2];
sx q[2];
rz(0.60962617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6440297) q[1];
sx q[1];
rz(-2.7179238) q[1];
sx q[1];
rz(-0.81987857) q[1];
rz(-pi) q[2];
rz(1.0172306) q[3];
sx q[3];
rz(-0.38700208) q[3];
sx q[3];
rz(-1.6978881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2018373) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(-1.6377829) q[2];
rz(2.2683897) q[3];
sx q[3];
rz(-0.88770294) q[3];
sx q[3];
rz(-2.8100815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37209508) q[0];
sx q[0];
rz(-2.7517419) q[0];
sx q[0];
rz(1.9918342) q[0];
rz(2.276543) q[1];
sx q[1];
rz(-1.1999835) q[1];
sx q[1];
rz(1.4335776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652611) q[0];
sx q[0];
rz(-1.0076242) q[0];
sx q[0];
rz(2.7050326) q[0];
rz(-0.63460372) q[2];
sx q[2];
rz(-0.57159492) q[2];
sx q[2];
rz(1.7224479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99231863) q[1];
sx q[1];
rz(-2.092835) q[1];
sx q[1];
rz(3.0401413) q[1];
rz(-0.20345466) q[3];
sx q[3];
rz(-1.9793812) q[3];
sx q[3];
rz(2.3060498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4992497) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(-0.7439822) q[2];
rz(2.3104525) q[3];
sx q[3];
rz(-1.7404514) q[3];
sx q[3];
rz(-2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0396073) q[0];
sx q[0];
rz(-0.87431854) q[0];
sx q[0];
rz(0.71514493) q[0];
rz(1.2932418) q[1];
sx q[1];
rz(-1.5903541) q[1];
sx q[1];
rz(0.23983461) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3823925) q[0];
sx q[0];
rz(-2.0682997) q[0];
sx q[0];
rz(-0.63775179) q[0];
x q[1];
rz(-1.5101899) q[2];
sx q[2];
rz(-2.1327634) q[2];
sx q[2];
rz(2.2417703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4308987) q[1];
sx q[1];
rz(-1.7393117) q[1];
sx q[1];
rz(0.88697834) q[1];
rz(1.6039944) q[3];
sx q[3];
rz(-2.0616459) q[3];
sx q[3];
rz(0.50118485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0750601) q[2];
sx q[2];
rz(-1.5873453) q[2];
sx q[2];
rz(1.4934322) q[2];
rz(3.1028124) q[3];
sx q[3];
rz(-1.6438234) q[3];
sx q[3];
rz(-1.405976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99603242) q[0];
sx q[0];
rz(-1.4936438) q[0];
sx q[0];
rz(-2.5160568) q[0];
rz(1.372555) q[1];
sx q[1];
rz(-2.1452429) q[1];
sx q[1];
rz(-2.5575584) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9720219) q[0];
sx q[0];
rz(-1.2002647) q[0];
sx q[0];
rz(1.2967923) q[0];
x q[1];
rz(0.66472424) q[2];
sx q[2];
rz(-1.4493086) q[2];
sx q[2];
rz(0.45035502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0528193) q[1];
sx q[1];
rz(-0.86630834) q[1];
sx q[1];
rz(1.2831844) q[1];
x q[2];
rz(0.68525641) q[3];
sx q[3];
rz(-1.7820089) q[3];
sx q[3];
rz(1.2880471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0293616) q[2];
sx q[2];
rz(-1.4697305) q[2];
sx q[2];
rz(-1.3545797) q[2];
rz(-2.3692047) q[3];
sx q[3];
rz(-1.200501) q[3];
sx q[3];
rz(1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1188485) q[0];
sx q[0];
rz(-2.4130555) q[0];
sx q[0];
rz(1.5424315) q[0];
rz(-2.6842885) q[1];
sx q[1];
rz(-0.87304579) q[1];
sx q[1];
rz(2.9541448) q[1];
rz(-3.0627261) q[2];
sx q[2];
rz(-1.3838242) q[2];
sx q[2];
rz(0.48505062) q[2];
rz(-2.0304092) q[3];
sx q[3];
rz(-1.2818973) q[3];
sx q[3];
rz(-3.1218798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
