OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(3.0732529) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(-3.1037722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28107444) q[0];
sx q[0];
rz(-1.3776508) q[0];
sx q[0];
rz(-3.0618969) q[0];
rz(-pi) q[1];
rz(-1.7224563) q[2];
sx q[2];
rz(-0.66705396) q[2];
sx q[2];
rz(2.5487713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0953857) q[1];
sx q[1];
rz(-2.3933105) q[1];
sx q[1];
rz(-0.074949646) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9904198) q[3];
sx q[3];
rz(-1.9133948) q[3];
sx q[3];
rz(-1.7455938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35090703) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(1.2834056) q[2];
rz(3.0154199) q[3];
sx q[3];
rz(-1.7254555) q[3];
sx q[3];
rz(-3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693817) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(2.5449261) q[0];
rz(1.5860575) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(1.7780875) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1648646) q[0];
sx q[0];
rz(-0.70978998) q[0];
sx q[0];
rz(-0.68517942) q[0];
rz(-3.0599669) q[2];
sx q[2];
rz(-2.584169) q[2];
sx q[2];
rz(-2.4068085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1494257) q[1];
sx q[1];
rz(-2.4352695) q[1];
sx q[1];
rz(-2.7041433) q[1];
x q[2];
rz(2.2927106) q[3];
sx q[3];
rz(-0.97182453) q[3];
sx q[3];
rz(2.2357383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3229225) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(2.4070516) q[2];
rz(2.5143886) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(0.27221361) q[0];
rz(-2.294337) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(2.8289657) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6486559) q[0];
sx q[0];
rz(-1.3515633) q[0];
sx q[0];
rz(1.5159025) q[0];
rz(-2.4467144) q[2];
sx q[2];
rz(-0.67145295) q[2];
sx q[2];
rz(1.7365255) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9885013) q[1];
sx q[1];
rz(-1.830849) q[1];
sx q[1];
rz(-0.071916332) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39194312) q[3];
sx q[3];
rz(-2.3392945) q[3];
sx q[3];
rz(2.9585569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.12038885) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(-2.1728544) q[3];
sx q[3];
rz(-2.0961943) q[3];
sx q[3];
rz(2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3880436) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(0.62336212) q[0];
rz(-0.81758824) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(3.0923016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.277963) q[0];
sx q[0];
rz(-0.81252126) q[0];
sx q[0];
rz(-0.066594007) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2055779) q[2];
sx q[2];
rz(-0.74539241) q[2];
sx q[2];
rz(0.12860194) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.149718) q[1];
sx q[1];
rz(-1.8135934) q[1];
sx q[1];
rz(-1.9096096) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0443346) q[3];
sx q[3];
rz(-1.9658078) q[3];
sx q[3];
rz(-0.21316866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-2.1566186) q[2];
rz(-0.90304053) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(-0.27339098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7966998) q[0];
sx q[0];
rz(-1.6051689) q[0];
sx q[0];
rz(0.087619089) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-0.55115288) q[1];
sx q[1];
rz(-0.87096754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.464251) q[0];
sx q[0];
rz(-2.8528677) q[0];
sx q[0];
rz(-3.0407453) q[0];
x q[1];
rz(1.9647314) q[2];
sx q[2];
rz(-1.7125868) q[2];
sx q[2];
rz(-2.1243492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37967967) q[1];
sx q[1];
rz(-1.2907791) q[1];
sx q[1];
rz(-0.67358394) q[1];
x q[2];
rz(-2.3351339) q[3];
sx q[3];
rz(-1.4440923) q[3];
sx q[3];
rz(-0.63092953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0481723) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(2.7887153) q[2];
rz(0.86587632) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.6784994) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(3.034806) q[0];
rz(1.9550025) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(2.1170763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2811919) q[0];
sx q[0];
rz(-1.8587347) q[0];
sx q[0];
rz(-0.47173758) q[0];
rz(0.60501955) q[2];
sx q[2];
rz(-1.9540678) q[2];
sx q[2];
rz(-0.59125102) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5974821) q[1];
sx q[1];
rz(-1.0162309) q[1];
sx q[1];
rz(2.5764562) q[1];
rz(-2.5173353) q[3];
sx q[3];
rz(-2.5578824) q[3];
sx q[3];
rz(0.39238413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.9275894) q[2];
sx q[2];
rz(-2.270703) q[2];
rz(1.332256) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-0.55091888) q[0];
rz(-2.6761966) q[1];
sx q[1];
rz(-0.67133633) q[1];
sx q[1];
rz(2.904772) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4846372) q[0];
sx q[0];
rz(-1.4716363) q[0];
sx q[0];
rz(-3.1037472) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4004164) q[2];
sx q[2];
rz(-2.9393662) q[2];
sx q[2];
rz(1.9066332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4983066) q[1];
sx q[1];
rz(-0.28043881) q[1];
sx q[1];
rz(-1.2966869) q[1];
x q[2];
rz(1.9067494) q[3];
sx q[3];
rz(-2.3785605) q[3];
sx q[3];
rz(1.2515765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44234309) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(-0.44011763) q[2];
rz(-1.0951428) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(-1.346689) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(-3.112088) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(-3.0227919) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122075) q[0];
sx q[0];
rz(-0.31353024) q[0];
sx q[0];
rz(-3.0010812) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88840862) q[2];
sx q[2];
rz(-2.4783122) q[2];
sx q[2];
rz(1.0379438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6606632) q[1];
sx q[1];
rz(-0.56023635) q[1];
sx q[1];
rz(-0.57628298) q[1];
rz(-pi) q[2];
rz(1.8700065) q[3];
sx q[3];
rz(-2.6282103) q[3];
sx q[3];
rz(-1.3083003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3999195) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(-1.5734394) q[2];
rz(-1.167477) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(-1.2742111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2385999) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(2.9883244) q[0];
rz(-1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.1538039) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62350863) q[0];
sx q[0];
rz(-2.6797446) q[0];
sx q[0];
rz(0.64103809) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93038721) q[2];
sx q[2];
rz(-0.79158917) q[2];
sx q[2];
rz(-1.1861578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24473083) q[1];
sx q[1];
rz(-1.3618999) q[1];
sx q[1];
rz(-0.056785866) q[1];
rz(-1.7180175) q[3];
sx q[3];
rz(-0.69268337) q[3];
sx q[3];
rz(-2.6976531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8447421) q[2];
sx q[2];
rz(-1.1721609) q[2];
sx q[2];
rz(-1.6142169) q[2];
rz(-0.99689233) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(2.2629288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8940354) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(2.9515008) q[0];
rz(2.4709573) q[1];
sx q[1];
rz(-1.1963444) q[1];
sx q[1];
rz(-1.488283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9207536) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(-0.12551813) q[0];
rz(-pi) q[1];
rz(-2.6121717) q[2];
sx q[2];
rz(-0.501907) q[2];
sx q[2];
rz(1.5159964) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6675122) q[1];
sx q[1];
rz(-2.0733359) q[1];
sx q[1];
rz(-0.71943347) q[1];
rz(-pi) q[2];
rz(-2.2404352) q[3];
sx q[3];
rz(-2.3271051) q[3];
sx q[3];
rz(-2.8709473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1378479) q[2];
sx q[2];
rz(-2.2403084) q[2];
sx q[2];
rz(-0.89912644) q[2];
rz(-2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0638194) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(-0.28941119) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(-2.4317447) q[2];
sx q[2];
rz(-1.3550497) q[2];
sx q[2];
rz(1.4646127) q[2];
rz(1.9990986) q[3];
sx q[3];
rz(-2.8430568) q[3];
sx q[3];
rz(-1.6381016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
