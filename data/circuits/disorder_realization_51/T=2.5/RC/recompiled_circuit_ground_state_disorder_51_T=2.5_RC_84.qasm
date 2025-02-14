OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2188323) q[0];
sx q[0];
rz(-0.063475944) q[0];
sx q[0];
rz(-1.3702962) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(-1.4501362) q[1];
sx q[1];
rz(-2.9210747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0974183) q[0];
sx q[0];
rz(-0.83651354) q[0];
sx q[0];
rz(0.091139779) q[0];
rz(-pi) q[1];
rz(2.0283255) q[2];
sx q[2];
rz(-1.2731247) q[2];
sx q[2];
rz(0.13099094) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.045956417) q[1];
sx q[1];
rz(-1.7929309) q[1];
sx q[1];
rz(1.5820811) q[1];
rz(-pi) q[2];
rz(1.9265367) q[3];
sx q[3];
rz(-1.3827207) q[3];
sx q[3];
rz(-2.637685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0567131) q[2];
sx q[2];
rz(-0.006338174) q[2];
sx q[2];
rz(-2.32178) q[2];
rz(0.020717185) q[3];
sx q[3];
rz(-0.31532225) q[3];
sx q[3];
rz(1.0516385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3725975) q[0];
sx q[0];
rz(-2.9561434) q[0];
sx q[0];
rz(0.73082596) q[0];
rz(-1.7164879) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(2.6740429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1913203) q[0];
sx q[0];
rz(-1.6015656) q[0];
sx q[0];
rz(-1.6732209) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.076084332) q[2];
sx q[2];
rz(-0.4108763) q[2];
sx q[2];
rz(2.8166688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4380355) q[1];
sx q[1];
rz(-2.272637) q[1];
sx q[1];
rz(1.3753533) q[1];
x q[2];
rz(-2.8465956) q[3];
sx q[3];
rz(-2.2824691) q[3];
sx q[3];
rz(1.0906578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1656437) q[2];
sx q[2];
rz(-3.1303945) q[2];
sx q[2];
rz(0.33640081) q[2];
rz(-2.8339556) q[3];
sx q[3];
rz(-3.1308789) q[3];
sx q[3];
rz(0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0129608) q[0];
sx q[0];
rz(-2.9427981) q[0];
sx q[0];
rz(0.13191731) q[0];
rz(1.4384653) q[1];
sx q[1];
rz(-1.7273644) q[1];
sx q[1];
rz(-1.6579113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8298873) q[0];
sx q[0];
rz(-2.1714806) q[0];
sx q[0];
rz(-1.2313103) q[0];
rz(-pi) q[1];
rz(2.5222579) q[2];
sx q[2];
rz(-0.029994596) q[2];
sx q[2];
rz(0.66823792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4705088) q[1];
sx q[1];
rz(-1.5189063) q[1];
sx q[1];
rz(1.6747649) q[1];
rz(-3.0555757) q[3];
sx q[3];
rz(-1.5002076) q[3];
sx q[3];
rz(0.020078192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5710473) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(1.4578311) q[2];
rz(-2.3633862) q[3];
sx q[3];
rz(-3.1359378) q[3];
sx q[3];
rz(0.2125423) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197765) q[0];
sx q[0];
rz(-1.3874522) q[0];
sx q[0];
rz(2.8507267) q[0];
rz(-1.5485171) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(-0.026550857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7598068) q[0];
sx q[0];
rz(-1.9142173) q[0];
sx q[0];
rz(-1.9476101) q[0];
x q[1];
rz(1.5174116) q[2];
sx q[2];
rz(-1.1204655) q[2];
sx q[2];
rz(-2.2105697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77912731) q[1];
sx q[1];
rz(-0.82632768) q[1];
sx q[1];
rz(0.67693784) q[1];
rz(-pi) q[2];
rz(-1.4064571) q[3];
sx q[3];
rz(-0.019329179) q[3];
sx q[3];
rz(-0.011474495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2619541) q[2];
sx q[2];
rz(-0.018435437) q[2];
sx q[2];
rz(2.7035942) q[2];
rz(-2.0336464) q[3];
sx q[3];
rz(-3.1236533) q[3];
sx q[3];
rz(-1.7855135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172025) q[0];
sx q[0];
rz(-2.6438535) q[0];
sx q[0];
rz(-2.9955731) q[0];
rz(3.0085425) q[1];
sx q[1];
rz(-1.0597884) q[1];
sx q[1];
rz(0.45566794) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9216292) q[0];
sx q[0];
rz(-1.5908003) q[0];
sx q[0];
rz(3.0476493) q[0];
rz(-pi) q[1];
rz(-0.44134042) q[2];
sx q[2];
rz(-2.5020863) q[2];
sx q[2];
rz(-2.9616982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2777777) q[1];
sx q[1];
rz(-1.3425273) q[1];
sx q[1];
rz(-0.14240188) q[1];
x q[2];
rz(2.9140329) q[3];
sx q[3];
rz(-1.5858082) q[3];
sx q[3];
rz(-0.14666569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2666152) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.2561426) q[2];
rz(-0.52136326) q[3];
sx q[3];
rz(-2.8837236) q[3];
sx q[3];
rz(2.7969587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94619036) q[0];
sx q[0];
rz(-0.38637105) q[0];
sx q[0];
rz(2.572701) q[0];
rz(2.9733114) q[1];
sx q[1];
rz(-1.5852837) q[1];
sx q[1];
rz(0.010244244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3061188) q[0];
sx q[0];
rz(-2.9635307) q[0];
sx q[0];
rz(2.8639069) q[0];
rz(-1.7513137) q[2];
sx q[2];
rz(-1.5517934) q[2];
sx q[2];
rz(0.190616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.467874) q[1];
sx q[1];
rz(-1.5281902) q[1];
sx q[1];
rz(-0.03454476) q[1];
rz(-0.041268392) q[3];
sx q[3];
rz(-2.5070243) q[3];
sx q[3];
rz(2.7260425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16113082) q[2];
sx q[2];
rz(-0.22393301) q[2];
sx q[2];
rz(1.0978318) q[2];
rz(2.791642) q[3];
sx q[3];
rz(-1.2772468) q[3];
sx q[3];
rz(2.2866975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.8769281) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(-0.40633416) q[0];
rz(1.8304652) q[1];
sx q[1];
rz(-2.6268112) q[1];
sx q[1];
rz(-0.11022551) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.742934) q[0];
sx q[0];
rz(-0.1029108) q[0];
sx q[0];
rz(-2.8714931) q[0];
rz(0.95944698) q[2];
sx q[2];
rz(-0.0065718647) q[2];
sx q[2];
rz(-0.98192865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.266763) q[1];
sx q[1];
rz(-1.5795261) q[1];
sx q[1];
rz(1.1432024) q[1];
x q[2];
rz(0.74177058) q[3];
sx q[3];
rz(-2.2549155) q[3];
sx q[3];
rz(-2.6009862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.280764) q[2];
sx q[2];
rz(-3.132756) q[2];
sx q[2];
rz(2.3542118) q[2];
rz(-2.8619838) q[3];
sx q[3];
rz(-3.0968554) q[3];
sx q[3];
rz(0.71981049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0019919458) q[0];
sx q[0];
rz(-2.3878492) q[0];
sx q[0];
rz(-1.4308223) q[0];
rz(2.9665973) q[1];
sx q[1];
rz(-0.7376968) q[1];
sx q[1];
rz(-1.7134679) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7044107) q[0];
sx q[0];
rz(-1.5488273) q[0];
sx q[0];
rz(0.0030513838) q[0];
rz(-pi) q[1];
rz(-1.0662717) q[2];
sx q[2];
rz(-3.1331535) q[2];
sx q[2];
rz(-1.5167459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0261532) q[1];
sx q[1];
rz(-2.8255531) q[1];
sx q[1];
rz(1.0277961) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.346502) q[3];
sx q[3];
rz(-1.845775) q[3];
sx q[3];
rz(-2.0338273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3565107) q[2];
sx q[2];
rz(-0.76866895) q[2];
sx q[2];
rz(-1.53995) q[2];
rz(-0.29940638) q[3];
sx q[3];
rz(-1.6573903) q[3];
sx q[3];
rz(-0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141075) q[0];
sx q[0];
rz(-3.1306559) q[0];
sx q[0];
rz(2.6543044) q[0];
rz(1.4393073) q[1];
sx q[1];
rz(-0.7936365) q[1];
sx q[1];
rz(-2.7557441) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0087465) q[0];
sx q[0];
rz(-1.9632959) q[0];
sx q[0];
rz(-1.5061492) q[0];
rz(-3.0294777) q[2];
sx q[2];
rz(-1.1476074) q[2];
sx q[2];
rz(1.8433169) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1362366) q[1];
sx q[1];
rz(-0.87974977) q[1];
sx q[1];
rz(-1.5278396) q[1];
rz(-1.5351546) q[3];
sx q[3];
rz(-1.9058133) q[3];
sx q[3];
rz(0.50594508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7172598) q[2];
sx q[2];
rz(-0.029742664) q[2];
sx q[2];
rz(2.7359803) q[2];
rz(0.82617104) q[3];
sx q[3];
rz(-0.043353733) q[3];
sx q[3];
rz(2.1521173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1694323) q[0];
sx q[0];
rz(-0.5652453) q[0];
sx q[0];
rz(1.7461079) q[0];
rz(1.6204429) q[1];
sx q[1];
rz(-2.4902159) q[1];
sx q[1];
rz(-1.7924538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682705) q[0];
sx q[0];
rz(-1.3600469) q[0];
sx q[0];
rz(0.83619946) q[0];
x q[1];
rz(-0.17041309) q[2];
sx q[2];
rz(-1.8080343) q[2];
sx q[2];
rz(1.3875675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1101035) q[1];
sx q[1];
rz(-1.5851136) q[1];
sx q[1];
rz(-1.6308466) q[1];
x q[2];
rz(3.1114499) q[3];
sx q[3];
rz(-1.5858012) q[3];
sx q[3];
rz(-0.74149132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32538432) q[2];
sx q[2];
rz(-0.00486972) q[2];
sx q[2];
rz(-1.7083141) q[2];
rz(-1.8520744) q[3];
sx q[3];
rz(-3.0472445) q[3];
sx q[3];
rz(-1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0636487) q[0];
sx q[0];
rz(-1.0366806) q[0];
sx q[0];
rz(0.59564577) q[0];
rz(3.0972277) q[1];
sx q[1];
rz(-2.8652419) q[1];
sx q[1];
rz(-1.2038632) q[1];
rz(1.1677713) q[2];
sx q[2];
rz(-2.6296305) q[2];
sx q[2];
rz(0.59811022) q[2];
rz(-1.1607004) q[3];
sx q[3];
rz(-1.5103419) q[3];
sx q[3];
rz(1.7132669) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
