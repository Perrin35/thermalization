OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.53192294) q[0];
sx q[0];
rz(1.6165531) q[0];
sx q[0];
rz(11.308148) q[0];
rz(-1.582107) q[1];
sx q[1];
rz(-2.6546302) q[1];
sx q[1];
rz(-0.43509405) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6995876) q[0];
sx q[0];
rz(-2.3908092) q[0];
sx q[0];
rz(-2.4155099) q[0];
rz(-1.2266394) q[2];
sx q[2];
rz(-2.0822869) q[2];
sx q[2];
rz(0.79424266) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3341325) q[1];
sx q[1];
rz(-1.301348) q[1];
sx q[1];
rz(-0.079748591) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55032309) q[3];
sx q[3];
rz(-1.1512412) q[3];
sx q[3];
rz(1.5997831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2121928) q[2];
sx q[2];
rz(-1.2779028) q[2];
sx q[2];
rz(2.0983992) q[2];
rz(-2.7405401) q[3];
sx q[3];
rz(-1.6845104) q[3];
sx q[3];
rz(2.5636087) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9143518) q[0];
sx q[0];
rz(-0.14509097) q[0];
sx q[0];
rz(2.5703854) q[0];
rz(-0.97194833) q[1];
sx q[1];
rz(-0.74058878) q[1];
sx q[1];
rz(2.0170225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034649523) q[0];
sx q[0];
rz(-0.9616344) q[0];
sx q[0];
rz(0.27584313) q[0];
rz(-pi) q[1];
rz(-0.45887453) q[2];
sx q[2];
rz(-0.70415184) q[2];
sx q[2];
rz(1.6353232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8001061) q[1];
sx q[1];
rz(-0.84034318) q[1];
sx q[1];
rz(0.12612242) q[1];
rz(-0.092512802) q[3];
sx q[3];
rz(-2.9399151) q[3];
sx q[3];
rz(2.8254333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9198415) q[2];
sx q[2];
rz(-1.5849179) q[2];
sx q[2];
rz(-2.705503) q[2];
rz(2.3927205) q[3];
sx q[3];
rz(-0.80513969) q[3];
sx q[3];
rz(2.3124636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0332758) q[0];
sx q[0];
rz(-1.746614) q[0];
sx q[0];
rz(3.0535611) q[0];
rz(0.85405675) q[1];
sx q[1];
rz(-2.4805534) q[1];
sx q[1];
rz(1.0850272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9097901) q[0];
sx q[0];
rz(-0.2806305) q[0];
sx q[0];
rz(-1.8156194) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7872115) q[2];
sx q[2];
rz(-1.4206352) q[2];
sx q[2];
rz(0.79126287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8768118) q[1];
sx q[1];
rz(-1.1327883) q[1];
sx q[1];
rz(-2.9609738) q[1];
rz(-pi) q[2];
rz(2.7522467) q[3];
sx q[3];
rz(-1.5028477) q[3];
sx q[3];
rz(0.24711497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4714841) q[2];
sx q[2];
rz(-1.7210311) q[2];
sx q[2];
rz(2.3954083) q[2];
rz(-0.21607312) q[3];
sx q[3];
rz(-1.8862628) q[3];
sx q[3];
rz(-0.20572534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3716607) q[0];
sx q[0];
rz(-0.07337229) q[0];
sx q[0];
rz(-1.368847) q[0];
rz(2.5313077) q[1];
sx q[1];
rz(-1.3657602) q[1];
sx q[1];
rz(-2.761421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1606522) q[0];
sx q[0];
rz(-1.9987172) q[0];
sx q[0];
rz(-1.227427) q[0];
rz(0.45079239) q[2];
sx q[2];
rz(-1.7480228) q[2];
sx q[2];
rz(1.4361567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3750864) q[1];
sx q[1];
rz(-1.0712873) q[1];
sx q[1];
rz(0.16721027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.065297619) q[3];
sx q[3];
rz(-1.9717798) q[3];
sx q[3];
rz(-2.9915265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2569106) q[2];
sx q[2];
rz(-2.1941049) q[2];
sx q[2];
rz(1.0379637) q[2];
rz(-0.3642309) q[3];
sx q[3];
rz(-2.3528152) q[3];
sx q[3];
rz(-2.4558333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.6895741) q[0];
sx q[0];
rz(-1.7114534) q[0];
sx q[0];
rz(-0.83056393) q[0];
rz(-3.0914302) q[1];
sx q[1];
rz(-0.12718931) q[1];
sx q[1];
rz(0.62279472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2723214) q[0];
sx q[0];
rz(-1.4928482) q[0];
sx q[0];
rz(2.9236631) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3741983) q[2];
sx q[2];
rz(-1.9834347) q[2];
sx q[2];
rz(-2.430344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.078256814) q[1];
sx q[1];
rz(-1.2163269) q[1];
sx q[1];
rz(0.88651231) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5707309) q[3];
sx q[3];
rz(-1.5483466) q[3];
sx q[3];
rz(-2.8714436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.043776) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(-2.3338976) q[2];
rz(1.6037327) q[3];
sx q[3];
rz(-1.5065498) q[3];
sx q[3];
rz(2.9174771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.9648801) q[0];
sx q[0];
rz(-1.3310615) q[0];
sx q[0];
rz(-1.4759195) q[0];
rz(-1.9077979) q[1];
sx q[1];
rz(-1.4101615) q[1];
sx q[1];
rz(-2.9880611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79685211) q[0];
sx q[0];
rz(-1.5670876) q[0];
sx q[0];
rz(3.0319503) q[0];
rz(-pi) q[1];
rz(2.6903908) q[2];
sx q[2];
rz(-1.1179194) q[2];
sx q[2];
rz(1.7141327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0316089) q[1];
sx q[1];
rz(-0.69797198) q[1];
sx q[1];
rz(2.4137761) q[1];
rz(-pi) q[2];
rz(0.12795398) q[3];
sx q[3];
rz(-2.4989486) q[3];
sx q[3];
rz(1.2817739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70225707) q[2];
sx q[2];
rz(-2.3372237) q[2];
sx q[2];
rz(0.55629936) q[2];
rz(0.016911658) q[3];
sx q[3];
rz(-0.97512475) q[3];
sx q[3];
rz(2.4771966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21338129) q[0];
sx q[0];
rz(-3.0688372) q[0];
sx q[0];
rz(2.4640006) q[0];
rz(1.3202336) q[1];
sx q[1];
rz(-2.0537387) q[1];
sx q[1];
rz(2.5755612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460577) q[0];
sx q[0];
rz(-1.9060152) q[0];
sx q[0];
rz(0.4636824) q[0];
rz(-pi) q[1];
rz(-1.2957057) q[2];
sx q[2];
rz(-0.054271532) q[2];
sx q[2];
rz(1.4473297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.54864751) q[1];
sx q[1];
rz(-0.8742395) q[1];
sx q[1];
rz(-0.19754438) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0637897) q[3];
sx q[3];
rz(-0.49178472) q[3];
sx q[3];
rz(3.0932757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2322106) q[2];
sx q[2];
rz(-0.38572329) q[2];
sx q[2];
rz(0.56987008) q[2];
rz(-2.0991523) q[3];
sx q[3];
rz(-1.9614377) q[3];
sx q[3];
rz(1.0038143) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7637699) q[0];
sx q[0];
rz(-2.7144987) q[0];
sx q[0];
rz(-2.1183993) q[0];
rz(2.2238253) q[1];
sx q[1];
rz(-2.387391) q[1];
sx q[1];
rz(2.2775211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26041481) q[0];
sx q[0];
rz(-0.27363187) q[0];
sx q[0];
rz(1.1009965) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0351428) q[2];
sx q[2];
rz(-1.2304145) q[2];
sx q[2];
rz(2.2156773) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.041124972) q[1];
sx q[1];
rz(-2.0534614) q[1];
sx q[1];
rz(-0.80759279) q[1];
rz(2.4141563) q[3];
sx q[3];
rz(-1.2890649) q[3];
sx q[3];
rz(1.7121332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0920022) q[2];
sx q[2];
rz(-1.8747753) q[2];
sx q[2];
rz(-0.76466307) q[2];
rz(0.90096724) q[3];
sx q[3];
rz(-0.67470208) q[3];
sx q[3];
rz(2.0425792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6242999) q[0];
sx q[0];
rz(-0.39532548) q[0];
sx q[0];
rz(0.98315352) q[0];
rz(0.82955018) q[1];
sx q[1];
rz(-1.7770276) q[1];
sx q[1];
rz(-2.5628536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5924517) q[0];
sx q[0];
rz(-1.5835174) q[0];
sx q[0];
rz(1.1770269) q[0];
rz(0.54429599) q[2];
sx q[2];
rz(-1.0981264) q[2];
sx q[2];
rz(-0.67677697) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86218802) q[1];
sx q[1];
rz(-1.1694238) q[1];
sx q[1];
rz(2.1251232) q[1];
rz(-pi) q[2];
x q[2];
rz(0.345295) q[3];
sx q[3];
rz(-1.5362979) q[3];
sx q[3];
rz(-1.9138699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4911554) q[2];
sx q[2];
rz(-1.2885685) q[2];
sx q[2];
rz(-0.11162652) q[2];
rz(-0.95587436) q[3];
sx q[3];
rz(-0.99146944) q[3];
sx q[3];
rz(-1.2607695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2535506) q[0];
sx q[0];
rz(-2.084806) q[0];
sx q[0];
rz(-0.014130935) q[0];
rz(-2.0290532) q[1];
sx q[1];
rz(-2.2281149) q[1];
sx q[1];
rz(-1.5412451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80373549) q[0];
sx q[0];
rz(-1.4232262) q[0];
sx q[0];
rz(1.7002559) q[0];
x q[1];
rz(-0.92629536) q[2];
sx q[2];
rz(-1.6716812) q[2];
sx q[2];
rz(0.48357329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75493357) q[1];
sx q[1];
rz(-1.5112097) q[1];
sx q[1];
rz(-0.3646551) q[1];
x q[2];
rz(-1.4061116) q[3];
sx q[3];
rz(-1.7385118) q[3];
sx q[3];
rz(-2.8693143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9666226) q[2];
sx q[2];
rz(-1.9113767) q[2];
sx q[2];
rz(-2.5959385) q[2];
rz(-0.086006554) q[3];
sx q[3];
rz(-1.514148) q[3];
sx q[3];
rz(0.29451323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.1957112) q[0];
sx q[0];
rz(-1.042689) q[0];
sx q[0];
rz(-1.9324017) q[0];
rz(0.60278268) q[1];
sx q[1];
rz(-2.5310015) q[1];
sx q[1];
rz(-2.3878154) q[1];
rz(2.2529765) q[2];
sx q[2];
rz(-1.7352413) q[2];
sx q[2];
rz(-2.2189463) q[2];
rz(-0.15074853) q[3];
sx q[3];
rz(-2.5133532) q[3];
sx q[3];
rz(-2.8379074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
