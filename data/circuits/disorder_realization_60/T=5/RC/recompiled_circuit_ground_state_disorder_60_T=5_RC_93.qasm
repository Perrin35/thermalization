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
rz(-0.54028264) q[1];
sx q[1];
rz(-1.5958495) q[1];
sx q[1];
rz(-0.63313142) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608434) q[0];
sx q[0];
rz(-2.9414163) q[0];
sx q[0];
rz(0.60524888) q[0];
x q[1];
rz(-1.6923793) q[2];
sx q[2];
rz(-2.0544397) q[2];
sx q[2];
rz(-1.759296) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3772845) q[1];
sx q[1];
rz(-1.8486684) q[1];
sx q[1];
rz(-1.4382526) q[1];
rz(-pi) q[2];
rz(-2.5217275) q[3];
sx q[3];
rz(-0.64398089) q[3];
sx q[3];
rz(-2.5247273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4151609) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(-0.64217448) q[3];
sx q[3];
rz(-0.29688501) q[3];
sx q[3];
rz(1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.17382962) q[0];
sx q[0];
rz(-2.1673295) q[0];
sx q[0];
rz(-2.1948788) q[0];
rz(0.058454839) q[1];
sx q[1];
rz(-0.36403257) q[1];
sx q[1];
rz(2.2296925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4218876) q[0];
sx q[0];
rz(-1.9870409) q[0];
sx q[0];
rz(3.1001607) q[0];
x q[1];
rz(-1.4824268) q[2];
sx q[2];
rz(-2.4659133) q[2];
sx q[2];
rz(2.8203204) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8337149) q[1];
sx q[1];
rz(-2.6243997) q[1];
sx q[1];
rz(2.6343995) q[1];
rz(1.5188317) q[3];
sx q[3];
rz(-1.1295737) q[3];
sx q[3];
rz(0.085746229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6902206) q[2];
sx q[2];
rz(-2.8896832) q[2];
sx q[2];
rz(-2.8042262) q[2];
rz(2.1515324) q[3];
sx q[3];
rz(-1.3288574) q[3];
sx q[3];
rz(-1.2216074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5628691) q[0];
sx q[0];
rz(-2.8674485) q[0];
sx q[0];
rz(1.3124056) q[0];
rz(2.9263272) q[1];
sx q[1];
rz(-1.5039517) q[1];
sx q[1];
rz(2.2832835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68033531) q[0];
sx q[0];
rz(-1.3051998) q[0];
sx q[0];
rz(3.0895825) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1444001) q[2];
sx q[2];
rz(-2.7421087) q[2];
sx q[2];
rz(-1.4956719) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.55187622) q[1];
sx q[1];
rz(-1.405542) q[1];
sx q[1];
rz(2.698353) q[1];
rz(-pi) q[2];
rz(-2.2792072) q[3];
sx q[3];
rz(-0.38194712) q[3];
sx q[3];
rz(2.2464584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12576088) q[2];
sx q[2];
rz(-0.55700004) q[2];
sx q[2];
rz(-2.1324943) q[2];
rz(0.5674181) q[3];
sx q[3];
rz(-1.6332473) q[3];
sx q[3];
rz(1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0030982) q[0];
sx q[0];
rz(-2.4479285) q[0];
sx q[0];
rz(0.19288119) q[0];
rz(-1.0640954) q[1];
sx q[1];
rz(-0.27609438) q[1];
sx q[1];
rz(2.1853133) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67502695) q[0];
sx q[0];
rz(-2.1888632) q[0];
sx q[0];
rz(1.3542372) q[0];
rz(-0.85727875) q[2];
sx q[2];
rz(-0.801221) q[2];
sx q[2];
rz(-0.73357936) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6234474) q[1];
sx q[1];
rz(-1.8817187) q[1];
sx q[1];
rz(0.53849935) q[1];
x q[2];
rz(0.90518454) q[3];
sx q[3];
rz(-1.815261) q[3];
sx q[3];
rz(2.8595231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3095779) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(0.94069707) q[2];
rz(-2.7637123) q[3];
sx q[3];
rz(-2.2418719) q[3];
sx q[3];
rz(-2.5009724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2864895) q[0];
sx q[0];
rz(-1.531895) q[0];
sx q[0];
rz(1.8163649) q[0];
rz(2.8638966) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(1.206254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0632478) q[0];
sx q[0];
rz(-1.6499015) q[0];
sx q[0];
rz(1.3782756) q[0];
rz(-2.5857193) q[2];
sx q[2];
rz(-0.68965675) q[2];
sx q[2];
rz(1.7574312) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6305651) q[1];
sx q[1];
rz(-1.2113809) q[1];
sx q[1];
rz(0.014822249) q[1];
rz(-0.60449497) q[3];
sx q[3];
rz(-2.8504728) q[3];
sx q[3];
rz(-1.1230942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18579927) q[2];
sx q[2];
rz(-2.9904371) q[2];
sx q[2];
rz(-1.7929057) q[2];
rz(-2.1336613) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0808829) q[0];
sx q[0];
rz(-0.01132975) q[0];
sx q[0];
rz(-2.9445373) q[0];
rz(-1.9673678) q[1];
sx q[1];
rz(-2.2815727) q[1];
sx q[1];
rz(-2.7161782) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9369318) q[0];
sx q[0];
rz(-0.46260329) q[0];
sx q[0];
rz(-3.0499056) q[0];
x q[1];
rz(-1.2666918) q[2];
sx q[2];
rz(-1.5291404) q[2];
sx q[2];
rz(1.1343985) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1716291) q[1];
sx q[1];
rz(-2.6416314) q[1];
sx q[1];
rz(0.5699711) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69611551) q[3];
sx q[3];
rz(-0.98721877) q[3];
sx q[3];
rz(2.7903874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5328488) q[2];
sx q[2];
rz(-3.0259735) q[2];
sx q[2];
rz(0.70333424) q[2];
rz(-2.9902847) q[3];
sx q[3];
rz(-1.1727099) q[3];
sx q[3];
rz(-2.5942588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587826) q[0];
sx q[0];
rz(-2.3334239) q[0];
sx q[0];
rz(2.1852344) q[0];
rz(-0.76902485) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(1.0395315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2295264) q[0];
sx q[0];
rz(-1.5472224) q[0];
sx q[0];
rz(2.0458771) q[0];
rz(0.86918594) q[2];
sx q[2];
rz(-2.1034046) q[2];
sx q[2];
rz(2.5857596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36827792) q[1];
sx q[1];
rz(-1.2864783) q[1];
sx q[1];
rz(-1.8892688) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9305228) q[3];
sx q[3];
rz(-1.8976334) q[3];
sx q[3];
rz(-2.0321996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9397554) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(1.5038097) q[2];
rz(-2.2683897) q[3];
sx q[3];
rz(-0.88770294) q[3];
sx q[3];
rz(2.8100815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7694976) q[0];
sx q[0];
rz(-0.3898507) q[0];
sx q[0];
rz(1.1497585) q[0];
rz(-2.276543) q[1];
sx q[1];
rz(-1.1999835) q[1];
sx q[1];
rz(1.7080151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85808774) q[0];
sx q[0];
rz(-2.4437987) q[0];
sx q[0];
rz(-2.1608864) q[0];
rz(-pi) q[1];
rz(1.2064798) q[2];
sx q[2];
rz(-1.120036) q[2];
sx q[2];
rz(-2.4414666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5123957) q[1];
sx q[1];
rz(-1.4828955) q[1];
sx q[1];
rz(1.0465266) q[1];
rz(-pi) q[2];
rz(1.9870342) q[3];
sx q[3];
rz(-1.7572973) q[3];
sx q[3];
rz(-2.4881252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4992497) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(-2.3976105) q[2];
rz(2.3104525) q[3];
sx q[3];
rz(-1.7404514) q[3];
sx q[3];
rz(-2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0396073) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(2.4264477) q[0];
rz(1.2932418) q[1];
sx q[1];
rz(-1.5512385) q[1];
sx q[1];
rz(2.901758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3837357) q[0];
sx q[0];
rz(-2.3546887) q[0];
sx q[0];
rz(-0.73946865) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5101899) q[2];
sx q[2];
rz(-1.0088293) q[2];
sx q[2];
rz(2.2417703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4308987) q[1];
sx q[1];
rz(-1.7393117) q[1];
sx q[1];
rz(-0.88697834) q[1];
rz(-pi) q[2];
rz(0.062021755) q[3];
sx q[3];
rz(-2.6497132) q[3];
sx q[3];
rz(2.5700701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0665325) q[2];
sx q[2];
rz(-1.5542474) q[2];
sx q[2];
rz(1.4934322) q[2];
rz(-3.1028124) q[3];
sx q[3];
rz(-1.4977692) q[3];
sx q[3];
rz(-1.405976) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1455602) q[0];
sx q[0];
rz(-1.4936438) q[0];
sx q[0];
rz(-2.5160568) q[0];
rz(-1.7690376) q[1];
sx q[1];
rz(-2.1452429) q[1];
sx q[1];
rz(-2.5575584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2997949) q[0];
sx q[0];
rz(-1.3158176) q[0];
sx q[0];
rz(0.38354446) q[0];
rz(-2.4768684) q[2];
sx q[2];
rz(-1.692284) q[2];
sx q[2];
rz(-0.45035502) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0528193) q[1];
sx q[1];
rz(-0.86630834) q[1];
sx q[1];
rz(-1.8584083) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3006421) q[3];
sx q[3];
rz(-2.2380201) q[3];
sx q[3];
rz(3.0285578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11223101) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(1.7870129) q[2];
rz(-2.3692047) q[3];
sx q[3];
rz(-1.9410917) q[3];
sx q[3];
rz(-1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1188485) q[0];
sx q[0];
rz(-0.72853715) q[0];
sx q[0];
rz(-1.5991612) q[0];
rz(2.6842885) q[1];
sx q[1];
rz(-2.2685469) q[1];
sx q[1];
rz(-0.1874478) q[1];
rz(-0.078866581) q[2];
sx q[2];
rz(-1.7577684) q[2];
sx q[2];
rz(-2.656542) q[2];
rz(-1.1111835) q[3];
sx q[3];
rz(-1.8596953) q[3];
sx q[3];
rz(0.019712899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
