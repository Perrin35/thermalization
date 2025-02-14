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
rz(7.0452375) q[0];
sx q[0];
rz(10.308916) q[0];
rz(2.60131) q[1];
sx q[1];
rz(1.5958495) q[1];
sx q[1];
rz(11.933239) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5954689) q[0];
sx q[0];
rz(-1.7350539) q[0];
sx q[0];
rz(-1.4558653) q[0];
rz(-0.22692899) q[2];
sx q[2];
rz(-0.49751505) q[2];
sx q[2];
rz(-2.016248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31189141) q[1];
sx q[1];
rz(-0.30712488) q[1];
sx q[1];
rz(2.7077371) q[1];
rz(-0.54852529) q[3];
sx q[3];
rz(-1.9270634) q[3];
sx q[3];
rz(-1.4725722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7264318) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(2.4994182) q[3];
sx q[3];
rz(-0.29688501) q[3];
sx q[3];
rz(1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.967763) q[0];
sx q[0];
rz(-2.1673295) q[0];
sx q[0];
rz(-0.94671384) q[0];
rz(3.0831378) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(2.2296925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0094441) q[0];
sx q[0];
rz(-1.5329038) q[0];
sx q[0];
rz(1.9873585) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2445685) q[2];
sx q[2];
rz(-1.5155715) q[2];
sx q[2];
rz(-1.8230452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2649041) q[1];
sx q[1];
rz(-1.1238681) q[1];
sx q[1];
rz(-1.8403711) q[1];
rz(-pi) q[2];
rz(1.5188317) q[3];
sx q[3];
rz(-1.1295737) q[3];
sx q[3];
rz(-3.0558464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6902206) q[2];
sx q[2];
rz(-2.8896832) q[2];
sx q[2];
rz(-2.8042262) q[2];
rz(-0.99006027) q[3];
sx q[3];
rz(-1.3288574) q[3];
sx q[3];
rz(-1.2216074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57872355) q[0];
sx q[0];
rz(-0.27414411) q[0];
sx q[0];
rz(-1.8291871) q[0];
rz(-2.9263272) q[1];
sx q[1];
rz(-1.637641) q[1];
sx q[1];
rz(-0.85830918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4612573) q[0];
sx q[0];
rz(-1.8363928) q[0];
sx q[0];
rz(-3.0895825) q[0];
rz(0.99719259) q[2];
sx q[2];
rz(-0.39948398) q[2];
sx q[2];
rz(1.4956719) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2006113) q[1];
sx q[1];
rz(-2.0075782) q[1];
sx q[1];
rz(-1.753356) q[1];
x q[2];
rz(-1.8668601) q[3];
sx q[3];
rz(-1.8157457) q[3];
sx q[3];
rz(1.7942269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.12576088) q[2];
sx q[2];
rz(-2.5845926) q[2];
sx q[2];
rz(-2.1324943) q[2];
rz(-0.5674181) q[3];
sx q[3];
rz(-1.5083454) q[3];
sx q[3];
rz(-1.6217888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0030982) q[0];
sx q[0];
rz(-2.4479285) q[0];
sx q[0];
rz(-2.9487115) q[0];
rz(2.0774972) q[1];
sx q[1];
rz(-0.27609438) q[1];
sx q[1];
rz(2.1853133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4665657) q[0];
sx q[0];
rz(-2.1888632) q[0];
sx q[0];
rz(1.3542372) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5474693) q[2];
sx q[2];
rz(-0.99677873) q[2];
sx q[2];
rz(-1.627315) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5259486) q[1];
sx q[1];
rz(-2.5275349) q[1];
sx q[1];
rz(2.581852) q[1];
x q[2];
rz(1.9547012) q[3];
sx q[3];
rz(-0.70261803) q[3];
sx q[3];
rz(1.5877569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8320147) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(0.94069707) q[2];
rz(-2.7637123) q[3];
sx q[3];
rz(-0.89972073) q[3];
sx q[3];
rz(-0.64062029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2864895) q[0];
sx q[0];
rz(-1.6096977) q[0];
sx q[0];
rz(1.3252277) q[0];
rz(-2.8638966) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(-1.206254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6645443) q[0];
sx q[0];
rz(-1.7627075) q[0];
sx q[0];
rz(0.080587793) q[0];
rz(2.5304753) q[2];
sx q[2];
rz(-1.9131994) q[2];
sx q[2];
rz(2.8813206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5884456) q[1];
sx q[1];
rz(-0.35970769) q[1];
sx q[1];
rz(-1.5313696) q[1];
x q[2];
rz(0.60449497) q[3];
sx q[3];
rz(-0.2911199) q[3];
sx q[3];
rz(-1.1230942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9557934) q[2];
sx q[2];
rz(-2.9904371) q[2];
sx q[2];
rz(1.7929057) q[2];
rz(1.0079314) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.86001992) q[1];
sx q[1];
rz(-0.42541447) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2046609) q[0];
sx q[0];
rz(-0.46260329) q[0];
sx q[0];
rz(3.0499056) q[0];
x q[1];
rz(-1.2666918) q[2];
sx q[2];
rz(-1.6124523) q[2];
sx q[2];
rz(-1.1343985) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0531333) q[1];
sx q[1];
rz(-1.8324553) q[1];
sx q[1];
rz(0.43105521) q[1];
rz(-pi) q[2];
rz(2.4454771) q[3];
sx q[3];
rz(-0.98721877) q[3];
sx q[3];
rz(-2.7903874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5328488) q[2];
sx q[2];
rz(-3.0259735) q[2];
sx q[2];
rz(2.4382584) q[2];
rz(-0.15130791) q[3];
sx q[3];
rz(-1.1727099) q[3];
sx q[3];
rz(2.5942588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.5587826) q[0];
sx q[0];
rz(-2.3334239) q[0];
sx q[0];
rz(2.1852344) q[0];
rz(-0.76902485) q[1];
sx q[1];
rz(-1.290657) q[1];
sx q[1];
rz(2.1020611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9120662) q[0];
sx q[0];
rz(-1.5472224) q[0];
sx q[0];
rz(1.0957155) q[0];
x q[1];
rz(0.65724397) q[2];
sx q[2];
rz(-0.98117706) q[2];
sx q[2];
rz(2.5319665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7733147) q[1];
sx q[1];
rz(-1.8551143) q[1];
sx q[1];
rz(-1.8892688) q[1];
rz(-1.9045181) q[3];
sx q[3];
rz(-1.3710554) q[3];
sx q[3];
rz(0.39272308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9397554) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(-1.5038097) q[2];
rz(2.2683897) q[3];
sx q[3];
rz(-0.88770294) q[3];
sx q[3];
rz(0.33151117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7694976) q[0];
sx q[0];
rz(-2.7517419) q[0];
sx q[0];
rz(-1.9918342) q[0];
rz(-2.276543) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(1.4335776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2835049) q[0];
sx q[0];
rz(-2.4437987) q[0];
sx q[0];
rz(2.1608864) q[0];
x q[1];
rz(2.6636567) q[2];
sx q[2];
rz(-1.897287) q[2];
sx q[2];
rz(-0.70604347) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3506602) q[1];
sx q[1];
rz(-2.6106842) q[1];
sx q[1];
rz(-1.3965308) q[1];
rz(-pi) q[2];
rz(-1.1341489) q[3];
sx q[3];
rz(-0.45387196) q[3];
sx q[3];
rz(-1.8271104) q[3];
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
rz(-2.3104525) q[3];
sx q[3];
rz(-1.7404514) q[3];
sx q[3];
rz(2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1019854) q[0];
sx q[0];
rz(-0.87431854) q[0];
sx q[0];
rz(2.4264477) q[0];
rz(1.2932418) q[1];
sx q[1];
rz(-1.5512385) q[1];
sx q[1];
rz(2.901758) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3837357) q[0];
sx q[0];
rz(-0.78690398) q[0];
sx q[0];
rz(-2.402124) q[0];
x q[1];
rz(1.6314028) q[2];
sx q[2];
rz(-1.0088293) q[2];
sx q[2];
rz(0.89982239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71069392) q[1];
sx q[1];
rz(-1.7393117) q[1];
sx q[1];
rz(-2.2546143) q[1];
rz(-pi) q[2];
rz(2.650514) q[3];
sx q[3];
rz(-1.5415191) q[3];
sx q[3];
rz(1.0539583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0750601) q[2];
sx q[2];
rz(-1.5542474) q[2];
sx q[2];
rz(-1.4934322) q[2];
rz(0.0387803) q[3];
sx q[3];
rz(-1.6438234) q[3];
sx q[3];
rz(1.405976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1455602) q[0];
sx q[0];
rz(-1.6479489) q[0];
sx q[0];
rz(2.5160568) q[0];
rz(1.372555) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(2.5575584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9720219) q[0];
sx q[0];
rz(-1.2002647) q[0];
sx q[0];
rz(-1.8448004) q[0];
rz(-pi) q[1];
rz(-1.4169078) q[2];
sx q[2];
rz(-2.2297572) q[2];
sx q[2];
rz(1.9264592) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4702813) q[1];
sx q[1];
rz(-1.7886506) q[1];
sx q[1];
rz(2.4163428) q[1];
rz(-pi) q[2];
rz(-2.8149393) q[3];
sx q[3];
rz(-2.4295837) q[3];
sx q[3];
rz(2.6077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0293616) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(1.7870129) q[2];
rz(0.77238798) q[3];
sx q[3];
rz(-1.200501) q[3];
sx q[3];
rz(-1.5310009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1188485) q[0];
sx q[0];
rz(-0.72853715) q[0];
sx q[0];
rz(-1.5991612) q[0];
rz(-0.45730418) q[1];
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
