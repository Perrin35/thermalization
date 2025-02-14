OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(0.37707075) q[1];
sx q[1];
rz(-1.7989676) q[1];
sx q[1];
rz(1.0599729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19883991) q[0];
sx q[0];
rz(-1.4497541) q[0];
sx q[0];
rz(-0.11902703) q[0];
x q[1];
rz(-3.0994514) q[2];
sx q[2];
rz(-2.4486447) q[2];
sx q[2];
rz(2.9646089) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2419338) q[1];
sx q[1];
rz(-0.32131689) q[1];
sx q[1];
rz(-2.4216273) q[1];
x q[2];
rz(2.9971231) q[3];
sx q[3];
rz(-1.472755) q[3];
sx q[3];
rz(1.6603744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9609191) q[2];
sx q[2];
rz(-0.91327614) q[2];
sx q[2];
rz(1.712435) q[2];
rz(2.5299431) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(3.070224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3058158) q[0];
sx q[0];
rz(-2.340305) q[0];
sx q[0];
rz(0.74547705) q[0];
rz(-1.9396793) q[1];
sx q[1];
rz(-1.3246091) q[1];
sx q[1];
rz(-2.1048996) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7734402) q[0];
sx q[0];
rz(-0.37640171) q[0];
sx q[0];
rz(-1.3163811) q[0];
rz(1.5444504) q[2];
sx q[2];
rz(-2.3374632) q[2];
sx q[2];
rz(-0.94601099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0526317) q[1];
sx q[1];
rz(-2.3792069) q[1];
sx q[1];
rz(0.41684581) q[1];
rz(-2.8193982) q[3];
sx q[3];
rz(-0.60880843) q[3];
sx q[3];
rz(-0.52458602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1408954) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(-2.1689283) q[2];
rz(0.85727143) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.641441) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(2.3386173) q[0];
rz(0.650644) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(0.21779901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686435) q[0];
sx q[0];
rz(-1.1945219) q[0];
sx q[0];
rz(-0.89366389) q[0];
rz(-pi) q[1];
rz(-0.62816633) q[2];
sx q[2];
rz(-1.7587314) q[2];
sx q[2];
rz(-2.4357901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.928339) q[1];
sx q[1];
rz(-1.6633185) q[1];
sx q[1];
rz(1.0800581) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5434615) q[3];
sx q[3];
rz(-2.7418828) q[3];
sx q[3];
rz(-2.3593115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95152068) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(1.5104843) q[2];
rz(-1.914628) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25201061) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(-2.4776283) q[0];
rz(-3.0162759) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(-2.1024316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1215724) q[0];
sx q[0];
rz(-1.5680321) q[0];
sx q[0];
rz(-3.1317668) q[0];
rz(-1.9661994) q[2];
sx q[2];
rz(-1.3229473) q[2];
sx q[2];
rz(1.445766) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0137456) q[1];
sx q[1];
rz(-1.7798073) q[1];
sx q[1];
rz(-2.0621939) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9431875) q[3];
sx q[3];
rz(-2.6909749) q[3];
sx q[3];
rz(3.0991447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(-0.077979716) q[2];
rz(2.0237427) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(-1.0361205) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9341105) q[0];
sx q[0];
rz(-2.9353862) q[0];
sx q[0];
rz(0.41314405) q[0];
rz(-1.235599) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(-1.8280169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.006047) q[0];
sx q[0];
rz(-2.1925547) q[0];
sx q[0];
rz(-2.9154587) q[0];
rz(-pi) q[1];
rz(-2.7817725) q[2];
sx q[2];
rz(-2.0096547) q[2];
sx q[2];
rz(-3.1391337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8660248) q[1];
sx q[1];
rz(-2.0179664) q[1];
sx q[1];
rz(2.7204896) q[1];
rz(-0.78727874) q[3];
sx q[3];
rz(-1.2378927) q[3];
sx q[3];
rz(1.6407788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1818587) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(-2.6971297) q[2];
rz(-1.2005165) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(-0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34894094) q[0];
sx q[0];
rz(-0.3949202) q[0];
sx q[0];
rz(-3.0643903) q[0];
rz(-0.96241799) q[1];
sx q[1];
rz(-1.9541157) q[1];
sx q[1];
rz(0.033800689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8706521) q[0];
sx q[0];
rz(-1.1198638) q[0];
sx q[0];
rz(-1.1876039) q[0];
rz(0.044174657) q[2];
sx q[2];
rz(-0.54810537) q[2];
sx q[2];
rz(-1.4352611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5742366) q[1];
sx q[1];
rz(-2.1891382) q[1];
sx q[1];
rz(0.38783698) q[1];
x q[2];
rz(1.46557) q[3];
sx q[3];
rz(-2.3447373) q[3];
sx q[3];
rz(-1.5060671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9683711) q[2];
sx q[2];
rz(-1.5366448) q[2];
sx q[2];
rz(0.25924337) q[2];
rz(0.94414583) q[3];
sx q[3];
rz(-2.9278432) q[3];
sx q[3];
rz(1.3313782) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060870085) q[0];
sx q[0];
rz(-2.935077) q[0];
sx q[0];
rz(-0.61332214) q[0];
rz(0.90336409) q[1];
sx q[1];
rz(-2.3009243) q[1];
sx q[1];
rz(-0.14796743) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8362691) q[0];
sx q[0];
rz(-0.9137872) q[0];
sx q[0];
rz(-0.61886529) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7025928) q[2];
sx q[2];
rz(-1.5996965) q[2];
sx q[2];
rz(2.7759107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92098713) q[1];
sx q[1];
rz(-1.9927532) q[1];
sx q[1];
rz(0.76442952) q[1];
x q[2];
rz(2.1359753) q[3];
sx q[3];
rz(-1.5907173) q[3];
sx q[3];
rz(1.7651778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2030187) q[2];
sx q[2];
rz(-1.4116986) q[2];
sx q[2];
rz(-0.99986783) q[2];
rz(-1.3153007) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(-0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95241791) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(2.1536105) q[0];
rz(1.8723764) q[1];
sx q[1];
rz(-2.3321584) q[1];
sx q[1];
rz(-0.16407897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7456147) q[0];
sx q[0];
rz(-1.4854637) q[0];
sx q[0];
rz(3.0239168) q[0];
rz(1.1804232) q[2];
sx q[2];
rz(-0.52902856) q[2];
sx q[2];
rz(1.0323553) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9050818) q[1];
sx q[1];
rz(-1.0781546) q[1];
sx q[1];
rz(0.0035159455) q[1];
rz(-pi) q[2];
rz(-1.3173728) q[3];
sx q[3];
rz(-2.4755423) q[3];
sx q[3];
rz(-2.1061153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.24171955) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(-1.0618173) q[2];
rz(-0.1344943) q[3];
sx q[3];
rz(-0.89956346) q[3];
sx q[3];
rz(-0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6580842) q[0];
sx q[0];
rz(-2.4299419) q[0];
sx q[0];
rz(-0.2051951) q[0];
rz(2.9070053) q[1];
sx q[1];
rz(-1.7154452) q[1];
sx q[1];
rz(-2.9451784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.097986624) q[0];
sx q[0];
rz(-2.4825486) q[0];
sx q[0];
rz(2.090467) q[0];
rz(1.4142198) q[2];
sx q[2];
rz(-1.051826) q[2];
sx q[2];
rz(1.9760946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8512021) q[1];
sx q[1];
rz(-2.7702077) q[1];
sx q[1];
rz(-1.0642278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.100575) q[3];
sx q[3];
rz(-0.8634206) q[3];
sx q[3];
rz(1.1288957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74497574) q[2];
sx q[2];
rz(-1.9605512) q[2];
sx q[2];
rz(-1.1783925) q[2];
rz(-0.34934238) q[3];
sx q[3];
rz(-2.0125466) q[3];
sx q[3];
rz(1.0497302) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0013393764) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(-0.70571357) q[0];
rz(-2.4440675) q[1];
sx q[1];
rz(-1.768528) q[1];
sx q[1];
rz(0.90550214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5030497) q[0];
sx q[0];
rz(-2.7075504) q[0];
sx q[0];
rz(-0.040157138) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.414213) q[2];
sx q[2];
rz(-2.670521) q[2];
sx q[2];
rz(-1.7076275) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6433555) q[1];
sx q[1];
rz(-1.61803) q[1];
sx q[1];
rz(-2.8308353) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8373413) q[3];
sx q[3];
rz(-0.47347927) q[3];
sx q[3];
rz(-2.480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2424348) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(-0.71851292) q[2];
rz(-2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8231507) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(-1.5009343) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(2.9983042) q[2];
sx q[2];
rz(-1.6874416) q[2];
sx q[2];
rz(-0.018358827) q[2];
rz(-0.80397687) q[3];
sx q[3];
rz(-2.0686276) q[3];
sx q[3];
rz(-0.008240464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
