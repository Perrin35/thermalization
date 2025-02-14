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
rz(1.9999003) q[0];
sx q[0];
rz(6.5394149) q[0];
rz(3.5186634) q[1];
sx q[1];
rz(4.9405603) q[1];
sx q[1];
rz(8.3648051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9427527) q[0];
sx q[0];
rz(-1.4497541) q[0];
sx q[0];
rz(-0.11902703) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.04214123) q[2];
sx q[2];
rz(-0.69294792) q[2];
sx q[2];
rz(-0.17698374) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6457451) q[1];
sx q[1];
rz(-1.331067) q[1];
sx q[1];
rz(-1.354753) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9971231) q[3];
sx q[3];
rz(-1.6688377) q[3];
sx q[3];
rz(1.6603744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18067351) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(-1.712435) q[2];
rz(-0.6116496) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(-0.071368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3058158) q[0];
sx q[0];
rz(-2.340305) q[0];
sx q[0];
rz(2.3961156) q[0];
rz(1.2019134) q[1];
sx q[1];
rz(-1.8169836) q[1];
sx q[1];
rz(-1.036693) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6408143) q[0];
sx q[0];
rz(-1.2070859) q[0];
sx q[0];
rz(3.0424434) q[0];
rz(0.027341893) q[2];
sx q[2];
rz(-2.3745656) q[2];
sx q[2];
rz(-0.98397827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0889609) q[1];
sx q[1];
rz(-0.76238576) q[1];
sx q[1];
rz(0.41684581) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5573331) q[3];
sx q[3];
rz(-1.7528894) q[3];
sx q[3];
rz(0.77891536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1408954) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(-0.97266436) q[2];
rz(-2.2843212) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.5001517) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(2.3386173) q[0];
rz(-2.4909486) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(0.21779901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9150216) q[0];
sx q[0];
rz(-0.75998291) q[0];
sx q[0];
rz(-1.0081916) q[0];
rz(2.5134263) q[2];
sx q[2];
rz(-1.3828613) q[2];
sx q[2];
rz(2.4357901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1863292) q[1];
sx q[1];
rz(-2.6429089) q[1];
sx q[1];
rz(1.7651943) q[1];
rz(-pi) q[2];
rz(2.8057116) q[3];
sx q[3];
rz(-1.7917197) q[3];
sx q[3];
rz(0.22798746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95152068) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(-1.5104843) q[2];
rz(1.914628) q[3];
sx q[3];
rz(-1.0534143) q[3];
sx q[3];
rz(-1.0043043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.889582) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(0.66396436) q[0];
rz(-0.12531677) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(-1.0391611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923415) q[0];
sx q[0];
rz(-1.5806222) q[0];
sx q[0];
rz(1.568032) q[0];
x q[1];
rz(0.2676312) q[2];
sx q[2];
rz(-1.1881141) q[2];
sx q[2];
rz(3.1185993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.127847) q[1];
sx q[1];
rz(-1.7798073) q[1];
sx q[1];
rz(1.0793988) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6658721) q[3];
sx q[3];
rz(-1.1296484) q[3];
sx q[3];
rz(-0.26218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30109721) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(3.0636129) q[2];
rz(-1.1178499) q[3];
sx q[3];
rz(-1.0406787) q[3];
sx q[3];
rz(1.0361205) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2074821) q[0];
sx q[0];
rz(-2.9353862) q[0];
sx q[0];
rz(2.7284486) q[0];
rz(1.235599) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(1.8280169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020371) q[0];
sx q[0];
rz(-1.7540723) q[0];
sx q[0];
rz(0.9367783) q[0];
x q[1];
rz(-0.35982015) q[2];
sx q[2];
rz(-1.1319379) q[2];
sx q[2];
rz(-3.1391337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.48653719) q[1];
sx q[1];
rz(-1.1933206) q[1];
sx q[1];
rz(1.0869763) q[1];
rz(-pi) q[2];
rz(0.78727874) q[3];
sx q[3];
rz(-1.2378927) q[3];
sx q[3];
rz(1.5008139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95973394) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(2.6971297) q[2];
rz(-1.9410761) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926517) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(0.077202395) q[0];
rz(-2.1791747) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(-3.107792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12594906) q[0];
sx q[0];
rz(-1.2275877) q[0];
sx q[0];
rz(2.6604466) q[0];
rz(-pi) q[1];
rz(-0.044174657) q[2];
sx q[2];
rz(-2.5934873) q[2];
sx q[2];
rz(1.7063315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2359743) q[1];
sx q[1];
rz(-1.8840569) q[1];
sx q[1];
rz(-2.226023) q[1];
rz(-pi) q[2];
rz(1.6760227) q[3];
sx q[3];
rz(-0.79685539) q[3];
sx q[3];
rz(1.6355255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.17322156) q[2];
sx q[2];
rz(-1.5366448) q[2];
sx q[2];
rz(2.8823493) q[2];
rz(-2.1974468) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(-1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-3.0807226) q[0];
sx q[0];
rz(-2.935077) q[0];
sx q[0];
rz(2.5282705) q[0];
rz(-2.2382286) q[1];
sx q[1];
rz(-2.3009243) q[1];
sx q[1];
rz(-0.14796743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3053235) q[0];
sx q[0];
rz(-0.9137872) q[0];
sx q[0];
rz(-0.61886529) q[0];
rz(-pi) q[1];
x q[1];
rz(0.067909165) q[2];
sx q[2];
rz(-0.43988827) q[2];
sx q[2];
rz(-1.875017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1175673) q[1];
sx q[1];
rz(-2.2541775) q[1];
sx q[1];
rz(1.014381) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6079801) q[3];
sx q[3];
rz(-2.5761009) q[3];
sx q[3];
rz(0.16298207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2030187) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(-2.1417248) q[2];
rz(1.8262919) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(-0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.1891747) q[0];
sx q[0];
rz(-0.84680951) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(-0.16407897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9768499) q[0];
sx q[0];
rz(-1.4535507) q[0];
sx q[0];
rz(1.4848723) q[0];
x q[1];
rz(-1.0751749) q[2];
sx q[2];
rz(-1.3775423) q[2];
sx q[2];
rz(0.19710625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9050818) q[1];
sx q[1];
rz(-2.063438) q[1];
sx q[1];
rz(3.1380767) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92045356) q[3];
sx q[3];
rz(-1.4152539) q[3];
sx q[3];
rz(2.807164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8998731) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(-1.0618173) q[2];
rz(-0.1344943) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(0.2051951) q[0];
rz(-2.9070053) q[1];
sx q[1];
rz(-1.7154452) q[1];
sx q[1];
rz(2.9451784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72442833) q[0];
sx q[0];
rz(-1.0104033) q[0];
sx q[0];
rz(-2.7743894) q[0];
rz(-1.7273728) q[2];
sx q[2];
rz(-2.0897667) q[2];
sx q[2];
rz(1.165498) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80322335) q[1];
sx q[1];
rz(-1.3937989) q[1];
sx q[1];
rz(-1.8990252) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78077448) q[3];
sx q[3];
rz(-1.9650243) q[3];
sx q[3];
rz(-0.0782644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74497574) q[2];
sx q[2];
rz(-1.1810415) q[2];
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
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1402533) q[0];
sx q[0];
rz(-2.0297191) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(0.69752518) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(-0.90550214) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0374091) q[0];
sx q[0];
rz(-1.5876803) q[0];
sx q[0];
rz(2.707858) q[0];
rz(1.2442676) q[2];
sx q[2];
rz(-1.2249607) q[2];
sx q[2];
rz(-2.2188733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9230629) q[1];
sx q[1];
rz(-0.31421146) q[1];
sx q[1];
rz(-2.988222) q[1];
rz(-pi) q[2];
rz(0.3304146) q[3];
sx q[3];
rz(-1.9163696) q[3];
sx q[3];
rz(-1.4525849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2424348) q[2];
sx q[2];
rz(-2.0826715) q[2];
sx q[2];
rz(-0.71851292) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.9971137) q[3];
sx q[3];
rz(-3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.318442) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(1.5009343) q[1];
sx q[1];
rz(-1.2702912) q[1];
sx q[1];
rz(2.695695) q[1];
rz(0.6871625) q[2];
sx q[2];
rz(-0.18451118) q[2];
sx q[2];
rz(-0.9103734) q[2];
rz(-2.4950776) q[3];
sx q[3];
rz(-0.91520354) q[3];
sx q[3];
rz(-1.147816) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
