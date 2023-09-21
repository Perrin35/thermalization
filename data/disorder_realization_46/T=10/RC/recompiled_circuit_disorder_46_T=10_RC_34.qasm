OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(3.1449218) q[0];
sx q[0];
rz(9.7700906) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(-2.8621434) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011443519) q[0];
sx q[0];
rz(-0.37405095) q[0];
sx q[0];
rz(0.98056294) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3596256) q[2];
sx q[2];
rz(-1.3227533) q[2];
sx q[2];
rz(-2.9458407) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2024723) q[1];
sx q[1];
rz(-1.7608374) q[1];
sx q[1];
rz(-3.128016) q[1];
rz(0.64198288) q[3];
sx q[3];
rz(-1.6025474) q[3];
sx q[3];
rz(0.25105219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44148579) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(0.95735615) q[2];
rz(-2.8422614) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(-2.7295952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608202) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(2.7278996) q[0];
rz(-1.3445688) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(-0.63562524) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5168415) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(-1.2463039) q[0];
rz(-1.7027249) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(-0.012243587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2555274) q[1];
sx q[1];
rz(-0.35554245) q[1];
sx q[1];
rz(-0.76053263) q[1];
rz(-2.440968) q[3];
sx q[3];
rz(-1.0430338) q[3];
sx q[3];
rz(-0.81567314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3893434) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(-2.6129369) q[2];
rz(1.2403437) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(-2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813331) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4814574) q[0];
sx q[0];
rz(-2.0148206) q[0];
sx q[0];
rz(-0.40159479) q[0];
rz(-1.4591818) q[2];
sx q[2];
rz(-0.53225213) q[2];
sx q[2];
rz(-2.3864631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79954051) q[1];
sx q[1];
rz(-2.6375348) q[1];
sx q[1];
rz(-0.44517681) q[1];
x q[2];
rz(0.86359777) q[3];
sx q[3];
rz(-2.2798385) q[3];
sx q[3];
rz(1.9883224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(-0.93786401) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(2.8985033) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(3.0010624) q[0];
rz(-0.17164104) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(-0.26352873) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7292273) q[0];
sx q[0];
rz(-1.0050887) q[0];
sx q[0];
rz(-1.6914781) q[0];
rz(-pi) q[1];
rz(-2.0758923) q[2];
sx q[2];
rz(-2.1720338) q[2];
sx q[2];
rz(-0.87841735) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0861974) q[1];
sx q[1];
rz(-0.89790422) q[1];
sx q[1];
rz(1.124568) q[1];
rz(-pi) q[2];
rz(-1.2372381) q[3];
sx q[3];
rz(-1.8309438) q[3];
sx q[3];
rz(-3.1317657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.2496703) q[2];
rz(-1.9469056) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(0.67888129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-0.93976218) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(3.1125267) q[0];
rz(1.7395696) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(0.32863858) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6304566) q[0];
sx q[0];
rz(-1.4689313) q[0];
sx q[0];
rz(1.0826375) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26942307) q[2];
sx q[2];
rz(-1.4002443) q[2];
sx q[2];
rz(1.5829057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8481816) q[1];
sx q[1];
rz(-1.1383346) q[1];
sx q[1];
rz(-1.6407938) q[1];
rz(-pi) q[2];
rz(0.5248431) q[3];
sx q[3];
rz(-2.3224324) q[3];
sx q[3];
rz(-2.6231349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(2.9809791) q[2];
rz(1.1462071) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(-0.78654003) q[0];
rz(-0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(0.4424817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8724976) q[0];
sx q[0];
rz(-1.2051393) q[0];
sx q[0];
rz(-0.38537607) q[0];
rz(-pi) q[1];
rz(-0.45583506) q[2];
sx q[2];
rz(-0.6243394) q[2];
sx q[2];
rz(2.3221743) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.22873951) q[1];
sx q[1];
rz(-2.733426) q[1];
sx q[1];
rz(-0.91650448) q[1];
rz(-0.093591452) q[3];
sx q[3];
rz(-1.2564661) q[3];
sx q[3];
rz(0.67925727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(-0.9712514) q[2];
rz(-2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-2.1321645) q[0];
rz(-0.55074739) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(1.9783463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35818415) q[0];
sx q[0];
rz(-1.9005214) q[0];
sx q[0];
rz(2.9924336) q[0];
x q[1];
rz(-0.86261729) q[2];
sx q[2];
rz(-1.7628981) q[2];
sx q[2];
rz(-2.5196599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59730676) q[1];
sx q[1];
rz(-1.7339216) q[1];
sx q[1];
rz(0.2225999) q[1];
rz(-pi) q[2];
rz(-0.30927741) q[3];
sx q[3];
rz(-2.2429372) q[3];
sx q[3];
rz(-2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(-2.362137) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(-1.9688169) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9629795) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(-0.12284199) q[0];
rz(-0.12610647) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(1.925148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2467263) q[0];
sx q[0];
rz(-0.77573949) q[0];
sx q[0];
rz(2.5345483) q[0];
x q[1];
rz(0.33151303) q[2];
sx q[2];
rz(-0.47795948) q[2];
sx q[2];
rz(0.31994672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.1864099) q[1];
sx q[1];
rz(1.7840506) q[1];
rz(1.3835195) q[3];
sx q[3];
rz(-0.98687275) q[3];
sx q[3];
rz(2.5919979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-0.043126062) q[2];
rz(2.96636) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(-1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50424987) q[0];
sx q[0];
rz(-0.054333996) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(-1.6437221) q[1];
sx q[1];
rz(-1.4893963) q[1];
sx q[1];
rz(1.1245022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61125206) q[0];
sx q[0];
rz(-1.7065305) q[0];
sx q[0];
rz(1.432857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30377109) q[2];
sx q[2];
rz(-2.4035932) q[2];
sx q[2];
rz(0.27429013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3885755) q[1];
sx q[1];
rz(-1.5696226) q[1];
sx q[1];
rz(-1.4374103) q[1];
x q[2];
rz(2.2665958) q[3];
sx q[3];
rz(-1.9584624) q[3];
sx q[3];
rz(-2.8029203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0008529) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(1.595165) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(-2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7889325) q[0];
sx q[0];
rz(-1.9910318) q[0];
sx q[0];
rz(1.0349405) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4193383) q[0];
sx q[0];
rz(-2.1968578) q[0];
sx q[0];
rz(-1.1115848) q[0];
x q[1];
rz(2.855004) q[2];
sx q[2];
rz(-1.5865241) q[2];
sx q[2];
rz(0.77428267) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1356218) q[1];
sx q[1];
rz(-2.0704198) q[1];
sx q[1];
rz(-2.6702704) q[1];
rz(2.8372739) q[3];
sx q[3];
rz(-2.439194) q[3];
sx q[3];
rz(1.1549293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4118816) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(-0.40851545) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(0.48172054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9020486) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(-1.3760024) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(2.3464936) q[2];
sx q[2];
rz(-1.8021402) q[2];
sx q[2];
rz(-2.6593593) q[2];
rz(1.8134712) q[3];
sx q[3];
rz(-0.96488733) q[3];
sx q[3];
rz(0.78028954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];