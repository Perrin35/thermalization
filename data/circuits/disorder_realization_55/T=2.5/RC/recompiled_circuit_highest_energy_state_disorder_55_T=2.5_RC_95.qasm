OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7776547) q[0];
sx q[0];
rz(-0.41843709) q[0];
sx q[0];
rz(1.8399746) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(1.4459223) q[1];
sx q[1];
rz(12.583153) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.421464) q[0];
sx q[0];
rz(-1.6454433) q[0];
sx q[0];
rz(-1.2089085) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4799825) q[2];
sx q[2];
rz(-0.1802643) q[2];
sx q[2];
rz(1.4992876) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0019579828) q[1];
sx q[1];
rz(-3.1345001) q[1];
sx q[1];
rz(-2.3435739) q[1];
rz(2.7990325) q[3];
sx q[3];
rz(-2.9800219) q[3];
sx q[3];
rz(1.8752961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5590543) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(1.5793229) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(0.16099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6120537) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(1.8147234) q[0];
rz(0.57948411) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(-2.5025867) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20637437) q[0];
sx q[0];
rz(-2.2274096) q[0];
sx q[0];
rz(2.4103863) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5875568) q[2];
sx q[2];
rz(-1.4498561) q[2];
sx q[2];
rz(-0.018509381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.745805) q[1];
sx q[1];
rz(-1.5823872) q[1];
sx q[1];
rz(-3.1245924) q[1];
rz(0.78044807) q[3];
sx q[3];
rz(-1.6326346) q[3];
sx q[3];
rz(2.0831747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7063893) q[2];
sx q[2];
rz(-0.13627626) q[2];
sx q[2];
rz(1.5380247) q[2];
rz(-1.5806574) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(-0.030979009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2494217) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(-2.7601335) q[0];
rz(0.70746607) q[1];
sx q[1];
rz(-3.1222157) q[1];
sx q[1];
rz(-2.0170508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0369506) q[0];
sx q[0];
rz(-2.7890887) q[0];
sx q[0];
rz(2.332325) q[0];
rz(1.7922282) q[2];
sx q[2];
rz(-0.11781684) q[2];
sx q[2];
rz(-3.0750781) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8360236) q[1];
sx q[1];
rz(-1.5828805) q[1];
sx q[1];
rz(-3.0786242) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74060969) q[3];
sx q[3];
rz(-2.4428074) q[3];
sx q[3];
rz(1.6388338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4418929) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(0.049467889) q[2];
rz(2.5345645) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(1.2073257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599051) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(-3.1244151) q[0];
rz(2.848564) q[1];
sx q[1];
rz(-2.3511062) q[1];
sx q[1];
rz(1.5471829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.806694) q[0];
sx q[0];
rz(-2.3133754) q[0];
sx q[0];
rz(-0.81636274) q[0];
rz(-0.54173476) q[2];
sx q[2];
rz(-1.3770665) q[2];
sx q[2];
rz(-1.3044747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80611011) q[1];
sx q[1];
rz(-1.4467738) q[1];
sx q[1];
rz(-1.5107442) q[1];
rz(-pi) q[2];
rz(-0.64597102) q[3];
sx q[3];
rz(-1.5630018) q[3];
sx q[3];
rz(-2.2199059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8881417) q[2];
sx q[2];
rz(-0.47051045) q[2];
sx q[2];
rz(-2.4593501) q[2];
rz(-0.055179723) q[3];
sx q[3];
rz(-0.0076871593) q[3];
sx q[3];
rz(1.3050219) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76160112) q[0];
sx q[0];
rz(-2.70607) q[0];
sx q[0];
rz(-2.7165661) q[0];
rz(1.60166) q[1];
sx q[1];
rz(-2.6585177) q[1];
sx q[1];
rz(-2.3262598) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352933) q[0];
sx q[0];
rz(-1.5599129) q[0];
sx q[0];
rz(-0.061200415) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5940395) q[2];
sx q[2];
rz(-1.5675307) q[2];
sx q[2];
rz(0.68319418) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40774959) q[1];
sx q[1];
rz(-1.6899278) q[1];
sx q[1];
rz(-1.6581737) q[1];
rz(-pi) q[2];
rz(0.32834239) q[3];
sx q[3];
rz(-0.36188618) q[3];
sx q[3];
rz(-0.34235024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1347947) q[2];
sx q[2];
rz(-3.1289913) q[2];
sx q[2];
rz(1.6689782) q[2];
rz(2.2115479) q[3];
sx q[3];
rz(-0.01447066) q[3];
sx q[3];
rz(-0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3227661) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(-1.7148788) q[0];
rz(0.73000437) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(1.1013365) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28601413) q[0];
sx q[0];
rz(-1.8402303) q[0];
sx q[0];
rz(0.84330171) q[0];
x q[1];
rz(-2.1997994) q[2];
sx q[2];
rz(-2.859349) q[2];
sx q[2];
rz(1.4636702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7633865) q[1];
sx q[1];
rz(-0.15608938) q[1];
sx q[1];
rz(-0.915145) q[1];
x q[2];
rz(-1.8227449) q[3];
sx q[3];
rz(-3.0130606) q[3];
sx q[3];
rz(2.4737308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71658984) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(-1.3072183) q[2];
rz(-0.37846765) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(-2.4881261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8398447) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(2.2140455) q[0];
rz(1.7842267) q[1];
sx q[1];
rz(-2.3097242) q[1];
sx q[1];
rz(1.6102788) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0831628) q[0];
sx q[0];
rz(-1.6091804) q[0];
sx q[0];
rz(-3.1192664) q[0];
rz(-pi) q[1];
rz(-0.39674098) q[2];
sx q[2];
rz(-1.7974241) q[2];
sx q[2];
rz(2.6643348) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5824052) q[1];
sx q[1];
rz(-1.6998763) q[1];
sx q[1];
rz(-0.0027133769) q[1];
x q[2];
rz(-1.0590943) q[3];
sx q[3];
rz(-1.8476433) q[3];
sx q[3];
rz(-0.83864318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7808468) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(1.2537664) q[2];
rz(-0.69416657) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(-0.23301253) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53741443) q[0];
sx q[0];
rz(-1.0056714) q[0];
sx q[0];
rz(2.1042714) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(-1.6745837) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1572185) q[0];
sx q[0];
rz(-1.7658002) q[0];
sx q[0];
rz(-0.068679811) q[0];
x q[1];
rz(-0.080341415) q[2];
sx q[2];
rz(-0.27896491) q[2];
sx q[2];
rz(-0.34948784) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.212939) q[1];
sx q[1];
rz(-0.0012081971) q[1];
sx q[1];
rz(1.9429643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19123938) q[3];
sx q[3];
rz(-2.6789224) q[3];
sx q[3];
rz(1.7709875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1303225) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(3.0976963) q[2];
rz(2.6210426) q[3];
sx q[3];
rz(-0.0046516727) q[3];
sx q[3];
rz(-1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540102) q[0];
sx q[0];
rz(-0.0024604877) q[0];
sx q[0];
rz(1.3186697) q[0];
rz(1.7240546) q[1];
sx q[1];
rz(-0.28957614) q[1];
sx q[1];
rz(-1.5444548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8108687) q[0];
sx q[0];
rz(-1.4557342) q[0];
sx q[0];
rz(-2.6601391) q[0];
rz(-0.6575281) q[2];
sx q[2];
rz(-1.3597466) q[2];
sx q[2];
rz(-1.576265) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.006598) q[1];
sx q[1];
rz(-1.8812351) q[1];
sx q[1];
rz(0.18532345) q[1];
x q[2];
rz(-2.0943589) q[3];
sx q[3];
rz(-2.1890543) q[3];
sx q[3];
rz(2.0766192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80457193) q[2];
sx q[2];
rz(-1.2995517) q[2];
sx q[2];
rz(-2.9447832) q[2];
rz(1.9491516) q[3];
sx q[3];
rz(-2.9350023) q[3];
sx q[3];
rz(-2.9406252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8404959) q[0];
sx q[0];
rz(-1.3571285) q[0];
sx q[0];
rz(-1.1916196) q[0];
rz(1.5246897) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(1.5764538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7708262) q[0];
sx q[0];
rz(-1.2773371) q[0];
sx q[0];
rz(-1.7459041) q[0];
rz(1.1168408) q[2];
sx q[2];
rz(-1.364214) q[2];
sx q[2];
rz(1.9898741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8305739) q[1];
sx q[1];
rz(-1.5717447) q[1];
sx q[1];
rz(-1.5702973) q[1];
rz(-0.11348806) q[3];
sx q[3];
rz(-1.4295414) q[3];
sx q[3];
rz(-3.0761513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9578751) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(1.4584165) q[2];
rz(-0.030473907) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(-0.20235801) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6644345) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(1.5674113) q[1];
sx q[1];
rz(-1.3290783) q[1];
sx q[1];
rz(-3.0507416) q[1];
rz(0.0049736918) q[2];
sx q[2];
rz(-1.6463574) q[2];
sx q[2];
rz(-2.9815383) q[2];
rz(-1.0394161) q[3];
sx q[3];
rz(-2.3134091) q[3];
sx q[3];
rz(-0.29534657) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
