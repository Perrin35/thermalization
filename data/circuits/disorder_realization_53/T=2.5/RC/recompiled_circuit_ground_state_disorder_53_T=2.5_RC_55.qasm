OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98874918) q[0];
sx q[0];
rz(-0.56892836) q[0];
sx q[0];
rz(2.1898354) q[0];
rz(0.8079575) q[1];
sx q[1];
rz(-3.0264049) q[1];
sx q[1];
rz(0.99339956) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4848413) q[0];
sx q[0];
rz(-0.97545058) q[0];
sx q[0];
rz(-0.30466051) q[0];
x q[1];
rz(2.1615218) q[2];
sx q[2];
rz(-0.23071846) q[2];
sx q[2];
rz(-1.7144817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7639072) q[1];
sx q[1];
rz(-2.0983258) q[1];
sx q[1];
rz(-2.1930013) q[1];
rz(1.7073329) q[3];
sx q[3];
rz(-2.2280558) q[3];
sx q[3];
rz(0.11634377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6318165) q[2];
sx q[2];
rz(-0.6002554) q[2];
sx q[2];
rz(0.14524761) q[2];
rz(0.97233573) q[3];
sx q[3];
rz(-1.626868) q[3];
sx q[3];
rz(1.2783031) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7250605) q[0];
sx q[0];
rz(-2.5711377) q[0];
sx q[0];
rz(0.78713) q[0];
rz(2.9540673) q[1];
sx q[1];
rz(-1.3860605) q[1];
sx q[1];
rz(-3.131391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8387209) q[0];
sx q[0];
rz(-0.10940675) q[0];
sx q[0];
rz(1.495269) q[0];
rz(-pi) q[1];
rz(2.0554916) q[2];
sx q[2];
rz(-1.5124413) q[2];
sx q[2];
rz(-1.7444407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9157313) q[1];
sx q[1];
rz(-0.26090947) q[1];
sx q[1];
rz(-0.4930851) q[1];
rz(-1.418347) q[3];
sx q[3];
rz(-0.746519) q[3];
sx q[3];
rz(0.47990769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0804704) q[2];
sx q[2];
rz(-3.0219813) q[2];
sx q[2];
rz(-1.97869) q[2];
rz(-1.844678) q[3];
sx q[3];
rz(-1.7087405) q[3];
sx q[3];
rz(-0.0059303693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6892683) q[0];
sx q[0];
rz(-0.56832123) q[0];
sx q[0];
rz(-0.34533056) q[0];
rz(3.0863702) q[1];
sx q[1];
rz(-2.5352812) q[1];
sx q[1];
rz(2.9323554) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064406618) q[0];
sx q[0];
rz(-1.4923054) q[0];
sx q[0];
rz(-2.2839943) q[0];
rz(1.3232068) q[2];
sx q[2];
rz(-0.074562975) q[2];
sx q[2];
rz(0.32081315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.71007939) q[1];
sx q[1];
rz(-2.4071879) q[1];
sx q[1];
rz(2.9148352) q[1];
rz(-pi) q[2];
rz(1.8322631) q[3];
sx q[3];
rz(-0.76305721) q[3];
sx q[3];
rz(-2.7558568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.050134) q[2];
sx q[2];
rz(-1.4860934) q[2];
sx q[2];
rz(2.6996108) q[2];
rz(0.4778536) q[3];
sx q[3];
rz(-1.7171532) q[3];
sx q[3];
rz(1.2598239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.1314989) q[0];
sx q[0];
rz(-0.57259125) q[0];
sx q[0];
rz(1.2226489) q[0];
rz(1.6652416) q[1];
sx q[1];
rz(-1.0226701) q[1];
sx q[1];
rz(-2.8388035) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0500844) q[0];
sx q[0];
rz(-1.2909527) q[0];
sx q[0];
rz(-0.39892964) q[0];
x q[1];
rz(-1.6499182) q[2];
sx q[2];
rz(-1.5357568) q[2];
sx q[2];
rz(-1.5372045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5440793) q[1];
sx q[1];
rz(-1.5248393) q[1];
sx q[1];
rz(1.8848898) q[1];
rz(0.92118951) q[3];
sx q[3];
rz(-1.0471034) q[3];
sx q[3];
rz(0.37674784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0786232) q[2];
sx q[2];
rz(-1.6958232) q[2];
sx q[2];
rz(-1.7930188) q[2];
rz(0.013269987) q[3];
sx q[3];
rz(-2.477406) q[3];
sx q[3];
rz(1.2191023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.211798) q[0];
sx q[0];
rz(-3.0199265) q[0];
sx q[0];
rz(-1.9450872) q[0];
rz(0.78367805) q[1];
sx q[1];
rz(-1.0107026) q[1];
sx q[1];
rz(-2.3616621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62155277) q[0];
sx q[0];
rz(-0.85264684) q[0];
sx q[0];
rz(1.2214127) q[0];
rz(-2.7042208) q[2];
sx q[2];
rz(-1.6212045) q[2];
sx q[2];
rz(0.1121262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.096141013) q[1];
sx q[1];
rz(-0.91288954) q[1];
sx q[1];
rz(-1.0784574) q[1];
rz(-2.6013589) q[3];
sx q[3];
rz(-1.8461986) q[3];
sx q[3];
rz(-1.6364975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1448867) q[2];
sx q[2];
rz(-0.11652623) q[2];
sx q[2];
rz(0.04280002) q[2];
rz(-1.608611) q[3];
sx q[3];
rz(-1.3442842) q[3];
sx q[3];
rz(-0.90095055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3881989) q[0];
sx q[0];
rz(-2.2008984) q[0];
sx q[0];
rz(0.70145506) q[0];
rz(0.88297168) q[1];
sx q[1];
rz(-1.3621623) q[1];
sx q[1];
rz(2.0515474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7146637) q[0];
sx q[0];
rz(-1.8986456) q[0];
sx q[0];
rz(-1.1529684) q[0];
x q[1];
rz(2.4187614) q[2];
sx q[2];
rz(-1.9246274) q[2];
sx q[2];
rz(-1.3530255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0193041) q[1];
sx q[1];
rz(-1.5221223) q[1];
sx q[1];
rz(-1.4282872) q[1];
rz(1.8102856) q[3];
sx q[3];
rz(-1.7305995) q[3];
sx q[3];
rz(-1.0996475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28114721) q[2];
sx q[2];
rz(-0.88938418) q[2];
sx q[2];
rz(1.3775728) q[2];
rz(-3.0456165) q[3];
sx q[3];
rz(-2.4317661) q[3];
sx q[3];
rz(-0.53203741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1808566) q[0];
sx q[0];
rz(-2.2611698) q[0];
sx q[0];
rz(-1.2197422) q[0];
rz(-0.60918728) q[1];
sx q[1];
rz(-1.6222619) q[1];
sx q[1];
rz(-0.46527299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.935509) q[0];
sx q[0];
rz(-0.92496189) q[0];
sx q[0];
rz(-0.95316621) q[0];
rz(-0.96107595) q[2];
sx q[2];
rz(-1.3111918) q[2];
sx q[2];
rz(-2.0849092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9309931) q[1];
sx q[1];
rz(-0.44856724) q[1];
sx q[1];
rz(1.468352) q[1];
rz(1.3749397) q[3];
sx q[3];
rz(-2.3367685) q[3];
sx q[3];
rz(0.63279654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51218) q[2];
sx q[2];
rz(-0.30567726) q[2];
sx q[2];
rz(-2.7194887) q[2];
rz(-0.01853881) q[3];
sx q[3];
rz(-1.5944578) q[3];
sx q[3];
rz(-2.5277933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9372395) q[0];
sx q[0];
rz(-0.48870191) q[0];
sx q[0];
rz(-0.82116425) q[0];
rz(-0.69860727) q[1];
sx q[1];
rz(-2.3804074) q[1];
sx q[1];
rz(-0.0084812976) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9152731) q[0];
sx q[0];
rz(-2.2157359) q[0];
sx q[0];
rz(-1.1862832) q[0];
rz(-pi) q[1];
rz(1.0462748) q[2];
sx q[2];
rz(-1.2226186) q[2];
sx q[2];
rz(-0.55931567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5183309) q[1];
sx q[1];
rz(-2.3185638) q[1];
sx q[1];
rz(-2.7743913) q[1];
rz(-pi) q[2];
rz(-0.92681411) q[3];
sx q[3];
rz(-1.5379526) q[3];
sx q[3];
rz(-0.72073267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2681793) q[2];
sx q[2];
rz(-0.11205967) q[2];
sx q[2];
rz(2.6123135) q[2];
rz(1.4389634) q[3];
sx q[3];
rz(-1.3114503) q[3];
sx q[3];
rz(0.24229351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5731803) q[0];
sx q[0];
rz(-2.4127164) q[0];
sx q[0];
rz(-0.53701425) q[0];
rz(0.88045398) q[1];
sx q[1];
rz(-1.302779) q[1];
sx q[1];
rz(-2.3939078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76067257) q[0];
sx q[0];
rz(-0.88762807) q[0];
sx q[0];
rz(-1.1364514) q[0];
rz(-pi) q[1];
rz(2.2373392) q[2];
sx q[2];
rz(-1.8996933) q[2];
sx q[2];
rz(2.5118206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97148731) q[1];
sx q[1];
rz(-0.4071281) q[1];
sx q[1];
rz(-0.92252964) q[1];
rz(2.5292812) q[3];
sx q[3];
rz(-1.3638168) q[3];
sx q[3];
rz(1.6918382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.487454) q[2];
sx q[2];
rz(-2.4345001) q[2];
sx q[2];
rz(2.1716165) q[2];
rz(-1.6449432) q[3];
sx q[3];
rz(-2.1275438) q[3];
sx q[3];
rz(2.1192726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3269761) q[0];
sx q[0];
rz(-1.6642445) q[0];
sx q[0];
rz(0.60272637) q[0];
rz(3.1372435) q[1];
sx q[1];
rz(-1.3163722) q[1];
sx q[1];
rz(-2.0306921) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85767001) q[0];
sx q[0];
rz(-2.1956148) q[0];
sx q[0];
rz(0.83832534) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87238042) q[2];
sx q[2];
rz(-1.1495628) q[2];
sx q[2];
rz(-2.5861458) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3746158) q[1];
sx q[1];
rz(-2.5363508) q[1];
sx q[1];
rz(0.21506439) q[1];
x q[2];
rz(0.59278005) q[3];
sx q[3];
rz(-2.5912632) q[3];
sx q[3];
rz(1.5717446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5000308) q[2];
sx q[2];
rz(-1.9226473) q[2];
sx q[2];
rz(-0.59263372) q[2];
rz(-0.57989132) q[3];
sx q[3];
rz(-0.65120828) q[3];
sx q[3];
rz(1.4402703) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92689571) q[0];
sx q[0];
rz(-1.0746645) q[0];
sx q[0];
rz(-0.97434531) q[0];
rz(3.123507) q[1];
sx q[1];
rz(-0.34307243) q[1];
sx q[1];
rz(-3.051563) q[1];
rz(0.68721107) q[2];
sx q[2];
rz(-0.39068475) q[2];
sx q[2];
rz(3.0013553) q[2];
rz(2.6206489) q[3];
sx q[3];
rz(-0.81650618) q[3];
sx q[3];
rz(3.0405844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
