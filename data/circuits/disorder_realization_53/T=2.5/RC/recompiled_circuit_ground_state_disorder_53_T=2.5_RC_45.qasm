OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1528435) q[0];
sx q[0];
rz(-2.5726643) q[0];
sx q[0];
rz(-2.1898354) q[0];
rz(0.8079575) q[1];
sx q[1];
rz(-3.0264049) q[1];
sx q[1];
rz(0.99339956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4848413) q[0];
sx q[0];
rz(-2.1661421) q[0];
sx q[0];
rz(2.8369321) q[0];
rz(-pi) q[1];
rz(-1.7634704) q[2];
sx q[2];
rz(-1.698508) q[2];
sx q[2];
rz(-0.72205262) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.722734) q[1];
sx q[1];
rz(-0.79241792) q[1];
sx q[1];
rz(0.78559135) q[1];
x q[2];
rz(2.967013) q[3];
sx q[3];
rz(-0.66922656) q[3];
sx q[3];
rz(2.804047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50977612) q[2];
sx q[2];
rz(-0.6002554) q[2];
sx q[2];
rz(-2.996345) q[2];
rz(-0.97233573) q[3];
sx q[3];
rz(-1.626868) q[3];
sx q[3];
rz(1.8632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7250605) q[0];
sx q[0];
rz(-0.57045492) q[0];
sx q[0];
rz(0.78713) q[0];
rz(2.9540673) q[1];
sx q[1];
rz(-1.3860605) q[1];
sx q[1];
rz(0.010201605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30287179) q[0];
sx q[0];
rz(-0.10940675) q[0];
sx q[0];
rz(-1.6463237) q[0];
rz(-pi) q[1];
rz(-1.446063) q[2];
sx q[2];
rz(-0.48791781) q[2];
sx q[2];
rz(-0.063340576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4233077) q[1];
sx q[1];
rz(-1.8000291) q[1];
sx q[1];
rz(1.4450816) q[1];
x q[2];
rz(-1.7232456) q[3];
sx q[3];
rz(-2.3950737) q[3];
sx q[3];
rz(-2.661685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.06112222) q[2];
sx q[2];
rz(-3.0219813) q[2];
sx q[2];
rz(-1.1629026) q[2];
rz(-1.844678) q[3];
sx q[3];
rz(-1.7087405) q[3];
sx q[3];
rz(-0.0059303693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6892683) q[0];
sx q[0];
rz(-0.56832123) q[0];
sx q[0];
rz(-2.7962621) q[0];
rz(0.055222424) q[1];
sx q[1];
rz(-0.60631141) q[1];
sx q[1];
rz(-0.20923722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.077186) q[0];
sx q[0];
rz(-1.4923054) q[0];
sx q[0];
rz(-0.85759832) q[0];
x q[1];
rz(-0.018304869) q[2];
sx q[2];
rz(-1.6430815) q[2];
sx q[2];
rz(-2.5725281) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40869409) q[1];
sx q[1];
rz(-2.2823221) q[1];
sx q[1];
rz(1.3705322) q[1];
rz(-1.3093295) q[3];
sx q[3];
rz(-0.76305721) q[3];
sx q[3];
rz(0.3857358) q[3];
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
rz(-2.6637391) q[3];
sx q[3];
rz(-1.7171532) q[3];
sx q[3];
rz(-1.8817687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(2.0100937) q[0];
sx q[0];
rz(-0.57259125) q[0];
sx q[0];
rz(-1.9189438) q[0];
rz(1.6652416) q[1];
sx q[1];
rz(-2.1189225) q[1];
sx q[1];
rz(-0.30278912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0593392) q[0];
sx q[0];
rz(-2.6586464) q[0];
sx q[0];
rz(2.5046136) q[0];
x q[1];
rz(3.1064433) q[2];
sx q[2];
rz(-1.4917231) q[2];
sx q[2];
rz(3.1052232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.099953) q[1];
sx q[1];
rz(-1.8845468) q[1];
sx q[1];
rz(0.048317475) q[1];
x q[2];
rz(-2.2204031) q[3];
sx q[3];
rz(-2.0944893) q[3];
sx q[3];
rz(-0.37674784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0786232) q[2];
sx q[2];
rz(-1.4457694) q[2];
sx q[2];
rz(-1.7930188) q[2];
rz(3.1283227) q[3];
sx q[3];
rz(-2.477406) q[3];
sx q[3];
rz(-1.2191023) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9297946) q[0];
sx q[0];
rz(-0.12166611) q[0];
sx q[0];
rz(-1.9450872) q[0];
rz(-0.78367805) q[1];
sx q[1];
rz(-2.13089) q[1];
sx q[1];
rz(0.7799305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570628) q[0];
sx q[0];
rz(-1.3100782) q[0];
sx q[0];
rz(0.74913673) q[0];
rz(-pi) q[1];
rz(3.0230396) q[2];
sx q[2];
rz(-2.701512) q[2];
sx q[2];
rz(1.5660182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1577252) q[1];
sx q[1];
rz(-1.1874501) q[1];
sx q[1];
rz(-2.4216837) q[1];
rz(-pi) q[2];
rz(0.50239222) q[3];
sx q[3];
rz(-2.5414782) q[3];
sx q[3];
rz(2.6502018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1448867) q[2];
sx q[2];
rz(-3.0250664) q[2];
sx q[2];
rz(-3.0987926) q[2];
rz(-1.608611) q[3];
sx q[3];
rz(-1.3442842) q[3];
sx q[3];
rz(-0.90095055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75339371) q[0];
sx q[0];
rz(-2.2008984) q[0];
sx q[0];
rz(2.4401376) q[0];
rz(2.258621) q[1];
sx q[1];
rz(-1.3621623) q[1];
sx q[1];
rz(1.0900452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7713666) q[0];
sx q[0];
rz(-0.52506444) q[0];
sx q[0];
rz(0.87319002) q[0];
rz(-2.4187614) q[2];
sx q[2];
rz(-1.2169653) q[2];
sx q[2];
rz(-1.3530255) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5970814) q[1];
sx q[1];
rz(-1.4284572) q[1];
sx q[1];
rz(3.0924209) q[1];
rz(-1.3313071) q[3];
sx q[3];
rz(-1.4109932) q[3];
sx q[3];
rz(-2.0419451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8604454) q[2];
sx q[2];
rz(-2.2522085) q[2];
sx q[2];
rz(1.3775728) q[2];
rz(-0.095976202) q[3];
sx q[3];
rz(-2.4317661) q[3];
sx q[3];
rz(-2.6095552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9607361) q[0];
sx q[0];
rz(-0.88042283) q[0];
sx q[0];
rz(-1.2197422) q[0];
rz(2.5324054) q[1];
sx q[1];
rz(-1.5193308) q[1];
sx q[1];
rz(-2.6763197) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03929654) q[0];
sx q[0];
rz(-1.0900153) q[0];
sx q[0];
rz(-2.395438) q[0];
x q[1];
rz(-0.31330534) q[2];
sx q[2];
rz(-2.1572972) q[2];
sx q[2];
rz(-2.4500606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0969995) q[1];
sx q[1];
rz(-2.0168418) q[1];
sx q[1];
rz(0.049179463) q[1];
x q[2];
rz(-0.77561982) q[3];
sx q[3];
rz(-1.4300775) q[3];
sx q[3];
rz(-1.0746806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6294127) q[2];
sx q[2];
rz(-2.8359154) q[2];
sx q[2];
rz(2.7194887) q[2];
rz(-0.01853881) q[3];
sx q[3];
rz(-1.5471349) q[3];
sx q[3];
rz(-0.6137994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2043532) q[0];
sx q[0];
rz(-0.48870191) q[0];
sx q[0];
rz(2.3204284) q[0];
rz(2.4429854) q[1];
sx q[1];
rz(-0.76118529) q[1];
sx q[1];
rz(-3.1331114) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5584921) q[0];
sx q[0];
rz(-1.8752397) q[0];
sx q[0];
rz(0.68171708) q[0];
rz(2.1979546) q[2];
sx q[2];
rz(-2.5211589) q[2];
sx q[2];
rz(0.47867423) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20336452) q[1];
sx q[1];
rz(-1.3044323) q[1];
sx q[1];
rz(0.78856923) q[1];
rz(-pi) q[2];
x q[2];
rz(0.04106122) q[3];
sx q[3];
rz(-0.92721894) q[3];
sx q[3];
rz(2.3161771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2681793) q[2];
sx q[2];
rz(-3.029533) q[2];
sx q[2];
rz(-2.6123135) q[2];
rz(-1.4389634) q[3];
sx q[3];
rz(-1.8301423) q[3];
sx q[3];
rz(0.24229351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5731803) q[0];
sx q[0];
rz(-2.4127164) q[0];
sx q[0];
rz(-0.53701425) q[0];
rz(2.2611387) q[1];
sx q[1];
rz(-1.8388137) q[1];
sx q[1];
rz(-2.3939078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76067257) q[0];
sx q[0];
rz(-0.88762807) q[0];
sx q[0];
rz(-1.1364514) q[0];
rz(-pi) q[1];
rz(-0.90425348) q[2];
sx q[2];
rz(-1.8996933) q[2];
sx q[2];
rz(2.5118206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1701053) q[1];
sx q[1];
rz(-0.4071281) q[1];
sx q[1];
rz(-0.92252964) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3196094) q[3];
sx q[3];
rz(-2.1681941) q[3];
sx q[3];
rz(-0.26439127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.487454) q[2];
sx q[2];
rz(-2.4345001) q[2];
sx q[2];
rz(2.1716165) q[2];
rz(-1.4966494) q[3];
sx q[3];
rz(-1.0140489) q[3];
sx q[3];
rz(-1.0223201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8146166) q[0];
sx q[0];
rz(-1.6642445) q[0];
sx q[0];
rz(2.5388663) q[0];
rz(-3.1372435) q[1];
sx q[1];
rz(-1.8252204) q[1];
sx q[1];
rz(1.1109005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13647381) q[0];
sx q[0];
rz(-2.2181151) q[0];
sx q[0];
rz(2.3939483) q[0];
rz(-pi) q[1];
rz(-0.87238042) q[2];
sx q[2];
rz(-1.9920298) q[2];
sx q[2];
rz(-2.5861458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0266089) q[1];
sx q[1];
rz(-0.98140684) q[1];
sx q[1];
rz(1.4242054) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5488126) q[3];
sx q[3];
rz(-0.55032941) q[3];
sx q[3];
rz(-1.5717446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6415619) q[2];
sx q[2];
rz(-1.9226473) q[2];
sx q[2];
rz(0.59263372) q[2];
rz(-2.5617013) q[3];
sx q[3];
rz(-0.65120828) q[3];
sx q[3];
rz(-1.4402703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(1.3152348) q[2];
sx q[2];
rz(-1.2719874) q[2];
sx q[2];
rz(-0.86612305) q[2];
rz(0.74541253) q[3];
sx q[3];
rz(-1.941962) q[3];
sx q[3];
rz(1.844187) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
