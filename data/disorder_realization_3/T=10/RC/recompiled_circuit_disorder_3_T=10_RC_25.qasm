OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(3.1728035) q[0];
sx q[0];
rz(6.7682545) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6611377) q[0];
sx q[0];
rz(-1.9248795) q[0];
sx q[0];
rz(0.35868355) q[0];
x q[1];
rz(1.7180874) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(1.6342083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5536082) q[1];
sx q[1];
rz(-1.506664) q[1];
sx q[1];
rz(-1.9967805) q[1];
x q[2];
rz(2.8768086) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(-0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-2.6020715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(-0.051785843) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9794481) q[2];
sx q[2];
rz(-2.2505629) q[2];
sx q[2];
rz(2.0603927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0961699) q[1];
sx q[1];
rz(-1.759937) q[1];
sx q[1];
rz(2.0778836) q[1];
rz(-2.2803454) q[3];
sx q[3];
rz(-2.0112787) q[3];
sx q[3];
rz(1.4351821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66723055) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(0.24307069) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.8977785) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(2.6932122) q[0];
rz(1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(0.2562491) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1218131) q[0];
sx q[0];
rz(-0.67987961) q[0];
sx q[0];
rz(0.42827423) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44552866) q[2];
sx q[2];
rz(-1.1330714) q[2];
sx q[2];
rz(-2.90403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5970522) q[1];
sx q[1];
rz(-1.2433194) q[1];
sx q[1];
rz(-2.271133) q[1];
rz(-pi) q[2];
rz(-2.9647397) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(2.0166486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-2.4096699) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5314732) q[0];
sx q[0];
rz(-2.0391132) q[0];
sx q[0];
rz(-2.5584695) q[0];
rz(-2.4243083) q[2];
sx q[2];
rz(-1.6589763) q[2];
sx q[2];
rz(2.7862273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1720097) q[1];
sx q[1];
rz(-1.3130377) q[1];
sx q[1];
rz(1.3700563) q[1];
rz(1.6963523) q[3];
sx q[3];
rz(-1.7896277) q[3];
sx q[3];
rz(0.20727508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(2.990492) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.2623825) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(2.343822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7724458) q[0];
sx q[0];
rz(-1.4506842) q[0];
sx q[0];
rz(3.1390879) q[0];
rz(-pi) q[1];
rz(-0.30349489) q[2];
sx q[2];
rz(-1.7804838) q[2];
sx q[2];
rz(1.4022624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0610173) q[1];
sx q[1];
rz(-1.6643545) q[1];
sx q[1];
rz(0.94228014) q[1];
rz(1.9005152) q[3];
sx q[3];
rz(-1.1493249) q[3];
sx q[3];
rz(-0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(2.8395555) q[2];
rz(1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-2.5938477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.1557895) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(0.070080431) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018774059) q[0];
sx q[0];
rz(-0.2304603) q[0];
sx q[0];
rz(-2.3019058) q[0];
rz(-pi) q[1];
rz(2.1224535) q[2];
sx q[2];
rz(-2.2746804) q[2];
sx q[2];
rz(-2.4920419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49680432) q[1];
sx q[1];
rz(-1.4092688) q[1];
sx q[1];
rz(-2.5700388) q[1];
rz(-pi) q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.6453214) q[3];
sx q[3];
rz(2.7217334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(0.62136674) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(-2.5682209) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-0.85987464) q[0];
rz(-1.9372008) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(-3.133657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.313109) q[0];
sx q[0];
rz(-0.61921739) q[0];
sx q[0];
rz(0.45066582) q[0];
x q[1];
rz(-0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(2.0932587) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4850033) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(0.90799241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31422024) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4456711) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(1.3683866) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(2.2369475) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96173441) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(-0.94003135) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(-1.6392802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84425981) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(2.8914333) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0566696) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(-0.67509292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2316206) q[1];
sx q[1];
rz(-2.2044704) q[1];
sx q[1];
rz(2.2909067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7765462) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(2.0765182) q[2];
rz(-0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(2.705943) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57396736) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(-0.91659878) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4531104) q[2];
sx q[2];
rz(-2.3496029) q[2];
sx q[2];
rz(2.288523) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0028487) q[1];
sx q[1];
rz(-1.4321623) q[1];
sx q[1];
rz(1.0671875) q[1];
x q[2];
rz(-0.14108373) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(-1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(0.22582516) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-0.58661714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.3226002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029862558) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(-1.8882621) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0820504) q[2];
sx q[2];
rz(-2.6419123) q[2];
sx q[2];
rz(-0.79364712) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0364914) q[1];
sx q[1];
rz(-1.0579915) q[1];
sx q[1];
rz(-2.3829616) q[1];
rz(1.7694468) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(-2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14324698) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-3.1238363) q[2];
sx q[2];
rz(-2.6323071) q[2];
sx q[2];
rz(1.4321362) q[2];
rz(-3.0512814) q[3];
sx q[3];
rz(-1.3406546) q[3];
sx q[3];
rz(-1.8484074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];