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
rz(2.3447073) q[0];
sx q[0];
rz(-1.7303884) q[0];
sx q[0];
rz(2.222173) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(3.0923831) q[1];
sx q[1];
rz(10.13848) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714234) q[0];
sx q[0];
rz(-1.7180301) q[0];
sx q[0];
rz(1.2933613) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4422779) q[2];
sx q[2];
rz(-1.6331731) q[2];
sx q[2];
rz(2.2888052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9992396) q[1];
sx q[1];
rz(-2.0217784) q[1];
sx q[1];
rz(-1.5723448) q[1];
rz(0.034293745) q[3];
sx q[3];
rz(-2.0418797) q[3];
sx q[3];
rz(0.33609875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0539703) q[2];
sx q[2];
rz(-0.38603187) q[2];
sx q[2];
rz(2.7719882) q[2];
rz(1.9965648) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(-2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10577781) q[0];
sx q[0];
rz(-1.4666297) q[0];
sx q[0];
rz(0.57304397) q[0];
rz(2.0447958) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(-2.6699452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99804634) q[0];
sx q[0];
rz(-0.61206383) q[0];
sx q[0];
rz(1.7630811) q[0];
rz(-0.71895968) q[2];
sx q[2];
rz(-2.3078354) q[2];
sx q[2];
rz(-1.3477088) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0309033) q[1];
sx q[1];
rz(-2.2784) q[1];
sx q[1];
rz(-0.75548197) q[1];
rz(-pi) q[2];
rz(0.40741323) q[3];
sx q[3];
rz(-1.0828567) q[3];
sx q[3];
rz(0.98516243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(-2.0527077) q[2];
rz(-2.0770843) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(0.3256807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.0786781) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(-0.56627083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1120934) q[0];
sx q[0];
rz(-1.5595961) q[0];
sx q[0];
rz(2.311934) q[0];
x q[1];
rz(-0.60795075) q[2];
sx q[2];
rz(-0.81038108) q[2];
sx q[2];
rz(-1.2255526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6614187) q[1];
sx q[1];
rz(-1.4858477) q[1];
sx q[1];
rz(-0.18710356) q[1];
x q[2];
rz(-1.503946) q[3];
sx q[3];
rz(-1.306157) q[3];
sx q[3];
rz(2.8422249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8038586) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(-0.0090948661) q[2];
rz(-3.005262) q[3];
sx q[3];
rz(-2.3970042) q[3];
sx q[3];
rz(-1.2135308) q[3];
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
rz(0.71753865) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(-2.9754382) q[0];
rz(1.1135788) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(-0.16214935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0320764) q[0];
sx q[0];
rz(-1.6639532) q[0];
sx q[0];
rz(-1.7514125) q[0];
x q[1];
rz(-3.1012898) q[2];
sx q[2];
rz(-2.5922814) q[2];
sx q[2];
rz(1.2072762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0905793) q[1];
sx q[1];
rz(-2.3115578) q[1];
sx q[1];
rz(0.49511893) q[1];
rz(-1.1433853) q[3];
sx q[3];
rz(-2.1899438) q[3];
sx q[3];
rz(-2.038508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1589511) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(-0.52784935) q[2];
rz(-2.2900901) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-0.270179) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(-2.340509) q[0];
rz(-0.78394765) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(0.98091006) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1005412) q[0];
sx q[0];
rz(-1.1642191) q[0];
sx q[0];
rz(1.340614) q[0];
rz(-2.0774242) q[2];
sx q[2];
rz(-1.9031823) q[2];
sx q[2];
rz(-1.6695963) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0438761) q[1];
sx q[1];
rz(-1.4773158) q[1];
sx q[1];
rz(3.082117) q[1];
rz(1.0193056) q[3];
sx q[3];
rz(-1.3634342) q[3];
sx q[3];
rz(-0.99121782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(2.2583466) q[2];
rz(-0.20032459) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(-1.0909572) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9969295) q[0];
sx q[0];
rz(-2.0893593) q[0];
sx q[0];
rz(2.1494179) q[0];
rz(-2.3400173) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(-1.7209524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039889) q[0];
sx q[0];
rz(-2.1891356) q[0];
sx q[0];
rz(1.9644587) q[0];
rz(-1.0384485) q[2];
sx q[2];
rz(-1.4803998) q[2];
sx q[2];
rz(2.9351075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3657377) q[1];
sx q[1];
rz(-2.4147662) q[1];
sx q[1];
rz(-0.50623843) q[1];
rz(-pi) q[2];
rz(1.4735953) q[3];
sx q[3];
rz(-1.6080813) q[3];
sx q[3];
rz(0.71747045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8617323) q[2];
sx q[2];
rz(-1.7143152) q[2];
sx q[2];
rz(-0.53708616) q[2];
rz(-0.01072695) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.895973) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(1.8396359) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(1.1713015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7599277) q[0];
sx q[0];
rz(-1.6070843) q[0];
sx q[0];
rz(0.29706232) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76979678) q[2];
sx q[2];
rz(-0.30210051) q[2];
sx q[2];
rz(0.090902791) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.724546) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(1.6032277) q[1];
rz(1.7460489) q[3];
sx q[3];
rz(-0.95357663) q[3];
sx q[3];
rz(-0.72886906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12477144) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-2.4701414) q[2];
rz(-2.6050383) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(-2.0898537) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-2.0874513) q[0];
sx q[0];
rz(-0.89624727) q[0];
rz(-0.25262901) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(-1.0505229) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4481215) q[0];
sx q[0];
rz(-2.6127609) q[0];
sx q[0];
rz(0.099683925) q[0];
rz(-2.975198) q[2];
sx q[2];
rz(-0.67182589) q[2];
sx q[2];
rz(-2.2937798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.412789) q[1];
sx q[1];
rz(-1.1096242) q[1];
sx q[1];
rz(0.16766602) q[1];
rz(-pi) q[2];
x q[2];
rz(1.518431) q[3];
sx q[3];
rz(-0.81365055) q[3];
sx q[3];
rz(2.6770075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5184021) q[2];
sx q[2];
rz(-2.9727327) q[2];
sx q[2];
rz(-1.5625578) q[2];
rz(1.7744428) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(0.76213837) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(0.38994625) q[0];
rz(-2.3155616) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(-2.0894076) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1601801) q[0];
sx q[0];
rz(-2.3773068) q[0];
sx q[0];
rz(2.0212964) q[0];
x q[1];
rz(0.39689831) q[2];
sx q[2];
rz(-1.7328615) q[2];
sx q[2];
rz(-1.6610749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3996399) q[1];
sx q[1];
rz(-1.9583079) q[1];
sx q[1];
rz(-0.31419973) q[1];
x q[2];
rz(2.1972606) q[3];
sx q[3];
rz(-1.9667224) q[3];
sx q[3];
rz(2.348071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0612001) q[2];
sx q[2];
rz(-1.281851) q[2];
sx q[2];
rz(0.52919394) q[2];
rz(-2.3906294) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(-2.9474337) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879409) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(1.1645114) q[0];
rz(1.8544082) q[1];
sx q[1];
rz(-2.4559805) q[1];
sx q[1];
rz(0.50416344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8541478) q[0];
sx q[0];
rz(-1.3519671) q[0];
sx q[0];
rz(1.8356505) q[0];
rz(-1.590056) q[2];
sx q[2];
rz(-1.1941205) q[2];
sx q[2];
rz(-2.6835359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2626942) q[1];
sx q[1];
rz(-2.8126908) q[1];
sx q[1];
rz(0.94543381) q[1];
rz(-pi) q[2];
rz(-2.9243117) q[3];
sx q[3];
rz(-0.8978399) q[3];
sx q[3];
rz(2.6126044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9265499) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(-0.034493383) q[2];
rz(-1.6132332) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(0.31149402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16019776) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(1.002671) q[1];
sx q[1];
rz(-1.4613338) q[1];
sx q[1];
rz(2.7636539) q[1];
rz(0.45838452) q[2];
sx q[2];
rz(-1.3946563) q[2];
sx q[2];
rz(-1.3763225) q[2];
rz(1.6335842) q[3];
sx q[3];
rz(-1.8973402) q[3];
sx q[3];
rz(-0.8473581) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
