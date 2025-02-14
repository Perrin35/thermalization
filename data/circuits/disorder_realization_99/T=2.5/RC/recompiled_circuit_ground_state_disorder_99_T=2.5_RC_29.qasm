OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34243256) q[0];
sx q[0];
rz(2.7437796) q[0];
sx q[0];
rz(7.5510511) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(2.3550912) q[1];
sx q[1];
rz(13.785706) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.08377) q[0];
sx q[0];
rz(-1.738213) q[0];
sx q[0];
rz(0.98982248) q[0];
rz(-pi) q[1];
rz(-1.7389347) q[2];
sx q[2];
rz(-1.0285707) q[2];
sx q[2];
rz(-0.10054345) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3883512) q[1];
sx q[1];
rz(-1.2920391) q[1];
sx q[1];
rz(0.21215794) q[1];
rz(-pi) q[2];
rz(2.9079535) q[3];
sx q[3];
rz(-0.92687449) q[3];
sx q[3];
rz(2.0104017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5775602) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.7731898) q[2];
rz(-2.2300301) q[3];
sx q[3];
rz(-1.7715958) q[3];
sx q[3];
rz(3.1172359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1034705) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(-0.11904112) q[0];
rz(1.1301522) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(1.4281323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85066909) q[0];
sx q[0];
rz(-2.4570298) q[0];
sx q[0];
rz(-0.67770358) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4006673) q[2];
sx q[2];
rz(-1.180449) q[2];
sx q[2];
rz(-2.0290749) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2547639) q[1];
sx q[1];
rz(-0.48476754) q[1];
sx q[1];
rz(-1.280275) q[1];
x q[2];
rz(2.8359706) q[3];
sx q[3];
rz(-1.247331) q[3];
sx q[3];
rz(-2.4572069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52593645) q[2];
sx q[2];
rz(-1.6776513) q[2];
sx q[2];
rz(-2.3770135) q[2];
rz(2.6584451) q[3];
sx q[3];
rz(-1.316148) q[3];
sx q[3];
rz(0.84883261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371542) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(1.6695439) q[0];
rz(0.56023359) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(1.4405506) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.154523) q[0];
sx q[0];
rz(-1.8049585) q[0];
sx q[0];
rz(-0.20407853) q[0];
rz(-pi) q[1];
rz(2.4666489) q[2];
sx q[2];
rz(-1.5663354) q[2];
sx q[2];
rz(-2.5738112) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.24748188) q[1];
sx q[1];
rz(-1.2318582) q[1];
sx q[1];
rz(2.1268658) q[1];
rz(-pi) q[2];
rz(2.3769242) q[3];
sx q[3];
rz(-2.8199785) q[3];
sx q[3];
rz(-2.4177616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(2.7991925) q[2];
rz(-1.418669) q[3];
sx q[3];
rz(-0.72903052) q[3];
sx q[3];
rz(-2.5820144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320084) q[0];
sx q[0];
rz(-2.7689458) q[0];
sx q[0];
rz(1.5268071) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5680983) q[1];
sx q[1];
rz(1.2015013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1781834) q[0];
sx q[0];
rz(-0.82287517) q[0];
sx q[0];
rz(1.9848518) q[0];
x q[1];
rz(-2.2204578) q[2];
sx q[2];
rz(-2.1138722) q[2];
sx q[2];
rz(1.7004101) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0240942) q[1];
sx q[1];
rz(-1.1765685) q[1];
sx q[1];
rz(0.82347639) q[1];
x q[2];
rz(-1.6225927) q[3];
sx q[3];
rz(-1.3148309) q[3];
sx q[3];
rz(-1.3802176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3379007) q[2];
sx q[2];
rz(-2.1521229) q[2];
sx q[2];
rz(1.2376415) q[2];
rz(1.8303653) q[3];
sx q[3];
rz(-1.5179736) q[3];
sx q[3];
rz(-1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(-0.9088687) q[0];
rz(2.2591649) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(-0.73208255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013618795) q[0];
sx q[0];
rz(-2.2165856) q[0];
sx q[0];
rz(0.60198204) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3734748) q[2];
sx q[2];
rz(-1.3663732) q[2];
sx q[2];
rz(-0.26917514) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4849347) q[1];
sx q[1];
rz(-1.1371433) q[1];
sx q[1];
rz(-2.0051433) q[1];
x q[2];
rz(-2.3957) q[3];
sx q[3];
rz(-0.98222662) q[3];
sx q[3];
rz(2.4323127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0714134) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(1.2893691) q[2];
rz(1.6728801) q[3];
sx q[3];
rz(-1.9332935) q[3];
sx q[3];
rz(2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656723) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(-2.0400203) q[0];
rz(0.22483243) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(0.98519957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3168047) q[0];
sx q[0];
rz(-1.390762) q[0];
sx q[0];
rz(2.2145674) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2217789) q[2];
sx q[2];
rz(-2.1563081) q[2];
sx q[2];
rz(1.3427918) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0538865) q[1];
sx q[1];
rz(-0.83982491) q[1];
sx q[1];
rz(2.5165416) q[1];
x q[2];
rz(1.5456301) q[3];
sx q[3];
rz(-1.0428793) q[3];
sx q[3];
rz(0.86564579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3605986) q[2];
sx q[2];
rz(-0.66086078) q[2];
sx q[2];
rz(1.194427) q[2];
rz(-3.0418975) q[3];
sx q[3];
rz(-1.3705148) q[3];
sx q[3];
rz(2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7128971) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(2.0275443) q[0];
rz(1.3975337) q[1];
sx q[1];
rz(-1.7618529) q[1];
sx q[1];
rz(-0.31375113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0359441) q[0];
sx q[0];
rz(-1.2782845) q[0];
sx q[0];
rz(-3.0894795) q[0];
x q[1];
rz(-2.9767562) q[2];
sx q[2];
rz(-1.4408752) q[2];
sx q[2];
rz(0.31229737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0278247) q[1];
sx q[1];
rz(-1.271444) q[1];
sx q[1];
rz(-2.37948) q[1];
rz(-2.4558732) q[3];
sx q[3];
rz(-1.7758533) q[3];
sx q[3];
rz(3.0577537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.067001192) q[2];
sx q[2];
rz(-0.40412298) q[2];
sx q[2];
rz(-0.60603777) q[2];
rz(-1.0780942) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(1.3585453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17230497) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(1.1267927) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-0.4387478) q[1];
sx q[1];
rz(2.9235358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0261966) q[0];
sx q[0];
rz(-0.65803981) q[0];
sx q[0];
rz(-0.81320073) q[0];
rz(2.3534691) q[2];
sx q[2];
rz(-1.0255073) q[2];
sx q[2];
rz(2.9719549) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3266284) q[1];
sx q[1];
rz(-1.4451298) q[1];
sx q[1];
rz(2.2180422) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0239437) q[3];
sx q[3];
rz(-1.434064) q[3];
sx q[3];
rz(-2.4916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-0.67535526) q[2];
sx q[2];
rz(-2.8524354) q[2];
rz(-2.1469927) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(0.4579671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7116123) q[0];
sx q[0];
rz(-2.0531605) q[0];
sx q[0];
rz(0.76643884) q[0];
rz(-1.6038766) q[1];
sx q[1];
rz(-1.3130554) q[1];
sx q[1];
rz(-2.2235353) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72836711) q[0];
sx q[0];
rz(-2.0371975) q[0];
sx q[0];
rz(-0.8829929) q[0];
x q[1];
rz(-1.7154022) q[2];
sx q[2];
rz(-1.7376126) q[2];
sx q[2];
rz(1.5125991) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.854771) q[1];
sx q[1];
rz(-2.3769551) q[1];
sx q[1];
rz(2.7809793) q[1];
rz(-pi) q[2];
rz(1.9647801) q[3];
sx q[3];
rz(-0.26991649) q[3];
sx q[3];
rz(-0.22701015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90091577) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(0.051699836) q[2];
rz(-2.2895571) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(-0.26028546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32605115) q[0];
sx q[0];
rz(-1.6764078) q[0];
sx q[0];
rz(3.0503804) q[0];
rz(-0.7971898) q[1];
sx q[1];
rz(-1.7557095) q[1];
sx q[1];
rz(-0.43112722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.5231966) q[0];
sx q[0];
rz(-3.0831227) q[0];
x q[1];
rz(-3.0431137) q[2];
sx q[2];
rz(-0.92442552) q[2];
sx q[2];
rz(-1.3502094) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7550274) q[1];
sx q[1];
rz(-1.7885099) q[1];
sx q[1];
rz(-0.45611195) q[1];
rz(-pi) q[2];
rz(2.71118) q[3];
sx q[3];
rz(-2.1280043) q[3];
sx q[3];
rz(-0.9924953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.572523) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(1.265906) q[2];
rz(-1.8484533) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(-2.7354447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010392808) q[0];
sx q[0];
rz(-2.4507903) q[0];
sx q[0];
rz(-2.3660085) q[0];
rz(-2.5776183) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(1.3971267) q[2];
sx q[2];
rz(-0.43890719) q[2];
sx q[2];
rz(-0.70499805) q[2];
rz(1.0913783) q[3];
sx q[3];
rz(-1.2751083) q[3];
sx q[3];
rz(0.90838065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
