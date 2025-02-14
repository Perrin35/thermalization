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
rz(-0.67686358) q[0];
sx q[0];
rz(-1.6114177) q[0];
sx q[0];
rz(-2.8443008) q[0];
rz(-5.5041432) q[1];
sx q[1];
rz(6.1225523) q[1];
sx q[1];
rz(10.860722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4477618) q[0];
sx q[0];
rz(-1.0658456) q[0];
sx q[0];
rz(2.2044066) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8298684) q[2];
sx q[2];
rz(-1.6951778) q[2];
sx q[2];
rz(2.4476515) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85992009) q[1];
sx q[1];
rz(-1.8919196) q[1];
sx q[1];
rz(2.4008277) q[1];
rz(-pi) q[2];
rz(-0.5663185) q[3];
sx q[3];
rz(-0.36091033) q[3];
sx q[3];
rz(0.64841112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5982738) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(2.9738026) q[2];
rz(1.1281891) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096864916) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(-2.0665533) q[0];
rz(-1.7754405) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(-1.3137821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15720651) q[0];
sx q[0];
rz(-1.3493378) q[0];
sx q[0];
rz(-2.8048672) q[0];
rz(-pi) q[1];
rz(0.73697097) q[2];
sx q[2];
rz(-0.24925512) q[2];
sx q[2];
rz(2.1971306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7915276) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(1.4552874) q[1];
x q[2];
rz(3.0079682) q[3];
sx q[3];
rz(-1.2425955) q[3];
sx q[3];
rz(-1.6673054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65110597) q[2];
sx q[2];
rz(-1.3765114) q[2];
sx q[2];
rz(2.7172078) q[2];
rz(0.73244798) q[3];
sx q[3];
rz(-1.4832352) q[3];
sx q[3];
rz(1.463416) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80377793) q[0];
sx q[0];
rz(-0.71500272) q[0];
sx q[0];
rz(-0.5027698) q[0];
rz(-2.1979507) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(2.3263993) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.316534) q[0];
sx q[0];
rz(-1.7586756) q[0];
sx q[0];
rz(1.5231537) q[0];
rz(-pi) q[1];
rz(-3.0222702) q[2];
sx q[2];
rz(-1.4652958) q[2];
sx q[2];
rz(-1.6245934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5126172) q[1];
sx q[1];
rz(-1.830258) q[1];
sx q[1];
rz(-1.8420868) q[1];
x q[2];
rz(0.92658073) q[3];
sx q[3];
rz(-2.332862) q[3];
sx q[3];
rz(2.1789428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.076685585) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(1.1265075) q[2];
rz(-1.3569776) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(-0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.5821238) q[0];
sx q[0];
rz(-0.72002763) q[0];
sx q[0];
rz(-0.35991392) q[0];
rz(-2.9008046) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(-0.37318939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2853965) q[0];
sx q[0];
rz(-1.4715949) q[0];
sx q[0];
rz(-1.7147934) q[0];
rz(0.56314205) q[2];
sx q[2];
rz(-1.7703319) q[2];
sx q[2];
rz(0.89479252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.14228257) q[1];
sx q[1];
rz(-1.9110084) q[1];
sx q[1];
rz(0.55026023) q[1];
x q[2];
rz(-2.6945962) q[3];
sx q[3];
rz(-1.742954) q[3];
sx q[3];
rz(-0.76154532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.084426247) q[2];
sx q[2];
rz(-1.4320222) q[2];
sx q[2];
rz(-0.90844321) q[2];
rz(2.0716136) q[3];
sx q[3];
rz(-1.9570276) q[3];
sx q[3];
rz(-2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1835566) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(2.2285484) q[0];
rz(-2.4519582) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(-2.3675809) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2406695) q[0];
sx q[0];
rz(-2.0271447) q[0];
sx q[0];
rz(-2.4805446) q[0];
rz(-pi) q[1];
rz(-0.031096259) q[2];
sx q[2];
rz(-2.3457639) q[2];
sx q[2];
rz(2.9253256) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8929114) q[1];
sx q[1];
rz(-1.9982583) q[1];
sx q[1];
rz(-0.21007725) q[1];
x q[2];
rz(-3.0930107) q[3];
sx q[3];
rz(-2.5932199) q[3];
sx q[3];
rz(-2.2480132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.245605) q[2];
sx q[2];
rz(-1.7481123) q[2];
sx q[2];
rz(-2.4020933) q[2];
rz(-1.5064404) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5189811) q[0];
sx q[0];
rz(-2.3464572) q[0];
sx q[0];
rz(2.0163037) q[0];
rz(0.13825026) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(-1.5743871) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.137321) q[0];
sx q[0];
rz(-1.9372276) q[0];
sx q[0];
rz(1.9204813) q[0];
x q[1];
rz(0.60092314) q[2];
sx q[2];
rz(-0.4888566) q[2];
sx q[2];
rz(-0.86684858) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.241227) q[1];
sx q[1];
rz(-1.8292184) q[1];
sx q[1];
rz(1.4421806) q[1];
x q[2];
rz(2.8254824) q[3];
sx q[3];
rz(-2.5694642) q[3];
sx q[3];
rz(2.8993949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3681727) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(-1.0991905) q[2];
rz(2.1558971) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(2.709008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2987591) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(-2.3327935) q[0];
rz(0.66566268) q[1];
sx q[1];
rz(-0.97949615) q[1];
sx q[1];
rz(-2.1601802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.736826) q[0];
sx q[0];
rz(-0.5739218) q[0];
sx q[0];
rz(1.6762275) q[0];
rz(-1.3150042) q[2];
sx q[2];
rz(-1.0793387) q[2];
sx q[2];
rz(2.0602399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1011103) q[1];
sx q[1];
rz(-0.48619929) q[1];
sx q[1];
rz(-1.8644488) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0043169) q[3];
sx q[3];
rz(-2.2386754) q[3];
sx q[3];
rz(1.560451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.064528) q[2];
sx q[2];
rz(-1.4746102) q[2];
sx q[2];
rz(0.36901739) q[2];
rz(-1.4024233) q[3];
sx q[3];
rz(-0.88715059) q[3];
sx q[3];
rz(-1.3212475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8333261) q[0];
sx q[0];
rz(-0.82887355) q[0];
sx q[0];
rz(-2.8114124) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-2.7062682) q[1];
sx q[1];
rz(2.5433069) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.79222) q[0];
sx q[0];
rz(-1.1354007) q[0];
sx q[0];
rz(2.4382082) q[0];
rz(-pi) q[1];
rz(-1.4315375) q[2];
sx q[2];
rz(-0.47310795) q[2];
sx q[2];
rz(2.3654761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51914267) q[1];
sx q[1];
rz(-0.40928091) q[1];
sx q[1];
rz(-2.5342788) q[1];
x q[2];
rz(-2.1586447) q[3];
sx q[3];
rz(-2.1390669) q[3];
sx q[3];
rz(-0.98230386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8699708) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(0.064621933) q[2];
rz(1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(-1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794063) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(-2.620328) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.9797378) q[1];
sx q[1];
rz(2.1102139) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9654664) q[0];
sx q[0];
rz(-0.83399189) q[0];
sx q[0];
rz(1.470572) q[0];
rz(-pi) q[1];
rz(0.29287405) q[2];
sx q[2];
rz(-1.4559064) q[2];
sx q[2];
rz(3.1174768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0887275) q[1];
sx q[1];
rz(-1.2825614) q[1];
sx q[1];
rz(-1.0239787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6691339) q[3];
sx q[3];
rz(-2.0458442) q[3];
sx q[3];
rz(-2.6779384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3970268) q[2];
sx q[2];
rz(-2.4137745) q[2];
sx q[2];
rz(-0.087372027) q[2];
rz(-1.9258063) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-0.89101321) q[0];
sx q[0];
rz(-2.0024306) q[0];
rz(2.6423404) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(-1.3334691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3402417) q[0];
sx q[0];
rz(-1.4805859) q[0];
sx q[0];
rz(-3.0349681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.049801783) q[2];
sx q[2];
rz(-1.8254455) q[2];
sx q[2];
rz(-1.7083502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6302609) q[1];
sx q[1];
rz(-1.3553265) q[1];
sx q[1];
rz(-0.25863077) q[1];
rz(0.74204294) q[3];
sx q[3];
rz(-2.2931406) q[3];
sx q[3];
rz(-0.53398856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6405876) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(-2.8583543) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-0.675942) q[3];
sx q[3];
rz(2.0961608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(0.82861154) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(3.037187) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(-0.5055867) q[2];
sx q[2];
rz(-1.7442295) q[2];
sx q[2];
rz(-2.8019047) q[2];
rz(0.66815175) q[3];
sx q[3];
rz(-1.2378319) q[3];
sx q[3];
rz(1.1938865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
