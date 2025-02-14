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
rz(2.4647291) q[0];
sx q[0];
rz(-1.530175) q[0];
sx q[0];
rz(-0.29729182) q[0];
rz(0.77904207) q[1];
sx q[1];
rz(-0.160633) q[1];
sx q[1];
rz(1.4359441) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53548649) q[0];
sx q[0];
rz(-1.0260884) q[0];
sx q[0];
rz(-2.5404055) q[0];
x q[1];
rz(-3.0129635) q[2];
sx q[2];
rz(-1.3137711) q[2];
sx q[2];
rz(2.2976053) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85992009) q[1];
sx q[1];
rz(-1.8919196) q[1];
sx q[1];
rz(-2.4008277) q[1];
rz(-0.5663185) q[3];
sx q[3];
rz(-2.7806823) q[3];
sx q[3];
rz(2.4931815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5433189) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(-0.16779009) q[2];
rz(-1.1281891) q[3];
sx q[3];
rz(-1.5725243) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
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
rz(-0.096864916) q[0];
sx q[0];
rz(-0.6701349) q[0];
sx q[0];
rz(1.0750394) q[0];
rz(1.3661522) q[1];
sx q[1];
rz(-1.7551883) q[1];
sx q[1];
rz(-1.8278106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2883817) q[0];
sx q[0];
rz(-2.7409005) q[0];
sx q[0];
rz(0.59817065) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4046217) q[2];
sx q[2];
rz(-2.8923375) q[2];
sx q[2];
rz(-2.1971306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.19590727) q[1];
sx q[1];
rz(-0.84896174) q[1];
sx q[1];
rz(3.0392749) q[1];
rz(-pi) q[2];
rz(1.9017392) q[3];
sx q[3];
rz(-1.4443436) q[3];
sx q[3];
rz(3.0883873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4904867) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(-2.7172078) q[2];
rz(0.73244798) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80377793) q[0];
sx q[0];
rz(-0.71500272) q[0];
sx q[0];
rz(-0.5027698) q[0];
rz(-0.94364199) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(0.81519333) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.40476) q[0];
sx q[0];
rz(-1.6175999) q[0];
sx q[0];
rz(2.9535049) q[0];
x q[1];
rz(0.11932245) q[2];
sx q[2];
rz(-1.4652958) q[2];
sx q[2];
rz(-1.6245934) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.80322641) q[1];
sx q[1];
rz(-0.37316457) q[1];
sx q[1];
rz(2.3514521) q[1];
rz(-pi) q[2];
x q[2];
rz(2.268154) q[3];
sx q[3];
rz(-2.0202352) q[3];
sx q[3];
rz(3.0118503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0649071) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(-1.1265075) q[2];
rz(-1.7846151) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(-3.0016541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55946881) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(-2.7816787) q[0];
rz(0.24078807) q[1];
sx q[1];
rz(-1.1477995) q[1];
sx q[1];
rz(-2.7684033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68615896) q[0];
sx q[0];
rz(-0.17466521) q[0];
sx q[0];
rz(0.96439488) q[0];
x q[1];
rz(2.5784506) q[2];
sx q[2];
rz(-1.7703319) q[2];
sx q[2];
rz(-0.89479252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2265985) q[1];
sx q[1];
rz(-2.0862596) q[1];
sx q[1];
rz(-1.1771918) q[1];
x q[2];
rz(-0.38244989) q[3];
sx q[3];
rz(-2.6646864) q[3];
sx q[3];
rz(-0.46602369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0571664) q[2];
sx q[2];
rz(-1.4320222) q[2];
sx q[2];
rz(2.2331494) q[2];
rz(-2.0716136) q[3];
sx q[3];
rz(-1.9570276) q[3];
sx q[3];
rz(2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1835566) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(-0.68963447) q[1];
sx q[1];
rz(-0.90877405) q[1];
sx q[1];
rz(0.77401179) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.000074) q[0];
sx q[0];
rz(-2.1545) q[0];
sx q[0];
rz(2.1271749) q[0];
rz(0.795587) q[2];
sx q[2];
rz(-1.5485816) q[2];
sx q[2];
rz(-1.8088248) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8929114) q[1];
sx q[1];
rz(-1.1433344) q[1];
sx q[1];
rz(-2.9315154) q[1];
x q[2];
rz(1.5411395) q[3];
sx q[3];
rz(-1.0231442) q[3];
sx q[3];
rz(-0.95049196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.245605) q[2];
sx q[2];
rz(-1.7481123) q[2];
sx q[2];
rz(-0.73949933) q[2];
rz(1.5064404) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(-0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5189811) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(-3.0033424) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(-1.5743871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7985021) q[0];
sx q[0];
rz(-2.6406086) q[0];
sx q[0];
rz(2.4128017) q[0];
x q[1];
rz(0.60092314) q[2];
sx q[2];
rz(-2.6527361) q[2];
sx q[2];
rz(-2.2747441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7789845) q[1];
sx q[1];
rz(-1.6951188) q[1];
sx q[1];
rz(2.881114) q[1];
rz(-2.592347) q[3];
sx q[3];
rz(-1.4016782) q[3];
sx q[3];
rz(-2.0813519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77341998) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(-1.0991905) q[2];
rz(-0.98569551) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(2.709008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2987591) q[0];
sx q[0];
rz(-0.99501139) q[0];
sx q[0];
rz(2.3327935) q[0];
rz(0.66566268) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(2.1601802) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.736826) q[0];
sx q[0];
rz(-0.5739218) q[0];
sx q[0];
rz(1.6762275) q[0];
rz(-0.44156011) q[2];
sx q[2];
rz(-0.54916635) q[2];
sx q[2];
rz(1.5541981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3700404) q[1];
sx q[1];
rz(-2.0345033) q[1];
sx q[1];
rz(0.15180219) q[1];
rz(-0.1372758) q[3];
sx q[3];
rz(-2.2386754) q[3];
sx q[3];
rz(-1.5811416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.064528) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(0.36901739) q[2];
rz(1.7391694) q[3];
sx q[3];
rz(-0.88715059) q[3];
sx q[3];
rz(1.8203452) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8333261) q[0];
sx q[0];
rz(-0.82887355) q[0];
sx q[0];
rz(-2.8114124) q[0];
rz(-2.6847243) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(-2.5433069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.79222) q[0];
sx q[0];
rz(-2.0061919) q[0];
sx q[0];
rz(2.4382082) q[0];
x q[1];
rz(2.0399698) q[2];
sx q[2];
rz(-1.6340877) q[2];
sx q[2];
rz(2.4710412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48409325) q[1];
sx q[1];
rz(-1.7998905) q[1];
sx q[1];
rz(0.34219663) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4264494) q[3];
sx q[3];
rz(-2.3480881) q[3];
sx q[3];
rz(1.8737829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2716219) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(0.064621933) q[2];
rz(-1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(-1.5456642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794063) q[0];
sx q[0];
rz(-0.33410826) q[0];
sx q[0];
rz(-0.52126467) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.9797378) q[1];
sx q[1];
rz(-1.0313787) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17612621) q[0];
sx q[0];
rz(-2.3076008) q[0];
sx q[0];
rz(1.470572) q[0];
rz(-2.7613369) q[2];
sx q[2];
rz(-0.31399841) q[2];
sx q[2];
rz(-1.2316201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0528652) q[1];
sx q[1];
rz(-1.2825614) q[1];
sx q[1];
rz(2.117614) q[1];
rz(-pi) q[2];
rz(1.6691339) q[3];
sx q[3];
rz(-2.0458442) q[3];
sx q[3];
rz(-0.46365427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3970268) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(-3.0542206) q[2];
rz(-1.2157863) q[3];
sx q[3];
rz(-1.6814503) q[3];
sx q[3];
rz(0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66607296) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(1.1391621) q[0];
rz(-2.6423404) q[1];
sx q[1];
rz(-1.4491932) q[1];
sx q[1];
rz(1.8081236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0718057) q[0];
sx q[0];
rz(-3.0020368) q[0];
sx q[0];
rz(-0.7044756) q[0];
rz(-pi) q[1];
rz(0.049801783) q[2];
sx q[2];
rz(-1.8254455) q[2];
sx q[2];
rz(1.7083502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4024248) q[1];
sx q[1];
rz(-0.3350733) q[1];
sx q[1];
rz(-0.70783028) q[1];
rz(-pi) q[2];
rz(2.2250149) q[3];
sx q[3];
rz(-0.98482705) q[3];
sx q[3];
rz(0.41205349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6405876) q[2];
sx q[2];
rz(-1.2747526) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129811) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(0.10440566) q[1];
sx q[1];
rz(-1.5900292) q[1];
sx q[1];
rz(1.1463696) q[1];
rz(0.5055867) q[2];
sx q[2];
rz(-1.3973631) q[2];
sx q[2];
rz(0.33968795) q[2];
rz(-0.5091359) q[3];
sx q[3];
rz(-2.4066299) q[3];
sx q[3];
rz(0.015711333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
