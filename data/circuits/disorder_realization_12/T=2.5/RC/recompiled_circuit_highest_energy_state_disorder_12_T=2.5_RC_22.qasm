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
rz(-1.9528376) q[0];
sx q[0];
rz(2.3030757) q[0];
sx q[0];
rz(7.6960201) q[0];
rz(0.38263327) q[1];
sx q[1];
rz(-1.9889979) q[1];
sx q[1];
rz(-1.9994073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.45158) q[0];
sx q[0];
rz(-0.40209189) q[0];
sx q[0];
rz(1.7507589) q[0];
rz(1.0866124) q[2];
sx q[2];
rz(-1.118045) q[2];
sx q[2];
rz(0.42175671) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9653017) q[1];
sx q[1];
rz(-1.7994731) q[1];
sx q[1];
rz(-0.92767529) q[1];
rz(-pi) q[2];
rz(-2.9441686) q[3];
sx q[3];
rz(-2.2490778) q[3];
sx q[3];
rz(-1.2048282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6483868) q[2];
sx q[2];
rz(-1.9923261) q[2];
sx q[2];
rz(0.0041848103) q[2];
rz(0.74797136) q[3];
sx q[3];
rz(-0.82704058) q[3];
sx q[3];
rz(0.23676693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0325322) q[0];
sx q[0];
rz(-2.2219658) q[0];
sx q[0];
rz(-0.064706651) q[0];
rz(-2.0789355) q[1];
sx q[1];
rz(-2.2149142) q[1];
sx q[1];
rz(2.0978755) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4499265) q[0];
sx q[0];
rz(-1.1083034) q[0];
sx q[0];
rz(-1.3426379) q[0];
rz(1.4725141) q[2];
sx q[2];
rz(-1.543902) q[2];
sx q[2];
rz(-1.8307476) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1295075) q[1];
sx q[1];
rz(-1.7237067) q[1];
sx q[1];
rz(-0.4985041) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51185913) q[3];
sx q[3];
rz(-2.8843237) q[3];
sx q[3];
rz(0.23990897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26754293) q[2];
sx q[2];
rz(-1.1288613) q[2];
sx q[2];
rz(-1.6740602) q[2];
rz(-2.2360146) q[3];
sx q[3];
rz(-2.0445243) q[3];
sx q[3];
rz(0.21152285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31731376) q[0];
sx q[0];
rz(-0.86238328) q[0];
sx q[0];
rz(-0.30340075) q[0];
rz(-0.41269451) q[1];
sx q[1];
rz(-1.3458601) q[1];
sx q[1];
rz(-2.5491098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6838786) q[0];
sx q[0];
rz(-1.3540097) q[0];
sx q[0];
rz(0.44332645) q[0];
rz(-pi) q[1];
rz(2.8635295) q[2];
sx q[2];
rz(-2.5603271) q[2];
sx q[2];
rz(-0.9424302) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5183181) q[1];
sx q[1];
rz(-0.94097301) q[1];
sx q[1];
rz(1.4131143) q[1];
x q[2];
rz(-2.0140225) q[3];
sx q[3];
rz(-0.89689287) q[3];
sx q[3];
rz(-1.8473089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6233643) q[2];
sx q[2];
rz(-2.6609504) q[2];
sx q[2];
rz(0.19409689) q[2];
rz(2.569681) q[3];
sx q[3];
rz(-1.5827936) q[3];
sx q[3];
rz(-1.5552049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4647575) q[0];
sx q[0];
rz(-2.0763626) q[0];
sx q[0];
rz(1.7179426) q[0];
rz(2.8061197) q[1];
sx q[1];
rz(-0.60233855) q[1];
sx q[1];
rz(0.64340341) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5470872) q[0];
sx q[0];
rz(-2.3427561) q[0];
sx q[0];
rz(-0.79669768) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57482608) q[2];
sx q[2];
rz(-1.4919623) q[2];
sx q[2];
rz(-3.0835033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4130385) q[1];
sx q[1];
rz(-0.87552363) q[1];
sx q[1];
rz(0.9630345) q[1];
x q[2];
rz(-1.2570862) q[3];
sx q[3];
rz(-0.47043741) q[3];
sx q[3];
rz(1.1251118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8947577) q[2];
sx q[2];
rz(-1.094123) q[2];
sx q[2];
rz(2.3168054) q[2];
rz(0.60683933) q[3];
sx q[3];
rz(-1.2068799) q[3];
sx q[3];
rz(1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7428699) q[0];
sx q[0];
rz(-2.5720808) q[0];
sx q[0];
rz(0.2670162) q[0];
rz(0.46785242) q[1];
sx q[1];
rz(-1.1111958) q[1];
sx q[1];
rz(-1.1763447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29821324) q[0];
sx q[0];
rz(-3.0867379) q[0];
sx q[0];
rz(-1.1712267) q[0];
x q[1];
rz(-0.77645923) q[2];
sx q[2];
rz(-2.0909703) q[2];
sx q[2];
rz(2.7732244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5084978) q[1];
sx q[1];
rz(-2.0208686) q[1];
sx q[1];
rz(0.061636713) q[1];
rz(0.72132206) q[3];
sx q[3];
rz(-1.112794) q[3];
sx q[3];
rz(1.7798701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0331369) q[2];
sx q[2];
rz(-1.5507853) q[2];
sx q[2];
rz(0.94672686) q[2];
rz(1.6999792) q[3];
sx q[3];
rz(-2.1082892) q[3];
sx q[3];
rz(-1.3942963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93447584) q[0];
sx q[0];
rz(-1.8659135) q[0];
sx q[0];
rz(-1.7913272) q[0];
rz(2.4840202) q[1];
sx q[1];
rz(-1.5627728) q[1];
sx q[1];
rz(-1.2616166) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3917106) q[0];
sx q[0];
rz(-0.80750033) q[0];
sx q[0];
rz(-0.9901643) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4665537) q[2];
sx q[2];
rz(-1.0770105) q[2];
sx q[2];
rz(1.5774278) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84307901) q[1];
sx q[1];
rz(-1.5061146) q[1];
sx q[1];
rz(-2.9119125) q[1];
rz(-0.17619074) q[3];
sx q[3];
rz(-1.3619845) q[3];
sx q[3];
rz(0.01419078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51800805) q[2];
sx q[2];
rz(-1.1058747) q[2];
sx q[2];
rz(0.27274954) q[2];
rz(-2.2128211) q[3];
sx q[3];
rz(-0.85799587) q[3];
sx q[3];
rz(-1.6509008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51468325) q[0];
sx q[0];
rz(-2.1105483) q[0];
sx q[0];
rz(2.2169901) q[0];
rz(-2.5968016) q[1];
sx q[1];
rz(-1.7926615) q[1];
sx q[1];
rz(-3.0371688) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2685933) q[0];
sx q[0];
rz(-1.6571123) q[0];
sx q[0];
rz(1.6593133) q[0];
rz(2.7265276) q[2];
sx q[2];
rz(-1.7156148) q[2];
sx q[2];
rz(-2.1803149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4778116) q[1];
sx q[1];
rz(-0.38214499) q[1];
sx q[1];
rz(-0.93375979) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35108836) q[3];
sx q[3];
rz(-1.9410053) q[3];
sx q[3];
rz(2.5552487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8203848) q[2];
sx q[2];
rz(-1.9309923) q[2];
sx q[2];
rz(-1.4804117) q[2];
rz(3.1373451) q[3];
sx q[3];
rz(-0.85631266) q[3];
sx q[3];
rz(-0.64928865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55411196) q[0];
sx q[0];
rz(-3.0453747) q[0];
sx q[0];
rz(-1.2531248) q[0];
rz(0.15086497) q[1];
sx q[1];
rz(-0.8395218) q[1];
sx q[1];
rz(1.8419267) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4010581) q[0];
sx q[0];
rz(-0.63176934) q[0];
sx q[0];
rz(-0.44332544) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9076882) q[2];
sx q[2];
rz(-2.680372) q[2];
sx q[2];
rz(0.55800822) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38176259) q[1];
sx q[1];
rz(-2.2564133) q[1];
sx q[1];
rz(-1.5900702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9392936) q[3];
sx q[3];
rz(-1.9194366) q[3];
sx q[3];
rz(1.543247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8396987) q[2];
sx q[2];
rz(-1.7599186) q[2];
sx q[2];
rz(-1.1559486) q[2];
rz(-2.4902952) q[3];
sx q[3];
rz(-2.7498701) q[3];
sx q[3];
rz(2.7294066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0363409) q[0];
sx q[0];
rz(-1.8008494) q[0];
sx q[0];
rz(0.80129188) q[0];
rz(-2.2835412) q[1];
sx q[1];
rz(-0.67887226) q[1];
sx q[1];
rz(2.7616995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7004922) q[0];
sx q[0];
rz(-2.0394562) q[0];
sx q[0];
rz(2.9392713) q[0];
x q[1];
rz(-0.11961898) q[2];
sx q[2];
rz(-1.3966171) q[2];
sx q[2];
rz(-1.2816789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84706789) q[1];
sx q[1];
rz(-0.89133584) q[1];
sx q[1];
rz(0.39107283) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5606468) q[3];
sx q[3];
rz(-2.0119609) q[3];
sx q[3];
rz(0.084189296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6259367) q[2];
sx q[2];
rz(-2.8069324) q[2];
sx q[2];
rz(-0.77896172) q[2];
rz(0.37911478) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(-0.085722119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15445408) q[0];
sx q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(-2.6692303) q[0];
rz(2.8624599) q[1];
sx q[1];
rz(-1.0529073) q[1];
sx q[1];
rz(-2.7739024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147509) q[0];
sx q[0];
rz(-1.923133) q[0];
sx q[0];
rz(-2.8295838) q[0];
rz(-pi) q[1];
rz(-1.3490178) q[2];
sx q[2];
rz(-1.1797172) q[2];
sx q[2];
rz(-1.329601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4039129) q[1];
sx q[1];
rz(-1.332799) q[1];
sx q[1];
rz(-0.42003553) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2982453) q[3];
sx q[3];
rz(-0.9002004) q[3];
sx q[3];
rz(-1.6141817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8184066) q[2];
sx q[2];
rz(-1.0419934) q[2];
sx q[2];
rz(2.6673178) q[2];
rz(-2.2222774) q[3];
sx q[3];
rz(-1.5654469) q[3];
sx q[3];
rz(-1.3924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267088) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(-0.23450163) q[1];
sx q[1];
rz(-2.0461743) q[1];
sx q[1];
rz(2.8428427) q[1];
rz(-0.089800553) q[2];
sx q[2];
rz(-1.0013781) q[2];
sx q[2];
rz(2.7120874) q[2];
rz(-1.6203703) q[3];
sx q[3];
rz(-1.5885316) q[3];
sx q[3];
rz(1.9902609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
