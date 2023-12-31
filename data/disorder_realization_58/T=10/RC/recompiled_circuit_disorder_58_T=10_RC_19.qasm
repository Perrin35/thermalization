OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(-2.3109205) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75738534) q[0];
sx q[0];
rz(-1.877458) q[0];
sx q[0];
rz(0.39461179) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0482424) q[2];
sx q[2];
rz(-2.2508143) q[2];
sx q[2];
rz(-0.29104656) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4315011) q[1];
sx q[1];
rz(-0.7845062) q[1];
sx q[1];
rz(2.0246519) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44736638) q[3];
sx q[3];
rz(-1.9379741) q[3];
sx q[3];
rz(0.51608738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.1581356) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-2.2297915) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(0.3266913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625898) q[0];
sx q[0];
rz(-2.0741182) q[0];
sx q[0];
rz(-2.3110564) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2454883) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(-0.59782366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3152299) q[1];
sx q[1];
rz(-0.76347199) q[1];
sx q[1];
rz(2.1748494) q[1];
x q[2];
rz(1.6750402) q[3];
sx q[3];
rz(-0.6904656) q[3];
sx q[3];
rz(0.0052099293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54863769) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(1.057391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1304504) q[0];
sx q[0];
rz(-1.5913977) q[0];
sx q[0];
rz(1.2304473) q[0];
x q[1];
rz(-0.51809394) q[2];
sx q[2];
rz(-0.38803852) q[2];
sx q[2];
rz(-2.5990017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2705546) q[1];
sx q[1];
rz(-1.241031) q[1];
sx q[1];
rz(2.2689181) q[1];
rz(-0.4226513) q[3];
sx q[3];
rz(-1.7461516) q[3];
sx q[3];
rz(2.5734176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574922) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388158) q[0];
sx q[0];
rz(-1.446615) q[0];
sx q[0];
rz(1.3846272) q[0];
rz(2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(1.4571112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0488102) q[1];
sx q[1];
rz(-2.4895992) q[1];
sx q[1];
rz(-2.584143) q[1];
rz(-pi) q[2];
rz(-2.1256251) q[3];
sx q[3];
rz(-1.3362243) q[3];
sx q[3];
rz(-2.6382584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(2.0139587) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(2.7868295) q[0];
rz(-1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-2.9096471) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1629099) q[0];
sx q[0];
rz(-1.9918348) q[0];
sx q[0];
rz(-2.1091503) q[0];
rz(1.4300214) q[2];
sx q[2];
rz(-0.88015926) q[2];
sx q[2];
rz(1.8292793) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2105337) q[1];
sx q[1];
rz(-1.7909044) q[1];
sx q[1];
rz(-1.348043) q[1];
x q[2];
rz(-2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(-0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.140124) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(1.9011964) q[2];
rz(2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(0.20275673) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(0.99745497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049012262) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(2.9655365) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73156725) q[2];
sx q[2];
rz(-1.8533857) q[2];
sx q[2];
rz(2.3716795) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.62085405) q[1];
sx q[1];
rz(-1.3315017) q[1];
sx q[1];
rz(-2.005307) q[1];
x q[2];
rz(-2.698425) q[3];
sx q[3];
rz(-0.89649761) q[3];
sx q[3];
rz(-0.5815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(-2.690199) q[2];
rz(-2.732892) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-0.61378941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95425883) q[0];
sx q[0];
rz(-1.3849392) q[0];
sx q[0];
rz(0.82877393) q[0];
rz(1.5067528) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(1.7664906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2144199) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(0.10417948) q[1];
x q[2];
rz(-1.1055787) q[3];
sx q[3];
rz(-1.4659766) q[3];
sx q[3];
rz(3.1302111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.9667352) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.90010086) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(-0.98446313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41868608) q[0];
sx q[0];
rz(-1.6081928) q[0];
sx q[0];
rz(0.54625578) q[0];
rz(2.369957) q[2];
sx q[2];
rz(-1.8688335) q[2];
sx q[2];
rz(0.76921295) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7998432) q[1];
sx q[1];
rz(-1.5647596) q[1];
sx q[1];
rz(-0.36532613) q[1];
rz(-pi) q[2];
rz(-0.55926178) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-0.53517503) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5110916) q[0];
sx q[0];
rz(-1.4303659) q[0];
sx q[0];
rz(1.3569843) q[0];
x q[1];
rz(2.4326153) q[2];
sx q[2];
rz(-0.93009863) q[2];
sx q[2];
rz(2.7880653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3013819) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(-2.0135897) q[1];
rz(-pi) q[2];
rz(2.9374398) q[3];
sx q[3];
rz(-1.6870058) q[3];
sx q[3];
rz(-1.0128563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(2.6055028) q[2];
rz(0.4195956) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-2.1283456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17908827) q[0];
sx q[0];
rz(-2.8694186) q[0];
sx q[0];
rz(-2.1129235) q[0];
rz(1.6147862) q[2];
sx q[2];
rz(-1.9242052) q[2];
sx q[2];
rz(0.24214889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71868616) q[1];
sx q[1];
rz(-1.6675073) q[1];
sx q[1];
rz(2.7194517) q[1];
x q[2];
rz(2.4471531) q[3];
sx q[3];
rz(-2.246292) q[3];
sx q[3];
rz(1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(-1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(-0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(2.5683174) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(3.0782386) q[3];
sx q[3];
rz(-2.2676716) q[3];
sx q[3];
rz(0.36361658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
