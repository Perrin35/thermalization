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
rz(-1.5440829) q[0];
sx q[0];
rz(-1.7680327) q[0];
sx q[0];
rz(-1.0184259) q[0];
rz(-0.51813689) q[1];
sx q[1];
rz(-2.3845446) q[1];
sx q[1];
rz(0.63176027) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951183) q[0];
sx q[0];
rz(-2.8993239) q[0];
sx q[0];
rz(2.1184068) q[0];
x q[1];
rz(-1.873513) q[2];
sx q[2];
rz(-1.5470501) q[2];
sx q[2];
rz(2.077092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.077191) q[1];
sx q[1];
rz(-2.4418412) q[1];
sx q[1];
rz(-2.0136858) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8381836) q[3];
sx q[3];
rz(-2.8041556) q[3];
sx q[3];
rz(-0.22663675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53039256) q[2];
sx q[2];
rz(-0.35553122) q[2];
sx q[2];
rz(1.4872888) q[2];
rz(-0.90855956) q[3];
sx q[3];
rz(-2.9049951) q[3];
sx q[3];
rz(-1.0049741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0205883) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(-2.9793136) q[0];
rz(2.8970215) q[1];
sx q[1];
rz(-1.9385612) q[1];
sx q[1];
rz(-0.29168209) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8957386) q[0];
sx q[0];
rz(-1.7632005) q[0];
sx q[0];
rz(-0.28783513) q[0];
rz(-1.8542833) q[2];
sx q[2];
rz(-1.3559196) q[2];
sx q[2];
rz(2.7708997) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6704935) q[1];
sx q[1];
rz(-0.80599497) q[1];
sx q[1];
rz(-0.32260311) q[1];
rz(-0.79370768) q[3];
sx q[3];
rz(-2.5849403) q[3];
sx q[3];
rz(-1.2964378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67819277) q[2];
sx q[2];
rz(-2.2726111) q[2];
sx q[2];
rz(-1.0758859) q[2];
rz(1.2470657) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4397044) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(-0.18371789) q[0];
rz(1.6366929) q[1];
sx q[1];
rz(-2.3972062) q[1];
sx q[1];
rz(-2.9768129) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713584) q[0];
sx q[0];
rz(-0.58588282) q[0];
sx q[0];
rz(-1.5536932) q[0];
rz(-pi) q[1];
rz(1.25827) q[2];
sx q[2];
rz(-2.222568) q[2];
sx q[2];
rz(1.5628017) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0244813) q[1];
sx q[1];
rz(-2.4035807) q[1];
sx q[1];
rz(2.5548425) q[1];
rz(-1.0298877) q[3];
sx q[3];
rz(-0.60324429) q[3];
sx q[3];
rz(-1.3689976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.160459) q[2];
sx q[2];
rz(-2.0351694) q[2];
rz(-2.4404081) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(-2.6760694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3727386) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(-0.94451529) q[0];
rz(-1.6230029) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(1.5400344) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9732653) q[0];
sx q[0];
rz(-2.1337318) q[0];
sx q[0];
rz(2.299978) q[0];
rz(-0.77727274) q[2];
sx q[2];
rz(-0.16675719) q[2];
sx q[2];
rz(2.6646032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9303927) q[1];
sx q[1];
rz(-0.89387776) q[1];
sx q[1];
rz(-2.6571214) q[1];
rz(0.047961162) q[3];
sx q[3];
rz(-0.75100079) q[3];
sx q[3];
rz(-0.50868552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.010178415) q[2];
sx q[2];
rz(-2.5203036) q[2];
sx q[2];
rz(2.9739001) q[2];
rz(3.0883279) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(1.3767327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662956) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(0.29933023) q[0];
rz(0.37295595) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(-1.6711055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9540967) q[0];
sx q[0];
rz(-2.4135655) q[0];
sx q[0];
rz(-2.4690829) q[0];
rz(-pi) q[1];
rz(0.46385455) q[2];
sx q[2];
rz(-1.6783444) q[2];
sx q[2];
rz(0.35075089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3186444) q[1];
sx q[1];
rz(-0.80893436) q[1];
sx q[1];
rz(-0.72973324) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3277131) q[3];
sx q[3];
rz(-2.5057) q[3];
sx q[3];
rz(-2.2896374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2028929) q[2];
sx q[2];
rz(-0.6627658) q[2];
sx q[2];
rz(-1.1104442) q[2];
rz(-2.2796196) q[3];
sx q[3];
rz(-1.6317261) q[3];
sx q[3];
rz(-1.8814258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2170169) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(0.4050912) q[0];
rz(-2.8054667) q[1];
sx q[1];
rz(-1.7205709) q[1];
sx q[1];
rz(2.1902671) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4753805) q[0];
sx q[0];
rz(-2.6586091) q[0];
sx q[0];
rz(-2.5673812) q[0];
rz(-pi) q[1];
rz(-2.4053965) q[2];
sx q[2];
rz(-1.412782) q[2];
sx q[2];
rz(2.406213) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0266708) q[1];
sx q[1];
rz(-1.7975866) q[1];
sx q[1];
rz(2.5610152) q[1];
rz(-pi) q[2];
rz(-1.3780932) q[3];
sx q[3];
rz(-1.2964222) q[3];
sx q[3];
rz(1.5128795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.090791) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(-0.23078272) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(-0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35293216) q[0];
sx q[0];
rz(-3.1020628) q[0];
sx q[0];
rz(2.2094862) q[0];
rz(-0.63356361) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(-1.3379785) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6659426) q[0];
sx q[0];
rz(-1.4650657) q[0];
sx q[0];
rz(-0.08203489) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7510178) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(-0.24402555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2436314) q[1];
sx q[1];
rz(-0.96968953) q[1];
sx q[1];
rz(-2.6578085) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5389493) q[3];
sx q[3];
rz(-1.0726811) q[3];
sx q[3];
rz(2.826988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46447095) q[2];
sx q[2];
rz(-1.9147583) q[2];
sx q[2];
rz(1.202549) q[2];
rz(-0.36457148) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(0.83474368) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637591) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(2.9649576) q[0];
rz(0.44003507) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(-0.96493351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4266708) q[0];
sx q[0];
rz(-1.3765125) q[0];
sx q[0];
rz(-1.705709) q[0];
x q[1];
rz(-2.4201323) q[2];
sx q[2];
rz(-0.74154343) q[2];
sx q[2];
rz(-2.3539345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.392004) q[1];
sx q[1];
rz(-2.4676305) q[1];
sx q[1];
rz(2.1753009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8252402) q[3];
sx q[3];
rz(-1.4973876) q[3];
sx q[3];
rz(-1.0496828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9289916) q[2];
sx q[2];
rz(-1.5233728) q[2];
sx q[2];
rz(-0.0099445899) q[2];
rz(-0.11659226) q[3];
sx q[3];
rz(-0.3392342) q[3];
sx q[3];
rz(2.5998083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-0.56124878) q[0];
sx q[0];
rz(-1.2224226) q[0];
sx q[0];
rz(0.62193459) q[0];
rz(-1.4382582) q[1];
sx q[1];
rz(-0.57241264) q[1];
sx q[1];
rz(-2.4780746) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2885482) q[0];
sx q[0];
rz(-2.2877734) q[0];
sx q[0];
rz(-1.2794446) q[0];
rz(-pi) q[1];
rz(1.8542669) q[2];
sx q[2];
rz(-1.6890235) q[2];
sx q[2];
rz(2.9663939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92163699) q[1];
sx q[1];
rz(-1.2531279) q[1];
sx q[1];
rz(-2.2457473) q[1];
rz(-2.8264753) q[3];
sx q[3];
rz(-0.2956008) q[3];
sx q[3];
rz(2.78418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62349391) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(-2.7790879) q[2];
rz(0.47719657) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(-1.1227192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057673205) q[0];
sx q[0];
rz(-0.73129439) q[0];
sx q[0];
rz(0.90173632) q[0];
rz(0.67507356) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(-0.14428446) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77567277) q[0];
sx q[0];
rz(-2.7808041) q[0];
sx q[0];
rz(1.3350639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7117908) q[2];
sx q[2];
rz(-0.24082213) q[2];
sx q[2];
rz(0.69320971) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54507885) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(-1.2596115) q[1];
rz(0.38548174) q[3];
sx q[3];
rz(-2.6691438) q[3];
sx q[3];
rz(0.89490283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6361864) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(2.9409161) q[2];
rz(-1.1096654) q[3];
sx q[3];
rz(-0.4626914) q[3];
sx q[3];
rz(2.5698575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5889482) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(0.48925346) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(0.092838661) q[2];
sx q[2];
rz(-1.1576817) q[2];
sx q[2];
rz(0.67410034) q[2];
rz(3.0808385) q[3];
sx q[3];
rz(-2.698632) q[3];
sx q[3];
rz(-2.5198577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
